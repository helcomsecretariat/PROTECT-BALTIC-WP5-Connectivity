library(terra)
library(igraph)

# Source and sink strength per species (run_log do_run == 1, passive archetype).
# Multi-session: species claimed in _run_log/strength/started, marked in finished.
source("_config.R")
source("_functions/species_archetype.R")
source("_functions/strength.R")

spec_df   <- read.csv("WP5_species_list.csv")
spec_arch <- species_archetype()
run_log   <- read.csv("run_log.csv")
species   <- run_log$scientific_name[run_log$do_run == 1]

started_dir  <- "_run_log/strength/started"
finished_dir <- "_run_log/strength/finished"
dir.create(started_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(finished_dir, recursive = TRUE, showWarnings = FALSE)

# Static inputs
grid_orig  <- rast(grid_file)
hex        <- vect(file.path(hex_grid_dir(), "hex_grid.gpkg"))
hex$hexID  <- as.integer(hex$hexID)   # integer, so 100000 isn't "1e+05"
hex_xy     <- crds(centroids(hex))
hex_coords <- data.frame(hexID = hex$hexID, x = hex_xy[, 1], y = hex_xy[, 2])
enough_particles <- hex$frac_of_max >= min_particle_frac
sea_g <- if (distance_method == "sea") sea_graph() else NULL
cat("Distance method:", distance_method,
    if (!is.null(sea_g)) paste0(" (", length(V(sea_g)), " hex nodes, ",
                                length(E(sea_g)), " edges)") else "",
    "| strength distance factor:", strength_distance, "mean distance\n")

# Hexagons with any habitat cell in a 1/NA raster
habitat_hexes <- function(r)
  is.finite(extract(r, hex, fun = max, na.rm = TRUE, exact = FALSE, ID = FALSE)[[1]])

# Strength map PNG (a constant field is drawn as one class)
plot_strength_map <- function(hex_bg, hex_val, field, breaks, title, out_png, pal = "gyr") {
  ncol <- if (length(breaks) == 1) breaks else length(breaks) - 1
  vals <- terra::values(hex_val)[[field]]
  png(out_png, width = 3000, height = 3000)
  on.exit(dev.off())
  plot(hex_bg, col = "grey92", border = NA, legend = FALSE, axes = FALSE, main = title)
  if (length(unique(vals[is.finite(vals)])) < 2) {
    u <- unique(vals[is.finite(vals)])
    plot(hex_val, col = map.pal(pal, 1), border = NA, add = TRUE)
    legend("topright", bty = "n", cex = 2,
           legend = paste0(field, " = ", if (length(u)) format(u) else "NA", " (constant)"))
  } else {
    plot(hex_val, field, type = "interval", breaks = breaks, col = map.pal(pal, ncol),
         border = NA, add = TRUE)
  }
}

for (k in seq_along(species)) {
  spec <- species[k]
  cat("\n[", k, "/", length(species), "] ", spec, "\n", sep = "")

  # ---- Skip checks ----
  if (file.exists(file.path(finished_dir, spec))) {
    cat(spec, "- already finished, skipping\n"); next
  }
  sa <- spec_arch[spec_arch$scientific_name == spec, ]
  if (nrow(sa) == 0 || is.na(sa$archetype_matrix[1])) {
    cat(spec, "- no passive archetype matrix, skipping\n"); next
  }
  pld  <- sa$pld[1]
  arch <- sa$archetype_matrix[1]
  pairs_file <- paste0(archetype_dir(pld, arch), "hex_connectivity_pairs.csv")
  if (!file.exists(pairs_file)) {
    cat(spec, "- archetype matrix", arch, "not built yet, skipping (", pairs_file, ")\n"); next
  }
  spec_group  <- spec_df$species_group[spec_df$scientific_name == spec]
  origin_file <- paste0("inputs/sdms_PBS_EFH/", spec, ".tif")
  dest_file   <- paste0("inputs/sdms/", spec_group, "/", sdm_type, "_tif/", spec, ".tif")
  has_origin_override <- file.exists(origin_file)
  if (!file.exists(dest_file) && !has_origin_override) {
    cat(spec, "- no SDM found, skipping\n"); next
  }

  # ---- Claim (atomic) ----
  claim <- file.path(started_dir, spec)
  if (!dir.create(claim, showWarnings = FALSE)) {
    cat(spec, "- claimed by another session, skipping\n"); next
  }

  ok <- tryCatch({
    # ---- Habitat: origin = PBS_EFH if present, destination = full SDM (== 4) ----
    if (has_origin_override)
      origin_hab <- ifel(rast(origin_file) == 1, 1, NA)
    if (file.exists(dest_file)) {
      dest_hab <- ifel(rast(dest_file) == 4, 1, 0) * grid_orig
      dest_hab[dest_hab == 0] <- NA
    } else {
      cat(spec, "- no binary_confidence SDM; using the PBS_EFH map for destination too\n")
      dest_hab <- origin_hab
    }
    cat(spec, if (has_origin_override) "- origin = PBS_EFH spawning, destination = full SDM\n"
              else "- origin = destination = full SDM\n")

    dest_has   <- habitat_hexes(dest_hab)
    origin_has <- if (has_origin_override) habitat_hexes(origin_hab) else dest_has
    dest_hex_ids   <- hex$hexID[dest_has & enough_particles]
    origin_hex_ids <- intersect(hex$hexID[origin_has & enough_particles], dest_hex_ids)
    cat("Destination (full) hexagons:", length(dest_hex_ids),
        " | origin (spawning) hexagons:", length(origin_hex_ids), "\n")

    # ---- Edges: origin -> destination, excluded hexagons dropped ----
    pairs <- data.table::fread(pairs_file, select = c("origin_hex", "dest_hex", "n"))
    pairs[, `:=`(origin_hex = as.integer(origin_hex),
                 dest_hex   = as.integer(dest_hex))]
    n_all <- nrow(pairs)
    pairs <- pairs[origin_hex %in% origin_hex_ids & dest_hex %in% dest_hex_ids &
                   !(origin_hex %in% exclude_hexagons) & !(dest_hex %in% exclude_hexagons)]
    cat("Edges before habitat filter:", n_all, " after:", nrow(pairs), "\n")
    data.table::setnames(pairs, "n", "weight")
    invisible(gc())

    # ---- Source and sink strength (nodes = destination hexagons) ----
    strength_df <- compute_strength(pairs, dest_hex_ids, hex_coords, sea_g,
                                    origin_ids = origin_hex_ids)
    rm(pairs); invisible(gc())
    for (metric in c("source_strength", "sink_strength")) {
      cat(metric, "summary (amount x variety-fraction x mean distance):\n")
      print(summary(strength_df[[metric]]))
    }

    # ---- Tukey filter: values above Q3 + k x IQR (of nonzero values) replaced ----
    for (metric in c("source_strength", "sink_strength")) {
      x     <- strength_df[[metric]]
      pos   <- x[x > 0]
      q     <- quantile(pos, c(0.25, 0.75), names = FALSE)
      fence <- if (length(pos)) q[2] + outlier_iqr_mult * (q[2] - q[1]) else Inf
      out   <- x > fence
      repl  <- switch(outlier_replace,
                      zero = 0,
                      max  = max(x[!out]),
                      stop("outlier_replace must be \"zero\" or \"max\""))
      strength_df[[paste0(metric, "_tukey")]] <- ifelse(out, repl, x)
      cat(metric, " outlier fence (Q3 + ", outlier_iqr_mult, " x IQR): ", round(fence), " - ",
          sum(out), " hexagon(s) set to ", round(repl), " (", outlier_replace, ") in ", metric, "_tukey",
          if (any(out)) paste0(": hexID ", paste(strength_df$hexID[out], collapse = ", ")), "\n", sep = "")
    }
    first <- c("hexID", "source_strength", "source_strength_tukey", "sink_strength", "sink_strength_tukey")
    strength_df <- strength_df[, c(first, setdiff(names(strength_df), first))]

    # ---- Outputs: csv, gpkg, source/sink maps ----
    out_dir <- strength_dir(spec, pld)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_stem <- paste0(out_dir, spec, "_strength")
    write.csv(strength_df, paste0(out_stem, ".csv"), row.names = FALSE)
    hex_out <- merge(hex, strength_df, by = "hexID", all.x = FALSE)
    writeVector(hex_out, paste0(out_stem, ".gpkg"), overwrite = TRUE)
    for (metric in c("source", "sink")) {
      plot_strength_map(hex, hex_out, paste0(metric, "_strength_tukey"), 10,
                        paste0(spec, " - ", metric, " strength, outliers filtered (",
                               if (metric == "source") "export" else "import",
                               " x variety-fraction x mean distance)\n",
                               "PLD=", pld, "d, ", hex_cellsize / 1000, "km hexagons"),
                        paste0(out_dir, spec, "_", metric, "_strength_map.png"))
    }
    cat("wrote", out_dir, paste0(spec, "_strength.{csv,gpkg} and _{source,sink}_strength_map.png\n"))
    TRUE
  }, error = function(e) {
    cat(spec, "- ERROR:", conditionMessage(e), "\n")
    FALSE
  })

  if (ok) file.create(file.path(finished_dir, spec))   # mark done
  else    unlink(claim, recursive = TRUE)                # release claim for retry
}
cat("\ndone\n")
