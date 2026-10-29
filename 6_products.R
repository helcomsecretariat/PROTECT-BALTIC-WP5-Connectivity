# Per-species products: connectivity gpkg and source/sink hex map PNGs
library(terra)
source("_config.R")
source("_functions/species_archetype.R")

metrics  <- c("source_strength", "sink_strength")
# Connectivity gpkg fields, in order
gpkg_fields <- c("hexID",
                 "source_strength", "source_strength_tukey", "source_settlers",
                 "source_n_dest", "source_mean_dist",
                 "sink_strength", "sink_strength_tukey", "sink_settlers",
                 "sink_n_origin", "sink_mean_dist")
prod_dir <- "outputs/products"
for (d in c(paste0(metrics, "_hex_png"), "connectivity_gpkg"))
  dir.create(file.path(prod_dir, d), recursive = TRUE, showWarnings = FALSE)
prod_file <- function(m, kind, spec, ext)
  file.path(prod_dir, paste0(m, "_", kind), paste0(spec, ".", ext))

# Species with a passive archetype
spec_arch <- species_archetype()
spec_arch <- spec_arch[!is.na(spec_arch$pld), ]

n_out <- setNames(rep(0L, 2), c("gpkg", "hex_png"))
for (k in seq_len(nrow(spec_arch))) {
  spec    <- spec_arch$scientific_name[k]
  pld     <- spec_arch$pld[k]
  net_dir <- strength_dir(spec, pld)

  gpkg <- paste0(net_dir, spec, "_strength.gpkg")
  if (!file.exists(gpkg)) next
  v <- vect(gpkg)
  miss <- setdiff(gpkg_fields, names(v))
  if (length(miss)) {
    cat(spec, "- missing field(s)", paste(miss, collapse = ", "),
        "(older step-5 output; re-run step 5), skipping\n")
    next
  }

  writeVector(v[, gpkg_fields], prod_file("connectivity", "gpkg", spec, "gpkg"), overwrite = TRUE)
  n_out["gpkg"] <- n_out["gpkg"] + 1L

  # Copy the step-5 maps
  for (m in metrics) {
    hex_png <- paste0(net_dir, spec, "_", m, "_map.png")
    if (file.exists(hex_png)) {
      file.copy(hex_png, prod_file(m, "hex_png", spec, "png"), overwrite = TRUE)
      n_out["hex_png"] <- n_out["hex_png"] + 1L
    }
  }
  cat(spec, "- done\n")
}

cat("\nWrote", n_out["gpkg"], "connectivity gpkg(s) and copied", n_out["hex_png"],
    "hex png(s) (source + sink) to", prod_dir, "\n")
