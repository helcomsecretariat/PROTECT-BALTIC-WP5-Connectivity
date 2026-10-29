library(terra)
library(dplyr)
library(arrow)
library(gtools)
library(tidyr)
library(data.table)

source("_config.R")

# Passive archetypes
arc <- read.csv("WP5_archetypes.csv")
arc <- filter(arc, dispersal_mode != "active")

# One matrix per unique (drift duration, season)
arc_distinct <-
  arc |>
  select(drift_duration_average, season) |>
  distinct()

arc_distinct$archetype_matrix <- paste0("ARCH-", 1:nrow(arc_distinct))
arc <- left_join(arc, arc_distinct, by = join_by(drift_duration_average, season))

arc_names <- paste0("ARCH-", 1:nrow(arc_distinct))
# Species <-> archetype assignments (delete the file to refresh)
if (!file.exists("WP5_archetype_assignments.csv"))
  write.csv(arc, "WP5_archetype_assignments.csv", row.names = FALSE)

# Multi-session claims
started_dir <- "_run_log/archetypes/started"
dir.create(started_dir, recursive = TRUE, showWarnings = FALSE)

if (length(exclude_dates))
  cat("Excluding", length(exclude_dates), "flagged date(s):",
      paste(exclude_dates, collapse = ", "), "\n")

years <- c("2023", "2024", "2025")

# One date's pairs table (NULL if absent)
read_pairs <- function(p, dt) {
  f <- paste0(matrix_dir_at(p, dt), "hex_connectivity_pairs.csv")
  if (file.exists(f)) as.data.frame(read_csv_arrow(f)) else NULL
}

# Average the per-date matrices: mean over years within each slot, then over slots
for (i in 1:nrow(arc_distinct)) {
  p  <- arc_distinct$drift_duration_average[i]
  out_dir <- archetype_dir(p, arc_names[i])

  # Skip if already built (.rds written last)
  if (file.exists(paste0(out_dir, "hex_connectivity_matrix.rds"))) {
    cat(arc_names[i], " (PLD ", p, "): already produced - skipping\n", sep = "")
    next
  }

  # Claim (atomic)
  claim <- file.path(started_dir, arc_names[i])
  if (!dir.create(claim, showWarnings = FALSE)) {
    cat(arc_names[i], " (PLD ", p, "): claimed by another session - skipping\n", sep = "")
    next
  }

  mi <- as.integer(unlist(strsplit(arc_distinct$season[i], ",")))
  slots <- paste0(rep(sprintf("%02d", mi), each = 2), "_", c("01", "15"))  # MM_DD

  # Usable dates per slot: present and not excluded
  usable <- lapply(slots, function(s) {
    dts <- paste0(years, "_", s)
    present <- file.exists(paste0(matrix_dir_at(p, dts), "hex_connectivity_pairs.csv"))
    list(ok = dts[present & !(dts %in% exclude_dates)],
         missing_notexcl = dts[!present & !(dts %in% exclude_dates)])
  })

  # Skip if any date is missing and not excluded
  miss <- unlist(lapply(usable, `[[`, "missing_notexcl"))
  if (length(miss)) {
    cat(arc_names[i], " (PLD ", p, "): ", length(miss),
        " matrix/matrices not produced and NOT flagged: ", paste(miss, collapse = ", "),
        " - skipping (generate them, or add to exclude_dates if bad)\n", sep = "")
    unlink(claim, recursive = TRUE); next
  }

  # Skip if any slot has no usable date
  n_ok <- vapply(usable, function(x) length(x$ok), integer(1))
  if (any(n_ok == 0)) {
    cat(arc_names[i], " (PLD ", p, "): season slot(s) with no usable matrices (",
        paste(slots[n_ok == 0], collapse = ", "), ") - skipping\n", sep = "")
    unlink(claim, recursive = TRUE); next
  }

  ok <- tryCatch({
    # Weighted sum: each date weighted 1 / (years in slot x number of slots)
    combined <- NULL
    t0 <- Sys.time()
    cat(arc_names[i], " (PLD ", p, "): averaging ", sum(n_ok), " date-matrices over ",
        length(usable), " slots\n", sep = "")
    for (si in seq_along(usable)) {
      x    <- usable[[si]]
      w    <- 1 / (length(x$ok) * length(usable))
      slot <- rbindlist(lapply(x$ok, function(dt) setDT(read_pairs(p, dt))))
      slot <- slot[, .(n = sum(n) * w), by = .(origin_hex, dest_hex)]
      combined <- if (is.null(combined)) slot else
        rbindlist(list(combined, slot))[, .(n = sum(n)), by = .(origin_hex, dest_hex)]
      rm(slot); invisible(gc())
      cat(sprintf("  %s slot %d/%d (%s, %d yr) - %s pairs so far - %.1f min\n",
                  arc_names[i], si, length(usable), slots[si], length(x$ok),
                  format(nrow(combined), big.mark = ","),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      flush.console()
    }

    # Pairs CSV (hexIDs as integer, so fwrite doesn't write "1e+05")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    combined[, `:=`(origin_hex = as.integer(origin_hex), dest_hex = as.integer(dest_hex))]
    fwrite(combined, paste0(out_dir, "hex_connectivity_pairs.csv"))

    # Sparse matrix .rds, written last (marks the archetype finished)
    hex_ids <- sort(unique(c(combined$origin_hex, combined$dest_hex)))
    mat <- Matrix::sparseMatrix(
      i = match(combined$origin_hex, hex_ids),
      j = match(combined$dest_hex, hex_ids),
      x = combined$n,
      dims = c(length(hex_ids), length(hex_ids)),
      dimnames = list(as.character(hex_ids), as.character(hex_ids)))
    saveRDS(mat, paste0(out_dir, "hex_connectivity_matrix.rds"))

    cat(arc_names[i], "(PLD", p, "): averaged", length(slots), "season slots (",
        sum(n_ok), "date-matrices ),", nrow(combined), "pairs, dim", length(hex_ids), "\n")
    TRUE
  }, error = function(e) {
    cat(arc_names[i], " (PLD ", p, "): ERROR - ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  if (!ok) unlink(claim, recursive = TRUE)   # release claim for retry
}

