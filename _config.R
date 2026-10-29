# Shared settings and path helpers, sourced by all pipeline scripts.

# OpenDrift parquet runs (run3 re-runs override run2 per block)
parquet_dir_run2 <- "D:/pb_opendrift_backup/outputs/run2_simtime_60days_calctime_60min/"
parquet_dir_run3 <- "D:/pb_opendrift_backup/outputs/run3_simtime_60days_calctime_60min/"
grid_file    <- "inputs/HELCOM_WP3_grid.tif"  # 250 m marine grid (EPSG:3035)
land_mask_file <- "inputs/WP3_land_mask.tif"  # 0 = sea, 1 = land
hex_cellsize <- 3000          # hexagon size (m)
#start_date   <- "2024_08_01"
pld_days     <- 40            # drift duration (days)
# Settlement window, first non-NA mode:
pld_buffer   <- NA            # fraction of drift: [1 - b, 1 + b] x drift
window_days  <- NA            # +/- days around drift
                              # both NA: drift_duration_min/max from WP5_archetypes.csv
sim_days     <- 60            # simulation length (days)
output_step_hours <- 3        # parquet output interval (h)

# Species habitat
sdm_type          <- "binary_confidence"   # inputs/sdms/<group>/<sdm_type>_tif/<species>.tif
min_particle_frac <- 0.0                  # drop hexagons below this fraction of max particles

# Distance between hexagons: "euclidean" or "sea" (shortest path via hex adjacency)
distance_method <- "sea"
stopifnot(distance_method %in% c("euclidean", "sea"))

# Strength distance factor: "weighted" (settler-weighted mean) or "plain" (unweighted mean)
strength_distance <- "weighted"
stopifnot(strength_distance %in% c("weighted", "plain"))

# Release dates excluded from the archetype averages, e.g. "2023_07_01"
exclude_dates <- sort(unique(c(

)))

# Hexagons excluded from the strength metrics (NA = none)

exclude_hexagons <- NA

# Tukey outlier filter: fence = Q3 + outlier_iqr_mult x IQR of nonzero values
outlier_iqr_mult <- 3
outlier_replace  <- "max"   # "zero" or "max" (highest non-outlier value)

# Hex grid folder
hex_grid_dir <- function() paste0("inputs/hex_grid/", hex_cellsize / 1000, "km")

# Sea-travel igraph (hexID vertex names, dist_km edge weights)
sea_graph <- function() {
  f <- file.path(hex_grid_dir(), "hex_adjacency.csv")
  if (!file.exists(f))
    stop("Sea adjacency network not found: ", f, " - (re)run 1_hex_grid.R.")
  adj <- data.table::fread(f, data.table = FALSE)
  adj$from_hex <- as.integer(adj$from_hex)   # integer, so 100000 isn't "1e+05"
  adj$to_hex   <- as.integer(adj$to_hex)
  igraph::graph_from_data_frame(adj, directed = FALSE)
}

# Parquet block paths for a date, run3 preferred per block
parquet_blocks <- function(date) {
  r2 <- file.path(parquet_dir_run2, date, "outputs_txt")
  r3 <- file.path(parquet_dir_run3, date, "outputs_txt")
  b2 <- if (dir.exists(r2)) list.files(r2, pattern = "parquet$") else character(0)
  b3 <- if (dir.exists(r3)) list.files(r3, pattern = "parquet$") else character(0)
  blocks <- gtools::mixedsort(unique(c(b3, b2)))
  file.path(ifelse(blocks %in% b3, r3, r2), blocks)
}

# All 72 release dates (1st and 15th of each month, 2023-2025)
all_release_dates <- function() paste0(
  rep(c("2023", "2024", "2025"), each = 24), "_",
  rep(sprintf("%02d", 1:12), each = 2), "_",
  rep(c("01", "15"), 36))

# Release dates for a season string of months, e.g. "4,5,6"
season_dates <- function(season) {
  mi <- as.integer(unlist(strsplit(season, ",")))
  mm <- rep(sprintf("%02d", mi), each = 2)
  dd <- rep(c("01", "15"), length(mi))
  paste0(rep(c("2023", "2024", "2025"), each = length(dd)), "_", rep(mm, 3), "_", rep(dd, 3))
}

# drift_duration_min/max (days) for a drift duration, from WP5_archetypes.csv
.drift_minmax_cache <- NULL
drift_minmax <- function(pld) {
  if (is.null(.drift_minmax_cache)) {
    a <- read.csv("WP5_archetypes.csv")
    a <- unique(a[!is.na(a$drift_duration_average),
                  c("drift_duration_average", "drift_duration_min", "drift_duration_max")])
    .drift_minmax_cache <<- a
  }
  row <- .drift_minmax_cache[.drift_minmax_cache$drift_duration_average == pld, ]
  if (nrow(row) == 0)
    stop("No drift_duration_average == ", pld, " in WP5_archetypes.csv")
  if (nrow(unique(row[, c("drift_duration_min", "drift_duration_max")])) > 1)
    stop("Ambiguous drift_duration_min/max for drift ", pld, " in WP5_archetypes.csv")
  c(min = row$drift_duration_min[1], max = row$drift_duration_max[1])
}

# Settlement window [lower, upper] in hours since release
settle_window_hours <- function(pld) {
  if (!is.na(pld_buffer))      c(pld * (1 - pld_buffer), pld * (1 + pld_buffer)) * 24
  else if (!is.na(window_days)) c(pld - window_days,      pld + window_days)      * 24
  else                         { mm <- drift_minmax(pld); c(mm[["min"]], mm[["max"]]) * 24 }
}

# Output path tag: pld40d_buffer20pct / pld40d_win5d / pld40d (min/max)
run_tag <- function(pld) {
  if (!is.na(pld_buffer))       paste0("pld", pld, "d_buffer", pld_buffer * 100, "pct")
  else if (!is.na(window_days)) paste0("pld", pld, "d_win", window_days, "d")
  else                          paste0("pld", pld, "d")
}

# Per-date matrix folder(s)
matrix_dir_at <- function(pld, date) paste0("outputs/matrices/", hex_cellsize / 1000, "km/",
                                            run_tag(pld), "/", date, "/")

matrix_dir <- function() matrix_dir_at(pld_days, start_date)

# Archetype matrix folder, e.g. arch = "ARCH-1"
archetype_dir <- function(pld, arch) paste0("outputs/archetype_matrices/", hex_cellsize / 1000, "km/",
                                            run_tag(pld), "/", arch, "/")

# Per-species strength output folder
strength_dir <- function(spec, pld) paste0("outputs/strength/", hex_cellsize / 1000, "km/", spec, "/",
                                           run_tag(pld), "/")
