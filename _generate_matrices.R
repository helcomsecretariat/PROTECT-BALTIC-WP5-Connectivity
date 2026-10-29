library(terra)
library(dplyr)
library(arrow)
library(gtools)

# Builds the connectivity matrix for the next unclaimed release date
source("_config.R")

dates <- all_release_dates()

if(length(list.files("_run_log/matrices/finished")) == length(dates)){stop("All matrices finished with current configuration")}

# Claim the first unstarted date
l <- list.files("_run_log/matrices/started")
date_miss <- first(dates[dates %in% l == FALSE])
if (is.na(date_miss)) stop("All matrices already started (none left to claim)")

write.csv(NA, paste0("_run_log/matrices/started/", date_miss))
start_date <- date_miss

stopifnot("pld_days must be shorter than the simulation length" = pld_days < sim_days)

# Settlement window (hours since release)
w <- settle_window_hours(pld_days)
window_lower <- w[1]
window_upper <- w[2]
stopifnot("window_lower must be >= 0 (check the window settings in _config.R)" = window_lower >= 0)
mode_lbl <- if (!is.na(pld_buffer)) sprintf("buffer = %g", pld_buffer) else
            if (!is.na(window_days)) sprintf("+/-%g d window", window_days) else "drift min/max window"
cat(sprintf("Drift = %g d, %s. Landing window = %.0f-%.0f h (%.1f-%.1f d)\n",
            pld_days, mode_lbl, window_lower, window_upper, window_lower/24, window_upper/24))

# Parquet time starts at 1 (= 0 h)
pt_lower <- window_lower + 1
pt_upper <- window_upper + 1

l <- parquet_blocks(start_date)
if (length(l) == 0) {
  stop("No parquet blocks for date '", start_date, "' in run2 (", parquet_dir_run2,
       ") or run3 (", parquet_dir_run3, "). Check the working directory (getwd() is '",
       getwd(), "') and the parquet_dir_run2/run3 paths in _config.R.")
}

# hexID raster values in memory, for fast lookup
hex_dir <- hex_grid_dir()
hex_id <- rast(file.path(hex_dir, "hex_id_wgs84.tif"))
hex_vals <- values(hex_id)[, 1]

out_dir <- matrix_dir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pair_list <- vector("list", length(l))

for (i in seq_along(l)) {

  f <- l[i]

  # Release rows and settlement-window rows only
  df <- open_dataset(f, format = "parquet") %>%
    filter(time == 1 | (time >= pt_lower & time <= pt_upper)) %>%
    select(particle, time, lat, lon) %>%
    collect()

  df$hours <- df$time - 1
  df <- filter(df, !is.na(lat), !is.na(lon))
  df_origin <- filter(df, hours == 0)
  df_window <- filter(df, hours >= window_lower & hours <= window_upper)

  # Origin and destination hexagons
  df_origin$origin_hex <- hex_vals[cellFromXY(hex_id, cbind(df_origin$lon, df_origin$lat))]
  df_origin <- filter(df_origin, !is.na(origin_hex))

  df_window$dest_hex <- hex_vals[cellFromXY(hex_id, cbind(df_window$lon, df_window$lat))]
  df_window <- filter(df_window, !is.na(dest_hex))
  # Each particle counts once per destination hexagon
  df_window <- distinct(df_window, particle, dest_hex)

  pairs <- inner_join(df_origin[, c("particle", "origin_hex")],
                      df_window[, c("particle", "dest_hex")],
                      by = "particle")

  pair_list[[i]] <- pairs %>% count(origin_hex, dest_hex, name = "n")

  cat("Iteration", i, "complete -", nrow(pair_list[[i]]), "pairs -")

  rm(df, df_origin, df_window, pairs)
  gc()

}
cat("\n")

all_pairs <- bind_rows(pair_list) %>%
  group_by(origin_hex, dest_hex) %>%
  summarise(n = sum(n), .groups = "drop")

write.csv(all_pairs, paste0(out_dir, "hex_connectivity_pairs.csv"), row.names = FALSE)

# Sparse origin x destination matrix
hex_ids <- sort(unique(c(all_pairs$origin_hex, all_pairs$dest_hex)))
mat <- Matrix::sparseMatrix(
  i = match(all_pairs$origin_hex, hex_ids),
  j = match(all_pairs$dest_hex, hex_ids),
  x = all_pairs$n,
  dims = c(length(hex_ids), length(hex_ids)),
  dimnames = list(as.character(hex_ids), as.character(hex_ids)))
saveRDS(mat, paste0(out_dir, "hex_connectivity_matrix.rds"))

cat("Done. Matrix dim:", nrow(mat), "x", ncol(mat),
    " total particle-pairs:", sum(mat), "\n")

# Mark date finished
write.csv(NA, paste0("_run_log/matrices/finished/", start_date))
