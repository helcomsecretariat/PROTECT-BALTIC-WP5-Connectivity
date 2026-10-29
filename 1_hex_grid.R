library(terra)
library(sf)

source("_config.R")
if(dir.exists(hex_grid_dir())){stop("Hexagon file of that size exists already")}

grid <- rast(grid_file)

# Hexagons over the grid extent (st_bbox on the raster keeps its CRS)
grid_hex <- st_make_grid(
  st_as_sfc(st_bbox(grid)),
  cellsize = hex_cellsize,
  square = FALSE) |>
  st_sf() |>
  vect()

grid_hex$hexID <- seq_len(nrow(grid_hex))

# Keep hexagons with at least one marine cell (max of all-NA is -Inf)
prop <- extract(grid, grid_hex, fun = max, na.rm = TRUE, exact = FALSE, ID = FALSE)
grid_hex$prop_grids <- prop[[1]]
grid_hex$valid <- is.finite(grid_hex$prop_grids)

cat("Total hexagons:", nrow(grid_hex), "\n")
cat("Valid (marine) hexagons:", sum(grid_hex$valid), "\n")

grid_hex_valid <- grid_hex[grid_hex$valid, ]

out_dir <- hex_grid_dir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# hexID raster in EPSG:4326, for fast particle -> hexagon lookup
hex_id_raster <- rasterize(grid_hex_valid, grid, field = "hexID")
hex_id_wgs84 <- project(hex_id_raster, "EPSG:4326", method = "near")
writeRaster(hex_id_wgs84, file.path(out_dir, "hex_id_wgs84.tif"), overwrite = TRUE)

# Particle release points per hexagon
pts <- read.csv("inputs/points_all.csv")
pts_v <- vect(pts, geom = c("lon", "lat"), crs = "EPSG:4326")
pts_hex <- terra::extract(hex_id_wgs84, pts_v)[, 2]
counts <- table(as.integer(pts_hex))   # integer, so 100000 isn't named "1e+05"
grid_hex_valid$n_particles <- as.integer(counts[as.character(grid_hex_valid$hexID)])
grid_hex_valid$n_particles[is.na(grid_hex_valid$n_particles)] <- 0L

# Fraction of the maximum particles (one per 250 m cell) a hexagon could hold
hex_area_m2 <- expanse(grid_hex_valid, unit = "m")[1]
theoretical_max_particles <- hex_area_m2 / (250 * 250)
grid_hex_valid$frac_of_max <- grid_hex_valid$n_particles / theoretical_max_particles

cat("Hexagon area:", hex_area_m2, "m^2 -> theoretical max particles per hexagon:",
    round(theoretical_max_particles, 1), "\n")
cat("n_particles summary:\n")
print(summary(grid_hex_valid$n_particles))
saveRDS(theoretical_max_particles, file.path(out_dir, "theoretical_max_particles.rds"))

writeVector(grid_hex_valid, file.path(out_dir, "hex_grid.gpkg"), overwrite = TRUE)

# ---- Sea-travel network: touching marine hexagons, centroid distance (km) ----
nb  <- relate(grid_hex_valid, relation = "touches", pairs = TRUE)
nb  <- nb[nb[, 1] < nb[, 2], , drop = FALSE]   # each pair once
xy  <- crds(centroids(grid_hex_valid))
adj <- data.frame(
  from_hex = grid_hex_valid$hexID[nb[, 1]],
  to_hex   = grid_hex_valid$hexID[nb[, 2]],
  dist_km  = sqrt(rowSums((xy[nb[, 1], , drop = FALSE] - xy[nb[, 2], , drop = FALSE])^2)) / 1000)
write.csv(adj, file.path(out_dir, "hex_adjacency.csv"), row.names = FALSE)
cat("Sea-travel network:", nrow(adj), "neighbour edges (",
    round(2 * nrow(adj) / nrow(grid_hex_valid), 1), "per hexagon )\n")

plot(grid_hex_valid, "prop_grids")
