# Source and sink strength per hexagon:
#   source_strength = settlers sent     x (n_dest / (N - 1))  x mean distance
#   sink_strength   = settlers received x (n_origin / max_in) x mean distance
# Self-recruitment excluded; distances in km; max_in = origins other than the hexagon.
# Mean distance is settler-weighted or plain (config strength_distance).

# Sea distance (km) per pair via shortest path, straight-line where disconnected
sea_distance <- function(sea_g, from_hex, to_hex, fx, fy, tx, ty, chunk = 500) {
  out   <- sqrt((fx - tx)^2 + (fy - ty)^2) / 1000
  verts <- as.integer(names(igraph::V(sea_g)))
  fv <- match(from_hex, verts); tv <- match(to_hex, verts)
  ok <- which(!is.na(fv) & !is.na(tv))
  if (length(ok)) {
    uo <- unique(fv[ok])
    grp <- split(ok, (match(fv[ok], uo) - 1L) %/% chunk)   # chunks of origins, to limit memory
    w   <- igraph::E(sea_g)$dist_km
    for (idx in grp) {
      co <- unique(fv[idx]); cd <- unique(tv[idx])
      D  <- igraph::distances(sea_g, v = co, to = cd, weights = w)
      sp <- D[cbind(match(fv[idx], co), match(tv[idx], cd))]
      fin <- is.finite(sp)
      out[idx[fin]] <- sp[fin]
    }
  }
  out
}

# pairs: origin_hex, dest_hex, weight; node_ids: all nodes; coords: hexID, x, y (m)
compute_strength <- function(pairs, node_ids, coords, sea_g = NULL, origin_ids = node_ids) {
  ids  <- as.integer(node_ids)
  N    <- length(ids)
  orig <- intersect(as.integer(origin_ids), ids)

  # Off-diagonal, positive edges from origins to nodes
  e <- data.table::data.table(o = as.integer(pairs$origin_hex), d = as.integer(pairs$dest_hex),
                              w = pairs$weight)
  e <- e[o != d & w > 0 & o %in% orig & d %in% ids]

  # Drop excluded hexagons; sum duplicate pairs
  if (exists("exclude_hexagons") && length(exclude_hexagons))
    e <- e[!(o %in% exclude_hexagons) & !(d %in% exclude_hexagons)]
  e <- e[, .(w = sum(w)), by = .(o, d)]

  # Edge distance (km)
  ci <- match(e$o, coords$hexID); cj <- match(e$d, coords$hexID)
  if (exists("distance_method") && identical(distance_method, "sea") && !is.null(sea_g)) {
    e[, dkm := sea_distance(sea_g, o, d, coords$x[ci], coords$y[ci], coords$x[cj], coords$y[cj])]
  } else {
    e[, dkm := sqrt((coords$x[ci] - coords$x[cj])^2 + (coords$y[ci] - coords$y[cj])^2) / 1000]
  }
  rm(ci, cj)

  # One metric from outgoing ("o") or incoming ("d") edges
  strength <- function(by_col, max_variety) {
    agg <- e[, .(settlers = sum(w),
                 fluxdist = sum(w * dkm),
                 variety  = .N,
                 distsum  = sum(dkm)),
             by = by_col]
    m <- match(ids, agg[[by_col]])
    pick <- function(x) { v <- x[m]; v[is.na(v)] <- 0; v }
    settlers <- pick(agg$settlers); fluxdist <- pick(agg$fluxdist)
    variety  <- pick(agg$variety);  distsum  <- pick(agg$distsum)
    mean_dist  <- ifelse(settlers > 0, fluxdist / settlers, 0)  # settler-weighted
    dist_plain <- ifelse(variety  > 0, distsum  / variety,  0)  # unweighted
    prop_variety <- variety / max_variety
    prop_variety[max_variety <= 0] <- 0
    dist_factor <- if (exists("strength_distance") && identical(strength_distance, "plain"))
                     dist_plain else mean_dist
    list(strength = settlers * prop_variety * dist_factor,
         settlers = settlers, variety = as.integer(variety),
         mean_dist = mean_dist, dist_plain = dist_plain)
  }
  src <- strength("o", N - 1)
  snk <- strength("d", length(orig) - (ids %in% orig))

  data.frame(
    hexID             = ids,
    source_strength   = src$strength,
    sink_strength     = snk$strength,
    source_settlers   = src$settlers,
    source_n_dest     = src$variety,
    source_mean_dist  = src$mean_dist,
    source_dist_plain = src$dist_plain,
    sink_settlers     = snk$settlers,
    sink_n_origin     = snk$variety,
    sink_mean_dist    = snk$mean_dist,
    sink_dist_plain   = snk$dist_plain
  )
}
