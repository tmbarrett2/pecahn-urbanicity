#' Boundary-distance travel-time estimate (internal)
#'
#' @description
#' Fallback used by the progressive search functions when no feature of the
#' requested type is found within the maximum search area. Estimates a
#' travel time by measuring least-cost travel from the community center to the
#' four cardinal edge-midpoints of the final (largest) search bounding box and
#' averaging the finite results. This yields a conservative "at least this far"
#' figure instead of a bare `NA`.
#'
#' @param final_bbox Numeric vector `c(left, bottom, right, top)` of the final
#'   (largest) search bbox actually reached by the search loop.
#' @param center_point An `sf` point (WGS84) for the community center.
#' @param tr_corrected A geo-corrected `gdistance` TransitionLayer, or `NULL`.
#' @param raster_crs CRS of the friction raster.
#'
#' @return Numeric scalar mean travel time (minutes) to the finite edge
#'   midpoints, or `NA_real_` if `tr_corrected` is `NULL` or none are finite.
#'
#' @details
#' The progressive search can expand `final_bbox` far beyond the friction raster,
#' which is only fetched for the regional bbox plus a modest buffer. Edge midpoints
#' that fall outside the raster are unreachable and return non-finite travel times,
#' which would collapse the estimate to `NA` for exactly the very remote communities
#' this fallback exists to serve. The edges are therefore clamped to the friction
#' raster extent (derived from `tr_corrected`) so the probe points always sit on
#' real, reachable cells, yielding a conservative "the nearest feature is at least
#' as far as the edge of our coverage" estimate.
#'
#' For the same reason the *origin* is clamped too. The MAP Global Friction Surface
#' is masked, not truly global, so the community center of a very remote site (deep
#' in GHS-SMOD class 11) can itself fall outside the raster. When that happens
#' `gdistance::costDistance()` has no valid origin cell, every edge probe returns a
#' non-finite time, and the estimate again collapses to `NA` — silently dropping the
#' very point this fallback exists to serve. A clamped *copy* of `center_point` is
#' therefore used as the cost-distance origin (only when the extent looks geographic,
#' mirroring the edge clamping). This is local to the probes inside this function;
#' the caller's `center_point` is never mutated.
#'
#' @keywords internal
#' @noRd
.boundary_distance_estimate <- function(final_bbox, center_point,
                                        tr_corrected, raster_crs) {
  if (is.null(tr_corrected)) return(NA_real_)

  left   <- as.numeric(final_bbox[1])
  bottom <- as.numeric(final_bbox[2])
  right  <- as.numeric(final_bbox[3])
  top    <- as.numeric(final_bbox[4])

  # Clamp the bbox edges to the friction raster extent. The friction surface is
  # built in EPSG:4326 (lon/lat), matching the probe points constructed below, so
  # only clamp when the extent looks geographic; otherwise fall back to the raw
  # (unclamped) bbox to avoid mixing projected units with lon/lat.
  ext <- tryCatch(raster::extent(tr_corrected), error = function(e) NULL)
  looks_geographic <- !is.null(ext) &&
    all(is.finite(c(ext@xmin, ext@xmax, ext@ymin, ext@ymax))) &&
    ext@xmin >= -180 && ext@xmax <= 180 &&
    ext@ymin >= -90  && ext@ymax <= 90
  # The origin defaults to the caller's center; only a clamped *copy* is used below
  # (never the caller's object). For non-geographic extents nothing is clamped, so
  # the origin and edges keep their pre-existing (unclamped) behavior.
  origin_point <- center_point

  if (looks_geographic) {
    # Inset by a tiny fraction so probe points sit strictly inside the grid.
    pad_x <- (ext@xmax - ext@xmin) * 1e-6
    pad_y <- (ext@ymax - ext@ymin) * 1e-6
    left   <- max(left,   ext@xmin + pad_x)
    right  <- min(right,  ext@xmax - pad_x)
    bottom <- max(bottom, ext@ymin + pad_y)
    top    <- min(top,    ext@ymax - pad_y)

    # Clamp a copy of the origin onto the same extent so cost-distance always has a
    # valid origin cell even when the community center lies outside the (masked)
    # friction surface. For an in-bounds center this is a no-op, so the estimate is
    # unchanged; for an off-raster center it recovers a conservative value instead
    # of collapsing to NA. center_point itself is left untouched.
    cxy   <- sf::st_coordinates(sf::st_transform(center_point, crs = 4326))[1, 1:2]
    c_lon <- min(max(cxy[1], ext@xmin + pad_x), ext@xmax - pad_x)
    c_lat <- min(max(cxy[2], ext@ymin + pad_y), ext@ymax - pad_y)
    origin_point <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(c_lon, c_lat)), crs = 4326))

    # If the clamped origin still lands on a masked / disconnected friction cell
    # (a nodata pixel INSIDE the extent -- e.g. a water mask in the MAP surface),
    # snap it to the nearest reachable cell so cost-distance has a valid start.
    # No-op when the origin already sits on a valid cell.
    origin_point <- .nearest_valid_origin(origin_point, tr_corrected)
  }

  center_lon <- (left + right) / 2
  center_lat <- (top + bottom) / 2

  # Midpoints of the four edges: North, East, South, West.
  edge_pts <- list(
    c(center_lon, top),
    c(right,      center_lat),
    c(center_lon, bottom),
    c(left,       center_lat)
  )

  times <- vapply(edge_pts, function(xy) {
    pt <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(xy), crs = 4326))
    tryCatch(
      calculate_travel_time(origin_point, pt, tr_corrected, raster_crs),
      error = function(e) NA_real_
    )
  }, numeric(1))

  finite_times <- times[is.finite(times)]
  if (length(finite_times) == 0) return(NA_real_)
  mean(finite_times)
}

#' Snap a point onto the nearest reachable friction cell (internal)
#'
#' @description
#' When the (extent-clamped) origin still lands on a masked / nodata friction
#' cell that lies *inside* the raster extent, `gdistance::costDistance()` has no
#' valid origin node and every probe collapses to `NA`. This finds the nearest
#' cell that is actually connected in the transition layer (non-zero conductance)
#' and returns an `sf` point on it. It is a no-op when the origin already sits on
#' a reachable cell, so a valid origin is returned unchanged -- only a masked
#' origin is moved, and only to the nearest real, reachable cell.
#'
#' @param point An `sf` point (WGS84) to snap.
#' @param tr_corrected A geo-corrected `gdistance` TransitionLayer.
#'
#' @return An `sf` point (WGS84): `point` unchanged if it is already on a valid
#'   cell (or if the grid/transition cannot be read), otherwise a point on the
#'   nearest reachable cell.
#'
#' @keywords internal
#' @noRd
.nearest_valid_origin <- function(point, tr_corrected) {
  grid <- tryCatch(raster::raster(tr_corrected), error = function(e) NULL)
  if (is.null(grid)) return(point)

  p     <- sf::st_coordinates(sf::st_transform(point, crs = 4326))[1, 1:2]
  ocell <- raster::cellFromXY(grid, matrix(p, ncol = 2))

  # Reachable cells = transition-layer nodes with any non-zero conductance.
  # Masked friction cells are disconnected, so their row/col sums are zero.
  tm    <- gdistance::transitionMatrix(tr_corrected)
  valid <- which(Matrix::rowSums(tm) > 0 | Matrix::colSums(tm) > 0)

  if (!is.na(ocell) && ocell %in% valid) return(point)  # origin already reachable
  if (length(valid) == 0) return(point)                 # nothing reachable -> give up

  xy <- raster::xyFromCell(grid, valid)
  i  <- which.min((xy[, 1] - p[1])^2 + (xy[, 2] - p[2])^2)
  sf::st_sf(geometry = sf::st_sfc(sf::st_point(xy[i, ]), crs = 4326))
}
