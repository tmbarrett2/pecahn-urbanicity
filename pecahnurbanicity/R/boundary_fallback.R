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
  if (looks_geographic) {
    # Inset by a tiny fraction so probe points sit strictly inside the grid.
    pad_x <- (ext@xmax - ext@xmin) * 1e-6
    pad_y <- (ext@ymax - ext@ymin) * 1e-6
    left   <- max(left,   ext@xmin + pad_x)
    right  <- min(right,  ext@xmax - pad_x)
    bottom <- max(bottom, ext@ymin + pad_y)
    top    <- min(top,    ext@ymax - pad_y)
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
      calculate_travel_time(center_point, pt, tr_corrected, raster_crs),
      error = function(e) NA_real_
    )
  }, numeric(1))

  finite_times <- times[is.finite(times)]
  if (length(finite_times) == 0) return(NA_real_)
  mean(finite_times)
}
