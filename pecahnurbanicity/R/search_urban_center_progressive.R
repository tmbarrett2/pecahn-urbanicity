#' Search progressively for the nearest urban center (GHS-SMOD)
#'
#' @description
#' Expands the search area incrementally to find the nearest urban center using the
#' GHSL Degree of Urbanisation settlement model (GHS-SMOD), retrieved from Google
#' Earth Engine, then computes least-cost travel time to it using a friction surface.
#' Replaces the previous OpenStreetMap `place = city` search.
#'
#' @param center_point An `sf` point (WGS84) representing the community center.
#' @param bbox Numeric vector defining the initial bounding box (`c(left, bottom, right, top)`).
#' @param friction_extent Optional `raster::extent` of the friction raster.
#' @param tr_corrected Optional `gdistance` TransitionLayer (geo-corrected) for travel time estimation.
#' @param raster_crs CRS of the friction raster.
#' @param max_search_multiplier Integer; maximum expansion factor for the bounding box (default = 50).
#' @param search_increment Integer; increment of expansion in each step (default = 2).
#' @param urban_center_codes Integer vector of GHS-SMOD `smod_code` values treated as an urban
#'   center. Defaults to `c(23, 30)` (Dense Urban Cluster + Urban Center). Use `30` for strict
#'   Urban Center only.
#' @param ghsl_year Integer; GHS-SMOD epoch to use (5-year epochs from 1975 to 2030, default = 2020).
#'   Non-epoch years are snapped to the nearest available epoch.
#'
#' @return
#' A named list with elements `value` (numeric; minimum estimated travel time in minutes to the nearest
#' urban center, or `NA`) and `is_boundary_est` (logical; `TRUE` when no urban center was found within the
#' maximum search area and a boundary-distance estimate was returned instead, `FALSE` otherwise).
#'
#' @details
#' GHS-SMOD `smod_code` uses Degree-of-Urbanisation level-2 codes: 30 = Urban Center,
#' 23 = Dense Urban Cluster, 22 = Semi-dense Urban Cluster, 21 = Suburban/Peri-urban,
#' 13 = Rural Cluster, 12 = Low Density Rural, 11 = Very Low Density Rural, 10 = Water.
#'
#' The search starts from `bbox` and expands by `search_increment` each step until urban-center pixels
#' (those whose `smod_code` is in `urban_center_codes`) are found within the extent, or
#' `max_search_multiplier` is reached. Travel time is then computed to the nearest reachable urban-center
#' centroid. If none are found — or if those found are all unreachable on the friction surface (beyond the
#' raster's coverage, so the travel time is `NA`) — a boundary-distance estimate is returned (the mean
#' least-cost travel time to the four cardinal edge-midpoints of the final search bounding box) with
#' `is_boundary_est = TRUE`.
#'
#' @examples
#' \dontrun{
#' # Earth Engine must be authenticated first: rgee::ee_Initialize()
#' tt_city <- search_urban_center_progressive(
#'   center_point, bbox,
#'   tr_corrected = tr, raster_crs = raster::crs(friction)
#' )
#' }
#'
#' @keywords internal
#' @importFrom sf st_bbox
#' @export
search_urban_center_progressive <- function(center_point, bbox,
                                            friction_extent = NULL,
                                            tr_corrected = NULL,
                                            raster_crs = NULL,
                                            max_search_multiplier = 50,
                                            search_increment = 2,
                                            urban_center_codes = c(23, 30),
                                            ghsl_year = 2020) {
  # Snap the requested year to the nearest available GHS-SMOD epoch (messages once).
  ghsl_year <- .snap_ghsl_year(ghsl_year)

  # Initialize
  search_multiplier <- 0
  found_centers <- NULL

  while (is.null(found_centers) && search_multiplier <= max_search_multiplier) {
    # Define search bbox
    if (search_multiplier == 0) {
      search_bbox <- c(bbox["left"], bbox["bottom"], bbox["right"], bbox["top"])
    } else {
      width  <- bbox["right"] - bbox["left"]
      height <- bbox["top"] - bbox["bottom"]
      search_bbox <- c(
        bbox["left"] - width * search_multiplier,
        bbox["bottom"] - height * search_multiplier,
        bbox["right"] + width * search_multiplier,
        bbox["top"] + height * search_multiplier
      )
    }

    # Pull urban-center pixel centroids from GHS-SMOD (Earth Engine, cached).
    centers <- ee_get_ghsl_urban_centers(search_bbox, ghsl_year, urban_center_codes)
    if (!is.null(centers) && nrow(centers) > 0) {
      found_centers <- centers
    }

    # Expand search if none found
    if (is.null(found_centers)) {
      search_multiplier <- search_multiplier + search_increment
      message(sprintf("  No urban center found, expanding search (multiplier: %d)...", search_multiplier))
    }
  }

  # Compute travel time if urban centers were found and a transition layer exists.
  if (!is.null(found_centers) && nrow(found_centers) > 0 && !is.null(tr_corrected)) {
    if (!is.null(friction_extent)) {
      fac_ext <- sf::st_bbox(found_centers)
      if (fac_ext["xmin"] < friction_extent@xmin ||
          fac_ext["xmax"] > friction_extent@xmax ||
          fac_ext["ymin"] < friction_extent@ymin ||
          fac_ext["ymax"] > friction_extent@ymax) {
        message("  Found urban center outside friction raster extent. Returning NA.")
        return(list(value = NA_real_, is_boundary_est = FALSE))
      }
    }

    # Pass all candidate centroids; calculate_travel_time() returns the minimum,
    # i.e. travel time to the nearest reachable urban center.
    tt <- calculate_travel_time(center_point, found_centers, tr_corrected, raster_crs)
    if (is.finite(tt)) {
      return(list(value = tt, is_boundary_est = FALSE))
    }
    # Urban centers were found but none is reachable on the friction surface (all
    # lie beyond the raster's coverage), so calculate_travel_time() returned NA.
    # Fall through to the boundary-distance estimate rather than reporting a bare NA.
    message("  Nearest urban center unreachable on friction surface - using boundary-distance estimate.")
  }

  # Search exhausted, found-but-unreachable, or no friction surface: boundary-distance estimate.
  message("  No reachable urban center found - returning boundary-distance estimate.")
  est <- .boundary_distance_estimate(search_bbox, center_point, tr_corrected, raster_crs)
  return(list(value = est, is_boundary_est = TRUE))
}
