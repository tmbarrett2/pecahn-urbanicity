#' Calculate minimum travel time between points using a friction surface
#'
#' @description
#' Computes the minimum travel time (in minutes) from a single origin point to one or
#' more destination points using a precomputed, geo-corrected transition layer from
#' a friction surface.
#'
#' @param from_point An `sf` point object representing the origin (typically community center).
#' @param to_points An `sf` object containing one or more destination points.
#' @param tr_corrected A `gdistance::TransitionLayer` object that has been geo-corrected.
#' @param raster_crs The coordinate reference system (CRS) of the friction surface raster.
#'
#' @return
#' A numeric value giving the minimum estimated travel time (in minutes) between the
#' origin and destinations, or `NA` if unavailable.
#'
#' @details
#' This function transforms both origin and destination coordinates to the CRS of the
#' friction surface, computes least-cost path distances using `gdistance::costDistance()`,
#' and returns the minimum finite value.
#'
#' @examples
#' \dontrun{
#' time_min <- calculate_travel_time(from_point = origin_sf,
#'                                   to_points = facilities_sf,
#'                                   tr_corrected = transition_layer,
#'                                   raster_crs = raster::crs(friction_raster))
#' }
#'
#' @importFrom sf st_transform st_coordinates
#' @importFrom gdistance costDistance
#' @export
calculate_travel_time <- function(from_point, to_points, tr_corrected, raster_crs) {
  # Return NA if destination points are missing
  if (is.null(to_points) || nrow(to_points) == 0) return(NA_real_)
  
  # Transform coordinates to friction surface CRS
  from_sp <- sf::st_transform(from_point, crs = raster_crs)
  to_sp   <- sf::st_transform(to_points,   crs = raster_crs)
  
  # Extract coordinates
  from_coords <- sf::st_coordinates(from_sp)
  to_coords   <- sf::st_coordinates(to_sp)
  
  # Compute cost-distance
  cost_dist <- gdistance::costDistance(tr_corrected, from_coords[1, ], to_coords)
  min_time  <- suppressWarnings(min(cost_dist, na.rm = TRUE))
  
  # Return minimum time if finite, otherwise NA
  if (is.finite(min_time)) min_time else NA_real_
}
