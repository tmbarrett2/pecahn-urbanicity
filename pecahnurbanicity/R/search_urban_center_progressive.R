#' Search progressively for nearest urban center (place = city)
#'
#' @description
#' Expands the search area incrementally to find the nearest city in OpenStreetMap (`place = city`),
#' optionally calculating least-cost travel time using a friction surface.
#'
#' @param center_point An `sf` point (WGS84) representing the community center.
#' @param bbox Numeric vector defining the initial bounding box (`c(left, bottom, right, top)`).
#' @param friction_extent Optional `raster::extent` of the friction raster.
#' @param tr_corrected Optional `gdistance` TransitionLayer (geo-corrected) for travel time estimation.
#' @param raster_crs CRS of the friction raster.
#' @param max_search_multiplier Integer; maximum expansion factor for the bounding box (default = 50).
#' @param search_increment Integer; increment of expansion in each step (default = 2).
#'
#' @return
#' Numeric; minimum estimated travel time (in minutes) to the nearest city, or `NA` if no city is found.
#'
#' @examples
#' \dontrun{
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
                                            search_increment = 2) {
  # Initialize
  search_multiplier <- 0
  found_cities <- NULL
  
  while (is.null(found_cities) && search_multiplier <= max_search_multiplier) {
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
    
    # Query OSM for cities
    city_data <- get_osm_data_cached(search_bbox, key = "place", value = "city")
    if (!is.null(city_data$osm_points) && nrow(city_data$osm_points) > 0) {
      found_cities <- city_data$osm_points
    }
    
    # Expand search if none found
    if (is.null(found_cities)) {
      search_multiplier <- search_multiplier + search_increment
      message(sprintf("  No cities found, expanding search (multiplier: %d)...", search_multiplier))
    }
  }
  
  # Compute travel time if cities were found and transition layer is available
  if (!is.null(found_cities) && nrow(found_cities) > 0 && !is.null(tr_corrected)) {
    if (!is.null(friction_extent)) {
      fac_ext <- sf::st_bbox(found_cities)
      if (fac_ext["xmin"] < friction_extent@xmin ||
          fac_ext["xmax"] > friction_extent@xmax ||
          fac_ext["ymin"] < friction_extent@ymin ||
          fac_ext["ymax"] > friction_extent@ymax) {
        message("  Found cities outside friction raster extent. Returning NA.")
        return(NA)
      }
    }
    
    return(calculate_travel_time(center_point, found_cities, tr_corrected, raster_crs))
  }
  
  message("  No cities found within maximum search extent.")
  return(NA)
}
