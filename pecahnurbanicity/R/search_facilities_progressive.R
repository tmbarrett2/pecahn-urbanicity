#' Progressively search for facilities (roads, schools, healthcare)
#'
#' @description
#' Expands the bounding box around a community incrementally until at least one facility of the specified type
#' is found in OpenStreetMap data. Optionally computes travel time to the nearest facility using a friction surface.
#'
#' @param center_point An `sf` point (WGS84) representing the community center.
#' @param bbox Numeric vector defining the initial bounding box (`c(left, bottom, right, top)`).
#' @param facility_type Character string; one of `"paved_road"`, `"healthcare"`, or `"school"`.
#' @param friction_extent Optional `raster::extent` defining the spatial coverage of the friction raster.
#' @param tr_corrected Optional `gdistance` TransitionLayer (geo-corrected) for travel time computation.
#' @param raster_crs CRS of the friction raster.
#' @param max_search_multiplier Integer; maximum expansion factor for the bounding box (default = 20).
#' @param search_increment Integer; amount to increase the expansion multiplier in each iteration (default = 2).
#'
#' @return
#' Numeric; minimum walking travel time (in minutes) from the community center to the nearest facility, or `NA`
#' if no facilities are found within the expanded area.
#'
#' @examples
#' \dontrun{
#' tt_school <- search_facilities_progressive(
#'   center_point, bbox, "school",
#'   tr_corrected = tr, raster_crs = raster::crs(friction)
#' )
#' }
#'
#' @keywords internal
#' @importFrom sf st_geometry st_centroid st_bbox st_cast st_sf
#' @export
search_facilities_progressive <- function(center_point, bbox, facility_type,
                                          friction_extent = NULL,
                                          tr_corrected = NULL,
                                          raster_crs = NULL,
                                          max_search_multiplier = 20,
                                          search_increment = 2) {
  # Initialize
  search_multiplier <- 0
  found_facilities <- NULL
  
  while (is.null(found_facilities) && search_multiplier <= max_search_multiplier) {
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
    
    # Query OSM by type (cached)
    if (facility_type == "paved_road") {
      roads_data <- get_osm_data_cached(search_bbox, key = "highway", value = NULL)
      if (!is.null(roads_data$osm_lines) && nrow(roads_data$osm_lines) > 0) {
        roads_sf <- roads_data$osm_lines
        if (!"surface" %in% colnames(roads_sf)) roads_sf$surface <- NA
        paved_surfaces <- c("paved", "asphalt", "concrete")
        paved_roads <- roads_sf[!is.na(roads_sf$surface) & roads_sf$surface %in% paved_surfaces, ]
        if (nrow(paved_roads) > 0) found_facilities <- paved_roads
      }
    } else if (facility_type == "healthcare") {
      health_data <- get_osm_data_cached(search_bbox, key = "amenity",
                                         value = c("hospital", "clinic", "doctors"))
      if (!is.null(health_data$osm_points) && nrow(health_data$osm_points) > 0) {
        found_facilities <- health_data$osm_points
      }
    } else if (facility_type == "school") {
      school_data <- get_osm_data_cached(search_bbox, key = "amenity", value = "school")
      school_points <- NULL
      if (!is.null(school_data$osm_points) && nrow(school_data$osm_points) > 0) {
        school_points <- sf::st_sf(geometry = sf::st_geometry(school_data$osm_points))
      }
      if (!is.null(school_data$osm_polygons) && nrow(school_data$osm_polygons) > 0) {
        school_centroids <- sf::st_sf(geometry = sf::st_centroid(sf::st_geometry(school_data$osm_polygons)))
        if (is.null(school_points)) {
          school_points <- school_centroids
        } else {
          school_points <- rbind(school_points, school_centroids)
        }
      }
      found_facilities <- school_points
    }
    
    # Expand search if none found
    if (is.null(found_facilities)) {
      search_multiplier <- search_multiplier + search_increment
      message(sprintf("  No %s found, expanding search area (multiplier: %d)...",
                      facility_type, search_multiplier))
    }
  }
  
  # Compute travel time (bounds check + convert lines → points)
  if (!is.null(found_facilities) && nrow(found_facilities) > 0 && !is.null(tr_corrected)) {
    if (!is.null(friction_extent)) {
      fac_ext <- sf::st_bbox(found_facilities)
      if (fac_ext["xmin"] < friction_extent@xmin ||
          fac_ext["xmax"] > friction_extent@xmax ||
          fac_ext["ymin"] < friction_extent@ymin ||
          fac_ext["ymax"] > friction_extent@ymax) {
        message(sprintf("  Found %s outside friction raster extent. Returning NA.", facility_type))
        return(NA)
      }
    }
    
    if (facility_type == "paved_road") {
      suppressWarnings({
        facility_points <- sf::st_sf(geometry = sf::st_cast(sf::st_geometry(found_facilities), "POINT"))
      })
    } else {
      facility_points <- found_facilities
    }
    
    return(calculate_travel_time(center_point, facility_points, tr_corrected, raster_crs))
  }
  
  message(sprintf("  No %s found within maximum search area.", facility_type))
  return(NA)
}
