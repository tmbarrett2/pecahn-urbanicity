#' Create bounding boxes for community locations
#'
#' @description
#' Generates bounding boxes for each community in a user-provided dataset of latitude and longitude coordinates.
#' If specific corner coordinates are provided in `custom_bbox`, the function uses those directly; otherwise,
#' it computes a bounding box that extends a fixed distance (in kilometers) around each community's central point.
#'
#' @param gps_df A data frame containing the following columns:
#'   \describe{
#'     \item{project}{Character; project or dataset name.}
#'     \item{ethnicity}{Character; population or ethnic group.}
#'     \item{community}{Character; community name.}
#'     \item{lat}{Numeric; community latitude in decimal degrees.}
#'     \item{lon}{Numeric; community longitude in decimal degrees.}
#'   }
#' @param distance_km Numeric; the radius (in kilometers) around each central point to define the bounding box.
#'   Defaults to 5 km.
#' @param custom_bbox Optional named list of bounding box coordinates, which may include either
#'   two corners (`nw_lat`, `nw_lon`, `se_lat`, `se_lon`) or four corners
#'   (`nw_lat`, `nw_lon`, `ne_lat`, `ne_lon`, `sw_lat`, `sw_lon`, `se_lat`, `se_lon`).
#'
#' @return
#' A named list of bounding box objects, each containing:
#'   \itemize{
#'     \item \code{bbox}: Numeric vector with named elements \code{left}, \code{bottom}, \code{right}, \code{top}.
#'     \item \code{project}, \code{ethnicity}, \code{community}: Metadata values from the input table.
#'   }
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   project = "CARE", ethnicity = "Qom/Toba", community = "Namqom",
#'   lat = -26.118, lon = -58.2254
#' )
#' bboxes <- create_bounding_boxes(data, distance_km = 10)
#' }
#'
#' @importFrom stats setNames
#' @export
create_bounding_boxes <- function(gps_df, distance_km = 5, custom_bbox = NULL) {
  # Validate columns and fill missing metadata with NA
  required_cols <- c("project", "ethnicity", "community", "lat", "lon")
  for (col in required_cols) {
    if (!col %in% names(gps_df)) {
      gps_df[[col]] <- NA
    }
  }
  
  # Ensure key metadata columns are character (avoid factor issues)
  gps_df$project   <- as.character(gps_df$project)
  gps_df$ethnicity <- as.character(gps_df$ethnicity)
  gps_df$community <- as.character(gps_df$community)
  
  # Initialize result list
  result <- vector("list", nrow(gps_df))
  
  # Process each community
  for (i in seq_len(nrow(gps_df))) {
    # Extract numeric lat/lon
    lat <- as.numeric(gps_df$lat[i])
    lon <- as.numeric(gps_df$lon[i])
    
    # Determine bounding box
    if (!is.null(custom_bbox) &&
        all(c("nw_lat", "nw_lon", "se_lat", "se_lon") %in% names(custom_bbox))) {
      # Use 2-corner custom box (NW–SE)
      bbox <- c(
        left   = custom_bbox$nw_lon,
        bottom = custom_bbox$se_lat,
        right  = custom_bbox$se_lon,
        top    = custom_bbox$nw_lat
      )
      
    } else if (!is.null(custom_bbox) &&
               all(c("nw_lat", "nw_lon", "ne_lat", "ne_lon",
                     "sw_lat", "sw_lon", "se_lat", "se_lon") %in% names(custom_bbox))) {
      # Use 4-corner custom box
      bbox <- c(
        left   = min(custom_bbox$nw_lon, custom_bbox$sw_lon),
        bottom = min(custom_bbox$sw_lat, custom_bbox$se_lat),
        right  = max(custom_bbox$ne_lon, custom_bbox$se_lon),
        top    = max(custom_bbox$nw_lat, custom_bbox$ne_lat)
      )
      
    } else {
      # Calculate offset-based box using distance_km
      lat_offset_deg <- distance_km / 110.57
      lon_km_per_degree <- 111.32 * cos(lat * pi / 180)
      lon_offset_deg <- distance_km / lon_km_per_degree
      
      bbox <- c(
        left   = lon - lon_offset_deg,
        bottom = lat - lat_offset_deg,
        right  = lon + lon_offset_deg,
        top    = lat + lat_offset_deg
      )
    }
    
    # Store bbox and associated metadata
    result[[i]] <- list(
      bbox      = bbox,
      project   = gps_df$project[i],
      ethnicity = gps_df$ethnicity[i],
      community = gps_df$community[i]
    )
  }
  
  # Name list elements by community for easier downstream reference
  names(result) <- gps_df$community
  
  # Return list of bounding boxes
  return(result)
}
