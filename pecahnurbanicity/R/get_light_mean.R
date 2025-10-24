#' Extract mean nighttime light intensity
#'
#' @description
#' Computes the mean nighttime light intensity within a given polygon using a VIIRS or other raster dataset.
#'
#' @param polygon_sf An `sf` polygon object defining the spatial boundary of interest.
#' @param raster_data A `RasterLayer` object (e.g., VIIRS raster of nighttime light intensity).
#'
#' @return
#' A numeric value giving the mean light intensity within the polygon, rounded to three decimal places.
#'
#' @details
#' The function converts the input `sf` polygon to a `Spatial` object and extracts the mean raster value,
#' ignoring missing data (`NA` values).
#'
#' @examples
#' \dontrun{
#' mean_light <- get_light_mean(polygon_sf, raster_data)
#' }
#'
#' @keywords internal
#' @importFrom methods as
#' @importFrom raster extract
#' @export
get_light_mean <- function(polygon_sf, raster_data) {
  # Convert sf polygon to Spatial for compatibility with raster::extract
  polygon_sp <- methods::as(polygon_sf, "Spatial")
  
  # Extract mean raster value within the polygon, ignoring NAs
  light_mean <- raster::extract(raster_data, polygon_sp, fun = mean, na.rm = TRUE)
  
  # Return rounded mean value (3 decimal places)
  return(round(light_mean, 3))
}
