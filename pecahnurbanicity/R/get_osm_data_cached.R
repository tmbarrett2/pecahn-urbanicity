#' Query and locally cache OpenStreetMap (OSM) data
#'
#' @description
#' Retrieves spatial data from OpenStreetMap for a specified bounding box and feature key/value pair.
#' Results are cached locally as `.rds` files to reduce redundant network calls in future runs.
#'
#' @param bbox Numeric vector defining the bounding box (`c(xmin, ymin, xmax, ymax)`).
#' @param key Character; the OSM key (e.g., `"highway"`, `"amenity"`).
#' @param value Optional character or vector; OSM value(s) (e.g., `"school"`, `"hospital"`).
#' @param cache_dir Directory where cached `.rds` files will be stored (default: `"osm_cache"`).
#'
#' @return
#' An `osmdata_sf` object (list with `$osm_points`, `$osm_lines`, `$osm_polygons`) containing OSM data.
#'
#' @details
#' Cached filenames are generated using a hashed digest of the input query parameters (`bbox`, `key`, and `value`).
#' This ensures unique caching for different feature types and regions while minimizing redundant API requests.
#'
#' @examples
#' \dontrun{
#' roads <- get_osm_data_cached(bbox = c(47.5, -14.3, 47.6, -14.2), key = "highway")
#' }
#'
#' @importFrom digest digest
#' @importFrom osmdata opq add_osm_feature osmdata_sf
#' @export
get_osm_data_cached <- function(bbox, key, value = NULL, cache_dir = "osm_cache") {
  # Ensure required package is available
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for caching.")
  }
  
  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Create unique cache identifier
  cache_id <- digest::digest(list(bbox = bbox, key = key, value = value))
  cache_path <- file.path(cache_dir, paste0(cache_id, ".rds"))
  
  # Return cached file if available
  if (file.exists(cache_path)) {
    cached_data <- readRDS(cache_path)
    return(cached_data)
  }
  
  # Otherwise query OSM and cache the results
  query <- osmdata::opq(bbox) |>
    osmdata::add_osm_feature(key = key, value = value)
  
  res <- osmdata::osmdata_sf(query)
  
  # Save to cache for future use
  saveRDS(res, cache_path)
  
  return(res)
}
