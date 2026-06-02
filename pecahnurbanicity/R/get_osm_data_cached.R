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
#' @param max_retries Integer; number of times to retry a failed Overpass query before giving up
#'   (default = 3, i.e. up to 4 attempts in total). A query that succeeds but finds no features is
#'   *not* an error and is never retried.
#' @param backoff Numeric; base wait in seconds between retries. The wait grows exponentially across
#'   attempts (`backoff`, `2 * backoff`, `4 * backoff`, ...) (default = 5).
#' @param timeout Integer; server-side query timeout in seconds, passed to `osmdata::opq()` (default = 60).
#'
#' @return
#' An `osmdata_sf` object (list with `$osm_points`, `$osm_lines`, `$osm_polygons`) containing OSM data.
#'
#' @details
#' Cached filenames are generated using a hashed digest of the input query parameters (`bbox`, `key`, and `value`).
#' This ensures unique caching for different feature types and regions while minimizing redundant API requests.
#'
#' The Overpass API frequently returns transient rate-limit (HTTP 429) or timeout/server (HTTP 504)
#' errors. To stop a single such failure from becoming a permanent `NA` in the caller's results, the
#' query is retried up to `max_retries` times with exponential backoff. Only successful results are
#' written to the cache (failures are never cached), so re-running a job retries only the cells that
#' previously failed. If every attempt fails the final error is re-raised, so the calling metric still
#' records `NA` for that query rather than silently masking a persistent outage.
#'
#' @examples
#' \dontrun{
#' roads <- get_osm_data_cached(bbox = c(47.5, -14.3, 47.6, -14.2), key = "highway")
#' }
#'
#' @importFrom digest digest
#' @importFrom osmdata opq add_osm_feature osmdata_sf
#' @export
get_osm_data_cached <- function(bbox, key, value = NULL, cache_dir = "osm_cache",
                                max_retries = 3, backoff = 5, timeout = 60) {
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

  # Otherwise query OSM, retrying transient (rate-limit / timeout / server) failures
  # with exponential backoff. A query that simply finds no features still returns a
  # valid osmdata_sf object, so only genuine errors trigger a retry.
  query <- osmdata::opq(bbox, timeout = timeout) |>
    osmdata::add_osm_feature(key = key, value = value)

  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    res <- tryCatch(osmdata::osmdata_sf(query), error = function(e) e)

    if (!inherits(res, "error")) break

    if (attempt > max_retries) {
      stop(sprintf("OSM query failed after %d attempt(s) (key = '%s'): %s",
                   attempt, key, conditionMessage(res)))
    }

    wait <- backoff * 2^(attempt - 1L)
    message(sprintf("  OSM query failed (attempt %d of %d), retrying in %gs: %s",
                    attempt, max_retries + 1L, wait, conditionMessage(res)))
    Sys.sleep(wait)
  }

  # Save successful result to cache for future use
  saveRDS(res, cache_path)

  return(res)
}
