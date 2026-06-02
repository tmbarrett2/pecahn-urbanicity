# -----------------------------------------------------------------------------
# Google Earth Engine raster extraction helpers (internal).
#
# These functions replace the local-file raster loads previously used by
# compute_urbanicity(). Each clips/reduces on the Earth Engine side before
# pulling data into R, and caches results to a local `gee_cache/` directory
# (mirroring get_osm_data_cached()) keyed on a digest of (asset, region, year).
#
# Default asset IDs were verified against the live Earth Engine catalog. Override
# any of them via the *_asset arguments on compute_urbanicity() /
# compute_urbanicity_iterative() if the catalog changes.
# -----------------------------------------------------------------------------

# Friction surface: MAP Global Friction Surface 2019 (walking-only).
#   ee.Image, single band 'friction', units minutes/meter, native scale ~927.67 m.
#   Feeds gdistance::transition() exactly as the previous local raster did.
.PECAHN_FRICTION_ASSET <- "projects/malariaatlasproject/assets/accessibility/friction_surface/2019_v5_1_walking_only"
.PECAHN_FRICTION_BAND  <- "friction"
.PECAHN_FRICTION_SCALE <- 927.67

# Population density: WorldPop next-generation Global Project Population (~100 m).
#   ImageCollection 'projects/sat-io/open-datasets/WORLDPOP/pop' from the
#   awesome-gee-community-catalog: band 'population' (estimated count per pixel),
#   annual coverage 2015-2030 (extends well past the 2021 limit of the official
#   'WorldPop/GP/100m/pop'). Image property 'year' (DOUBLE); we filterBounds() +
#   filter(year) and mosaic() to cover arbitrary sites. Because this is a count
#   raster, reduceRegion(sum) is run at the image's native nominalScale() (derived
#   per call) so totals stay correct regardless of the exact grid; .PECAHN_POP_SCALE
#   is only a fallback. For pre-2015 years, override population_asset with the
#   official 'WorldPop/GP/100m/pop' (band 'population', 2000-2021).
.PECAHN_POP_ASSET <- "projects/sat-io/open-datasets/WORLDPOP/pop"
.PECAHN_POP_BAND  <- "population"
.PECAHN_POP_SCALE <- 92.77

# Nighttime lights: VIIRS DNB annual composites.
#   ImageCollection 'NOAA/VIIRS/DNB/ANNUAL_V22' (covers 2012-2024), band
#   'average_masked' matches the masked annual product the package used locally.
#   Native scale ~463.83 m. Falls back to ANNUAL_V21 for out-of-range years.
.PECAHN_NTL_ASSET          <- "NOAA/VIIRS/DNB/ANNUAL_V22"
.PECAHN_NTL_ASSET_FALLBACK <- "NOAA/VIIRS/DNB/ANNUAL_V21"
.PECAHN_NTL_BAND           <- "average_masked"
.PECAHN_NTL_SCALE          <- 463.83

# Null-coalescing helper (internal).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Normalize a bounding box to an unnamed c(xmin, ymin, xmax, ymax) vector.
# Accepts the named left/bottom/right/top form used throughout the package as
# well as a plain xmin/ymin/xmax/ymax vector.
.bbox_xy <- function(bbox) {
  if (!is.null(names(bbox)) &&
      all(c("left", "bottom", "right", "top") %in% names(bbox))) {
    out <- c(bbox[["left"]], bbox[["bottom"]], bbox[["right"]], bbox[["top"]])
  } else {
    out <- as.numeric(bbox)[1:4]
  }
  as.numeric(unname(out))
}

# Build a path inside the cache directory from a digest of `parts`.
.gee_cache_path <- function(parts, ext = ".rds", cache_dir = "gee_cache") {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for caching.", call. = FALSE)
  }
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
  cache_id <- digest::digest(parts)
  file.path(cache_dir, paste0(cache_id, ext))
}

# Construct an EE rectangle geometry from c(xmin, ymin, xmax, ymax).
.ee_rectangle <- function(xy) {
  rgee::ee$Geometry$Rectangle(
    coords    = c(xy[1], xy[2], xy[3], xy[4]),
    proj      = "EPSG:4326",
    geodesic  = FALSE
  )
}

#' Fetch the friction surface for a bounding box from Earth Engine (internal)
#'
#' Returns a `RasterLayer` of the MAP walking-only friction surface, clipped on
#' the Earth Engine side to `bbox` expanded by `buffer_deg` degrees, ready for
#' `gdistance::transition()`. The cropped raster is written to (and re-read from)
#' a GeoTIFF in `gee_cache/`.
#'
#' @param bbox Bounding box (named left/bottom/right/top or xmin/ymin/xmax/ymax).
#' @param buffer_deg Numeric; degrees to expand the bbox before clipping.
#' @param asset Optional character override for the friction asset ID.
#' @param cache_dir Directory for cached GeoTIFFs (default `"gee_cache"`).
#'
#' @return A `RasterLayer`.
#' @keywords internal
#' @noRd
ee_get_friction_raster <- function(bbox, buffer_deg = 1, asset = NULL,
                                   cache_dir = "gee_cache") {
  asset <- asset %||% .PECAHN_FRICTION_ASSET
  band  <- .PECAHN_FRICTION_BAND

  xy     <- .bbox_xy(bbox)
  region <- c(xy[1] - buffer_deg, xy[2] - buffer_deg,
              xy[3] + buffer_deg, xy[4] + buffer_deg)

  cache_path <- .gee_cache_path(
    list(type = "friction", asset = asset, band = band, region = region),
    ext = ".tif", cache_dir = cache_dir
  )

  # Cache hit: re-open the stored GeoTIFF.
  if (file.exists(cache_path)) {
    r <- tryCatch(raster::raster(cache_path), error = function(e) NULL)
    if (!is.null(r)) return(r)
  }

  geom <- .ee_rectangle(region)
  img  <- rgee::ee$Image(asset)$select(band)

  # Pull the clipped raster synchronously via a download URL. This avoids
  # rgee::ee_as_raster() (deprecated, and its 'via' only supports the
  # Drive/GCS export path). The friction surface is single-band and, once
  # clipped to a site bbox at ~927 m, is well within getDownloadURL size limits.
  url <- img$getDownloadURL(list(
    region      = geom,
    scale       = .PECAHN_FRICTION_SCALE,
    crs         = "EPSG:4326",
    format      = "GEO_TIFF",
    filePerBand = FALSE
  ))

  tmp <- tempfile(fileext = ".dl")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)

  # GEO_TIFF normally returns a single .tif; some endpoints return a zip.
  magic <- readBin(tmp, "raw", n = 2L)
  if (identical(magic, as.raw(c(0x50, 0x4b)))) {          # "PK" -> zip archive
    files <- utils::unzip(tmp, exdir = tempdir())
    tif   <- files[grepl("\\.tif$", files, ignore.case = TRUE)][1]
    file.copy(tif, cache_path, overwrite = TRUE)
  } else {
    file.copy(tmp, cache_path, overwrite = TRUE)
  }

  raster::raster(cache_path)
}

#' Total population within a bounding box from Earth Engine (internal)
#'
#' Sums WorldPop population counts within `bbox` for `year`. The caller divides
#' the returned total by the bbox area to obtain density, exactly as before.
#'
#' @param bbox Bounding box (named left/bottom/right/top or xmin/ymin/xmax/ymax).
#' @param year Integer year (default WorldPop next-gen coverage 2015-2030).
#' @param asset Optional character override for the population asset ID.
#' @param cache_dir Directory for cached values (default `"gee_cache"`).
#'
#' @return A numeric scalar total population, or `NA_real_` on failure.
#' @keywords internal
#' @noRd
ee_get_population <- function(bbox, year, asset = NULL, cache_dir = "gee_cache") {
  asset <- asset %||% .PECAHN_POP_ASSET
  band  <- .PECAHN_POP_BAND
  year  <- as.integer(year)

  xy <- .bbox_xy(bbox)

  cache_path <- .gee_cache_path(
    list(type = "pop", asset = asset, band = band, region = xy, year = year),
    ext = ".rds", cache_dir = cache_dir
  )
  if (file.exists(cache_path)) {
    v <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(v)) return(v)
  }

  val <- tryCatch({
    geom <- .ee_rectangle(xy)
    # The next-gen sat-io WorldPop stores 'year' as a STRING ("2020"); the
    # official 'WorldPop/GP/100m/pop' stores it as a number. Match either.
    ic <- rgee::ee$ImageCollection(asset)$
      filter(rgee::ee$Filter$Or(
        rgee::ee$Filter$eq("year", year),
        rgee::ee$Filter$eq("year", as.character(year))
      ))$
      filterBounds(geom)$
      select(band)
    # Population is a count raster, so sum at the grid's native resolution.
    # mosaic() drops the native projection, so read nominalScale() from a source
    # image first; fall back to the constant if the collection is empty.
    native_scale <- tryCatch(
      ic$first()$projection()$nominalScale()$getInfo(),
      error = function(e) .PECAHN_POP_SCALE
    )
    if (is.null(native_scale) || !is.finite(native_scale)) {
      native_scale <- .PECAHN_POP_SCALE
    }
    img  <- ic$mosaic()
    stat <- img$reduceRegion(
      reducer   = rgee::ee$Reducer$sum(),
      geometry  = geom,
      scale     = native_scale,
      maxPixels = 1e13,
      bestEffort = TRUE
    )
    out <- stat$get(band)$getInfo()
    if (is.null(out)) NA_real_ else as.numeric(out)
  }, error = function(e) {
    warning("ee_get_population failed for year ", year, ": ",
            conditionMessage(e), call. = FALSE)
    NA_real_
  })

  # Only cache successes — never poison the cache with a transient NA.
  if (!is.na(val)) saveRDS(val, cache_path)
  val
}

#' Mean nighttime light radiance within a bounding box from Earth Engine (internal)
#'
#' Returns the mean VIIRS annual radiance within `bbox` for `year`, using
#' ANNUAL_V22 and falling back to ANNUAL_V21 when the primary collection has no
#' image for the requested year.
#'
#' @param bbox Bounding box (named left/bottom/right/top or xmin/ymin/xmax/ymax).
#' @param year Integer year.
#' @param asset Optional character override for the primary nighttime-light asset.
#' @param asset_fallback Optional character override for the fallback asset.
#' @param cache_dir Directory for cached values (default `"gee_cache"`).
#'
#' @return A numeric scalar mean radiance, or `NA_real_` on failure.
#' @keywords internal
#' @noRd
ee_get_nighttime_light <- function(bbox, year, asset = NULL,
                                   asset_fallback = NULL, cache_dir = "gee_cache") {
  asset          <- asset %||% .PECAHN_NTL_ASSET
  asset_fallback <- asset_fallback %||% .PECAHN_NTL_ASSET_FALLBACK
  band           <- .PECAHN_NTL_BAND
  year           <- as.integer(year)

  xy <- .bbox_xy(bbox)

  cache_path <- .gee_cache_path(
    list(type = "ntl", asset = asset, band = band, region = xy, year = year),
    ext = ".rds", cache_dir = cache_dir
  )
  if (file.exists(cache_path)) {
    v <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(v)) return(v)
  }

  mean_for_asset <- function(asset_id) {
    geom <- .ee_rectangle(xy)
    ic <- rgee::ee$ImageCollection(asset_id)$
      filterDate(paste0(year, "-01-01"), paste0(year + 1L, "-01-01"))$
      select(band)
    n <- ic$size()$getInfo()
    if (is.null(n) || n == 0) return(NA_real_)
    img  <- ic$first()
    stat <- img$reduceRegion(
      reducer   = rgee::ee$Reducer$mean(),
      geometry  = geom,
      scale     = .PECAHN_NTL_SCALE,
      maxPixels = 1e13,
      bestEffort = TRUE
    )
    out <- stat$get(band)$getInfo()
    if (is.null(out)) NA_real_ else as.numeric(out)
  }

  # V22 covers recent years only (2022+); V21 covers 2013-2021. An empty
  # collection returns NA (not an error), so the fallback fills earlier years.
  val <- tryCatch(mean_for_asset(asset), error = function(e) {
    warning("ee_get_nighttime_light (", asset, ") failed for year ", year, ": ",
            conditionMessage(e), call. = FALSE)
    NA_real_
  })
  if (is.na(val)) {
    val <- tryCatch(mean_for_asset(asset_fallback), error = function(e) {
      warning("ee_get_nighttime_light (", asset_fallback, ") failed for year ",
              year, ": ", conditionMessage(e), call. = FALSE)
      NA_real_
    })
  }

  # Only cache successes — never poison the cache with a transient NA.
  if (!is.na(val)) saveRDS(val, cache_path)
  val
}

# -----------------------------------------------------------------------------
# GHSL Degree-of-Urbanisation (GHS-SMOD) urban-center extraction (internal).
#
# Used by search_urban_center_progressive() to locate the nearest urban center.
# Asset 'JRC/GHSL/P2023A/GHS_SMOD_V2-0' is an ImageCollection whose per-epoch
# images are addressed directly as '.../{YEAR}' (e.g. '.../2020'); band
# 'smod_code', 1000 m resolution. Degree-of-Urbanisation level-2 codes:
# 30 = Urban Center, 23 = Dense Urban Cluster, 22 = Semi-dense Urban Cluster,
# 21 = Suburban/Peri-urban, 13 = Rural Cluster, 12 = Low Density Rural,
# 11 = Very Low Density Rural, 10 = Water. (Replaces the deprecated
# 'JRC/GHSL/P2023A/GHS_SMOD'.)
# -----------------------------------------------------------------------------
.PECAHN_GHSL_ASSET <- "JRC/GHSL/P2023A/GHS_SMOD_V2-0"
.PECAHN_GHSL_SCALE <- 1000
.GHSL_EPOCHS       <- seq(1975L, 2030L, by = 5L)

# Snap an arbitrary year to the nearest available GHS-SMOD epoch (messages if it changes).
.snap_ghsl_year <- function(year) {
  year <- as.integer(round(as.numeric(year)))
  if (year %in% .GHSL_EPOCHS) return(year)
  snapped <- as.integer(.GHSL_EPOCHS[which.min(abs(.GHSL_EPOCHS - year))])
  message(sprintf("  GHS-SMOD has no %d epoch; snapping to nearest available epoch %d.", year, snapped))
  snapped
}

#' Urban-center pixel centroids for a bounding box from GHS-SMOD (internal)
#'
#' Returns an `sf` object of point centroids of contiguous GHS-SMOD urban-center
#' regions (pixels whose `smod_code` is in `codes`) within `bbox` for `ghsl_year`,
#' or `NULL` if none are present. Non-empty results are cached in `gee_cache/`.
#'
#' @param bbox Bounding box (named left/bottom/right/top or xmin/ymin/xmax/ymax).
#' @param ghsl_year Integer GHS-SMOD epoch (snapped to the nearest valid epoch).
#' @param codes Integer vector of `smod_code` values defining "urban center".
#' @param asset Optional character override for the GHS-SMOD asset ID.
#' @param cache_dir Directory for cached results (default `"gee_cache"`).
#'
#' @return An `sf` POINT object (WGS84) of urban-center centroids, or `NULL`.
#' @keywords internal
#' @noRd
ee_get_ghsl_urban_centers <- function(bbox, ghsl_year, codes = c(23, 30),
                                      asset = NULL, cache_dir = "gee_cache") {
  asset <- asset %||% .PECAHN_GHSL_ASSET
  codes <- as.integer(codes)
  year  <- .snap_ghsl_year(ghsl_year)
  xy    <- .bbox_xy(bbox)

  cache_path <- .gee_cache_path(
    list(type = "ghsl", asset = asset, year = year, codes = sort(codes), region = xy),
    ext = ".rds", cache_dir = cache_dir
  )
  if (file.exists(cache_path)) {
    v <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(v)) return(v)
  }

  result <- tryCatch({
    geom <- .ee_rectangle(xy)
    smod <- rgee::ee$Image(paste0(asset, "/", year))$select("smod_code")

    # Build a mask for the requested urban classes (smod_code %in% codes).
    mask <- NULL
    for (cd in codes) {
      m <- smod$eq(cd)
      mask <- if (is.null(mask)) m else mask$Or(m)
    }
    urban <- smod$updateMask(mask)

    # Vectorise contiguous urban regions to their centroids.
    fc <- urban$reduceToVectors(
      geometry       = geom,
      scale          = .PECAHN_GHSL_SCALE,
      geometryType   = "centroid",
      eightConnected = TRUE,
      maxPixels      = 1e10,
      bestEffort     = TRUE
    )

    n <- fc$size()$getInfo()
    if (is.null(n) || n == 0) {
      NULL
    } else {
      # Pull the centroids via getInfo() and build the sf directly. This avoids
      # rgee::ee_as_sf(), which depends on the rgee session-info file and the
      # Drive/GCS export machinery.
      gj     <- fc$getInfo()
      coords <- lapply(gj$features, function(f) f$geometry$coordinates)
      coords <- Filter(function(z) length(z) >= 2, coords)
      if (length(coords) == 0) {
        NULL
      } else {
        lon <- vapply(coords, function(z) as.numeric(z[[1]]), numeric(1))
        lat <- vapply(coords, function(z) as.numeric(z[[2]]), numeric(1))
        sf::st_as_sf(data.frame(lon = lon, lat = lat),
                     coords = c("lon", "lat"), crs = 4326)
      }
    }
  }, error = function(e) {
    warning("ee_get_ghsl_urban_centers failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })

  if (!is.null(result)) saveRDS(result, cache_path)
  result
}
