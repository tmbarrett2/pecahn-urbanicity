#' Compute urbanicity metrics for a single community
#'
#' @description
#' Computes multiple geospatial and infrastructural indicators of urbanicity for a single community,
#' including road access, facilities, travel times, building density, nighttime light intensity, and
#' population density. The function uses OpenStreetMap data combined with raster-based surfaces for
#' travel time and population estimates. Supports multi-year analysis for population density (WorldPop)
#' and nighttime lights (Earth Observation Group). Supports separate bounding boxes for local measures
#' (population, infrastructure) and regional measures (travel times to facilities).
#'
#' @param community_data A list containing community metadata (`bbox`, `project`, `ethnicity`, `community`).
#'   This bbox is used as the local_bbox for population density, nighttime lights, and infrastructure counts.
#' @param regional_bbox Optional list containing a larger bounding box for distance/travel time calculations.
#'   If NULL, uses the bbox from community_data for all measures.
#' @param name Optional character string naming the community.
#' @param roads,shops,healthcare,transport,financial,schools,urban_center,cell_towers,buildings,nighttime_light,population Logical flags controlling which metrics are computed.
#' @param search_buffer Numeric; degrees to expand the bounding box when cropping friction surfaces (default = 1).
#' @param friction_surface_path Path to the friction surface raster used for travel time calculations.
#' @param friction_surface_global Optional pre-loaded global friction surface raster. If provided, this will be used
#'   instead of loading from friction_surface_path. This parameter is typically used internally by
#'   [compute_urbanicity_iterative()] to avoid reloading the same raster for multiple communities.
#' @param population_raster_paths Character vector of paths to population density rasters (e.g., WorldPop) for multiple years.
#' @param nighttime_light_paths Character vector of paths to nighttime light rasters (e.g., VIIRS) for multiple years.
#' @param verbose Logical; if `TRUE`, prints progress messages (default = `FALSE`).
#'
#' @return
#' A one-row `data.frame` containing numeric and categorical urbanicity metrics for the specified community.
#' Population density columns are named `pop_density_YYYY` and nighttime light columns are named `nighttime_light_YYYY`.
#'
#' @details
#' Travel time calculations rely on `gdistance::costDistance()` and a precomputed transition matrix derived
#' from the friction surface raster. OSM queries are automatically cached to minimize repeated downloads.
#' Years are automatically extracted from raster filenames using the pattern `_YYYY_`.
#' 
#' The local_bbox (from community_data) is used for: population density, nighttime lights, building density,
#' and counts of shops, financial services, transport stops, and cell towers.
#' 
#' The regional_bbox is used for: travel time to healthcare, schools, paved roads, and urban centers.
#'
#' @examples
#' \dontrun{
#' # Create both local (1km) and regional (5km) bounding boxes
#' bbox_1km <- create_bounding_boxes(test_data, distance_km = 1)
#' bbox_5km <- create_bounding_boxes(test_data, distance_km = 5)
#' 
#' results <- compute_urbanicity(
#'   community_data = bbox_1km[["Mandena"]],
#'   regional_bbox = bbox_5km[["Mandena"]],
#'   friction_surface_path = "friction_surface_walking.geotiff",
#'   population_raster_paths = c("worldpop_2015.tif", "worldpop_2020.tif"),
#'   nighttime_light_paths = c("viirs_2015.tif", "viirs_2020.tif")
#' )
#' }
#'
#' @importFrom sf st_polygon st_sfc st_point st_crs st_transform st_area st_sf
#' @importFrom raster raster extent crop crs extract projectRaster
#' @importFrom gdistance transition geoCorrection
#' @importFrom stats aggregate
#' @importFrom methods as
#' @export
compute_urbanicity <- function(community_data, 
                               regional_bbox = NULL,
                               name = NULL,
                               roads = TRUE,
                               shops = TRUE,
                               healthcare = TRUE,
                               transport = TRUE,
                               financial = TRUE,
                               schools = TRUE,
                               urban_center = TRUE,
                               cell_towers = TRUE,
                               buildings = TRUE,
                               nighttime_light = TRUE,
                               population = TRUE,
                               search_buffer = 1,  # degrees (~100km)
                               friction_surface_path = NULL,
                               friction_surface_global = NULL,  # NEW: pre-loaded global raster
                               population_raster_paths = NULL,
                               nighttime_light_paths = NULL,
                               verbose = FALSE) {
  
  # Redirect output if not verbose
  if (!verbose) {
    sink(tempfile())
    on.exit(sink(), add = TRUE)
  }
  
  # Ensure required packages are installed
  required_packages <- c("osmdata", "sf", "raster", "gdistance", "httr", "terra", "digest")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "needed for compute_urbanicity()."))
    }
  }
  
  # Extract local bounding box (for population, nighttime light, infrastructure)
  local_bbox <- community_data$bbox
  
  # Use regional bbox for distance measures if provided, otherwise use local
  regional_bbox_use <- if (!is.null(regional_bbox)) regional_bbox$bbox else local_bbox
  
  if (is.null(name)) {
    name <- paste(community_data$project, community_data$ethnicity, community_data$community, sep = " - ")
  }
  if (verbose) cat("Processing community:", name, "\n")
  
  # OSM bbox for local measures (infrastructure counts)
  osm_bbox_local   <- c(local_bbox["left"], local_bbox["bottom"], local_bbox["right"], local_bbox["top"])
  
  # OSM bbox for regional measures (distance searches)
  osm_bbox_regional <- c(regional_bbox_use["left"], regional_bbox_use["bottom"], 
                         regional_bbox_use["right"], regional_bbox_use["top"])
  
  # Center point (use local bbox)
  center_lon <- (local_bbox["left"] + local_bbox["right"]) / 2
  center_lat <- (local_bbox["bottom"] + local_bbox["top"]) / 2
  utm_zone   <- floor((center_lon + 180) / 6) + 1
  hemisphere <- ifelse(center_lat < 0, "+south", "")
  utm_crs    <- paste0("+proj=utm +zone=", utm_zone, " ", hemisphere, " +datum=WGS84")
  if (verbose) cat("  Location UTM zone:", utm_zone, hemisphere, "\n")
  
  # Polygon for local bbox (population, nighttime light, buildings)
  poly_local <- sf::st_polygon(list(matrix(
    c(local_bbox["left"], local_bbox["bottom"],
      local_bbox["left"], local_bbox["top"],
      local_bbox["right"], local_bbox["top"],
      local_bbox["right"], local_bbox["bottom"],
      local_bbox["left"], local_bbox["bottom"]),
    ncol = 2, byrow = TRUE
  ))) |> sf::st_sfc(crs = 4326)
  
  center_point <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
  
  # Initialize result container
  results <- list(project = community_data$project,
                  ethnicity = community_data$ethnicity,
                  community = community_data$community)
  
  # Load friction raster and transition (use regional bbox for extent)
  friction_raster <- NULL
  tr_corrected <- NULL
  raster_crs <- NULL
  
  if (healthcare || schools || roads || urban_center) {
    tryCatch({
      if (verbose) cat("  Loading friction surface data...\n")
      
      # Use pre-loaded global raster if provided, otherwise load it
      if (!is.null(friction_surface_global)) {
        friction_global <- friction_surface_global
        if (verbose) cat("  Using pre-loaded friction surface\n")
      } else if (!is.null(friction_surface_path) && file.exists(friction_surface_path)) {
        friction_global <- raster::raster(friction_surface_path)
        if (verbose) cat("  Loaded friction surface from file\n")
      } else {
        stop("Local friction surface file required but not found.")
      }
      
      extent_buffered <- raster::extent(
        regional_bbox_use["left"] - search_buffer,
        regional_bbox_use["right"] + search_buffer,
        regional_bbox_use["bottom"] - search_buffer,
        regional_bbox_use["top"] + search_buffer
      )
      friction_raster <- raster::crop(friction_global, extent_buffered)
      raster_crs <- raster::crs(friction_raster)
      tr <- gdistance::transition(friction_raster, function(x) 1 / mean(x), directions = 8)
      tr_corrected <- gdistance::geoCorrection(tr, type = "c")
      if (verbose) cat("  Friction surface loaded and transition built\n")
    }, error = function(e) {
      if (verbose) cat("  Error loading friction surface:", e$message, "\n")
      tr_corrected <- NULL
    })
  }
  
  # Roads (ratios + travel time) - uses regional bbox for distance search
  if (roads) {
    tryCatch({
      if (verbose) cat("  Analyzing roads...\n")
      roads_data <- get_osm_data_cached(osm_bbox_local, key = "highway", value = NULL)
      paved_surfaces   <- c("paved", "asphalt", "concrete")
      unpaved_surfaces <- c("unpaved", "dirt", "gravel", "fine_gravel", "sand", "grass", 
                            "ground", "earth", "mud", "compacted", "unclassified")
      
      if (!is.null(roads_data$osm_lines) && nrow(roads_data$osm_lines) > 0) {
        roads_sf <- roads_data$osm_lines
        if (!"surface" %in% colnames(roads_sf)) roads_sf$surface <- NA
        roads_sf$road_type <- "other"
        if (!all(is.na(roads_sf$surface))) {
          paved_mask   <- !is.na(roads_sf$surface) & roads_sf$surface %in% paved_surfaces
          unpaved_mask <- !is.na(roads_sf$surface) & roads_sf$surface %in% unpaved_surfaces
          if (any(paved_mask))   roads_sf$road_type[paved_mask]   <- "paved"
          if (any(unpaved_mask)) roads_sf$road_type[unpaved_mask] <- "unpaved"
        }
        if ("highway" %in% colnames(roads_sf)) {
          unpaved_highway_types <- c("track", "path")
          hmask <- !is.na(roads_sf$highway) & roads_sf$highway %in% unpaved_highway_types
          if (any(hmask)) roads_sf$road_type[hmask] <- "unpaved"
        }
        roads_sf <- sf::st_transform(roads_sf, crs = sf::st_crs(utm_crs))
        roads_sf$length <- as.numeric(sf::st_length(roads_sf))
        summary <- stats::aggregate(roads_sf$length, by = list(type = roads_sf$road_type), FUN = sum)
        total <- sum(summary$x, na.rm = TRUE)
        paved <- sum(summary$x[summary$type == "paved"], na.rm = TRUE)
        unpaved <- sum(summary$x[summary$type == "unpaved"], na.rm = TRUE)
        results$paved_to_unpaved_ratio <- if (unpaved > 0) paved / unpaved else NA
        results$pct_paved_roads <- if (total > 0) 100 * (paved / total) else NA
        
        paved_roads <- roads_sf[roads_sf$road_type == "paved", ]
        if (nrow(paved_roads) > 0 && !is.null(tr_corrected)) {
          paved_points <- suppressWarnings(sf::st_cast(sf::st_geometry(paved_roads), "POINT") |> sf::st_sf())
          results$travel_time_paved_road_min <- calculate_travel_time(center_point, paved_points, tr_corrected, raster_crs)
        } else {
          results$travel_time_paved_road_min <- search_facilities_progressive(
            center_point, regional_bbox_use, "paved_road", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 20
          )
        }
      } else {
        results$paved_to_unpaved_ratio <- NA
        results$pct_paved_roads <- NA
        results$travel_time_paved_road_min <- search_facilities_progressive(
          center_point, regional_bbox_use, "paved_road", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 20
        )
      }
    }, error = function(e) {
      if (verbose) cat("  Error processing roads:", e$message, "\n")
      results$paved_to_unpaved_ratio <- NA
      results$pct_paved_roads <- NA
      results$travel_time_paved_road_min <- NA
    })
  }
  
  # Shops - uses local bbox
  if (shops) {
    tryCatch({
      if (verbose) cat("  Analyzing shops and markets...\n")
      shops_data <- get_osm_data_cached(osm_bbox_local, key = "shop", value = NULL)
      results$n_shops <- if (!is.null(shops_data$osm_points) && nrow(shops_data$osm_points) > 0) nrow(shops_data$osm_points) else 0
    }, error = function(e) { if (verbose) cat("  Error processing shops:", e$message, "\n"); results$n_shops <- NA })
  }
  
  # Transport - uses local bbox
  if (transport) {
    tryCatch({
      if (verbose) cat("  Analyzing transport infrastructure...\n")
      transport_count <- 0
      t1 <- get_osm_data_cached(osm_bbox_local, key = "public_transport", value = NULL)
      if (!is.null(t1$osm_points)) transport_count <- transport_count + nrow(t1$osm_points)
      t2 <- get_osm_data_cached(osm_bbox_local, key = "highway", value = "bus_stop")
      if (!is.null(t2$osm_points)) transport_count <- transport_count + nrow(t2$osm_points)
      t3 <- get_osm_data_cached(osm_bbox_local, key = "railway", value = c("station", "halt"))
      if (!is.null(t3$osm_points))   transport_count <- transport_count + nrow(t3$osm_points)
      if (!is.null(t3$osm_polygons)) transport_count <- transport_count + nrow(t3$osm_polygons)
      results$n_transport_stops <- transport_count
    }, error = function(e) { if (verbose) cat("  Error processing transport infrastructure:", e$message, "\n"); results$n_transport_stops <- NA })
  }
  
  # Healthcare - uses regional bbox for search
  if (healthcare) {
    tryCatch({
      if (verbose) cat("  Analyzing healthcare facilities...\n")
      results$travel_time_healthcare_min <- search_facilities_progressive(
        center_point, regional_bbox_use, "healthcare", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 50
      )
    }, error = function(e) { if (verbose) cat("  Error processing healthcare facilities:", e$message, "\n"); results$travel_time_healthcare_min <- NA })
  }
  
  # Financial - uses local bbox
  if (financial) {
    tryCatch({
      if (verbose) cat("  Analyzing financial services...\n")
      fin_data <- get_osm_data_cached(osm_bbox_local, key = "amenity", value = c("bank", "atm"))
      results$n_financial <- if (!is.null(fin_data$osm_points) && nrow(fin_data$osm_points) > 0) nrow(fin_data$osm_points) else 0
    }, error = function(e) { if (verbose) cat("  Error processing financial services:", e$message, "\n"); results$n_financial <- NA })
  }
  
  # Schools - uses regional bbox for search
  if (schools) {
    tryCatch({
      if (verbose) cat("  Analyzing schools...\n")
      results$travel_time_school_min <- search_facilities_progressive(
        center_point, regional_bbox_use, "school", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 50
      )
    }, error = function(e) { if (verbose) cat("  Error processing schools:", e$message, "\n"); results$travel_time_school_min <- NA })
  }
  
  # Urban centers - uses regional bbox for search
  if (urban_center) {
    tryCatch({
      if (verbose) cat("  Analyzing travel time to nearest urban center (place = city)...\n")
      results$travel_time_urban_center_min <- search_urban_center_progressive(
        center_point, regional_bbox_use, tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 50
      )
    }, error = function(e) { if (verbose) cat("  Error processing urban centers:", e$message, "\n"); results$travel_time_urban_center_min <- NA })
  }
  
  # Cell towers - uses local bbox
  if (cell_towers) {
    tryCatch({
      if (verbose) cat("  Analyzing cell towers...\n")
      cell_data <- get_osm_data_cached(osm_bbox_local, key = "communication", value = "mobile_phone")
      results$n_cell_towers <- if (!is.null(cell_data$osm_points) && nrow(cell_data$osm_points) > 0) nrow(cell_data$osm_points) else 0
    }, error = function(e) { if (verbose) cat("  Error processing cell towers:", e$message, "\n"); results$n_cell_towers <- NA })
  }
  
  # Buildings - uses local bbox
  if (buildings) {
    tryCatch({
      if (verbose) cat("  Analyzing building footprints...\n")
      buildings_data <- get_osm_data_cached(osm_bbox_local, key = "building", value = NULL)
      if (!is.null(buildings_data$osm_polygons) && nrow(buildings_data$osm_polygons) > 0) {
        build_proj <- sf::st_transform(buildings_data$osm_polygons, crs = sf::st_crs(utm_crs))
        build_proj$area <- sf::st_area(build_proj)
        total_building_area <- sum(as.numeric(build_proj$area), na.rm = TRUE)
        poly_projected <- sf::st_transform(poly_local, crs = sf::st_crs(utm_crs))
        bbox_area <- as.numeric(sf::st_area(poly_projected))
        if (bbox_area == 0 || !is.finite(total_building_area) || !is.finite(bbox_area)) {
          results$building_density <- NA
        } else {
          bden <- (total_building_area / bbox_area) * 100
          results$building_density <- if (is.finite(bden) && bden >= 0 && bden <= 100) as.numeric(bden) else NA
        }
      } else results$building_density <- 0
    }, error = function(e) { if (verbose) cat("  Error processing buildings:", e$message, "\n"); results$building_density <- NA })
  }
  
  # Nighttime light (multi-year) - uses local bbox
  if (nighttime_light && !is.null(nighttime_light_paths)) {
    if (verbose) cat("  Analyzing nighttime light (multi-year)...\n")
    poly_for_nl <- sf::st_sf(geometry = poly_local); poly_for_nl$id <- 1
    
    for (nl_path in nighttime_light_paths) {
      tryCatch({
        # Extract year from filename
        year_match <- regmatches(nl_path, regexpr("_[0-9]{4}_", nl_path))
        if (length(year_match) > 0) {
          year <- gsub("_", "", year_match[1])
          col_name <- paste0("nighttime_light_", year)
          
          if (file.exists(nl_path)) {
            viirs_raster <- raster::raster(nl_path)
            nl <- get_light_mean(poly_for_nl, viirs_raster)[1]
            results[[col_name]] <- as.numeric(nl)
            if (verbose) cat("    Processed nighttime light for year", year, "\n")
          } else {
            results[[col_name]] <- NA
            if (verbose) cat("    File not found for year", year, "\n")
          }
        } else {
          if (verbose) cat("    Could not extract year from filename:", nl_path, "\n")
        }
      }, error = function(e) {
        if (verbose) cat("    Error processing nighttime light file:", nl_path, "-", e$message, "\n")
      })
    }
  }
  
  # Population (multi-year) - uses local bbox
  if (population && !is.null(population_raster_paths)) {
    if (verbose) cat("  Calculating population density (multi-year)...\n")
    poly_projected <- sf::st_transform(poly_local, crs = sf::st_crs(utm_crs))
    poly_sf <- sf::st_sf(geometry = poly_projected)
    area_km2 <- as.numeric(sf::st_area(poly_projected)) / 1e6
    
    for (pop_path in population_raster_paths) {
      tryCatch({
        # Extract year from filename
        year_match <- regmatches(pop_path, regexpr("_[0-9]{4}_", pop_path))
        if (length(year_match) > 0) {
          year <- gsub("_", "", year_match[1])
          col_name <- paste0("pop_density_", year)
          
          if (file.exists(pop_path)) {
            pop_global <- raster::raster(pop_path)
            extent_comm <- raster::extent(local_bbox["left"], local_bbox["right"], local_bbox["bottom"], local_bbox["top"])
            pop_r <- raster::crop(pop_global, extent_comm)
            pop_proj <- raster::projectRaster(pop_r, crs = raster::crs(poly_sf))
            pop_values <- raster::extract(pop_proj, poly_sf, fun = sum, na.rm = TRUE)
            total_pop <- sum(pop_values, na.rm = TRUE)
            results[[col_name]] <- as.numeric(total_pop / area_km2)
            if (verbose) cat("    Processed population for year", year, "\n")
          } else {
            results[[col_name]] <- NA
            if (verbose) cat("    File not found for year", year, "\n")
          }
        } else {
          if (verbose) cat("    Could not extract year from filename:", pop_path, "\n")
        }
      }, error = function(e) {
        if (verbose) cat("    Error processing population file:", pop_path, "-", e$message, "\n")
      })
    }
  }
  
  # Return (as 1-row data.frame)
  return(as.data.frame(lapply(results, function(x) if (length(x) == 0) NA else x)))
}