#' Compute urbanicity metrics for a single community
#'
#' @description
#' Computes multiple geospatial and infrastructural indicators of urbanicity for a single community,
#' including road access, facilities, travel times, building density, nighttime light intensity, and
#' population density. The function uses OpenStreetMap data combined with raster-based surfaces for
#' travel time and population estimates.
#'
#' @param community_data A list containing community metadata (`bbox`, `project`, `ethnicity`, `community`).
#' @param name Optional character string naming the community.
#' @param roads,shops,healthcare,transport,financial,schools,urban_center,cell_towers,buildings,nighttime_light,population Logical flags controlling which metrics are computed.
#' @param search_buffer Numeric; degrees to expand the bounding box when cropping friction surfaces (default = 1).
#' @param friction_surface_path Path to the friction surface raster used for travel time calculations.
#' @param population_raster_path Path to the population density raster (e.g., GPW).
#' @param nighttime_light_path Path to the nighttime light raster (e.g., VIIRS).
#' @param verbose Logical; if `TRUE`, prints progress messages (default = `FALSE`).
#'
#' @return
#' A one-row `data.frame` containing numeric and categorical urbanicity metrics for the specified community.
#'
#' @details
#' Travel time calculations rely on `gdistance::costDistance()` and a precomputed transition matrix derived
#' from the friction surface raster. OSM queries are automatically cached to minimize repeated downloads.
#'
#' @examples
#' \dontrun{
#' results <- compute_urbanicity(
#'   community_data = community_list[["Mandena"]],
#'   friction_surface_path = "friction_surface_walking.geotiff",
#'   population_raster_path = "pop_raster.tif",
#'   nighttime_light_path = "nighttime_lights.tif"
#' )
#' }
#'
#' @importFrom sf st_polygon st_sfc st_point st_crs st_transform st_area st_sf
#' @importFrom raster raster extent crop crs extract projectRaster
#' @importFrom gdistance transition geoCorrection
#' @importFrom stats aggregate
#' @importFrom methods as
#' @export
compute_urbanicity <- function(community_data, name = NULL,
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
                               population_raster_path = NULL,
                               nighttime_light_path = NULL,
                               verbose = FALSE) {
  
  # Redirect output if not verbose
  if (!verbose) {
    utils::sink(tempfile())
    on.exit(utils::sink(), add = TRUE)
  }
  
  # Ensure required packages are installed
  required_packages <- c("osmdata", "sf", "raster", "gdistance", "httr", "terra", "digest")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "needed for compute_urbanicity()."))
    }
  }
  
  # Extract bounding box and metadata
  bbox <- community_data$bbox
  if (is.null(name)) {
    name <- paste(community_data$project, community_data$ethnicity, community_data$community, sep = " - ")
  }
  if (verbose) cat("Processing community:", name, "\n")
  
  osm_bbox   <- c(bbox["left"], bbox["bottom"], bbox["right"], bbox["top"])
  center_lon <- (bbox["left"] + bbox["right"]) / 2
  center_lat <- (bbox["bottom"] + bbox["top"]) / 2
  utm_zone   <- floor((center_lon + 180) / 6) + 1
  hemisphere <- ifelse(center_lat < 0, "+south", "")
  utm_crs    <- paste0("+proj=utm +zone=", utm_zone, " ", hemisphere, " +datum=WGS84")
  if (verbose) cat("  Location UTM zone:", utm_zone, hemisphere, "\n")
  
  poly <- sf::st_polygon(list(matrix(
    c(bbox["left"], bbox["bottom"],
      bbox["left"], bbox["top"],
      bbox["right"], bbox["top"],
      bbox["right"], bbox["bottom"],
      bbox["left"], bbox["bottom"]),
    ncol = 2, byrow = TRUE
  ))) |> sf::st_sfc(crs = 4326)
  
  center_point <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
  
  # Initialize result container
  results <- list(project = community_data$project,
                  ethnicity = community_data$ethnicity,
                  community = community_data$community)
  
  # Load friction raster and transition
  friction_raster <- NULL
  tr_corrected <- NULL
  raster_crs <- NULL
  
  if (healthcare || schools || roads || urban_center) {
    tryCatch({
      if (verbose) cat("  Loading friction surface data...\n")
      if (is.null(friction_surface_path) || !file.exists(friction_surface_path))
        stop("Local friction surface file required but not found.")
      
      friction_global <- raster::raster(friction_surface_path)
      extent_buffered <- raster::extent(
        bbox["left"] - search_buffer,
        bbox["right"] + search_buffer,
        bbox["bottom"] - search_buffer,
        bbox["top"] + search_buffer
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
  
  # Roads (ratios + travel time)
  if (roads) {
    tryCatch({
      if (verbose) cat("  Analyzing roads...\n")
      roads_data <- get_osm_data_cached(osm_bbox, key = "highway", value = NULL)
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
            center_point, bbox, "paved_road", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 20
          )
        }
      } else {
        results$paved_to_unpaved_ratio <- NA
        results$pct_paved_roads <- NA
        results$travel_time_paved_road_min <- search_facilities_progressive(
          center_point, bbox, "paved_road", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 20
        )
      }
    }, error = function(e) {
      if (verbose) cat("  Error processing roads:", e$message, "\n")
      results$paved_to_unpaved_ratio <- NA
      results$pct_paved_roads <- NA
      results$travel_time_paved_road_min <- NA
    })
  }
  
  # Shops
  if (shops) {
    tryCatch({
      if (verbose) cat("  Analyzing shops and markets...\n")
      shops_data <- get_osm_data_cached(osm_bbox, key = "shop", value = NULL)
      results$n_shops <- if (!is.null(shops_data$osm_points) && nrow(shops_data$osm_points) > 0) nrow(shops_data$osm_points) else 0
    }, error = function(e) { if (verbose) cat("  Error processing shops:", e$message, "\n"); results$n_shops <- NA })
  }
  
  # Transport
  if (transport) {
    tryCatch({
      if (verbose) cat("  Analyzing transport infrastructure...\n")
      transport_count <- 0
      t1 <- get_osm_data_cached(osm_bbox, key = "public_transport", value = NULL)
      if (!is.null(t1$osm_points)) transport_count <- transport_count + nrow(t1$osm_points)
      t2 <- get_osm_data_cached(osm_bbox, key = "highway", value = "bus_stop")
      if (!is.null(t2$osm_points)) transport_count <- transport_count + nrow(t2$osm_points)
      t3 <- get_osm_data_cached(osm_bbox, key = "railway", value = c("station", "halt"))
      if (!is.null(t3$osm_points))   transport_count <- transport_count + nrow(t3$osm_points)
      if (!is.null(t3$osm_polygons)) transport_count <- transport_count + nrow(t3$osm_polygons)
      results$n_transport_stops <- transport_count
    }, error = function(e) { if (verbose) cat("  Error processing transport infrastructure:", e$message, "\n"); results$n_transport_stops <- NA })
  }
  
  # Healthcare
  if (healthcare) {
    tryCatch({
      if (verbose) cat("  Analyzing healthcare facilities...\n")
      results$travel_time_healthcare_min <- search_facilities_progressive(
        center_point, bbox, "healthcare", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 50
      )
    }, error = function(e) { if (verbose) cat("  Error processing healthcare facilities:", e$message, "\n"); results$travel_time_healthcare_min <- NA })
  }
  
  # Financial
  if (financial) {
    tryCatch({
      if (verbose) cat("  Analyzing financial services...\n")
      fin_data <- get_osm_data_cached(osm_bbox, key = "amenity", value = c("bank", "atm"))
      results$n_financial <- if (!is.null(fin_data$osm_points) && nrow(fin_data$osm_points) > 0) nrow(fin_data$osm_points) else 0
    }, error = function(e) { if (verbose) cat("  Error processing financial services:", e$message, "\n"); results$n_financial <- NA })
  }
  
  # Schools
  if (schools) {
    tryCatch({
      if (verbose) cat("  Analyzing schools...\n")
      results$travel_time_school_min <- search_facilities_progressive(
        center_point, bbox, "school", tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 50
      )
    }, error = function(e) { if (verbose) cat("  Error processing schools:", e$message, "\n"); results$travel_time_school_min <- NA })
  }
  
  # Urban centers
  if (urban_center) {
    tryCatch({
      if (verbose) cat("  Analyzing travel time to nearest urban center (place = city)...\n")
      results$travel_time_urban_center_min <- search_urban_center_progressive(
        center_point, bbox, tr_corrected = tr_corrected, raster_crs = raster_crs, max_search_multiplier = 50
      )
    }, error = function(e) { if (verbose) cat("  Error processing urban centers:", e$message, "\n"); results$travel_time_urban_center_min <- NA })
  }
  
  # Cell towers
  if (cell_towers) {
    tryCatch({
      if (verbose) cat("  Analyzing cell towers...\n")
      cell_data <- get_osm_data_cached(osm_bbox, key = "communication", value = "mobile_phone")
      results$n_cell_towers <- if (!is.null(cell_data$osm_points) && nrow(cell_data$osm_points) > 0) nrow(cell_data$osm_points) else 0
    }, error = function(e) { if (verbose) cat("  Error processing cell towers:", e$message, "\n"); results$n_cell_towers <- NA })
  }
  
  # Buildings
  if (buildings) {
    tryCatch({
      if (verbose) cat("  Analyzing building footprints...\n")
      buildings_data <- get_osm_data_cached(osm_bbox, key = "building", value = NULL)
      if (!is.null(buildings_data$osm_polygons) && nrow(buildings_data$osm_polygons) > 0) {
        build_proj <- sf::st_transform(buildings_data$osm_polygons, crs = sf::st_crs(utm_crs))
        build_proj$area <- sf::st_area(build_proj)
        total_building_area <- sum(as.numeric(build_proj$area), na.rm = TRUE)
        poly_projected <- sf::st_transform(poly, crs = sf::st_crs(utm_crs))
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
  
  # Nighttime light
  if (nighttime_light) {
    tryCatch({
      if (verbose) cat("  Analyzing nighttime light...\n")
      viirs_raster <- raster::raster(nighttime_light_path)
      poly_for_nl  <- sf::st_sf(geometry = poly); poly_for_nl$id <- 1
      nl <- get_light_mean(poly_for_nl, viirs_raster)[1]
      results$nighttime_light <- as.numeric(nl)
    }, error = function(e) { if (verbose) cat("  Error processing nighttime light:", e$message, "\n"); results$nighttime_light <- NA })
  }
  
  # Population
  if (population) {
    tryCatch({
      if (verbose) cat("  Calculating population density...\n")
      if (!is.null(population_raster_path) && file.exists(population_raster_path)) {
        pop_global <- raster::raster(population_raster_path)
        extent_comm <- raster::extent(bbox["left"], bbox["right"], bbox["bottom"], bbox["top"])
        pop_r <- raster::crop(pop_global, extent_comm)
        poly_projected <- sf::st_transform(poly, crs = sf::st_crs(utm_crs))
        poly_sf <- sf::st_sf(geometry = poly_projected)
        area_km2 <- as.numeric(sf::st_area(poly_projected)) / 1e6
        pop_proj <- raster::projectRaster(pop_r, crs = raster::crs(poly_sf))
        pop_values <- raster::extract(pop_proj, poly_sf, fun = sum, na.rm = TRUE)
        total_pop <- sum(pop_values, na.rm = TRUE)
        results$pop_density <- as.numeric(total_pop / area_km2)
      } else results$pop_density <- NA
    }, error = function(e) { if (verbose) cat("  Error calculating population density:", e$message, "\n"); results$pop_density <- NA })
  }
  
  # Return (as 1-row data.frame)
  return(as.data.frame(lapply(results, function(x) if (length(x) == 0) NA else x)))
}
