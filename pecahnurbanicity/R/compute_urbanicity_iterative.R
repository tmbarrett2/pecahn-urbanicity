#' Compute urbanicity metrics for multiple communities
#'
#' @description
#' Iteratively applies [compute_urbanicity()] across a list of communities, generating a combined data frame
#' of urbanicity indicators for each location. Displays a progress bar during execution. Supports multi-year
#' analysis for population density (WorldPop) and nighttime lights (Earth Observation Group). Supports
#' separate bounding boxes for local village-level measures and regional distance measures.
#'
#' @param local_bboxes A named list of community data (each containing `bbox`, `project`, `ethnicity`, and `community`).
#'   This represents the local bounding box for population density, nighttime lights, and infrastructure counts.
#' @param regional_bboxes Optional named list of community data with larger bounding boxes for distance/travel
#'   time calculations. If NULL, uses local_bboxes for all measures.
#' @param metrics Character vector specifying which indicators to compute (`"all"` or subset).
#' @param search_buffer Degrees to expand friction surface cropping area (default = 1).
#' @param friction_surface_path Path to local friction surface raster.
#' @param population_raster_paths Character vector of paths to population density rasters (e.g., WorldPop) for multiple years.
#' @param nighttime_light_paths Character vector of paths to nighttime light rasters (e.g., VIIRS) for multiple years.
#' @param verbose Logical; if `TRUE`, prints intermediate output from each community run (default = `FALSE`).
#'
#' @return
#' A combined `data.frame` with one row per community and columns for each calculated metric.
#' Population density columns are named `pop_density_YYYY` and nighttime light columns are named `nighttime_light_YYYY`.
#'
#' @details
#' Years are automatically extracted from raster filenames using the pattern `_YYYY_`.
#' 
#' The local bounding box (communities_list) is used for: population density, nighttime lights, building density,
#' and counts of shops, financial services, transport stops, and cell towers.
#' 
#' The regional bounding box (regional_communities_list) is used for: travel time to healthcare, schools, 
#' paved roads, and urban centers.
#'
#' @examples
#' \dontrun{
#' # Create both local (1km) and regional (5km) bounding boxes
#' bbox_1km <- create_bounding_boxes(test_data, distance_km = 1)
#' bbox_5km <- create_bounding_boxes(test_data, distance_km = 5)
#' 
#' results <- compute_urbanicity_iterative(
#'   local_bboxes = bbox_1km,
#'   regional_bboxes = bbox_5km,
#'   metrics = "all",
#'   friction_surface_path = "friction_surface_walking.geotiff",
#'   population_raster_paths = c(
#'     "worldpop/global_pop_2015_CN_1km_R2025A_UA_v1.tif",
#'     "worldpop/global_pop_2020_CN_1km_R2025A_UA_v1.tif",
#'     "worldpop/global_pop_2025_CN_1km_R2025A_UA_v1.tif"
#'   ),
#'   nighttime_light_paths = c(
#'     "viirs/VNL_npp_2015_global_vcmslcfg_v2.tif",
#'     "viirs/VNL_npp_2020_global_vcmslcfg_v2.tif",
#'     "viirs/VNL_npp_2024_global_vcmslcfg_v2.tif"
#'   )
#' )
#' }
#'
#' @export
compute_urbanicity_iterative <- function(local_bboxes,
                                         regional_bboxes = NULL,
                                         metrics = c("all"),
                                         search_buffer = 1,  # degrees (~100 km)
                                         friction_surface_path = NULL,
                                         population_raster_paths = NULL,
                                         nighttime_light_paths = NULL,
                                         verbose = FALSE) {
  
  # Determine which metrics to compute
  if (identical(metrics, "all") || "all" %in% metrics) {
    do_roads <- do_shops <- do_healthcare <- do_transport <- do_financial <- TRUE
    do_schools <- do_urban_center <- do_cell_towers <- do_buildings <- TRUE
    do_nighttime_light <- do_population <- TRUE
  } else {
    do_roads           <- "roads"           %in% metrics
    do_shops           <- "shops"           %in% metrics
    do_healthcare      <- "healthcare"      %in% metrics
    do_transport       <- "transport"       %in% metrics
    do_financial       <- "financial"       %in% metrics
    do_schools         <- "schools"         %in% metrics
    do_urban_center    <- "urban_center"    %in% metrics
    do_cell_towers     <- "cell_towers"     %in% metrics
    do_buildings       <- "buildings"       %in% metrics
    do_nighttime_light <- "nighttime_light" %in% metrics
    do_population      <- "population"      %in% metrics
  }
  
  # Load friction surface once if needed for distance calculations
  friction_global <- NULL
  if ((do_healthcare || do_schools || do_roads || do_urban_center) && 
      !is.null(friction_surface_path) && file.exists(friction_surface_path)) {
    if (verbose) cat("Loading global friction surface (one-time load)...\n")
    friction_global <- raster::raster(friction_surface_path)
    if (verbose) cat("Friction surface loaded\n")
  }
  
  # Initialize progress bar
  total <- length(local_bboxes)
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  
  # Run sequentially across all communities
  all_results <- vector("list", total)
  for (i in seq_along(local_bboxes)) {
    # Get regional bbox for this community if provided
    regional_data <- if (!is.null(regional_bboxes)) {
      regional_bboxes[[i]]
    } else {
      NULL
    }
    
    all_results[[i]] <- compute_urbanicity(
      local_bboxes[[i]],
      regional_bbox = regional_data,
      name = names(local_bboxes)[i],
      roads = do_roads,
      shops = do_shops,
      healthcare = do_healthcare,
      transport = do_transport,
      financial = do_financial,
      schools = do_schools,
      urban_center = do_urban_center,
      cell_towers = do_cell_towers,
      buildings = do_buildings,
      nighttime_light = do_nighttime_light,
      population = do_population,
      search_buffer = search_buffer,
      friction_surface_path = friction_surface_path,
      friction_surface_global = friction_global,  # Pass pre-loaded raster
      population_raster_paths = population_raster_paths,
      nighttime_light_paths = nighttime_light_paths,
      verbose = verbose
    )
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Combine results safely across all communities
  result_names <- unique(unlist(lapply(all_results, names)))
  all_results_df <- lapply(all_results, function(x) {
    x <- as.data.frame(x)
    for (col in result_names) if (!col %in% names(x)) x[[col]] <- NA
    x[, result_names, drop = FALSE]
  })
  combined_results <- do.call(rbind, all_results_df)
  
  return(combined_results)
}