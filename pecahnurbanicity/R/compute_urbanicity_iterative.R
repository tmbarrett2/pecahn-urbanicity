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
#' @param search_buffer Degrees to expand the friction surface fetch/crop area (default = 1).
#' @param ee_project Optional character; the Google Cloud project passed to `rgee::ee_Initialize()`. Earth Engine is
#'   initialized once at the start of the run; supply this if you have not already initialized Earth Engine in your
#'   session and have no default project configured.
#' @param population_years Integer vector of years for which to compute population density from WorldPop
#'   (e.g., `c(2015, 2020)`). If `NULL`, population density is skipped.
#' @param nighttime_light_years Integer vector of years for which to compute mean nighttime light from VIIRS.
#'   If `NULL`, nighttime light is skipped.
#' @param friction_asset,population_asset,nighttime_light_asset Optional character overrides for the Earth Engine
#'   asset IDs. Leave `NULL` to use the package defaults.
#' @param urban_center_codes Integer vector of GHS-SMOD `smod_code` values defining an "urban center".
#'   Defaults to `c(23, 30)` (Dense Urban Cluster + Urban Center); use `30` for strict Urban Center only.
#' @param ghsl_year Integer; GHS-SMOD epoch used for the urban-center search (5-year epochs 1975-2030,
#'   default = 2020). Non-epoch years are snapped to the nearest available epoch.
#' @param verbose Logical; if `TRUE`, prints intermediate output from each community run (default = `FALSE`).
#'
#' @return
#' A combined `data.frame` with one row per community and columns for each calculated metric.
#' Population density columns are named `pop_density_YYYY` and nighttime light columns are named `nighttime_light_YYYY`.
#' Each travel-time distance measure also has a companion `*_is_boundary_est` logical flag column (see
#' [compute_urbanicity()]).
#'
#' @details
#' Raster data (friction surface, population, nighttime lights) is retrieved from Google Earth Engine via `rgee`,
#' which must be authenticated beforehand. The friction surface is fetched per community (clipped on the Earth
#' Engine side) and cached locally in `gee_cache/`. Population and nighttime-light years are taken from
#' `population_years` and `nighttime_light_years`.
#'
#' The local bounding box (communities_list) is used for: population density, nighttime lights, building density,
#' and counts of shops, financial services, transport stops, and cell towers.
#' 
#' The regional bounding box (regional_communities_list) is used for: travel time to hospital, schools, 
#' paved roads, and urban centers.
#'
#' @examples
#' \dontrun{
#' # Earth Engine must be authenticated first: rgee::ee_Initialize()
#'
#' # Create both local (1km) and regional (5km) bounding boxes
#' bbox_1km <- create_bounding_boxes(test_data, distance_km = 1)
#' bbox_5km <- create_bounding_boxes(test_data, distance_km = 5)
#'
#' results <- compute_urbanicity_iterative(
#'   local_bboxes = bbox_1km,
#'   regional_bboxes = bbox_5km,
#'   metrics = "all",
#'   ee_project = "my-gcp-project",
#'   population_years = c(2015, 2020, 2021),
#'   nighttime_light_years = c(2015, 2020, 2024)
#' )
#' }
#'
#' @export
compute_urbanicity_iterative <- function(local_bboxes,
                                         regional_bboxes = NULL,
                                         metrics = c("all"),
                                         search_buffer = 1,  # degrees (~100 km)
                                         ee_project = NULL,
                                         population_years = NULL,
                                         nighttime_light_years = NULL,
                                         friction_asset = NULL,
                                         population_asset = NULL,
                                         nighttime_light_asset = NULL,
                                         urban_center_codes = c(23, 30),
                                         ghsl_year = 2020,
                                         verbose = FALSE) {
  
  # Determine which metrics to compute
  if (identical(metrics, "all") || "all" %in% metrics) {
    do_roads <- do_shops <- do_hospital <- do_transport <- do_financial <- TRUE
    do_schools <- do_urban_center <- do_cell_towers <- do_buildings <- TRUE
    do_nighttime_light <- do_population <- TRUE
  } else {
    do_roads           <- "roads"           %in% metrics
    do_shops           <- "shops"           %in% metrics
    do_hospital      <- "hospital"      %in% metrics
    do_transport       <- "transport"       %in% metrics
    do_financial       <- "financial"       %in% metrics
    do_schools         <- "schools"         %in% metrics
    do_urban_center    <- "urban_center"    %in% metrics
    do_cell_towers     <- "cell_towers"     %in% metrics
    do_buildings       <- "buildings"       %in% metrics
    do_nighttime_light <- "nighttime_light" %in% metrics
    do_population      <- "population"      %in% metrics
  }
  
  # Initialize Earth Engine once up front if any EE-backed metric is requested.
  # The friction surface is fetched per community (and cached in gee_cache/),
  # so no global raster is pre-loaded here.
  need_ee <- (do_hospital || do_schools || do_roads || do_urban_center) ||
    (do_population && !is.null(population_years)) ||
    (do_nighttime_light && !is.null(nighttime_light_years))
  if (need_ee) {
    if (verbose) cat("Initializing Earth Engine...\n")
    .ee_ensure_init(ee_project)
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
      hospital = do_hospital,
      transport = do_transport,
      financial = do_financial,
      schools = do_schools,
      urban_center = do_urban_center,
      cell_towers = do_cell_towers,
      buildings = do_buildings,
      nighttime_light = do_nighttime_light,
      population = do_population,
      search_buffer = search_buffer,
      ee_project = ee_project,
      population_years = population_years,
      nighttime_light_years = nighttime_light_years,
      friction_asset = friction_asset,
      population_asset = population_asset,
      nighttime_light_asset = nighttime_light_asset,
      urban_center_codes = urban_center_codes,
      ghsl_year = ghsl_year,
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