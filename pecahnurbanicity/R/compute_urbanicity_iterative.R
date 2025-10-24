#' Compute urbanicity metrics for multiple communities
#'
#' @description
#' Iteratively applies [compute_urbanicity()] across a list of communities, generating a combined data frame
#' of urbanicity indicators for each location. Displays a progress bar during execution.
#'
#' @param communities_list A named list of community data (each containing `bbox`, `project`, `ethnicity`, and `community`).
#' @param metrics Character vector specifying which indicators to compute (`"all"` or subset).
#' @param search_buffer Degrees to expand friction surface cropping area (default = 1).
#' @param friction_surface_path Path to local friction surface raster.
#' @param population_raster_path Path to local population density raster.
#' @param nighttime_light_path Path to local nighttime light raster.
#' @param verbose Logical; if `TRUE`, prints intermediate output from each community run (default = `FALSE`).
#'
#' @return
#' A combined `data.frame` with one row per community and columns for each calculated metric.
#'
#' @examples
#' \dontrun{
#' results <- compute_urbanicity_iterative(
#'   communities_list = bbox_5km,
#'   metrics = "all",
#'   friction_surface_path = "friction_surface_walking.geotiff",
#'   population_raster_path = "pop_raster.tif",
#'   nighttime_light_path = "nighttime_lights.tif"
#' )
#' }
#'
#' @export
compute_urbanicity_iterative <- function(communities_list,
                                         metrics = c("all"),
                                         search_buffer = 1,  # degrees (~100 km)
                                         friction_surface_path = NULL,
                                         population_raster_path = NULL,
                                         nighttime_light_path = NULL,
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
  
  # Initialize progress bar
  total <- length(communities_list)
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  
  # Run sequentially across all communities
  all_results <- vector("list", total)
  for (i in seq_along(communities_list)) {
    all_results[[i]] <- compute_urbanicity(
      communities_list[[i]],
      name = names(communities_list)[i],
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
      population_raster_path = population_raster_path,
      nighttime_light_path = nighttime_light_path,
      verbose = verbose
    )
    utils::setTxtProgressBar(pb, i)
  }
  utils::close(pb)
  
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
