#' Create temporal trend plots for multi-year variables
#'
#' @description
#' Generates line plots showing temporal trends for multi-year variables (population density
#' and nighttime light intensity) across communities. Automatically detects variables with
#' year suffixes (e.g., `pop_density_2015`, `nighttime_light_2020`) and creates separate
#' plots for each metric type.
#'
#' @param data A data frame output from [compute_urbanicity_iterative()] containing
#'   multi-year variables.
#'
#' @return
#' A named list of `ggplot` objects, one for each temporal metric type detected
#' (e.g., `pop_density`, `nighttime_light`). Returns an empty list if no temporal
#' variables are found.
#'
#' @details
#' The function identifies temporal variables by detecting patterns like `varname_YYYY`
#' where YYYY is a four-digit year. Each community is plotted as a separate line with
#' color determined by project affiliation. Points are added at each year to improve
#' readability.
#'
#' @examples
#' \dontrun{
#' # Compute urbanicity with multi-year data
#' results <- compute_urbanicity_iterative(
#'   bbox_5km,
#'   population_raster_paths = c("pop_2015.tif", "pop_2020.tif", "pop_2025.tif"),
#'   nighttime_light_paths = c("nl_2015.tif", "nl_2020.tif", "nl_2024.tif")
#' )
#' 
#' # Create temporal plots
#' temporal_plots <- create_temporal_plots(results)
#' temporal_plots$pop_density
#' temporal_plots$nighttime_light
#' }
#'
#' @importFrom dplyr select starts_with mutate group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal theme scale_color_viridis_d scale_x_continuous scale_y_continuous facet_wrap coord_cartesian expansion
#' @importFrom magrittr %>%
#' @export
create_temporal_plots <- function(data) {
  # Avoid R CMD check notes for NSE variables
  year <- value <- community_label <- mean_value <- NULL
  
  # Check for required columns
  if (!all(c("project", "community") %in% names(data))) {
    stop("Data must contain 'project' and 'community' columns")
  }
  
  # Identify temporal variable patterns (varname_YYYY)
  all_vars <- names(data)
  temporal_pattern <- "_[0-9]{4}$"
  temporal_vars <- grep(temporal_pattern, all_vars, value = TRUE)
  
  if (length(temporal_vars) == 0) {
    warning("No temporal variables detected. Variables should follow pattern: varname_YYYY")
    return(list())
  }
  
  # Extract base variable names (e.g., "pop_density" from "pop_density_2015")
  base_vars <- unique(gsub(temporal_pattern, "", temporal_vars))
  
  # Create plots for each base variable
  all_plots <- list()
  
  for (base_var in base_vars) {
    # Get all columns for this variable
    var_cols <- grep(paste0("^", base_var, "_[0-9]{4}$"), all_vars, value = TRUE)
    
    if (length(var_cols) == 0) next
    
    # Reshape data for plotting
    plot_data <- data %>%
      dplyr::select(project, community, dplyr::all_of(var_cols)) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with(base_var),
        names_to = "year",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        year = as.numeric(gsub(paste0(base_var, "_"), "", year)),
        community_label = paste(project, community, sep = ": ")
      ) %>%
      dplyr::filter(!is.na(value))
    
    # Skip if no valid data
    if (nrow(plot_data) == 0) next
    
    # Calculate project means
    project_means <- plot_data %>%
      dplyr::group_by(project, year) %>%
      dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Create formatted variable name for labels
    var_label <- gsub("_", " ", base_var)
    var_label <- tools::toTitleCase(var_label)
    
    # Add units to y-axis label
    if (base_var == "pop_density") {
      y_label <- paste0(var_label, " (people/km\u00B2)")
    } else if (base_var == "nighttime_light") {
      y_label <- paste0(var_label, " (nW/cm\u00B2/sr)")
    } else {
      y_label <- var_label
    }
    
    # Plot 1: Combined plot with project means (bold) and individual communities (faint)
    p_combined <- ggplot2::ggplot() +
      # Layer 1: Individual communities (faint)
      ggplot2::geom_line(
        data = plot_data,
        ggplot2::aes(x = year, y = value, color = project, group = community_label),
        linewidth = 0.5, 
        alpha = 0.3
      ) +
      # Layer 2: Project means (bold)
      ggplot2::geom_line(
        data = project_means,
        ggplot2::aes(x = year, y = mean_value, color = project, group = project),
        linewidth = 2,
        alpha = 1.0
      ) +
      ggplot2::scale_color_viridis_d(option = "turbo") +
      ggplot2::scale_x_continuous(breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
      ggplot2::labs(
        title = paste("Temporal Trends:", var_label, "(bold = project mean)"),
        x = "Year",
        y = y_label,
        color = "Project"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "top",
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(size = 11),
        axis.text = ggplot2::element_text(size = 10)
      )
    
    # Plot 2: Faceted by project
    p_faceted <- ggplot2::ggplot(plot_data, 
                                 ggplot2::aes(x = year, y = value, 
                                              color = community_label, 
                                              group = community_label)) +
      ggplot2::geom_line(linewidth = 1, alpha = 0.8) +
      ggplot2::facet_wrap(~project, scales = "free_y") +
      ggplot2::scale_color_viridis_d(option = "turbo") +
      ggplot2::scale_x_continuous(breaks = function(x) seq(ceiling(min(x)), floor(max(x)), by = 1)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::coord_cartesian(ylim = c(0, NA)) +
      ggplot2::labs(
        title = paste("Temporal Trends by Project:", var_label),
        x = "Year",
        y = y_label,
        color = "Community"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(size = 11),
        axis.text = ggplot2::element_text(size = 10),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        strip.text = ggplot2::element_text(size = 11, face = "bold")
      )
    
    # Add both plots to output list
    all_plots[[base_var]] <- p_combined
    all_plots[[paste0(base_var, "_faceted")]] <- p_faceted
  }
  
  return(all_plots)
}


#' Create all summary and temporal plots
#'
#' @description
#' Convenience wrapper that generates both cross-sectional summary plots (via
#' [create_summary_plots()]) and temporal trend plots (via [create_temporal_plots()])
#' for urbanicity data.
#'
#' @param data A data frame output from [compute_urbanicity_iterative()].
#'
#' @return
#' A named list containing two sub-lists:
#' \describe{
#'   \item{summary_plots}{Cross-sectional bar plots for each metric}
#'   \item{temporal_plots}{Line plots showing trends over time for multi-year variables}
#' }
#'
#' @details
#' This function provides a one-stop solution for visualizing urbanicity results,
#' automatically detecting which plots are appropriate for the data.
#'
#' @examples
#' \dontrun{
#' all_plots <- create_all_plots(urbanicity_results)
#' 
#' # Access summary plots
#' all_plots$summary_plots$pct_paved_roads
#' 
#' # Access temporal plots
#' all_plots$temporal_plots$pop_density
#' }
#'
#' @export
create_all_plots <- function(data) {
  summary_plots <- create_summary_plots(data)
  temporal_plots <- create_temporal_plots(data)
  
  return(list(
    summary_plots = summary_plots,
    temporal_plots = temporal_plots
  ))
}