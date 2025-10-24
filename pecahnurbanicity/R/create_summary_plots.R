#' Create numeric bar plots for individual variables
#'
#' @description
#' Internal helper function used by [create_summary_plots()] to generate
#' standardized horizontal bar plots for numeric variables in the urbanicity dataset.
#'
#' @param data A data frame containing numeric variables to be plotted.
#' @param var_name Character string specifying the variable name to plot.
#'
#' @return
#' A single `ggplot` object showing the distribution of the specified variable
#' across all communities.
#'
#' @details
#' Variables beginning with `"n_"` (e.g., counts of shops or facilities) are
#' labeled with whole numbers. All other variables are labeled with three-decimal precision.
#'
#' @examples
#' \dontrun{
#' p <- create_numeric_plot(urbanicity_results, "pct_paved_roads")
#' print(p)
#' }
#'
#' @keywords internal
#' @importFrom dplyr filter mutate arrange desc sym
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_minimal theme scale_y_continuous coord_flip scale_fill_viridis_d expansion
#' @importFrom magrittr %>%
#' @export
create_numeric_plot <- function(data, var_name) {
  # Filter data to non-missing values for this variable
  data_filtered <- dplyr::filter(data, !is.na(!!dplyr::sym(var_name)))
  
  # Choose label precision based on variable type
  format_string <- ifelse(startsWith(var_name, "n_"), "%.0f", "%.3f")
  
  # Build plot
  p <- data_filtered %>%
    dplyr::mutate(community = paste(project, community, sep = ": ")) %>%
    dplyr::arrange(dplyr::desc(!!dplyr::sym(var_name))) %>%
    dplyr::mutate(community = factor(community, levels = community)) %>%
    ggplot2::ggplot(ggplot2::aes(x = community, y = !!dplyr::sym(var_name), fill = project)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf(format_string, !!dplyr::sym(var_name))),
                       hjust = -0.2) +
    ggplot2::labs(title = paste("Summary of", var_name),
                  x = "Community",
                  y = var_name) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::coord_flip()
  
  return(p)
}


#' Create summary plots for all numeric metrics
#'
#' @description
#' Generates a collection of standardized summary bar plots for each numeric variable
#' in the results of [compute_urbanicity_iterative()]. Each plot shows values by community
#' and project, facilitating quick comparisons of urbanicity metrics.
#'
#' @param data A data frame output from [compute_urbanicity_iterative()].
#'
#' @return
#' A named list of `ggplot` objects, one for each numeric variable in the dataset.
#'
#' @details
#' Non-numeric columns such as `project` and `community` are automatically excluded.
#' Each plot is formatted with consistent color scales, labels, and axis expansion.
#'
#' @examples
#' \dontrun{
#' plots <- create_summary_plots(urbanicity_results)
#' plots$pct_paved_roads
#' plots$travel_time_school_min
#' }
#'
#' @importFrom dplyr select_if
#' @export
create_summary_plots <- function(data) {
  # Identify numeric columns (excluding metadata)
  all_vars <- names(data)
  all_vars <- all_vars[!all_vars %in% c("project", "community")]
  
  # Generate plots for numeric variables only
  all_plots <- list()
  for (var in all_vars) {
    if (is.numeric(data[[var]])) {
      all_plots[[var]] <- create_numeric_plot(data, var)
    }
  }
  
  return(all_plots)
}
