#' Compute composite urbanicity index from core features.
#'
#' @description
#' Computes a composite urbanicity index across a set of core features using
#' either z-score standardization or min-max normalization. The composite is
#' always computed as the row-wise mean of standardized feature values.
#' Travel time variables are reverse-coded so that higher scores
#' consistently reflect greater urbanicity (i.e., shorter travel times = higher
#' score). Multi-year variables (nighttime light, population density) are
#' averaged across available years before standardization when applicable.
#'
#' Core features included in the composite:
#' - `pct_paved_roads` -- percent of roads that are paved (higher = more urban)
#' - `paved_to_unpaved_ratio` -- ratio of paved to unpaved road (higher = more urban)
#' - `travel_time_paved_road_min` -- travel time to nearest paved road (lower = more urban)
#' - `travel_time_healthcare_min` -- travel time to nearest healthcare facility (lower = more urban)
#' - `travel_time_school_min` -- travel time to nearest school (lower = more urban)
#' - `travel_time_urban_center_min` -- travel time to nearest urban center (lower = more urban)
#' - `nighttime_light` -- mean nighttime light intensity, averaged across years if multiple years are included (higher = more urban)
#' - `pop_density` -- population density, averaged across years if multiple years are included (higher = more urban)
#'
#' @param data A `data.frame` with one row per community, as returned by
#'   [compute_urbanicity()] or [compute_urbanicity_iterative()]. Must contain
#'   at least some of the core feature columns listed above.
#'   Multi-year columns (e.g., `nighttime_light_2015`, `pop_density_2020`) are
#'   automatically detected and averaged prior to standardization.
#' @param features Optional character vector of column name stubs to include in
#'   the composite. Defaults to all eight core features. Use stubs (e.g.,
#'   `"nighttime_light"`) for multi-year variables -- all matching yearly columns
#'   will be averaged automatically.
#' @param method Character string specifying the standardization method.
#'   One of:
#'   \describe{
#'     \item{`"zscore"`}{(default) Each feature is z-score standardized
#'       (mean = 0, SD = 1). The composite is the row-wise mean of z-scores,
#'       interpretable as the average standardized deviation from the sample
#'       mean across features.}
#'     \item{`"minmax"`}{Each feature is min-max normalized to \[0, 1\].
#'       The composite is the row-wise mean of normalized values (0-1 scale).}
#'   }
#' @param min_features Integer; minimum number of non-missing features required
#'   to compute a composite score. Communities with fewer non-missing features
#'   than this threshold receive `NA` with a warning (default `6L`).
#' @param cor_warning Logical; if `TRUE`, prints a warning for any pair of
#'   features with Pearson |r| >= `cor_threshold`, which may indicate implicit
#'   double-weighting of a shared dimension (default `TRUE`).
#' @param cor_threshold Numeric; correlation threshold above which a warning
#'   is issued (default `0.7`). Only used when `cor_warning = TRUE`.
#' @param suffix Character string naming the output composite column
#'   (default `"urban_index"`).
#' @param keep_standardized Logical; if `TRUE`, the standardized version of
#'   each feature is appended to the returned data frame for inspection
#'   (default `FALSE`). Columns are suffixed `_z` (z-score) or `_norm`
#'   (min-max).
#' @param verbose Logical; if `TRUE`, prints progress and diagnostics
#'   (default `FALSE`).
#'
#' @return
#' The input `data.frame` with the following columns appended:
#' - `urban_index` -- composite urbanicity score computed as the row-wise mean
#'   of standardized feature values. For `method = "zscore"`, this is the mean
#'   z-score across features (unbounded, mean of 0 across communities, higher =
#'   more urban). For `method = "minmax"`, this is the mean normalized value
#'   (0-1, higher = more urban).
#' - `n_features_used` -- integer count of non-missing features that contributed
#'   to each community's score. Communities below `min_features` receive `NA`.
#' - Intermediate averaged columns (e.g., `nighttime_light_mean`,
#'   `pop_density_mean`) are added when multi-year data are detected.
#' - If `keep_standardized = TRUE`, standardized columns are appended with
#'   suffix `_z` or `_norm` depending on `method`.
#'
#' @details
#' **Z-score method (`method = "zscore"`):**
#' Each feature is standardized as:
#' \deqn{z = \frac{x - \bar{x}}{s}}
#' where \eqn{\bar{x}} is the sample mean and \eqn{s} is the sample standard
#' deviation. The composite is the row-wise mean of z-scores:
#' \deqn{\text{urban\_index}_i = \frac{1}{K_i} \sum_{k=1}^{K_i} z_{ik}}
#' where \eqn{K_i} is the number of non-missing features for community \eqn{i}.
#' Features with zero variance are excluded with a warning.
#'
#' **Min-max method (`method = "minmax"`):**
#' Each feature is normalized as:
#' \deqn{x_{norm} = \frac{x - \min(x)}{\max(x) - \min(x)}}
#' The composite is the row-wise mean of normalized values. Sample-dependent:
#' adding or removing communities changes all scores.
#'
#' **Direction:** Travel time variables are reverse-coded after standardization
#' so that higher scores uniformly indicate greater urbanicity. For z-scores,
#' the sign is negated (`-z`); for min-max, the value is flipped (`1 - x_norm`).
#'
#' **Missing data:** The composite is the mean of available standardized
#' features per community. Communities with fewer than `min_features` non-missing
#' features receive `NA`.
#'
#' **Correlation diagnostics:** When `cor_warning = TRUE`, feature pairs with
#' |r| >= `cor_threshold` are flagged. Highly correlated features effectively
#' upweight a shared latent dimension. Consider removing redundant features or
#' using PCA-based weighting if collinearity is severe.
#'
#' @examples
#' \dontrun{
#' results <- compute_urbanicity_iterative(
#'   communities             = community_list,
#'   friction_surface_path   = "friction_surface_walking.geotiff",
#'   population_raster_paths = c("worldpop_2015.tif", "worldpop_2020.tif"),
#'   nighttime_light_paths   = c("viirs_2015.tif",    "viirs_2020.tif")
#' )
#'
#' # Z-score composite (recommended)
#' results_indexed <- compute_urban_index(results)
#'
#' # Min-max composite
#' results_indexed <- compute_urban_index(results, method = "minmax")
#'
#' # Require at least 5 features for a valid score
#' results_indexed <- compute_urban_index(results, min_features = 5)
#'
#' # Inspect standardized components
#' results_indexed <- compute_urban_index(results, keep_standardized = TRUE)
#' }
#'
#' @seealso [compute_urbanicity()], [compute_urbanicity_iterative()]
#'
#' @importFrom stats cor sd aggregate
#' @export
compute_urbanicity_composite_index <- function(data,
                                features = c(
                                  "pct_paved_roads",
                                  "paved_to_unpaved_ratio",
                                  "travel_time_paved_road_min",
                                  "travel_time_healthcare_min",
                                  "travel_time_school_min",
                                  "travel_time_urban_center_min",
                                  "nighttime_light",
                                  "pop_density"
                                ),
                                method            = c("zscore", "minmax"),
                                min_features      = 6L,
                                cor_warning       = TRUE,
                                cor_threshold     = 0.7,
                                suffix            = "urban_index",
                                keep_standardized = FALSE,
                                verbose           = FALSE) {

  method <- match.arg(method)

  # --- Input validation -------------------------------------------------------
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (nrow(data) < 2) {
    warning(
      "Cannot standardize with fewer than 2 rows. `urban_index` will be NA. ",
      "Run compute_urbanicity_iterative() across multiple communities first."
    )
    data[[suffix]]            <- NA_real_
    data[["n_features_used"]] <- NA_integer_
    return(data)
  }

  reverse_vars    <- c(
    "travel_time_paved_road_min",
    "travel_time_healthcare_min",
    "travel_time_school_min",
    "travel_time_urban_center_min"
  )
  multiyear_stubs <- c("nighttime_light", "pop_density")
  std_suffix      <- if (method == "zscore") "_z" else "_norm"

  # --- Resolve columns --------------------------------------------------------
  resolved_cols <- character(0)

  for (feat in features) {
    if (feat %in% multiyear_stubs) {
      pattern   <- paste0("^", feat, "_[0-9]{4}$")
      year_cols <- grep(pattern, colnames(data), value = TRUE)

      if (length(year_cols) == 0) {
        if (feat %in% colnames(data)) {
          resolved_cols <- c(resolved_cols, feat)
          if (verbose) cat("  Using bare column:", feat, "\n")
        } else {
          warning("Feature '", feat, "' not found in data -- skipping.")
        }
        next
      }

      mean_col <- paste0(feat, "_mean")
      year_matrix <- as.matrix(data[, year_cols, drop = FALSE])
      storage.mode(year_matrix) <- "double"
      data[[mean_col]] <- rowMeans(year_matrix, na.rm = TRUE)
      data[[mean_col]][rowSums(!is.na(year_matrix)) == 0] <- NA_real_
      resolved_cols <- c(resolved_cols, mean_col)
      if (verbose) cat("  Averaged", length(year_cols), "year(s) for",
                       feat, "->", mean_col, "\n")
    } else {
      if (feat %in% colnames(data)) {
        resolved_cols <- c(resolved_cols, feat)
      } else {
        warning("Feature '", feat, "' not found in data -- skipping.")
      }
    }
  }

  if (length(resolved_cols) == 0)
    stop("No valid feature columns found in `data`.")

  # --- Correlation diagnostic -------------------------------------------------
  if (cor_warning && length(resolved_cols) >= 2) {
    feat_matrix <- as.matrix(data[, resolved_cols, drop = FALSE])
    storage.mode(feat_matrix) <- "double"
    cor_matrix  <- suppressWarnings(cor(feat_matrix, use = "pairwise.complete.obs"))
    pairs       <- which(abs(cor_matrix) >= cor_threshold & upper.tri(cor_matrix),
                         arr.ind = TRUE)
    if (nrow(pairs) > 0) {
      pair_strs <- apply(pairs, 1, function(idx) {
        sprintf("  %s & %s (r = %.2f)",
                resolved_cols[idx[1]], resolved_cols[idx[2]],
                cor_matrix[idx[1], idx[2]])
      })
      warning(
        "Highly correlated feature pairs (|r| >= ", cor_threshold, ") -- ",
        "may implicitly double-weight a shared dimension:\n",
        paste(pair_strs, collapse = "\n"), "\n",
        "Consider removing redundant features or using PCA-based weighting."
      )
    }
  }

  # --- Standardize features ---------------------------------------------------
  std_cols <- character(0)

  for (col in resolved_cols) {
    x     <- as.numeric(data[[col]])
    n_obs <- sum(!is.na(x))

    if (n_obs < 3) {
      warning("Feature '", col, "' has fewer than 3 non-missing values -- ",
              "excluding from composite.")
      next
    }

    std_col    <- paste0(col, std_suffix)
    is_reverse <- any(sapply(reverse_vars, function(rv) startsWith(col, rv)))

    if (method == "zscore") {
      x_mean <- mean(x, na.rm = TRUE)
      x_sd   <- sd(x,   na.rm = TRUE)
      if (x_sd == 0) {
        warning("Feature '", col, "' has zero variance -- excluding from composite.")
        next
      }
      z               <- (x - x_mean) / x_sd
      data[[std_col]] <- if (is_reverse) -z else z
      if (verbose) cat("  Z-scored:", col,
                       if (is_reverse) "(reverse-coded)\n" else "\n")

    } else {
      x_min <- min(x, na.rm = TRUE)
      x_max <- max(x, na.rm = TRUE)
      if ((x_max - x_min) == 0) {
        warning("Feature '", col, "' has zero range -- excluding from composite.")
        next
      }
      nm              <- (x - x_min) / (x_max - x_min)
      data[[std_col]] <- if (is_reverse) 1 - nm else nm
      if (verbose) cat("  Min-max normalized:", col,
                       if (is_reverse) "(reverse-coded)\n" else "\n")
    }

    std_cols <- c(std_cols, std_col)
  }

  if (length(std_cols) == 0)
    stop("No features could be standardized. Check that columns contain numeric data.")

  # --- Composite: row-wise mean -------------------------------------
  std_matrix  <- as.matrix(data[, std_cols, drop = FALSE])
  storage.mode(std_matrix) <- "double"

  n_available         <- rowSums(!is.na(std_matrix))
  data[[suffix]]      <- rowMeans(std_matrix, na.rm = TRUE)

  # Apply min_features threshold -- insufficient data -> NA
  below_threshold             <- n_available < min_features
  data[[suffix]][below_threshold] <- NA_real_
  if (any(below_threshold)) {
    warning(
      sum(below_threshold), " community/communities had fewer than ", min_features,
      " non-missing features and received NA. ",
      "Adjust `min_features` to change this threshold."
    )
  }

  # Communities missing everything -> NA
  data[[suffix]][n_available == 0] <- NA_real_
  data[["n_features_used"]]        <- as.integer(n_available)

  if (verbose) {
    cat("  Method: mean of", method, "standardized scores\n")
    cat("  Features (", length(std_cols), "):",
        paste(std_cols, collapse = ", "), "\n")
    cat("  Score range: [",
        round(min(data[[suffix]], na.rm = TRUE), 3), ",",
        round(max(data[[suffix]], na.rm = TRUE), 3), "]\n")
    cat("  Communities below min_features threshold:", sum(below_threshold), "\n")
  }

  # --- Optionally retain standardized columns --------------------------------
  if (!keep_standardized) data[, std_cols] <- NULL

  return(data)
}
