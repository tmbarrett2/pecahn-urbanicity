#' Compute composite urbanicity index from core features.
#'
#' @description
#' Computes a composite urbanicity index across a set of core features using
#' either z-score standardization or min-max normalization. Standardization
#' parameters are always derived from the PEcAHN global reference dataset
#' (1,400 globally stratified communities) pre-computed at package build time,
#' anchoring scores to a fixed global scale that is comparable across study
#' sites, datasets, and future data collection waves.
#'
#' The composite score is the row-wise mean of standardized feature values.
#' Travel time variables are reverse-coded so that higher scores consistently
#' reflect greater urbanicity (shorter travel times = higher score). For
#' multi-year variables (nighttime light, population density), the most recent
#' available year column is used automatically. Reference parameters for these
#' features are keyed by stub name (e.g., \code{"nighttime_light"}) rather than
#' by year, so input data from any year is correctly standardized
#' against the reference distribution regardless of which year the reference
#' was built from.
#'
#' Core features included in the composite:
#' \describe{
#'   \item{`pct_paved_roads`}{Percent of roads that are paved (higher = more urban)}
#'   \item{`travel_time_hospital_min`}{Travel time to nearest hospital in minutes (lower = more urban, reverse-coded)}
#'   \item{`travel_time_school_min`}{Travel time to nearest school in minutes (lower = more urban, reverse-coded)}
#'   \item{`travel_time_urban_center_min`}{Travel time to nearest urban center (lower = more urban, reverse-coded)}
#'   \item{`nighttime_light`}{Nighttime light intensity, most recent year available (higher = more urban)}
#'   \item{`pop_density`}{Population density, most recent year available (higher = more urban)}
#' }
#'
#' @param data A `data.frame` with one row per community, as returned by
#'   [compute_urbanicity()] or [compute_urbanicity_iterative()]. Must contain
#'   at least some of the core feature columns listed above. For multi-year
#'   variables, supply either a bare column (e.g., `nighttime_light`) or
#'   year-suffixed columns (e.g., `nighttime_light_2022`, `nighttime_light_2024`);
#'   the most recent year present is used automatically. The year in `data` does
#'   not need to match the year used when building the reference parameters.
#' @param features Optional character vector of column name stubs to include in
#'   the composite. Defaults to all six core features. Use stubs (e.g.,
#'   `"nighttime_light"`) for multi-year variables — the most recent available
#'   year column is selected automatically.
#' @param method Character string specifying the standardization method.
#'   One of:
#'   \describe{
#'     \item{`"zscore"`}{(default) Each feature is standardized using the
#'       reference mean and SD. The composite is the row-wise mean of z-scores,
#'       interpretable as standard deviations above/below the global reference
#'       mean.}
#'     \item{`"minmax"`}{Each feature is normalized to \[0, 1\] using the
#'       reference min and max. Input values outside the reference range are
#'       clamped to \[0, 1\] with a warning (see `clamp`).}
#'   }
#' @param clamp Logical; only relevant when `method = "minmax"`. If `TRUE`
#'   (default), input values outside the reference \[min, max\] range are
#'   clamped to \[0, 1\] with a warning. If `FALSE`, out-of-range values
#'   produce scores outside \[0, 1\], preserving their relative position on
#'   the global scale.
#' @param min_features Integer; minimum number of non-missing features required
#'   to compute a composite score. Communities below this threshold receive
#'   `NA` with a warning (default `6L`).
#' @param cor_warning Logical; if `TRUE`, warns when any feature pair has
#'   Pearson |r| >= `cor_threshold` in `data` (default `TRUE`).
#' @param cor_threshold Numeric; correlation threshold for the collinearity
#'   warning (default `0.7`). Only used when `cor_warning = TRUE`.
#' @param suffix Character string naming the output composite column
#'   (default `"urban_index"`).
#' @param keep_standardized Logical; if `TRUE`, the standardized version of
#'   each feature is appended to the returned data frame (default `FALSE`).
#'   Columns are suffixed `_z` (z-score) or `_norm` (min-max).
#' @param verbose Logical; if `TRUE`, prints per-feature diagnostics including
#'   the resolved column name, reference parameters used, and whether the
#'   feature was reverse-coded (default `FALSE`).
#'
#' @return
#' The input `data.frame` with the following columns appended:
#' \describe{
#'   \item{`urban_index`}{Composite urbanicity score (row-wise mean of
#'     standardized features). Higher values = more urban. For
#'     `method = "zscore"`, units are standard deviations from the global
#'     reference mean. For `method = "minmax"`, scores are on a \[0, 1\]
#'     scale anchored to the global reference range.}
#'   \item{`n_features_used`}{Integer count of non-missing features that
#'     contributed to each community's composite score. Communities below
#'     `min_features` receive `NA`.}
#' }
#' If `keep_standardized = TRUE`, per-feature standardized columns are
#' appended with suffix `_z` (z-score) or `_norm` (min-max).
#'
#' @details
#' ## Reference-anchored standardization
#'
#' Standardization parameters are pre-computed from the PEcAHN global
#' reference dataset at package build time and stored internally. This means:
#'
#' \enumerate{
#'   \item **Fixed scale.** The parameters never change based on the
#'     composition of `data`. Adding or removing communities from your study
#'     dataset does not shift scores.
#'   \item **Cross-site comparability.** Urbanicity scores computed for
#'     different field sites, at different times, or with different sample
#'     sizes are directly comparable.
#'   \item **Year-agnostic multi-year lookup.** Reference parameters for
#'     `nighttime_light` and `pop_density` are stored under their stub names,
#'     not the specific year used when building the reference. Input data
#'     from any collection year is correctly matched.
#' }
#'
#' To inspect the reference parameters (including which year was used to
#' build them), call [urbanicity_ref_params()].
#'
#' ## Z-score method (`method = "zscore"`)
#' \deqn{z_i = \frac{x_i - \mu_{ref}}{\sigma_{ref}}}
#' Composite: \eqn{\text{urban\_index}_i = \frac{1}{K_i} \sum_k z_{ik}}
#'
#' ## Min-max method (`method = "minmax"`)
#' \deqn{x_{norm,i} = \frac{x_i - \min_{ref}}{\max_{ref} - \min_{ref}}}
#'
#' ## Multi-year variables
#' For `nighttime_light` and `pop_density`, if year-suffixed columns are
#' present (e.g., `nighttime_light_2020`, `nighttime_light_2024`), the column
#' with the highest year suffix is selected. A bare column (e.g.,
#' `nighttime_light`) is used if no year-suffixed columns exist. The resolved
#' column is always looked up in the reference parameters under the stub name
#' (e.g., `"nighttime_light"`), so the year in `data` never needs to match
#' the year the reference was built from.
#'
#' @examples
#' \dontrun{
#' # Z-score composite (default, recommended)
#' results_indexed <- compute_urbanicity_composite_index(results)
#'
#' # Min-max composite with clamping
#' results_indexed <- compute_urbanicity_composite_index(results, method = "minmax", clamp = TRUE)
#'
#' # Inspect standardized feature scores
#' results_indexed <- compute_urbanicity_composite_index(results, keep_standardized = TRUE)
#'
#' # Check what reference parameters are being used
#' urbanicity_ref_params()
#' }
#'
#' @seealso [compute_urbanicity()], [compute_urbanicity_iterative()],
#'   [urbanicity_ref_params()]
#'
#' @importFrom stats cor sd
#' @export
compute_urbanicity_composite_index <- function(data,
                                features = c(
                                  "pct_paved_roads",
                                  "travel_time_hospital_min",
                                  "travel_time_school_min",
                                  "travel_time_urban_center_min",
                                  "nighttime_light",
                                  "pop_density"
                                ),
                                method            = c("zscore", "minmax"),
                                clamp             = TRUE,
                                min_features      = 6L,
                                cor_warning       = TRUE,
                                cor_threshold     = 0.7,
                                suffix            = "urban_index",
                                keep_standardized = FALSE,
                                verbose           = FALSE) {

  method <- match.arg(method)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (nrow(data) < 1) {
    warning("`data` has no rows. Returning unchanged.")
    data[[suffix]]            <- NA_real_
    data[["n_features_used"]] <- NA_integer_
    return(data)
  }

  # .urbanicity_ref_params is loaded from R/sysdata.rda at package load time.
  # Parameters for multi-year stubs are stored under the stub name
  # (e.g., "nighttime_light"), never under a year-specific column name.
  ref_params <- .urbanicity_ref_params[[method]]

  reverse_vars    <- c(
    "travel_time_paved_road_min",
    "travel_time_hospital_min",
    "travel_time_school_min",
    "travel_time_urban_center_min"
  )
  multiyear_stubs <- c("nighttime_light", "pop_density")
  std_suffix      <- if (method == "zscore") "_z" else "_norm"

  # ---------------------------------------------------------------------------
  # Resolve feature columns in `data`
  # ---------------------------------------------------------------------------
  resolved      <- .resolve_feature_cols(data, features, multiyear_stubs, verbose)
  data          <- resolved$data
  resolved_cols <- resolved$cols   # actual column names, e.g. "nighttime_light_2024"
  stub_map      <- resolved$stubs  # named vector: col -> stub, e.g. "nighttime_light_2024" -> "nighttime_light"

  if (length(resolved_cols) == 0)
    stop("No valid feature columns found in `data`.")

  # ---------------------------------------------------------------------------
  # Correlation diagnostic
  # ---------------------------------------------------------------------------
  if (cor_warning && length(resolved_cols) >= 2) {
    feat_matrix <- as.matrix(data[, resolved_cols, drop = FALSE])
    storage.mode(feat_matrix) <- "double"
    cor_matrix <- suppressWarnings(
      cor(feat_matrix, use = "pairwise.complete.obs")
    )
    pairs <- which(
      abs(cor_matrix) >= cor_threshold & upper.tri(cor_matrix),
      arr.ind = TRUE
    )
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

  # ---------------------------------------------------------------------------
  # Standardize features
  # ---------------------------------------------------------------------------
  std_cols <- character(0)

  for (col in resolved_cols) {
    x     <- as.numeric(data[[col]])
    n_obs <- sum(!is.na(x))

    if (n_obs < 3) {
      warning("Feature '", col, "' has fewer than 3 non-missing values -- ",
              "excluding from composite.")
      next
    }

    # Always look up via stub so year mismatches between input and reference
    # are never a problem (e.g. input has nighttime_light_2024, reference was
    # built from nighttime_light_2025 -- both resolve to stub "nighttime_light")
    ref_key <- stub_map[[col]]

    if (!ref_key %in% names(ref_params)) {
      warning(
        "No reference parameters found for '", col, "' (stub: '", ref_key, "'). ",
        "Excluding from composite. Ensure the feature is present in the PEcAHN ",
        "reference dataset, then rebuild R/sysdata.rda."
      )
      next
    }

    center  <- ref_params[[ref_key]]$center
    scale   <- ref_params[[ref_key]]$scale
    std_col <- paste0(col, std_suffix)

    is_reverse <- any(sapply(reverse_vars, function(rv) startsWith(col, rv)))

    if (method == "zscore") {
      z               <- (x - center) / scale
      data[[std_col]] <- if (is_reverse) -z else z

    } else {
      nm <- (x - center) / scale
      if (clamp) {
        n_below <- sum(nm < 0, na.rm = TRUE)
        n_above <- sum(nm > 1, na.rm = TRUE)
        if (n_below + n_above > 0) {
          warning(
            "Feature '", col, "': ", n_below, " value(s) below reference min ",
            "and ", n_above, " value(s) above reference max. ",
            "Clamping to [0, 1]. Set clamp = FALSE to allow out-of-range scores."
          )
          nm <- pmin(pmax(nm, 0), 1)
        }
      }
      data[[std_col]] <- if (is_reverse) 1 - nm else nm
    }

    if (verbose) {
      cat(sprintf(
        "  %-40s  ref_key = %-30s  center = %8.3f  scale = %8.3f%s\n",
        col, ref_key, center, scale,
        if (is_reverse) "  (reverse-coded)" else ""
      ))
    }

    std_cols <- c(std_cols, std_col)
  }

  if (length(std_cols) == 0)
    stop("No features could be standardized. Check that columns exist and ",
         "contain numeric data matching the core feature set.")

  # ---------------------------------------------------------------------------
  # Composite: row-wise mean of standardized features
  # ---------------------------------------------------------------------------
  std_matrix            <- as.matrix(data[, std_cols, drop = FALSE])
  storage.mode(std_matrix) <- "double"

  n_available           <- rowSums(!is.na(std_matrix))
  data[[suffix]]        <- rowMeans(std_matrix, na.rm = TRUE)

  below_threshold                  <- n_available < min_features
  data[[suffix]][below_threshold]  <- NA_real_
  data[[suffix]][n_available == 0] <- NA_real_
  data[["n_features_used"]]        <- as.integer(n_available)

  if (any(below_threshold)) {
    warning(
      sum(below_threshold), " community/communities had fewer than ",
      min_features, " non-missing features and received NA. ",
      "Adjust `min_features` to change this threshold."
    )
  }

  if (verbose) {
    cat("  Method:", method, "| Standardization: PEcAHN reference (n =",
        .urbanicity_ref_params$meta$n_communities, ")\n")
    cat("  Features used (", length(std_cols), "):",
        paste(gsub(std_suffix, "", std_cols), collapse = ", "), "\n")
    cat("  Score range: [",
        round(min(data[[suffix]], na.rm = TRUE), 3), ",",
        round(max(data[[suffix]], na.rm = TRUE), 3), "]\n")
    cat("  Communities below min_features threshold:", sum(below_threshold), "\n")
  }

  if (!keep_standardized) data[, std_cols] <- NULL

  return(data)
}


# ==============================================================================
# urbanicity_ref_params() — provenance / inspection helper
# ==============================================================================

#' Inspect the pre-computed urbanicity reference parameters.
#'
#' @description
#' Returns and optionally prints the standardization parameters derived from
#' the PEcAHN global reference dataset that are used internally by
#' [compute_urbanicity_composite_index()]. Useful for documenting methods, checking which
#' year was used for multi-year variables, or verifying that the package was
#' built against the expected reference dataset.
#'
#' @param method One of `"zscore"` (default) or `"minmax"`.
#' @param print Logical; if `TRUE` (default), prints a formatted summary.
#'
#' @return Invisibly returns a list with elements `params` and `meta`.
#'
#' @examples
#' urbanicity_ref_params()
#' p <- urbanicity_ref_params(method = "minmax", print = FALSE)
#'
#' @seealso [compute_urbanicity_composite_index()]
#' @export
urbanicity_ref_params <- function(method = c("zscore", "minmax"),
                                  print  = TRUE) {
  method <- match.arg(method)
  params <- .urbanicity_ref_params[[method]]
  meta   <- .urbanicity_ref_params$meta

  if (print) {
    cat(sprintf(
      "PEcAHN urbanicity reference parameters  [method = %s]\n", method
    ))
    cat(sprintf(
      "Reference dataset: n = %d communities | Generated: %s\n",
      meta$n_communities,
      format(meta$generated_at, "%Y-%m-%d %H:%M %Z")
    ))
    cat(sprintf(
      "Multi-year source columns: %s\n\n",
      paste(grep("[0-9]{4}$", meta$resolved_cols, value = TRUE), collapse = ", ")
    ))
    label <- if (method == "zscore") c("mean", "SD") else c("min", "range")
    cat(sprintf("  %-35s  %10s  %10s\n", "Feature (ref key)", label[1], label[2]))
    cat(strrep("-", 59), "\n")
    for (nm in names(params)) {
      cat(sprintf("  %-35s  %10.4f  %10.4f\n",
                  nm, params[[nm]]$center, params[[nm]]$scale))
    }
  }

  invisible(list(params = params, meta = meta))
}


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Resolve feature stubs to actual column names; returns cols and stub mapping.
#'
#' For multi-year stubs, selects the most recent year column present in `data`.
#' Also returns a named character vector mapping each resolved column name back
#' to its canonical stub, which is the key used for reference param lookup.
#'
#' @return list(data = data, cols = character vector, stubs = named character vector)
#' @keywords internal
.resolve_feature_cols <- function(data, features, multiyear_stubs, verbose) {
  resolved_cols <- character(0)
  stub_map      <- character(0)  # col_name -> stub (or col_name itself)

  for (feat in features) {
    if (feat %in% multiyear_stubs) {
      pattern   <- paste0("^", feat, "_[0-9]{4}$")
      year_cols <- grep(pattern, colnames(data), value = TRUE)

      if (length(year_cols) > 0) {
        years  <- as.integer(regmatches(year_cols, regexpr("[0-9]{4}$", year_cols)))
        chosen <- year_cols[which.max(years)]
        resolved_cols        <- c(resolved_cols, chosen)
        stub_map[chosen]     <- feat   # "nighttime_light_2024" -> "nighttime_light"
        if (verbose) {
          cat(sprintf("  Multi-year stub '%s': using '%s' (year %d) -> ref key '%s'\n",
                      feat, chosen, max(years), feat))
        }
      } else if (feat %in% colnames(data)) {
        resolved_cols    <- c(resolved_cols, feat)
        stub_map[feat]   <- feat       # bare column; stub is itself
        if (verbose) cat("  Using bare column:", feat, "-> ref key:", feat, "\n")
      } else {
        warning("Feature '", feat, "' not found in data -- skipping.")
      }

    } else {
      if (feat %in% colnames(data)) {
        resolved_cols  <- c(resolved_cols, feat)
        stub_map[feat] <- feat         # non-multi-year; stub is itself
      } else {
        warning("Feature '", feat, "' not found in data -- skipping.")
      }
    }
  }

  list(data = data, cols = resolved_cols, stubs = stub_map)
}
