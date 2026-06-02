# Internal package state. Tracks whether Earth Engine has been initialized in
# this R session so that .ee_ensure_init() is a no-op on repeat calls.
.pkg_env <- new.env(parent = emptyenv())
.pkg_env$ee_initialized <- FALSE

#' Ensure Earth Engine is initialized (internal)
#'
#' @description
#' Initializes the `rgee` / Google Earth Engine session exactly once per R
#' session. Called automatically at the top of [compute_urbanicity_iterative()]
#' and [compute_urbanicity()] whenever Earth Engine data is required. Subsequent
#' calls are a no-op once initialization has succeeded.
#'
#' @param ee_project Optional character giving the Google Cloud project passed to
#'   `rgee::ee_Initialize(project = ...)`. If `NULL`, `ee_Initialize()` is called
#'   without an explicit project (relying on any default project configured for
#'   the authenticated user).
#'
#' @details
#' The user must have installed `rgee` and its Python dependencies
#' (`rgee::ee_install()`) and completed Earth Engine authentication beforehand.
#' If Earth Engine is already live in the session (the user ran
#' `rgee::ee_Initialize()` themselves), that session is reused as-is and no
#' re-initialization is attempted. Otherwise `rgee::ee_Initialize()` is called.
#' If initialization fails, an informative error is raised explaining how to
#' authenticate and/or supply a Cloud project.
#'
#' @return Invisibly `TRUE` once Earth Engine is initialized; otherwise stops
#'   with an informative error.
#'
#' @keywords internal
#' @noRd
.ee_ensure_init <- function(ee_project = NULL) {
  # No-op if we have already initialized in this session.
  if (isTRUE(.pkg_env$ee_initialized)) {
    return(invisible(TRUE))
  }

  if (!requireNamespace("rgee", quietly = TRUE)) {
    stop(
      "Package 'rgee' is required for Google Earth Engine access.\n",
      "Install it with install.packages('rgee'), then run rgee::ee_install() ",
      "and authenticate once with rgee::ee_Initialize() before using this package.",
      call. = FALSE
    )
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required by 'rgee' for Google Earth Engine access.\n",
      "Install it with install.packages('reticulate').",
      call. = FALSE
    )
  }

  # If Earth Engine is already live in this session (e.g. the user ran
  # rgee::ee_Initialize() themselves before calling this package), use it as-is.
  # Re-initializing can fail or clobber a perfectly working session, so a cheap
  # liveness probe avoids that.
  already_live <- tryCatch({
    invisible(rgee::ee$Number(1L)$getInfo())
    TRUE
  }, error = function(e) FALSE)
  if (isTRUE(already_live)) {
    .pkg_env$ee_initialized <- TRUE
    return(invisible(TRUE))
  }

  tryCatch({
    if (!is.null(ee_project)) {
      rgee::ee_Initialize(project = ee_project)
    } else {
      rgee::ee_Initialize()
    }
    .pkg_env$ee_initialized <- TRUE
  }, error = function(e) {
    # Non-intrusive diagnostic: report whether reticulate already has a Python
    # session, without forcing Python to initialize (initialize = FALSE).
    py_ready <- tryCatch(
      reticulate::py_available(initialize = FALSE),
      error = function(...) FALSE
    )
    stop(
      "Failed to initialize Earth Engine: ", conditionMessage(e), "\n",
      if (!isTRUE(py_ready)) {
        paste0(
          "No active Python session was found for reticulate; run ",
          "rgee::ee_install() to configure the Earth Engine Python environment. "
        )
      } else "",
      if (is.null(ee_project)) {
        paste0(
          "If you have not configured a default Google Cloud project, supply ",
          "one via the 'ee_project' argument. "
        )
      } else "",
      "Ensure you have run rgee::ee_install() and authenticated with ",
      "rgee::ee_Initialize().",
      call. = FALSE
    )
  })

  invisible(.pkg_env$ee_initialized)
}
