# 06_compute_study_sites.R
#
# Computes the urbanicity metrics for the PEcAHN study sites against the SAME
# self-hosted Overpass mirror (OSM planet 2026-06-08) and the same settings as
# the reference run (04_compute_urbanicity_reference.Rmd): 1 km local / 5 km
# regional boxes, metrics = "all", search_buffer = 1, years 2015-2025.
#
# Unattended-friendly:
#   * resumable  - batches already saved to OUT_DIR are skipped on re-run
#   * fault-tolerant - a failed batch is retried one site at a time, so one bad
#     site cannot sink its batch-mates; sites that still fail are logged
#   * isolated caches - runs with OUT_DIR as working dir, so osm_cache/ and
#     gee_cache/ are fresh (no reuse of results from the public Overpass API)
#
# Usage (Rscript is not on PATH on this machine):
#   & "C:\Program Files\R\R-4.4.3\bin\Rscript.exe" 06_compute_study_sites.R          # full run
#   & "C:\Program Files\R\R-4.4.3\bin\Rscript.exe" 06_compute_study_sites.R --check  # preflight only
#
# Output (in OUT_DIR):
#   study_sites_batch_NNN.rds          per-batch results
#   study_sites_failed.csv             sites that failed even individually (if any)
#   study_sites_urbanicity.rds / .csv  final bound table, one row per site

####################
##  CONFIGURATION ##
####################

	SITES_FILE   <- "C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn/manuscript_data/PEcAHN_2026_CommunityInfo_All_9July2026.txt"
	OUT_DIR      <- "C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn/study_sites_results"
	OVERPASS_URL <- "http://localhost:12345/api/interpreter"
	EE_PROJECT   <- "ee-tmbarrett2"
	YEARS        <- 2015:2025
	BATCH_SIZE   <- 5

	CHECK_ONLY <- "--check" %in% commandArgs(trailingOnly = TRUE)

#	Must be set BEFORE reticulate loads, or it grabs a Python without `ee`.
	Sys.setenv(RETICULATE_PYTHON = "C:/Users/tmbar/AppData/Local/r-miniconda/envs/rgee/python.exe")

	suppressPackageStartupMessages({
		library(pecahnurbanicity)
		library(dplyr)
		library(readr)
	})

	stamp <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), sprintf(...), "\n", sep = "")


#######################
##  LOAD STUDY SITES ##
#######################

	raw <- read_tsv(SITES_FILE, show_col_types = FALSE, trim_ws = TRUE)

	sites <- raw |>
		filter(!(is.na(Lat) & is.na(Long))) |>
		transmute(
			project   = trimws(Project),
			community = trimws(Sampling_Location),
			lat       = as.numeric(Lat),
			lon       = as.numeric(Long)
		)

	bad <- is.na(sites$lat) | is.na(sites$lon) | abs(sites$lat) > 90 | abs(sites$lon) > 180
	if (any(bad)) stop(sprintf("%d site(s) have missing/out-of-range coordinates: %s",
		sum(bad), paste(sites$community[bad], collapse = ", ")))
	if (anyDuplicated(sites$community)) stop("Duplicate community names: ",
		paste(unique(sites$community[duplicated(sites$community)]), collapse = ", "))

	stamp("Loaded %d study sites across %d projects.", nrow(sites), n_distinct(sites$project))


####################
##   PREFLIGHT    ##
####################

#	Overpass mirror must answer before we spend hours on Earth Engine calls.
#	Random node id: the server rejects a query identical to one still in flight
#	(OVERPASS_ALLOW_DUPLICATE_QUERIES=no), so a fixed ping can collide with itself.
	ping <- tryCatch({
		con <- url(sprintf("%s?data=%%5Bout%%3Ajson%%5D%%3Bnode%%28%d%%29%%3Bout%%3B",
			OVERPASS_URL, sample.int(1e7, 1)))
		on.exit(close(con), add = TRUE)
		paste(readLines(con, warn = FALSE), collapse = "")
	}, error = function(e) NA_character_)
	overpass_ok <- !is.na(ping) && grepl("\"elements\"", ping)
	stamp("Overpass mirror at %s: %s", OVERPASS_URL, if (overpass_ok) "OK" else "NOT RESPONDING")

	ee <- reticulate::import("ee")
	ee$Initialize(project = EE_PROJECT)
	stamp("Earth Engine initialised (project %s).", EE_PROJECT)

	if (CHECK_ONLY) {
		stamp("--check: site list, package, and Earth Engine OK.%s",
			if (overpass_ok) "" else " (Overpass not up yet - expected before the rebuild.)")
		quit(status = 0)
	}
	if (!overpass_ok) stop("Overpass mirror is not responding; refusing to run against the public API.")

	osmdata::set_overpass_url(OVERPASS_URL)

	dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
	setwd(OUT_DIR)   # keeps osm_cache/ and gee_cache/ local to this run


####################
##  BATCHED RUN   ##
####################

	run_sites <- function(df) {
		compute_urbanicity_iterative(
			local_bboxes          = create_bounding_boxes(df, distance_km = 1),
			regional_bboxes       = create_bounding_boxes(df, distance_km = 5),
			metrics               = "all",
			ee_project            = EE_PROJECT,
			search_buffer         = 1,
			population_years      = YEARS,
			nighttime_light_years = YEARS,
			verbose               = TRUE
		)
	}

	starts <- seq(1, nrow(sites), by = BATCH_SIZE)
	failed <- character(0)

	for (i in seq_along(starts)) {
		f   <- file.path(OUT_DIR, sprintf("study_sites_batch_%03d.rds", i))
		idx <- starts[i]:min(starts[i] + BATCH_SIZE - 1, nrow(sites))
		if (file.exists(f)) { stamp("Batch %d/%d already done - skipping.", i, length(starts)); next }

		stamp("Batch %d/%d: rows %d-%d (%s)", i, length(starts), min(idx), max(idx),
			paste(sites$community[idx], collapse = "; "))
		df <- sites[idx, ]

		res <- tryCatch(run_sites(df), error = function(e) {
			stamp("  batch error: %s -> retrying sites one at a time", conditionMessage(e)); NULL
		})

		if (is.null(res)) {
			parts <- lapply(seq_len(nrow(df)), function(j) {
				tryCatch(run_sites(df[j, ]), error = function(e) {
					stamp("  [FAILED] %s: %s", df$community[j], conditionMessage(e))
					failed <<- c(failed, df$community[j]); NULL
				})
			})
			parts <- Filter(Negate(is.null), parts)
			res <- if (length(parts)) bind_rows(parts) else NULL
		}

		if (!is.null(res)) { saveRDS(res, f); stamp("  saved %s", basename(f)) }
	}


####################
##   ASSEMBLE     ##
####################

	batch_files <- list.files(OUT_DIR, pattern = "^study_sites_batch_\\d{3}\\.rds$", full.names = TRUE)
	out <- bind_rows(lapply(batch_files, readRDS))

	missing <- setdiff(sites$community, out$community)
	if (length(missing)) {
		write_csv(filter(sites, community %in% missing), file.path(OUT_DIR, "study_sites_failed.csv"))
		stamp("%d site(s) missing from results - see study_sites_failed.csv. Re-run this script to retry them", length(missing))
		stamp("  (delete the batch .rds that should contain them first, or they will be skipped).")
	}

	saveRDS(out, file.path(OUT_DIR, "study_sites_urbanicity.rds"))
	write_csv(out, file.path(OUT_DIR, "study_sites_urbanicity.csv"))
	stamp("Done: %d / %d sites -> study_sites_urbanicity.rds/.csv in %s", nrow(out), nrow(sites), OUT_DIR)
