# 05_assemble_reference_dataset.R
#
# Assembles the urbanicity reference dataset: binds the per-class urbanicity
# batch results for the 1,400 reference points, joins each point's design weight
# and settlement class, and writes one tidy table (one row per reference point).
#
# Inputs:
#   - per-class urbanicity batch results .rds   (from 04_compute_urbanicity_reference.Rmd)
#   - ref_pts.csv                               (design weights w; this folder)
# Output:
#   - reference_data/urbanicity_reference_dataset.rds

####################
##    PACKAGES    ##
####################

library(dplyr)


########################
##    CONFIGURATION   ##
########################

#	Directory holding the per-class batch results from
#	04_compute_urbanicity_reference.Rmd, and the glob that matches them
#	(save_prefix = "<output_fp>subsetN", file = "<prefix>_batch_NNN.rds").
	RESULTS_DIR  <- "C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn"
	RESULTS_GLOB <- "ref_pts_resultssubset[0-9]*_batch_[0-9]*\\.rds$"

#	Canonical reference point table (carries the design weights w). Use the copy
#	shipped in this folder so weights track the committed sample exactly.
	REF_PTS_CSV  <- "ref_pts.csv"

#	Assembled dataset is written next to this script.
	OUTPUT_RDS   <- "urbanicity_reference_dataset.rds"

#	Provenance recorded onto the assembled object (as attributes).
	OSM_SNAPSHOT <- "2026-06-08"
	EXPECTED_N   <- 1400L
	EXPECTED_CLASSES <- c(11L, 12L, 13L, 21L, 22L, 23L, 30L)


###########################
##  LOAD BATCH RESULTS   ##
###########################

#	Glob every per-class batch .rds and bind into one long results table. Each
#	.rds is the bound output of up to five communities from the urbanicity run.
	rds_files <- list.files(RESULTS_DIR, pattern = RESULTS_GLOB, full.names = TRUE)
	if (length(rds_files) == 0) {
		stop(sprintf("No batch results matched '%s' in %s", RESULTS_GLOB, RESULTS_DIR))
	}

	results <- bind_rows(lapply(rds_files, readRDS))

	cat(sprintf("Loaded %d batch file(s) -> %d community rows.\n",
		length(rds_files), nrow(results)))

#	Fold in the targeted re-run of the off-grid / masked points
#	(rerun_problem_points.R). For any community present in the rerun file, drop the
#	stale (NA) original row and keep the recovered one. Skipped if no rerun file.
	RERUN_RDS <- file.path(RESULTS_DIR, "ref_pts_results_rerun.rds")
	if (file.exists(RERUN_RDS)) {
		rerun  <- readRDS(RERUN_RDS)
		n_repl <- sum(results$community %in% rerun$community)
		results <- bind_rows(results[!results$community %in% rerun$community, ], rerun)
		cat(sprintf("Merged rerun: replaced %d stale community row(s) -> %d total.\n",
			n_repl, nrow(results)))
	}

#	The urbanicity results store the point identity in `community` as "lat, lon"
#	and the class label in `project`. Warn (not stop) on a short/duplicated set.
	if (anyDuplicated(results$community)) {
		warning(sum(duplicated(results$community)),
			" duplicate community key(s) found in the batch results.")
	}
	if (nrow(results) != EXPECTED_N) {
		message(sprintf("Note: %d rows loaded (expected %d). Continuing.",
			nrow(results), EXPECTED_N))
	}


########################
##    JOIN WEIGHTS    ##
########################

#	Read ref_pts.csv exactly as 04_compute_urbanicity_reference.Rmd did (base
#	read.csv) and rebuild the same "lat, lon" key so the IEEE-754 doubles format
#	identically on both sides. Do NOT re-round lat/lon — that would break the join.
#	fileEncoding strips the UTF-8 BOM so the first column name stays "Classified"
#	on every platform; lat/lon parse to identical doubles either way, so the
#	"lat, lon" join key matches the urbanicity results regardless.
	ref_pts <- read.csv(REF_PTS_CSV, fileEncoding = "UTF-8-BOM")
	stopifnot(all(c("Classified", "w", "lat", "lon") %in% names(ref_pts)))

	ref_pts$community <- paste0(ref_pts$lat, ", ", ref_pts$lon)
	weights <- ref_pts[, c("community", "Classified", "w")]
	names(weights) <- c("community", "smod_class", "weight")

	dataset <- left_join(results, weights, by = "community")

	n_unmatched <- sum(is.na(dataset$weight))
	if (n_unmatched > 0) {
		stop(sprintf("%d community key(s) did not match ref_pts.csv — check lat/lon precision.",
			n_unmatched))
	}


#################################
##   ASSEMBLE & ORDER TABLE    ##
#################################

#	Move the identity / design columns to the front and keep every computed
#	urbanicity metric so users can summarise or weight the data however they like.
#	w = 0 points (uninhabited class-11 cells) are kept — they carry zero design
#	weight, so any weighted summary gives them no mass, but their metrics remain.
	front_cols <- c("community", "smod_class", "weight", "project")
	dataset <- dataset[, c(front_cols[front_cols %in% names(dataset)],
		setdiff(names(dataset), front_cols))]

	dataset <- dataset[order(dataset$smod_class, dataset$community), ]
	rownames(dataset) <- NULL

#	Record provenance as attributes so it travels with the .rds.
	attr(dataset, "osm_snapshot")    <- OSM_SNAPSHOT
	attr(dataset, "n_points")        <- nrow(dataset)
	attr(dataset, "package_version") <- tryCatch(
		as.character(utils::packageVersion("pecahnurbanicity")), error = function(e) NA_character_)
	attr(dataset, "built_on")        <- as.character(Sys.Date())

	cat(sprintf("Assembled reference dataset: %d points, %d classes.\n",
		nrow(dataset), length(unique(dataset$smod_class))))


#################################
##        WRITE DATASET        ##
#################################

#	Guard: only write when the full stratified set is present (all 1,400 points,
#	all seven classes). On a partial run this self-skips so a class-incomplete
#	dataset is never committed.
	classes_present <- sort(unique(dataset$smod_class))
	is_complete <- nrow(dataset) == EXPECTED_N &&
		all(EXPECTED_CLASSES %in% classes_present)

	if (is_complete) {
		saveRDS(dataset, OUTPUT_RDS)
		cat(sprintf("Wrote %s (%d points).\n", OUTPUT_RDS, nrow(dataset)))
	} else {
		message(sprintf(
			paste0("Reference set incomplete (%d rows; classes present: %s). ",
				"Skipping write. Re-run after all seven classes finish 04_compute_urbanicity_reference.Rmd."),
			nrow(dataset), paste(classes_present, collapse = ", ")))
	}
