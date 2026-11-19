# 02_generate_local_buffers.R
#
# Prepares the clip mask for the self-hosted Overpass backend. Builds the
# dissolved polygon used to clip the OSM "local 1 km" layers (buildings, shops,
# transport, financial, cell towers) BEFORE tag-filtering, so the extract holds
# only features near the reference points instead of the whole planet.
#
# Input:  ref_pts.csv            (1,400 reference points; cols incl. lon, lat)
#         + study-site CSV       (PEcAHN study sites; cols incl. lon, lat)
# Output: local_buffers.geojson  (dissolved point buffers, WGS84; osmium -p)
#
# The original 2026-06 build buffered only ref_pts.csv. Study sites were added
# 2026-11 so a single mirror serves both; reference-point coverage is unchanged.

####################
##    PACKAGES    ##
####################

library(sf)
library(readr)

########################
##    CONFIGURATION   ##
########################

#	Every CSV listed here is buffered into ONE clip mask. The mirror is shared:
#	reference points + study sites live in the same Overpass instance, so both
#	are measured against the same OSM snapshot. Each CSV needs lon and lat
#	columns. Paths can be overridden from the command line:
#	  Rscript 02_generate_local_buffers.R out.geojson in1.csv in2.csv ...
#	Keep OUTPUT in C:/osm/data so the Docker container sees it at /data.
	INPUT_CSVS     <- c(
		"C:/Users/tmbar/OneDrive/Documents/GitHub/pecahn-urbanicity/reference_data/ref_pts.csv",
		"C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn/manuscript_data/PEcAHN_2026_CommunityInfo_All_9July2026.txt"
	)
	OUTPUT_GEOJSON <- "C:/osm/data/local_buffers.geojson"

	cli_args <- commandArgs(trailingOnly = TRUE)
	if (length(cli_args) >= 2) {
		OUTPUT_GEOJSON <- cli_args[1]
		INPUT_CSVS     <- cli_args[-1]
	}

#	Buffer radius in metres. 3 km covers the 1 km analysis box (corner ~1.41 km)
#	with margin for distance distortion in the global equal-area projection.
	BUFFER_M       <- 3000
	EQUAL_AREA_CRS <- "ESRI:54009"   # World Mollweide (metres); buffer unit is metres
	EXPECTED_N     <- NA             # set to a count to get a warning on mismatch


############################
##    LOAD AND VALIDATE   ##
############################

#	readr strips the UTF-8 BOM that would otherwise corrupt the first column name.
#	The study-site master list is tab-delimited with Lat/Long headers; normalise.
	pts <- do.call(rbind, lapply(INPUT_CSVS, function(f) {
		d <- if (grepl("\\.txt$", f)) read_tsv(f, show_col_types = FALSE) else read_csv(f, show_col_types = FALSE)
		names(d)[tolower(names(d)) == "lat"]  <- "lat"
		names(d)[tolower(names(d)) %in% c("long", "lon")] <- "lon"
		d <- d[!(is.na(d$lat) & is.na(d$lon)), ]   # drop blank trailing rows
		stopifnot("lon and lat columns required" = all(c("lon", "lat") %in% names(d)))
		cat(sprintf("  %4d points from %s\n", nrow(d), basename(f)))
		data.frame(lon = as.numeric(d$lon), lat = as.numeric(d$lat))
	}))

#	Report the count; warn (do not stop) if it differs from the expected draw.
	if (!is.na(EXPECTED_N) && nrow(pts) != EXPECTED_N) {
		message(sprintf("Note: %d points loaded (expected %d). Continuing.",
			nrow(pts), EXPECTED_N))
	}

#	Reject missing or out-of-range coordinates before buffering.
	bad_coords <- is.na(pts$lon) | is.na(pts$lat) |
		pts$lon < -180 | pts$lon > 180 | pts$lat < -90 | pts$lat > 90
	if (any(bad_coords)) {
		stop(sprintf("%d row(s) have missing or out-of-range coordinates.",
			sum(bad_coords)))
	}

	cat(sprintf("Loaded %d points. lon [%.2f, %.2f], lat [%.2f, %.2f]\n",
		nrow(pts), min(pts$lon), max(pts$lon), min(pts$lat), max(pts$lat)))


##################################
##    BUILD DISSOLVED BUFFERS   ##
##################################
	
#	Project to metres, buffer, dissolve overlaps, return to lon/lat for osmium.
	pts_sf <- st_as_sf(pts, coords = c("lon", "lat"), crs = 4326)

	local_buf <- pts_sf |>
		st_transform(EQUAL_AREA_CRS) |>
		st_buffer(BUFFER_M) |>
		st_union() |>
		st_transform(4326)


#################################
##    WRITE AND SANITY CHECK   ##
#################################
	
	st_write(local_buf, OUTPUT_GEOJSON, delete_dsn = TRUE, quiet = TRUE)

#	Count disjoint buffer clusters and total footprint as a quick gut-check.
	n_parts  <- length(st_geometry(st_cast(local_buf, "POLYGON")))
	area_km2 <- as.numeric(st_area(st_transform(local_buf, EQUAL_AREA_CRS))) / 1e6

	cat(sprintf("Wrote %s\n", OUTPUT_GEOJSON))
	cat(sprintf("Buffered %d points at %d m -> %d dissolved polygon(s), total %.0f km^2\n",
		nrow(pts), BUFFER_M, n_parts, area_km2))
	cat("Buffer generation complete.\n")
