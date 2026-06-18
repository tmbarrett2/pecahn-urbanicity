# 02_generate_local_buffers.R
#
# Prepares the clip mask for the self-hosted Overpass backend. Builds the
# dissolved polygon used to clip the OSM "local 1 km" layers (buildings, shops,
# transport, financial, cell towers) BEFORE tag-filtering, so the extract holds
# only features near the reference points instead of the whole planet.
#
# Input:  ref_pts.csv            (1,400 reference points; cols incl. lon, lat)
# Output: local_buffers.geojson  (dissolved point buffers, WGS84; osmium -p)

####################
##    PACKAGES    ##
####################

library(sf)
library(readr)

########################
##    CONFIGURATION   ##
########################

#	Point INPUT_CSV at the file's real location; keep OUTPUT in C:/osm/data
#	so the Docker container can reach it at /data during the osmium step.
	INPUT_CSV      <- "C:/osm/data/ref_pts.csv"
	OUTPUT_GEOJSON <- "C:/osm/data/local_buffers.geojson"

#	Buffer radius in metres. 3 km covers the 1 km analysis box (corner ~1.41 km)
#	with margin for distance distortion in the global equal-area projection.
	BUFFER_M       <- 3000
	EQUAL_AREA_CRS <- "ESRI:54009"   # World Mollweide (metres); buffer unit is metres
	EXPECTED_N     <- 1400


############################
##    LOAD AND VALIDATE   ##
############################
	
#	readr strips the UTF-8 BOM that would otherwise corrupt the first column name.
	pts <- read_csv(INPUT_CSV, show_col_types = FALSE)

	stopifnot("lon and lat columns required" = all(c("lon", "lat") %in% names(pts)))

#	Report the count; warn (do not stop) if it differs from the expected draw.
	if (nrow(pts) != EXPECTED_N) {
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
