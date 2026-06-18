# Shared test fixture: a tiny geographic (EPSG:4326) friction transition layer,
# built exactly as compute_urbanicity() builds the real one (transition() +
# geoCorrection()). Covers [0,2] x [0,2] degrees with uniform friction. No Earth
# Engine / OSM required.
make_geo_transition <- function() {
  r <- raster::raster(nrows = 20, ncols = 20,
                      xmn = 0, xmx = 2, ymn = 0, ymx = 2,
                      crs = "EPSG:4326")
  raster::values(r) <- 0.01  # minutes/meter, uniform
  tr <- gdistance::transition(r, function(x) 1 / mean(x), directions = 8)
  tr <- gdistance::geoCorrection(tr, type = "c")
  list(tr = tr, crs = raster::crs(r))
}
