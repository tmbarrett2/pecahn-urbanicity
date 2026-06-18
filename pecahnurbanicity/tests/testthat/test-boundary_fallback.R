# Tests for .boundary_distance_estimate() and its boundary-estimate flag wiring.
#
# The bug these guard: for very remote communities the community center can fall
# outside the (masked, non-global) friction surface. The fallback clamps the four
# edge probes onto the raster but must also clamp a *copy* of the origin, otherwise
# gdistance::costDistance() has no valid origin cell, every probe is non-finite, and
# the estimate collapses to NA -- silently dropping the point. These tests use a
# tiny synthetic friction raster (no Earth Engine / OSM) built exactly as
# compute_urbanicity() builds the real one. The friction fixture make_geo_transition()
# lives in helper-friction.R.

test_that("an off-raster origin is clamped so the estimate stays finite (the bug)", {
  g <- make_geo_transition()
  # Community center well outside the [0,2] x [0,2] raster extent.
  center_out <- sf::st_sfc(sf::st_point(c(5, 5)), crs = 4326)
  # Final search bbox much larger than the raster; its edges get clamped in too.
  bbox <- c(left = -10, bottom = -10, right = 10, top = 10)

  est <- .boundary_distance_estimate(bbox, center_out, g$tr, g$crs)

  expect_true(is.finite(est))   # before the fix this was NA
  expect_gt(est, 0)
})

test_that("the recovered estimate is flagged is_boundary_est = TRUE by the caller", {
  g <- make_geo_transition()
  # Force the 'no facility ever found' path so the caller hits the fallback.
  testthat::local_mocked_bindings(
    get_osm_data_cached = function(bbox, key, value = NULL, ...) list()
  )
  center_out <- sf::st_sfc(sf::st_point(c(5, 5)), crs = 4326)  # off the raster
  bbox <- c(left = -10, bottom = -10, right = 10, top = 10)

  res <- suppressMessages(search_facilities_progressive(
    center_point = center_out, bbox = bbox, facility_type = "paved_road",
    tr_corrected = g$tr, raster_crs = g$crs, max_search_multiplier = 0
  ))

  expect_true(is.finite(res$value))
  expect_true(res$is_boundary_est)
})

test_that("an in-bounds center is unaffected by the origin clamp (regression)", {
  g <- make_geo_transition()
  center_in <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)  # well inside [0,2]^2
  bbox <- c(left = -10, bottom = -10, right = 10, top = 10)

  est <- .boundary_distance_estimate(bbox, center_in, g$tr, g$crs)

  # For an in-bounds center the origin clamp is a no-op, so the function must
  # return exactly the pre-fix value: the mean of the four clamped edge probes
  # computed from the *unclamped* center.
  ext   <- raster::extent(g$tr)
  pad_x <- (ext@xmax - ext@xmin) * 1e-6
  pad_y <- (ext@ymax - ext@ymin) * 1e-6
  left   <- max(-10, ext@xmin + pad_x)
  right  <- min( 10, ext@xmax - pad_x)
  bottom <- max(-10, ext@ymin + pad_y)
  top    <- min( 10, ext@ymax - pad_y)
  clon <- (left + right) / 2
  clat <- (top + bottom) / 2
  edges <- list(c(clon, top), c(right, clat), c(clon, bottom), c(left, clat))
  manual <- vapply(edges, function(xy) {
    pt <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(xy), crs = 4326))
    calculate_travel_time(center_in, pt, g$tr, g$crs)
  }, numeric(1))

  expect_equal(est, mean(manual[is.finite(manual)]))
})

test_that("a missing transition layer returns NA_real_ (unchanged)", {
  center <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)
  bbox <- c(left = -10, bottom = -10, right = 10, top = 10)
  expect_identical(.boundary_distance_estimate(bbox, center, NULL, NULL), NA_real_)
})

test_that("a non-geographic extent is not clamped and returns NA_real_ (unchanged)", {
  # Projected raster (UTM, units = m): the extent is outside lon/lat bounds, so
  # looks_geographic is FALSE and neither edges nor origin are clamped. The center
  # and edge midpoints transform off this grid, so every probe is non-finite ->
  # NA. This also guards against the origin clamp leaking into the projected path.
  r <- raster::raster(nrows = 20, ncols = 20,
                      xmn = 500000, xmx = 600000, ymn = 4000000, ymx = 4100000,
                      crs = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs")
  raster::values(r) <- 0.01
  tr <- gdistance::transition(r, function(x) 1 / mean(x), directions = 8)
  tr <- gdistance::geoCorrection(tr, type = "c")
  center <- sf::st_sfc(sf::st_point(c(5, 5)), crs = 4326)
  bbox <- c(left = -10, bottom = -10, right = 10, top = 10)

  # costDistance warns "some coordinates not found and omitted" for each off-grid
  # probe -- that is the expected mechanism here, so suppress it; we assert the NA.
  expect_identical(
    suppressWarnings(.boundary_distance_estimate(bbox, center, tr, raster::crs(r))),
    NA_real_
  )
})

test_that("an origin on a masked in-extent friction cell is snapped so the estimate stays finite", {
  # Friction raster with one interior cell masked (NA) -- a nodata pixel INSIDE
  # the extent, like a water/desert mask edge in the real MAP surface. Built as
  # compute_urbanicity() builds the real transition layer.
  r <- raster::raster(nrows = 20, ncols = 20,
                      xmn = 0, xmx = 2, ymn = 0, ymx = 2, crs = "EPSG:4326")
  raster::values(r) <- 0.01
  masked_cell <- raster::cellFromXY(r, cbind(1, 1))
  r[masked_cell] <- NA
  tr <- gdistance::transition(r, function(x) 1 / mean(x), directions = 8)
  tr <- gdistance::geoCorrection(tr, type = "c")

  # Center sits on the masked cell -> pre-snap cost-distance has no valid origin
  # node and the estimate would collapse to NA. The snap moves it to the nearest
  # reachable cell, recovering a finite (conservative) boundary estimate.
  center_masked <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)
  bbox <- c(left = -10, bottom = -10, right = 10, top = 10)

  est <- suppressWarnings(.boundary_distance_estimate(bbox, center_masked, tr, raster::crs(r)))
  expect_true(is.finite(est))
  expect_gt(est, 0)
})

test_that(".nearest_valid_origin is a no-op for an origin already on a valid cell", {
  g <- make_geo_transition()
  pt <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326))  # valid cell
  snapped <- .nearest_valid_origin(pt, g$tr)
  expect_equal(sf::st_coordinates(snapped)[1, 1:2],
               sf::st_coordinates(pt)[1, 1:2])
})
