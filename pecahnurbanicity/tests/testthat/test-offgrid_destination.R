# Tests for the off-grid-destination handling: calculate_travel_time() must not
# throw when the origin or all destinations fall outside the friction grid, and the
# progressive search functions must fall back to a boundary-distance estimate when
# the *found* feature is unreachable on the friction surface (rather than letting the
# error propagate and the metric collapse to a bare NA).
#
# Fixture make_geo_transition() (EPSG:4326, extent [0,2]x[0,2]) is in helper-friction.R.

test_that("calculate_travel_time returns NA (not an error) when all destinations are off-grid", {
  g <- make_geo_transition()
  origin <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)                              # on-grid
  dest   <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(50, 50)), crs = 4326))      # far off-grid

  expect_no_error(tt <- calculate_travel_time(origin, dest, g$tr, g$crs))
  expect_true(is.na(tt))
})

test_that("calculate_travel_time returns NA (not an error) when the origin is off-grid", {
  g <- make_geo_transition()
  origin <- sf::st_sfc(sf::st_point(c(50, 50)), crs = 4326)                            # off-grid
  dest   <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326))        # on-grid

  expect_no_error(tt <- calculate_travel_time(origin, dest, g$tr, g$crs))
  expect_true(is.na(tt))
})

test_that("calculate_travel_time still measures travel time to an on-grid destination", {
  g <- make_geo_transition()
  origin <- sf::st_sfc(sf::st_point(c(0.5, 0.5)), crs = 4326)
  dest   <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(1.5, 1.5)), crs = 4326))

  tt <- calculate_travel_time(origin, dest, g$tr, g$crs)
  expect_true(is.finite(tt))
  expect_gt(tt, 0)
})

test_that("search_facilities_progressive falls back to a boundary estimate when the found facility is off-grid", {
  g <- make_geo_transition()
  off_grid <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(50, 50)), crs = 4326))
  # Every OSM query returns a hospital point that sits off the friction grid.
  testthat::local_mocked_bindings(
    get_osm_data_cached = function(bbox, key, value = NULL, ...) list(osm_points = off_grid)
  )
  center <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)   # on-grid origin
  bbox <- c(left = 0.5, bottom = 0.5, right = 1.5, top = 1.5)

  res <- suppressMessages(search_facilities_progressive(
    center_point = center, bbox = bbox, facility_type = "hospital",
    tr_corrected = g$tr, raster_crs = g$crs, max_search_multiplier = 0
  ))

  expect_true(is.finite(res$value))     # not a bare NA
  expect_true(res$is_boundary_est)      # flagged as a boundary estimate
})

test_that("search_urban_center_progressive falls back to a boundary estimate when found centers are off-grid", {
  g <- make_geo_transition()
  off_grid <- sf::st_as_sf(data.frame(lon = 50, lat = 50),
                           coords = c("lon", "lat"), crs = 4326)
  testthat::local_mocked_bindings(
    ee_get_ghsl_urban_centers = function(bbox, ghsl_year, codes, ...) off_grid
  )
  center <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)
  bbox <- c(left = 0.5, bottom = 0.5, right = 1.5, top = 1.5)

  res <- suppressMessages(search_urban_center_progressive(
    center_point = center, bbox = bbox,
    tr_corrected = g$tr, raster_crs = g$crs, max_search_multiplier = 0
  ))

  expect_true(is.finite(res$value))
  expect_true(res$is_boundary_est)
})
