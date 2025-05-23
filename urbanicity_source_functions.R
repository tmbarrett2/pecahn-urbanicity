# Urbanicity Index for Anthropological Fieldsites - Source Functions
# Developed as Part of the Population Ecology, Aging, and Health Network (PEcAHN)

# Create Boundary for Each Community
# Provide specific boundary coordinates (custom_bbox) or use a pre-specified distance around a central point (distance_km).
# Function defaults to the custom coordinates if provided.
  create_bounding_boxes <- function(gps_df, distance_km = 5, custom_bbox = NULL) {
    
    # Initialize result list
      result <- list()
    
    # Process each point
      for (i in seq_len(nrow(gps_df))) {
        
        # Ensure lat/lon are numeric
          lat <- as.numeric(gps_df$lat[i])
          lon <- as.numeric(gps_df$lon[i])
        
        # Determine the bounding box
          if (!is.null(custom_bbox) && 
              all(c("nw_lat", "nw_lon", "se_lat", "se_lon") %in% names(custom_bbox))) {
            
            # Use the custom bounding box if valid coordinates are provided
              bbox <- c(
                left = custom_bbox$nw_lon,
                bottom = custom_bbox$se_lat,
                right = custom_bbox$se_lon,
                top = custom_bbox$nw_lat
              )
          } else if (!is.null(custom_bbox) && 
                     all(c("nw_lat", "nw_lon", "ne_lat", "ne_lon", "sw_lat", "sw_lon", "se_lat", "se_lon") %in% names(custom_bbox))) {
            
            # Use all four corners if provided
              bbox <- c(
                left = min(custom_bbox$nw_lon, custom_bbox$sw_lon),
                bottom = min(custom_bbox$sw_lat, custom_bbox$se_lat),
                right = max(custom_bbox$ne_lon, custom_bbox$se_lon),
                top = max(custom_bbox$nw_lat, custom_bbox$ne_lat)
              )
          } else {
            
            # Calculate offsets based on the point's latitude
              lat_offset_km <- distance_km / 110.57
            
            # Longitude degrees vary with latitude
              lon_km_per_degree <- 111.32 * cos(lat * pi/180)
              lon_offset_km <- distance_km / lon_km_per_degree
            
            # Create standard bounding box based on distance
              bbox <- c(
                left = lon - lon_offset_km, 
                bottom = lat - lat_offset_km,
                right = lon + lon_offset_km, 
                top = lat + lat_offset_km
              )
          }
        
        # Name the bbox coordinates
          names(bbox) <- c("left", "bottom", "right", "top")
        
        # Create bounding box with original metadata
          result[[i]] <- list(
            bbox = bbox,
            project = as.character(gps_df$project[i]),
            ethnicity = as.character(gps_df$ethnicity[i]),
            community = as.character(gps_df$community[i])
        )
      }
    
    # Add names to the list for easier reference
      names(result) <- paste0("community_", seq_along(result))
    
    return(result)
  }
  
# Helper function to calculate travel time using friction surface
  calculate_travel_time <- function(from_point, to_points, friction_raster, max_search_km = 50) {
    if (is.null(to_points) || nrow(to_points) == 0) {
      return(NA)
    }
    
    # Convert points to same CRS as raster
      from_sp <- sf::st_transform(from_point, crs = raster::crs(friction_raster))
      to_sp <- sf::st_transform(to_points, crs = raster::crs(friction_raster))
    
    # Extract coordinates
      from_coords <- sf::st_coordinates(from_sp)
      to_coords <- sf::st_coordinates(to_sp)
    
    # Create transition layer (inverse of friction for speed)
      tr <- gdistance::transition(friction_raster, 
                                  transitionFunction = function(x) 1/mean(x), 
                                  directions = 8)
      tr_corrected <- gdistance::geoCorrection(tr, type = "c")
    
    # Calculate cost distance from the origin point
      cost_dist <- gdistance::costDistance(tr_corrected, 
                                           from_coords[1,], 
                                           to_coords)
    
    # Return minimum travel time in minutes
      min_time <- min(cost_dist, na.rm = TRUE)
      return(ifelse(is.finite(min_time), min_time, NA))
  }
  
# Compute Urbanicity Metric For Single community
  compute_urbanicity <- function(community_data, name = NULL,
                                 roads = TRUE,
                                 shops = TRUE,
                                 healthcare = TRUE,
                                 transport = TRUE,
                                 financial = TRUE,
                                 schools = TRUE,
                                 cell_towers = TRUE,
                                 buildings = TRUE,
                                 friction_surface = "walking",
                                 friction_surface_path = NULL) {
    # Load required libraries
      required_packages <- c("osmdata", "sf", "raster", "gdistance", "httr", "terra")
      for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(paste("Package", pkg, "needed for this function to work."))
        }
      }
    
    # Import functions
      library(osmdata)
      library(sf)
      library(raster)
      library(gdistance)
      library(httr)
      library(terra)
    
    # Extract the bbox from the community_data
      bbox <- community_data$bbox
    
    # If name is not provided, construct it from metadata
      if (is.null(name)) {
        name <- paste(community_data$project, community_data$ethnicity, community_data$community, sep = " - ")
      }
    
    # Print processing status
      cat("Processing community:", name, "\n")
    
    # Create sf polygon directly from the bbox coordinates
      bbox_matrix <- matrix(
        c(bbox["left"], bbox["bottom"],
          bbox["left"], bbox["top"],
          bbox["right"], bbox["top"],
          bbox["right"], bbox["bottom"],
          bbox["left"], bbox["bottom"]),
        ncol = 2, byrow = TRUE
      )
    
      poly <- sf::st_polygon(list(bbox_matrix)) %>% sf::st_sfc(crs = 4326)
    
    # Initialize results with base information
      results <- list(
        project = community_data$project,
        ethnicity = community_data$ethnicity,
        community = community_data$community
      )
    
    # Define bounding box for osmdata in the correct format
      osm_bbox <- c(bbox["left"], bbox["bottom"], bbox["right"], bbox["top"])
      
    # Calculate center point from bbox
      center_lon <- (bbox["left"] + bbox["right"]) / 2
      center_lat <- (bbox["bottom"] + bbox["top"]) / 2
      
    # Create center point sf object
      center_point <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
      
    # Load friction surface (with local file support)
      friction_raster <- NULL
      if (healthcare || schools || roads) {
        tryCatch({
          cat("  Loading friction surface data...\n")
          
          # Check if a local friction surface path is provided
          if (!is.null(friction_surface_path) && file.exists(friction_surface_path)) {
            cat("  Using local friction surface file:", friction_surface_path, "\n")
            friction_global <- raster::raster(friction_surface_path)
            
            # Crop to a larger area around the community (for searching)
            search_buffer <- 0.5  # degrees (~50km)
            extent_buffered <- raster::extent(
              bbox["left"] - search_buffer,
              bbox["right"] + search_buffer,
              bbox["bottom"] - search_buffer,
              bbox["top"] + search_buffer
            )
            
            friction_raster <- raster::crop(friction_global, extent_buffered)
            cat("  Friction surface loaded successfully\n")
            
          } else {
            # Fallback to download approach (keeping for backwards compatibility)
            cat("  No local friction surface provided, attempting download...\n")
            
            # Define URLs for MAP friction surfaces
            friction_urls <- list(
              walking = "https://malariaatlas.org/geoserver/Accessibility/ows?service=WCS&version=2.0.1&request=GetCoverage&format=image/geotiff&coverageid=Accessibility:202001_Global_Walking_Only_Travel_Time_To_Healthcare",
              motorized = "https://malariaatlas.org/geoserver/Accessibility/ows?service=WCS&version=2.0.1&request=GetCoverage&format=image/geotiff&coverageid=Accessibility:202001_Global_Motorized_Travel_Time_To_Healthcare"
            )
            
            # Create cache directory
            cache_dir <- file.path(tempdir(), "friction_cache")
            dir.create(cache_dir, showWarnings = FALSE)
            cache_file <- file.path(cache_dir, paste0("friction_", friction_surface, ".tif"))
            
            # Download or load from cache
            if (!file.exists(cache_file)) {
              cat("  Downloading friction surface (this may take a while on first run)...\n")
              download.file(friction_urls[[friction_surface]], cache_file, mode = "wb", quiet = TRUE)
            }
            
            # Load the global friction raster
            friction_global <- raster::raster(cache_file)
            
            # Crop to study area
            search_buffer <- 0.5
            extent_buffered <- raster::extent(
              bbox["left"] - search_buffer,
              bbox["right"] + search_buffer,
              bbox["bottom"] - search_buffer,
              bbox["top"] + search_buffer
            )
            
            friction_raster <- raster::crop(friction_global, extent_buffered)
          }
          
        }, error = function(e) {
          cat("  Error loading friction surface. Using fallback straight-line distance method.\n")
          cat("  Error details:", e$message, "\n")
          friction_raster <- NULL
        })
      }
      
    # Create center point sf object
      center_point <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
      center_lat <- as.numeric(community_data$lat)
      center_lon <- as.numeric(community_data$lon)
    
    # Ratio of Paved to Unpaved Roads and Percent Paved Roads
      if (roads) {
        tryCatch({
          cat("  Analyzing roads...\n")
          roads_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "highway") %>%
            osmdata::osmdata_sf()
          
          # Define surface types
          paved_surfaces <- c("paved", "asphalt", "concrete")
          unpaved_surfaces <- c("unpaved", "dirt", "gravel", "fine_gravel", "sand", "grass", "ground", "earth", "mud", "compacted", "unclassified")
          
          if (!is.null(roads_data$osm_lines) && nrow(roads_data$osm_lines) > 0) {
            roads_sf <- roads_data$osm_lines
            
            # Check if surface column exists
            if (!"surface" %in% colnames(roads_sf)) {
              cat("  Surface data does not exist in the OSM data. Using highway types for classification.\n")
              roads_sf$surface <- NA
            }
            
            # Classify roads based on surface tag
            roads_sf$road_type <- "other"
            roads_sf$road_type[roads_sf$surface %in% paved_surfaces] <- "paved"
            roads_sf$road_type[roads_sf$surface %in% unpaved_surfaces] <- "unpaved"
            
            # Classify roads based on highway tag for tracks and paths
            # This helps when surface information is missing but highway type indicates unpaved
            if ("highway" %in% colnames(roads_sf)) {
              unpaved_highway_types <- c("track", "path")
              roads_sf$road_type[roads_sf$highway %in% unpaved_highway_types] <- "unpaved"
            }
            
            # Calculate road lengths
            roads_sf <- sf::st_transform(roads_sf, crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84"))
            roads_sf$length <- as.numeric(sf::st_length(roads_sf))
            
            # Summarize by road type
            road_summary <- aggregate(as.numeric(roads_sf$length), by = list(road_type = roads_sf$road_type), FUN = sum)
            
            # Calculate paved road ratio
            total_length <- sum(road_summary$x, na.rm = TRUE)
            paved_length <- sum(road_summary$x[road_summary$road_type == "paved"], na.rm = TRUE)
            unpaved_length <- sum(road_summary$x[road_summary$road_type == "unpaved"], na.rm = TRUE)
            
            # Calculate ratio
            if (unpaved_length > 0) {
              results$paved_to_unpaved_ratio <- as.numeric(paved_length / unpaved_length)
            } else {
              results$paved_to_unpaved_ratio <- Inf
            }
            
            results$pct_paved_roads <- 100 * as.numeric(paved_length / total_length)
            
            # Calculate travel time to nearest paved road
            paved_roads <- roads_sf[roads_sf$road_type == "paved", ]
            
            if (nrow(paved_roads) > 0) {
              if (!is.null(friction_raster)) {
                # Convert lines to points for travel time calculation
                suppressWarnings({
                  paved_points <- sf::st_cast(paved_roads, "POINT")
                })
                results$travel_time_paved_road_min <- calculate_travel_time(center_point, paved_points, friction_raster)
              } else {
                # Fallback: straight-line distance
                distances <- sf::st_distance(center_point, paved_roads)
                min_dist_m <- min(distances, na.rm = TRUE)
                # Assume walking speed of 5 km/h
                results$travel_time_paved_road_min <- as.numeric(min_dist_m) / 1000 / 5 * 60
              }
            } else {
              # Search in larger area
              cat("  No paved roads in immediate area, searching wider region...\n")
              search_multiplier <- 2
              extended_bbox <- c(
                bbox["left"] - (bbox["right"] - bbox["left"]) * search_multiplier,
                bbox["bottom"] - (bbox["top"] - bbox["bottom"]) * search_multiplier,
                bbox["right"] + (bbox["right"] - bbox["left"]) * search_multiplier,
                bbox["top"] + (bbox["top"] - bbox["bottom"]) * search_multiplier
              )
              
              roads_extended <- osmdata::opq(extended_bbox) %>%
                osmdata::add_osm_feature(key = "highway") %>%
                osmdata::add_osm_feature(key = "surface", value = paved_surfaces) %>%
                osmdata::osmdata_sf()
              
              if (!is.null(roads_extended$osm_lines) && nrow(roads_extended$osm_lines) > 0) {
                if (!is.null(friction_raster)) {
                  suppressWarnings({
                    paved_points_ext <- sf::st_cast(roads_extended$osm_lines, "POINT")
                  })
                  paved_points_ext <- sf::st_cast(roads_extended$osm_lines, "POINT")
                  results$travel_time_paved_road_min <- calculate_travel_time(center_point, paved_points_ext, friction_raster)
                } else {
                  distances <- sf::st_distance(center_point, roads_extended$osm_lines)
                  min_dist_m <- min(distances, na.rm = TRUE)
                  results$travel_time_paved_road_min <- as.numeric(min_dist_m) / 1000 / 5 * 60
                }
              } else {
                results$travel_time_paved_road_min <- NA
              }
            }
            
          } else {
            results$paved_to_unpaved_ratio <- NA
            results$pct_paved_roads <- NA
            results$travel_time_paved_road_min <- NA
            cat("  No roads found in the area.\n")
          }
        }, error = function(e) {  # FIXED: Added comma before 'error'
          cat("  Error processing roads:", e$message, "\n")
          results$paved_to_unpaved_ratio <- NA
          results$pct_paved_roads <- NA
          results$travel_time_paved_road_min <- NA
        })  # ADDED: Closing parenthesis for tryCatch
      }  # ADDED: Closing brace for if (roads)
               
    
    # Number of Formal Shops
      if (shops) {
        tryCatch({
          cat("  Analyzing shops and markets...\n")
          shops_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "shop") %>%
            osmdata::osmdata_sf()
          
          if (!is.null(shops_data$osm_points) && nrow(shops_data$osm_points) > 0) {
            
            # Count all shops
              total_shops_count <- nrow(shops_data$osm_points)
              results$n_shops <- total_shops_count
            } else {
              results$n_shops <- 0
            }
          }, error = function(e) {
            cat("  Error processing shops:", e$message, "\n")
            results$n_shops <- NA
         })
      }
    
    # Number of Healthcare Facilities and Distance to Nearest Healthcare Facility
    # NOTE: jUST yes / no for now, was getting oddly large values for some communities
            # Travel time to nearest Healthcare Facility
            if (healthcare) {
              tryCatch({
                cat("  Analyzing healthcare facilities...\n")
                healthcare_data <- osmdata::opq(osm_bbox) %>%
                  osmdata::add_osm_feature(key = "amenity", value = c("hospital", "clinic", "doctors")) %>%
                  osmdata::osmdata_sf()
                
                healthcare_points <- NULL
                if (!is.null(healthcare_data$osm_points) && nrow(healthcare_data$osm_points) > 0) {
                  healthcare_points <- healthcare_data$osm_points
                }
                
                if (!is.null(healthcare_points) && nrow(healthcare_points) > 0) {
                  if (!is.null(friction_raster)) {
                    results$travel_time_healthcare_min <- calculate_travel_time(center_point, healthcare_points, friction_raster)
                  } else {
                    # Fallback: straight-line distance
                    distances <- sf::st_distance(center_point, healthcare_points)
                    min_dist_m <- min(distances, na.rm = TRUE)
                    results$travel_time_healthcare_min <- as.numeric(min_dist_m) / 1000 / 5 * 60
                  }
                } else {
                  # Search in larger area
                  cat("  No healthcare facilities in immediate area, searching wider region...\n")
                  search_multiplier <- 2
                  extended_bbox <- c(
                    bbox["left"] - (bbox["right"] - bbox["left"]) * search_multiplier,
                    bbox["bottom"] - (bbox["top"] - bbox["bottom"]) * search_multiplier,
                    bbox["right"] + (bbox["right"] - bbox["left"]) * search_multiplier,
                    bbox["top"] + (bbox["top"] - bbox["bottom"]) * search_multiplier
                  )
                  
                  healthcare_extended <- osmdata::opq(extended_bbox) %>%
                    osmdata::add_osm_feature(key = "amenity", value = c("hospital", "clinic", "doctors")) %>%
                    osmdata::osmdata_sf()
                  
                  if (!is.null(healthcare_extended$osm_points) && nrow(healthcare_extended$osm_points) > 0) {
                    if (!is.null(friction_raster)) {
                      results$travel_time_healthcare_min <- calculate_travel_time(center_point, healthcare_extended$osm_points, friction_raster)
                    } else {
                      distances <- sf::st_distance(center_point, healthcare_extended$osm_points)
                      min_dist_m <- min(distances, na.rm = TRUE)
                      results$travel_time_healthcare_min <- as.numeric(min_dist_m) / 1000 / 5 * 60
                    }
                  } else {
                    results$travel_time_healthcare_min <- NA
                  }
                }
                
              }, error = function(e) {
                cat("  Error processing healthcare facilities:", e$message, "\n")
                results$travel_time_healthcare_min <- NA
              })
            }
    
    # Number of Financial Services
      if (financial) {
        tryCatch({
          cat("  Analyzing financial services...\n")
          financial_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "amenity", value = c("bank", "atm")) %>%
            osmdata::osmdata_sf()
          
          if (!is.null(financial_data$osm_points) && nrow(financial_data$osm_points) > 0) {
            results$n_financial <- nrow(financial_data$osm_points)
          } else {
            results$n_financial <- 0
          }
        }, error = function(e) {
          cat("  Error processing financial services:", e$message, "\n")
          results$n_financial <- NA
        })
      }
    
    # Number of Schools -- NOTE: jUST yes / no for now, was getting oddly large values for some communities
      # Travel time to nearest School
      if (schools) {
        tryCatch({
          cat("  Analyzing schools...\n")
          
          school_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "amenity", value = "school") %>%
            osmdata::osmdata_sf()
          
          school_points <- NULL
          
          # Collect points
          if (!is.null(school_data$osm_points) && nrow(school_data$osm_points) > 0) {
            # Keep only geometry for consistency
            school_points <- sf::st_geometry(school_data$osm_points) %>% 
              sf::st_sf(geometry = .)
          }
          
          # Collect polygon centroids
          if (!is.null(school_data$osm_polygons) && nrow(school_data$osm_polygons) > 0) {
            school_centroids <- sf::st_centroid(school_data$osm_polygons) %>%
              sf::st_geometry() %>%
              sf::st_sf(geometry = .)
            
            if (is.null(school_points)) {
              school_points <- school_centroids
            } else {
              # Combine geometries only
              school_points <- rbind(school_points, school_centroids)
            }
          }
          
          if (!is.null(school_points) && nrow(school_points) > 0) {
            if (!is.null(friction_raster)) {
              results$travel_time_school_min <- calculate_travel_time(center_point, school_points, friction_raster)
            } else {
              # Fallback: straight-line distance
              distances <- sf::st_distance(center_point, school_points)
              min_dist_m <- min(distances, na.rm = TRUE)
              results$travel_time_school_min <- as.numeric(min_dist_m) / 1000 / 5 * 60
            }
          } else {
            # Search in larger area
            cat("  No schools in immediate area, searching wider region...\n")
            search_multiplier <- 2
            extended_bbox <- c(
              bbox["left"] - (bbox["right"] - bbox["left"]) * search_multiplier,
              bbox["bottom"] - (bbox["top"] - bbox["bottom"]) * search_multiplier,
              bbox["right"] + (bbox["right"] - bbox["left"]) * search_multiplier,
              bbox["top"] + (bbox["top"] - bbox["bottom"]) * search_multiplier
            )
            
            school_extended <- osmdata::opq(extended_bbox) %>%
              osmdata::add_osm_feature(key = "amenity", value = "school") %>%
              osmdata::osmdata_sf()
            
            extended_points <- NULL
            
            # Collect extended points
            if (!is.null(school_extended$osm_points) && nrow(school_extended$osm_points) > 0) {
              extended_points <- sf::st_geometry(school_extended$osm_points) %>%
                sf::st_sf(geometry = .)
            }
            
            # Collect extended polygon centroids
            if (!is.null(school_extended$osm_polygons) && nrow(school_extended$osm_polygons) > 0) {
              extended_centroids <- sf::st_centroid(school_extended$osm_polygons) %>%
                sf::st_geometry() %>%
                sf::st_sf(geometry = .)
              
              if (is.null(extended_points)) {
                extended_points <- extended_centroids
              } else {
                extended_points <- rbind(extended_points, extended_centroids)
              }
            }
            
            if (!is.null(extended_points) && nrow(extended_points) > 0) {
              if (!is.null(friction_raster)) {
                results$travel_time_school_min <- calculate_travel_time(center_point, extended_points, friction_raster)
              } else {
                distances <- sf::st_distance(center_point, extended_points)
                min_dist_m <- min(distances, na.rm = TRUE)
                results$travel_time_school_min <- as.numeric(min_dist_m) / 1000 / 5 * 60
              }
            } else {
              results$travel_time_school_min <- NA
            }
          }
          
        }, error = function(e) {
          cat("  Error processing schools:", e$message, "\n")
          results$travel_time_school_min <- NA
        })
      }
    
    # Number of Cell Towers
      if (cell_towers) {
        tryCatch({
          cat("  Analyzing cell towers...\n")
          cell_towers_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "communication", value = "mobile_phone") %>%
            osmdata::osmdata_sf()
          
          if (!is.null(cell_towers_data$osm_points) && nrow(cell_towers_data$osm_points) > 0) {
            results$n_cell_towers <- nrow(cell_towers_data$osm_points)
          } else {
            results$n_cell_towers <- 0
          }
        }, error = function(e) {
          cat("  Error processing cell towers:", e$message, "\n")
          results$n_cell_towers <- NA
        })
      }
    
    # Building Density (% of area covered by buildings)
      if (buildings) {
        tryCatch({
          cat("  Analyzing building footprints...\n")
          buildings_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "building") %>%
            osmdata::osmdata_sf()
          
          if (!is.null(buildings_data$osm_polygons) && nrow(buildings_data$osm_polygons) > 0) {
            # Transform to a suitable projection for accurate area calculation
              buildings_projected <- sf::st_transform(
                buildings_data$osm_polygons, 
                crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84")
              )
            
            # Calculate area of each building
              buildings_projected$area <- sf::st_area(buildings_projected)
            
            # Calculate total building area
              total_building_area <- sum(as.numeric(buildings_projected$area))
            
            # Calculate bounding box area
              poly_projected <- sf::st_transform(
                poly, 
                crs = sf::st_crs("+proj=utm +zone=32 +datum=WGS84")
              )
              bbox_area <- as.numeric(sf::st_area(poly_projected))
            
            # Calculate building density (% of area covered by buildings)
              building_density <- (total_building_area / bbox_area) * 100
              
              results$building_density_pct <- as.numeric(building_density)
            } else {
              results$building_density_pct <- 0
            }
          }, error = function(e) {
            cat("  Error processing buildings:", e$message, "\n")
            results$building_density_pct <- NA
          })
      }
    
    cat("Processing complete for:", name, "\n\n")
    
    # Convert results list to a dataframe
      return(as.data.frame(lapply(results, function(x) if(length(x) == 0) NA else x)))
  }
 
# Apply compute_urbanicity Function to a List of Communities
# communities_list is the output from the bounding boxes function above.
# metrics is a vector of the measures you want to compute.
  compute_urbanicity_iterative <- function(communities_list, 
                                           metrics = c("all"),
                                           friction_surface = "walking",
                                           friction_surface_path = NULL) {  # ADD THIS PARAMETER
    
    # Set up which metrics to analyze based on user input
    do_roads <- "all" %in% metrics || "roads" %in% metrics
    do_shops <- "all" %in% metrics || "shops" %in% metrics
    do_healthcare <- "all" %in% metrics || "healthcare" %in% metrics
    do_transport <- "all" %in% metrics || "transport" %in% metrics
    do_financial <- "all" %in% metrics || "financial" %in% metrics
    do_schools <- "all" %in% metrics || "schools" %in% metrics
    do_cell_towers <- "all" %in% metrics || "cell_towers" %in% metrics
    do_buildings <- "all" %in% metrics || "buildings" %in% metrics
    
    # Initialize an empty list to store results
    all_results <- list()
    
    # Process each community
    for (i in seq_along(communities_list)) {
      community_name <- names(communities_list)[i]
      community_data <- communities_list[[i]]
      
      # Process the community with selected metrics
      result <- compute_urbanicity(
        community_data, 
        name = community_name,
        roads = do_roads,
        shops = do_shops,
        healthcare = do_healthcare,
        transport = do_transport,
        financial = do_financial,
        schools = do_schools,
        cell_towers = do_cell_towers,
        buildings = do_buildings,
        friction_surface = friction_surface,
        friction_surface_path = friction_surface_path
      )
      
      # Add to results list
      all_results[[i]] <- result
    }
    
    # Combine all results into a single data frame
    result_names <- unique(unlist(lapply(all_results, names)))
    combined_results <- data.frame(matrix(NA, nrow = length(all_results), ncol = length(result_names)))
    names(combined_results) <- result_names
    
    for (i in seq_along(all_results)) {
      for (col in names(all_results[[i]])) {
        combined_results[i, col] <- all_results[[i]][1, col]
      }
    }
    
    return(combined_results)
  }
  
# Plot Results
  # Bar Plot For Numeric Variables
    create_numeric_plot <- function(data, var_name) {
    
    # Load packages
      library(ggplot2)
      library(dplyr)
    
    # Handle NA values
      data_filtered <- data %>% 
        filter(!is.na(!!sym(var_name)))
    
    # Determine format based on variable name
    # Use whole numbers for variables starting with "n_", three decimals for others
      format_string <- ifelse(startsWith(var_name, "n_"), "%.0f", "%.3f")
    
    # Create plot
      p <- data_filtered %>%
        mutate(community = paste(project, community, sep=": ")) %>%
        arrange(desc(!!sym(var_name))) %>%
        mutate(community = factor(community, levels = community)) %>%
        ggplot(aes(x = community, y = !!sym(var_name), fill = project)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d() +
        geom_text(aes(label = sprintf(format_string, !!sym(var_name))), 
                  hjust = -0.2) +
        labs(title = paste("Summary of", var_name),
             x = "Community",
             y = var_name) +
        theme_minimal() +
        theme(legend.position = "top") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        coord_flip()
    
    return(p)
  }
  
  # Iterate Plotting Function Over All Measures (Continuous Measures Only)
    create_summary_plots <- function(data) {
    
    # Load packages
      library(ggplot2)
      library(dplyr)
    
    # Get all column names except ID columns
      all_vars <- names(data)
      all_vars <- all_vars[!all_vars %in% c("project", "community")]
    
    # Create a list to store all plots
      all_plots <- list()
    
    # Create plots for each variable
      for (var in all_vars) {
        
        # Check if variable is numeric (continuous)
        if (is.numeric(data[[var]])) {
          all_plots[[var]] <- create_numeric_plot(data, var)
        }
      }
    
    # Return the list of plots
      return(all_plots)
  }
  
  
  