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
  
# Helper Function to Compute Mean Nighttime Light Value
  get_light_mean <- function(polygon_sf, raster_data) {
      polygon_sp <- as(polygon_sf, "Spatial")
      light_mean <- raster::extract(raster_data, polygon_sp, fun = mean, na.rm = TRUE)
      return(round(light_mean, 3))
    }
  
# Helper Function to Calculate Travel Time Using Friction Surface Raster
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
  
# Helper Function to Extend Search for Very Remote Communities
  search_facilities_progressive <- function(center_point, bbox, facility_type, 
                                            friction_raster = NULL, 
                                            max_search_multiplier = 20,
                                            search_increment = 2) {
    
    search_multiplier <- 0
    found_facilities <- NULL
    
    while (is.null(found_facilities) && search_multiplier <= max_search_multiplier) {
      # Define search bbox
        if (search_multiplier == 0) {
          search_bbox <- c(bbox["left"], bbox["bottom"], bbox["right"], bbox["top"])
        } else {
          width <- bbox["right"] - bbox["left"]
          height <- bbox["top"] - bbox["bottom"]
          search_bbox <- c(
            bbox["left"] - width * search_multiplier,
            bbox["bottom"] - height * search_multiplier,
            bbox["right"] + width * search_multiplier,
            bbox["top"] + height * search_multiplier
          )
        }
      
      # Search for facilities based on type
        if (facility_type == "paved_road") {
         # Search for ALL roads first
          roads_data <- osmdata::opq(search_bbox) %>%
            osmdata::add_osm_feature(key = "highway") %>%
            osmdata::osmdata_sf()
          
          if (!is.null(roads_data$osm_lines) && nrow(roads_data$osm_lines) > 0) {
            roads_sf <- roads_data$osm_lines
            
          # Check if surface column exists
            if (!"surface" %in% colnames(roads_sf)) {
              roads_sf$surface <- NA
            }
            
          # Define surface types
            paved_surfaces <- c("paved", "asphalt", "concrete")
          
          # Filter for paved roads - only those with explicit surface tags
            paved_roads <- roads_sf[roads_sf$surface %in% paved_surfaces, ]
          
          if (nrow(paved_roads) > 0) {
            found_facilities <- paved_roads
          }
          # If no paved roads found, found_facilities remains NULL and search will continue
              }
            } else if (facility_type == "healthcare") {
              healthcare_data <- osmdata::opq(search_bbox) %>%
                osmdata::add_osm_feature(key = "amenity", value = c("hospital", "clinic", "doctors")) %>%
                osmdata::osmdata_sf()
              
              if (!is.null(healthcare_data$osm_points) && nrow(healthcare_data$osm_points) > 0) {
                found_facilities <- healthcare_data$osm_points
              }
            } else if (facility_type == "school") {
              school_data <- osmdata::opq(search_bbox) %>%
                osmdata::add_osm_feature(key = "amenity", value = "school") %>%
                osmdata::osmdata_sf()
              
        # Collect both points and polygons
          school_points <- NULL
          if (!is.null(school_data$osm_points) && nrow(school_data$osm_points) > 0) {
            school_points <- sf::st_geometry(school_data$osm_points) %>% 
              sf::st_sf(geometry = .)
          }
        
          if (!is.null(school_data$osm_polygons) && nrow(school_data$osm_polygons) > 0) {
            school_centroids <- sf::st_centroid(sf::st_geometry(school_data$osm_polygons)) %>%
              sf::st_sf(geometry = .)
            
            if (is.null(school_points)) {
              school_points <- school_centroids
            } else {
              school_points <- rbind(school_points, school_centroids)
            }
          }
        
        found_facilities <- school_points
      }
      
      # Increment search area
        if (is.null(found_facilities)) {
          search_multiplier <- search_multiplier + search_increment
          cat(sprintf("  No %s found, expanding search area (multiplier: %d)...\n", 
                      facility_type, search_multiplier))
        }
      }
    
    # Calculate travel time if facilities found
      if (!is.null(found_facilities) && nrow(found_facilities) > 0) {
        if (!is.null(friction_raster)) {
        # Check if we need to expand friction raster
          facilities_extent <- sf::st_bbox(found_facilities)
          friction_extent <- raster::extent(friction_raster)
          
        # If facilities are outside friction raster extent, return NA
          if (facilities_extent["xmin"] < friction_extent@xmin ||
              facilities_extent["xmax"] > friction_extent@xmax ||
              facilities_extent["ymin"] < friction_extent@ymin ||
              facilities_extent["ymax"] > friction_extent@ymax) {
            
            cat(sprintf("  Found %s outside friction raster extent. Returning NA.\n", facility_type))
            return(NA)
        }
        
        # For roads, convert to points
          if (facility_type == "paved_road") {
            suppressWarnings({
              facility_points <- sf::st_cast(sf::st_geometry(found_facilities), "POINT") %>%
                sf::st_sf(geometry = .)
            })
          } else {
            facility_points <- found_facilities
          }
          
          travel_time <- calculate_travel_time(center_point, facility_points, friction_raster)
          return(travel_time)
        } else {
        # If no friction raster provided, return NA
          cat("  No friction raster provided. Returning NA.\n")
          return(NA)
      }
    }
    
    # If still no facilities found after maximum search
      cat(sprintf("  No %s found within maximum search area.\n", facility_type))
      return(NA)
  }
  
# Compute Urbanicity Metric for Single Community
  compute_urbanicity <- function(community_data, name = NULL,
                                 roads = TRUE,
                                 shops = TRUE,
                                 healthcare = TRUE,
                                 transport = TRUE,
                                 financial = TRUE,
                                 schools = TRUE,
                                 cell_towers = TRUE,
                                 buildings = TRUE,
                                 nighttime_light = TRUE,
                                 population = TRUE,
                                 search_buffer = 1,  # degrees (100km)
                                 friction_surface_path = NULL,
                                 population_raster_path = NULL,
                                 nighttime_light_path = NULL) {
    # Load required libraries
      required_packages <- c("osmdata", "sf", "raster", "gdistance", "httr", "terra")
      for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(paste("Package", pkg, "needed for this function to work."))
        }
      }
    
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' needed for this function to work.")
    }
    
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package 'raster' needed for this function to work.")
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
      
    # Determine appropriate UTM zone for this location
      utm_zone <- floor((center_lon + 180) / 6) + 1
      hemisphere <- ifelse(center_lat < 0, "+south", "")
      utm_crs <- paste0("+proj=utm +zone=", utm_zone, " ", hemisphere, " +datum=WGS84")
      cat("  Location UTM zone:", utm_zone, hemisphere, "\n")
      
    # Create center point sf object
      center_point <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
      
    # Load friction surface (local files only)
      friction_raster <- NULL
      if (healthcare || schools || roads) {
        tryCatch({
          cat("  Loading friction surface data...\n")
          
          # Check if a local friction surface path is provided
            if (is.null(friction_surface_path) || !file.exists(friction_surface_path)) {
              stop("Error: Local friction surface file is required but not provided or does not exist.\n",
                   "Please provide a valid path to a local friction surface file using the 'friction_surface_path' parameter.\n",
                   "Expected file format: GeoTIFF (.tif or .tiff)")
            }
          
            cat("  Using local friction surface file:", friction_surface_path, "\n")
            friction_global <- raster::raster(friction_surface_path)
          
          # Crop to a larger area around the community (for searching)
            extent_buffered <- raster::extent(
              bbox["left"] - search_buffer,
              bbox["right"] + search_buffer,
              bbox["bottom"] - search_buffer,
              bbox["top"] + search_buffer
            )
            
            friction_raster <- raster::crop(friction_global, extent_buffered)
            cat("  Friction surface loaded successfully\n")
            
          }, error = function(e) {
            cat("  Error loading friction surface:", e$message, "\n")
            cat("  Please ensure:\n")
            cat("    1. The friction_surface_path parameter points to a valid GeoTIFF file\n")
            cat("    2. The file exists and is readable\n")
            cat("    3. The file contains valid raster data\n")
            friction_raster <- NULL
          })
        }
      
    # Create center point sf object
      center_point <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
      center_lat <- as.numeric(community_data$lat)
      center_lon <- as.numeric(community_data$lon)
    
    # Ratio of Paved to Unpaved Roads, Percent Paved Roads, and Travel Time to Nearest Paved Road
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
              
            # Initialize road_type column first
              roads_sf$road_type <- "other"
              
            # Classify roads based on surface tag
              if (!all(is.na(roads_sf$surface))) {
              # Only try to assign if there are non-NA surface values
                paved_mask <- !is.na(roads_sf$surface) & roads_sf$surface %in% paved_surfaces
                if (any(paved_mask)) {
                  roads_sf$road_type[paved_mask] <- "paved"
                }
                
                unpaved_mask <- !is.na(roads_sf$surface) & roads_sf$surface %in% unpaved_surfaces
                if (any(unpaved_mask)) {
                  roads_sf$road_type[unpaved_mask] <- "unpaved"
                }
              }
              
            # Classify roads based on highway tag for tracks and paths
              if ("highway" %in% colnames(roads_sf) && !all(is.na(roads_sf$highway))) {
                unpaved_highway_types <- c("track", "path")
                highway_mask <- !is.na(roads_sf$highway) & roads_sf$highway %in% unpaved_highway_types
                if (any(highway_mask)) {
                  roads_sf$road_type[highway_mask] <- "unpaved"
                }
              }
            
            # Calculate road lengths
              # Determine appropriate UTM zone based on longitude
                center_lon <- (bbox["left"] + bbox["right"]) / 2
                utm_zone <- floor((center_lon + 180) / 6) + 1
                hemisphere <- ifelse(center_lat < 0, "+south", "")
                utm_crs <- paste0("+proj=utm +zone=", utm_zone, " ", hemisphere, " +datum=WGS84")
                roads_sf <- sf::st_transform(roads_sf, crs = sf::st_crs(utm_crs))
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
                results$paved_to_unpaved_ratio <- NA
              }
              
            # Calculate percent paved
              if (unpaved_length > 0) {
                results$pct_paved_roads <- 100 * as.numeric(paved_length / total_length)
              } else {
                results$pct_paved_roads <- NA
              }
            
            # Calculate travel time to nearest paved road
              paved_roads <- roads_sf[roads_sf$road_type == "paved", ]
              
              if (nrow(paved_roads) > 0 && !is.null(friction_raster)) {
              # Convert lines to points for travel time calculation
                  suppressWarnings({
                    paved_points <- sf::st_cast(sf::st_geometry(paved_roads), "POINT") %>%
                      sf::st_sf(geometry = .)
                  })
                  results$travel_time_paved_road_min <- calculate_travel_time(center_point, paved_points, friction_raster)
                } else {
              # No paved roads locally, use progressive search
                results$travel_time_paved_road_min <- search_facilities_progressive(
                  center_point = center_point,
                  bbox = bbox,
                  facility_type = "paved_road",
                  friction_raster = friction_raster,
                  max_search_multiplier = 20
              )
            }
            
          } else {
            # No roads found in the initial area at all
              cat("  No roads found in the immediate area.\n")
              results$paved_to_unpaved_ratio <- NA
              results$pct_paved_roads <- NA
              
            # Still try progressive search for paved roads in wider area
              results$travel_time_paved_road_min <- search_facilities_progressive(
                center_point = center_point,
                bbox = bbox,
                facility_type = "paved_road",
                friction_raster = friction_raster,
                max_search_multiplier = 20
              )
            }
          }, error = function(e) {
            cat("  Error processing roads:", e$message, "\n")
            results$paved_to_unpaved_ratio <- NA
            results$pct_paved_roads <- NA
            results$travel_time_paved_road_min <- NA
          })
      }
               
    
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
    
    # Number of Public Transport Stops
      if (transport) {
        tryCatch({
          cat("  Analyzing transport infrastructure...\n")
          
          transport_count <- 0
          
        # Query public transport nodes
          transport_data1 <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "public_transport") %>%
            osmdata::osmdata_sf()
          
          if (!is.null(transport_data1$osm_points) && nrow(transport_data1$osm_points) > 0) {
            transport_count <- transport_count + nrow(transport_data1$osm_points)
          }
          
        # Query bus stops
          bus_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "highway", value = "bus_stop") %>%
            osmdata::osmdata_sf()
          
          if (!is.null(bus_data$osm_points) && nrow(bus_data$osm_points) > 0) {
            transport_count <- transport_count + nrow(bus_data$osm_points)
          }
          
        # Query railway stations
          rail_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "railway", value = c("station", "halt")) %>%
            osmdata::osmdata_sf()
          
          if (!is.null(rail_data$osm_points) && nrow(rail_data$osm_points) > 0) {
            transport_count <- transport_count + nrow(rail_data$osm_points)
          }
          
        # Check for polygons too (some stations are mapped as buildings)
          if (!is.null(rail_data$osm_polygons) && nrow(rail_data$osm_polygons) > 0) {
            transport_count <- transport_count + nrow(rail_data$osm_polygons)
          }
          
          results$n_transport_stops <- transport_count
          
        }, error = function(e) {
          cat("  Error processing transport infrastructure:", e$message, "\n")
          results$n_transport_stops <- NA
        })
      }
      
    # Travel time to nearest Healthcare Facility
      if (healthcare) {
        tryCatch({
          cat("  Analyzing healthcare facilities...\n")
                
          results$travel_time_healthcare_min <- search_facilities_progressive(
            center_point = center_point,
            bbox = bbox,
            facility_type = "healthcare",
            friction_raster = friction_raster,
            max_search_multiplier = 50
            )
                
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
    
    # Travel time to nearest School
      if (schools) {
        tryCatch({
        cat("  Analyzing schools...\n")
                
        results$travel_time_school_min <- search_facilities_progressive(
          center_point = center_point,
          bbox = bbox,
          facility_type = "school",
          friction_raster = friction_raster,
          max_search_multiplier = 50
        )
                
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
          # Determine appropriate UTM zone based on longitude
            center_lon <- (bbox["left"] + bbox["right"]) / 2
            center_lat <- (bbox["bottom"] + bbox["top"]) / 2
            utm_zone <- floor((center_lon + 180) / 6) + 1
            hemisphere <- ifelse(center_lat < 0, "+south", "")
            utm_crs <- paste0("+proj=utm +zone=", utm_zone, " ", hemisphere, " +datum=WGS84")
            
            cat("  Using UTM zone", utm_zone, hemisphere, "for area calculations\n")
            
          # Transform to appropriate projection for accurate area calculation
            buildings_projected <- sf::st_transform(
              buildings_data$osm_polygons, 
              crs = sf::st_crs(utm_crs)
            )
            
          # Calculate area of each building
            buildings_projected$area <- sf::st_area(buildings_projected)
            
          # Calculate total building area
            total_building_area <- sum(as.numeric(buildings_projected$area), na.rm = TRUE)
            
          # Calculate bounding box area
            poly_projected <- sf::st_transform(
              poly, 
              crs = sf::st_crs(utm_crs)  # Use the same UTM zone!
            )
            bbox_area <- as.numeric(sf::st_area(poly_projected))
            
            cat("  Total building area:", round(total_building_area), "m²\n")
            cat("  Bounding box area:", round(bbox_area), "m²\n")
            
          # Check for potential issues
            if (bbox_area == 0) {
              cat("  WARNING: Bounding box area is 0!\n")
              results$building_density <- NA
            } else if (!is.finite(total_building_area) || !is.finite(bbox_area)) {
              cat("  WARNING: Non-finite values in area calculation\n")
              results$building_density <- NA
            } else {
            # Calculate building density (% of area covered by buildings)
              building_density <- (total_building_area / bbox_area) * 100
              
             # Check if result is valid
              if (is.finite(building_density) && building_density >= 0 && building_density <= 100) {
                results$building_density <- as.numeric(building_density)
              } else {
                cat("  WARNING: Building density outside valid range (0-100%): ", building_density, "\n")
                results$building_density <- NA
              }
            }
          } else {
            results$building_density <- 0
          }
        }, error = function(e) {
          cat("  Error processing buildings:", e$message, "\n")
          results$building_density <- NA
        })
      }
      
      # Nighttime Light (from VIIRS)
        if (nighttime_light) {
          tryCatch({
            cat("  Analyzing nighttime light...\n")
            
            viirs_raster <- raster(nighttime_light_path)
            
            # Convert to sf object and add id column
              poly_for_nighttime_light <- sf::st_sf(geometry = poly)
              poly_for_nighttime_light$id <- 1
              
              nighttime_light <- get_light_mean(poly_for_nighttime_light, viirs_raster)[1]
              
              results$nighttime_light <- as.numeric(nighttime_light)
            
          }, error = function(e) {
            cat("  Error processing nighttime light:", e$message, "\n")
            results$nighttime_light <- NA
          })
        }
    
    # Population Density from SEDAC/NASA Gridded Population of the World
      if (population) {
        tryCatch({
          cat("  Calculating population density...\n")
          
          if (!is.null(population_raster_path) && file.exists(population_raster_path)) {
            cat("  Using local population raster file:", population_raster_path, "\n")
            
          # Load the population raster
            population_global <- raster::raster(population_raster_path)
            
          # Crop to community area
            extent_community <- raster::extent(
              bbox["left"],
              bbox["right"],
              bbox["bottom"],
              bbox["top"]
            )
            
            population_raster <- raster::crop(population_global, extent_community)
            cat("  Population raster loaded and cropped successfully\n")
            
          # Calculate population density
          # Use appropriate UTM projection for accurate area calculation
            center_lon <- (bbox["left"] + bbox["right"]) / 2
            center_lat <- (bbox["bottom"] + bbox["top"]) / 2
            utm_zone <- floor((center_lon + 180) / 6) + 1
            hemisphere <- ifelse(center_lat < 0, "+south", "")
            utm_crs <- paste0("+proj=utm +zone=", utm_zone, " ", hemisphere, " +datum=WGS84")
            
            poly_projected <- sf::st_transform(poly, crs = sf::st_crs(utm_crs))
            
          # Convert sfc to sf object for extract function
            poly_projected_sf <- sf::st_sf(geometry = poly_projected)
            
          # Calculate area in square kilometers
            area_m2 <- as.numeric(sf::st_area(poly_projected))
            area_km2 <- area_m2 / 1000000
            
          # Extract population values within the community boundary
          # Reproject raster to match the projected polygon
            population_raster_proj <- raster::projectRaster(
              population_raster, 
              crs = raster::crs(poly_projected_sf)
            )
            
          # Extract values
            pop_values <- raster::extract(
              population_raster_proj, 
              poly_projected_sf, 
              fun = sum, 
              na.rm = TRUE
            )
            
            total_population <- sum(pop_values, na.rm = TRUE)
            
          # Calculate population density (people per km²)
            results$pop_density <- as.numeric(total_population / area_km2)
            
            cat("  Total population:", round(total_population), "\n")
            cat("  Area (km²):", round(area_km2, 2), "\n")
            cat("  Population density (per km²):", round(results$pop_density, 2), "\n")
            
          } else {
            cat("  No population raster file provided or file not found.\n")
            results$pop_density <- NA
          }
          
        }, error = function(e) {
          cat("  Error calculating population density:", e$message, "\n")
          results$pop_density <- NA
        })
      }
        
    # Convert results list to a dataframe
      return(as.data.frame(lapply(results, function(x) if(length(x) == 0) NA else x)))
  }
 
# Apply compute_urbanicity Function to a List of Communities
# communities_list is the output from the bounding boxes function above.
# metrics is a vector of the measures you want to compute.
  compute_urbanicity_iterative <- function(communities_list, 
                                           metrics = c("all"),
                                           search_buffer = 1,  # degrees (100km)
                                           friction_surface_path = NULL,
                                           population_raster_path = NULL,
                                           nighttime_light_path = NULL) {
    
    # Set up which metrics to analyze based on user input
      do_roads <- "all" %in% metrics || "roads" %in% metrics
      do_shops <- "all" %in% metrics || "shops" %in% metrics
      do_healthcare <- "all" %in% metrics || "healthcare" %in% metrics
      do_transport <- "all" %in% metrics || "transport" %in% metrics
      do_financial <- "all" %in% metrics || "financial" %in% metrics
      do_schools <- "all" %in% metrics || "schools" %in% metrics
      do_cell_towers <- "all" %in% metrics || "cell_towers" %in% metrics
      do_buildings <- "all" %in% metrics || "buildings" %in% metrics
      do_nighttime_light <- "all" %in% metrics || "nighttime_light" %in% metrics
      do_population <- "all" %in% metrics || "population" %in% metrics
    
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
            nighttime_light = do_nighttime_light,
            population = do_population,
            friction_surface_path = friction_surface_path,
            population_raster_path = population_raster_path,
            nighttime_light_path = nighttime_light_path
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
  
  
  