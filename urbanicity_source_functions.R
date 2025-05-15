# Urbanicity Index for Anthropological Fieldsites - Source Functions
# Developed as Part of the Population Ecology, Aging, and Health Network (PEcAHN)

# Create Boundary for Each Community
# Provide specific boundary coordinates (custom_bbox) or use a pre-specified distance around a central location (distance_km).
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
  
# Compute Urbanicity Metric For Single Location
  compute_urbanicity <- function(location_data, name = NULL,
                                 roads = TRUE,
                                 shops = TRUE,
                                 healthcare = TRUE,
                                 transport = TRUE,
                                 financial = TRUE,
                                 schools = TRUE,
                                 cell_towers = TRUE,
                                 buildings = TRUE) {
    # Load required libraries
      if (!requireNamespace("osmdata", quietly = TRUE)) {
        stop("Package 'osmdata' needed for this function to work.")
      }
      if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Package 'sf' needed for this function to work.")
      }
    
    # Import functions
      library(osmdata)
      library(sf)
    
    # Extract the bbox from the location_data
      bbox <- location_data$bbox
    
    # If name is not provided, construct it from metadata
      if (is.null(name)) {
        name <- paste(location_data$project, location_data$ethnicity, location_data$community, sep = " - ")
      }
    
    # Print processing status
      cat("Processing location:", name, "\n")
    
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
    
    # Center point for calculations
      center_lat <- location_data$lat
      center_lon <- location_data$lon
    
    # Initialize results with base information
      results <- list(
        project = location_data$project,
        ethnicity = location_data$ethnicity,
        community = location_data$community
      )
    
    # Define bounding box for osmdata in the correct format
      osm_bbox <- c(bbox["left"], bbox["bottom"], bbox["right"], bbox["top"])
    
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
            } else {
              results$paved_to_unpaved_ratio <- NA
              results$pct_paved_roads <- NA
              cat("  No roads found in the area.\n")
            }
            }, error = function(e) {
              cat("  Error processing roads:", e$message, "\n")
              results$paved_to_unpaved_ratio <- NA
              results$pct_paved_roads <- NA
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
            # Count formal shops
            formal_shops <- c("supermarket", "convenience", "department_store", "mall", "grocery", "general")
            formal_count <- sum(shops_data$osm_points$shop %in% formal_shops)
            results$n_formal_shops <- formal_count
          } else {
            results$n_formal_shops <- 0
          }
        }, error = function(e) {
          cat("  Error processing shops:", e$message, "\n")
          results$n_formal_shops <- NA
        })
      }
    
    # Number of Healthcare Facilities and Distance to Nearest Healthcare Facility
    # NOTE: jUST yes / no for now, was getting oddly large values for some communities
      if (healthcare) {
        tryCatch({
          cat("  Analyzing healthcare facilities...\n")
          healthcare_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "amenity", value = c("hospital", "clinic", "doctors")) %>%
            osmdata::osmdata_sf()
          
          # Check if healthcare facilities exist (yes/no)
          results$healthcare <- "No"
          if (!is.null(healthcare_data$osm_points) && nrow(healthcare_data$osm_points) > 0) {
            results$healthcare <- "Yes"
          }
          
        }, error = function(e) {
          cat("  Error processing healthcare facilities:", e$message, "\n")
          results$healthcare <- NA
        })
      }
    
    # Number of Public Transportation Stops
      if (transport) {
        tryCatch({
          cat("  Analyzing public transportation...\n")
          
          # Combined query for all transport stops
          transport_queries <- list(
            osmdata::opq(osm_bbox) %>% osmdata::add_osm_feature(key = "public_transport", value = c("stop_position", "station")),
            osmdata::opq(osm_bbox) %>% osmdata::add_osm_feature(key = "railway", value = "station"),
            osmdata::opq(osm_bbox) %>% osmdata::add_osm_feature(key = "highway", value = "bus_stop"),
            osmdata::opq(osm_bbox) %>% osmdata::add_osm_feature(key = "amenity", value = "bus_station")
          )
          
          # Count total stops
          total_stops <- 0
          
          # Process each query
          for (query in transport_queries) {
            data <- osmdata::osmdata_sf(query)
            if (!is.null(data$osm_points) && nrow(data$osm_points) > 0) {
              total_stops <- total_stops + nrow(data$osm_points)
            }
          }
          
          results$n_transport_stops <- total_stops
          
        }, error = function(e) {
          cat("  Error processing public transport:", e$message, "\n")
          results$n_transport_stops <- NA
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
      if (schools) {
        tryCatch({
          cat("  Analyzing schools...\n")
          
          # Query for schools using amenity=school
          school_data <- osmdata::opq(osm_bbox) %>%
            osmdata::add_osm_feature(key = "amenity", value = "school") %>%
            osmdata::osmdata_sf()
          
          # Check if any schools exist (either as points or polygons)
          has_school <- FALSE
          
          if ((!is.null(school_data$osm_points) && nrow(school_data$osm_points) > 0) || 
              (!is.null(school_data$osm_polygons) && nrow(school_data$osm_polygons) > 0)) {
            has_school <- TRUE
          }
          
          # Return "yes" or "no"
          results$school <- ifelse(has_school, "yes", "no")
          
        }, error = function(e) {
          cat("  Error processing schools:", e$message, "\n")
          results$school <- NA
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
 
# Apply compute_urbanicity Function to a List of Locations
# locations_list is the output from the bounding boxes function above.
# metrics is a vector of the measures you want to compute.
  compute_urbanicity_iterative <- function(locations_list, metrics = c("all")) {
    
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
    
    # Process each location
      for (i in seq_along(locations_list)) {
        location_name <- names(locations_list)[i]
        location_data <- locations_list[[i]]
        
        # Process the location with selected metrics
          result <- compute_urbanicity(
            location_data, 
            name = location_name,
            roads = do_roads,
            shops = do_shops,
            healthcare = do_healthcare,
            transport = do_transport,
            financial = do_financial,
            schools = do_schools,
            cell_towers = do_cell_towers,
            buildings = do_buildings
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
  
  
  
  
  