# PEcAHN Exploratory Analysis
# Created by Tyler Barrett on December 23, 2025

#################
#   PREAMBLE    #
#################

# Load Packages
library(tidyverse)
library(viridis)
library(gridExtra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(osmdata)

# Set File Path
  fp <- "C:/Users/tmb80/OneDrive - Duke University/Active Projects/pecahn/test_data"

# Load Data
  community_coords <- read_csv(paste0(fp, "/community_coordinates_5km.csv"))
  world <- ne_countries(scale = "medium", returnclass = "sf")

# Split Projects Covering a Large Geographic Space to Minimize Download Size  
  # Split REACH Into Separate Single-Point Projects
    community_coords <- community_coords %>%
      mutate(project = case_when(
        project == "REACH" & community == "Mississippi" ~ "REACH - Mississippi",
        project == "REACH" & community == "SW Illinois" ~ "REACH - Illinois",
        TRUE ~ project
      ))
    
  # Split Orang Asli Into Three Regions
    community_coords <- community_coords %>%
      mutate(project = case_when(
        project == "Orang Asli" & lat >= 4.8 ~ "Orang Asli - North",
        project == "Orang Asli" & lat >= 3.6 & lat < 4.8 ~ "Orang Asli - Central",
        project == "Orang Asli" & lat < 3.6 ~ "Orang Asli - South",
        TRUE ~ project
      ))

  # Split Turkana Into Four Regions
    community_coords <- community_coords %>%
      mutate(project = case_when(
        project == "Turkana" & lat >= 3.1 ~ "Turkana - North",
        project == "Turkana" & lat >= 1.5 & lat < 3.1 ~ "Turkana - Central North",
        project == "Turkana" & lat >= 0 & lat < 1.5 ~ "Turkana - Central South",
        project == "Turkana" & lat < 0 ~ "Turkana - South",
        TRUE ~ project
      ))

# Create Cache Directory to Save OSM Data
  cache_dir <- paste0(fp, "/osm_cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir)
    cat("Created cache directory:", cache_dir, "\n")
  }

# Create Maps Directory to Save Output
  maps_dir <- paste0(fp, "/maps")
  if (!dir.exists(maps_dir)) {
    dir.create(maps_dir)
    cat("Created maps directory:", maps_dir, "\n")
  }


#######################
#   MAP COMMUNITIES   #
#######################

# Global Map of Study Communities
  global_map <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "grey60", size = 0.3) +
    geom_point(data = community_coords,
               aes(x = lon, y = lat, fill = project),
               size = 1,
               alpha = 0.9,
               shape = 21,
               color = "black",
               stroke = 0.3) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.key.size = unit(0.8, "lines")
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Project",
      title = "Community Locations by Project"
    ) +
    scale_fill_viridis_d(option = "turbo")
  global_map

# Save Global Map
  ggsave(paste0(fp, "/maps/global_map.png"), plot = global_map, 
         width = 12, height = 7, dpi = 300, bg = "white")

# Individual Project Maps with OSM Data
  # Common Theme For Project Maps
    common_theme <- theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey90", size = 0.2),
        panel.background = element_rect(fill = "aliceblue", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 9)
      )

  # Get List of Unique Projects
    projects <- unique(community_coords$project)

  # Create Map For Each Project
    for (proj in projects) {
      cat("Creating map for:", proj, "\n")
      
      #	Filter data for this project
        project_data <- community_coords %>% filter(project == proj)
      
      #	Skip if no data
        if (nrow(project_data) == 0) {
          cat("  No data for project:", proj, "\n\n")
          next
      }
  
    #	Calculate bounding box with buffer
      lon_range <- range(project_data$lon)
      lat_range <- range(project_data$lat)
    
    #	Use 10km buffer for single communities (~0.09 degrees), larger for multi-community
      min_buffer <- if (nrow(project_data) == 1) 0.09 else 0.5
    
      lon_buffer <- max((lon_range[2] - lon_range[1]) * 0.3, min_buffer)
      lat_buffer <- max((lat_range[2] - lat_range[1]) * 0.3, min_buffer)
      lon_lim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
      lat_lim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
      bbox <- c(lon_lim[1], lat_lim[1], lon_lim[2], lat_lim[2])
    
    #	Check cache or download OSM data (all features in one query)
      cache_file <- file.path(cache_dir, paste0(proj, "_osm_v2.rds"))
    
      if (file.exists(cache_file)) {
        cat("  Loading cached OSM data\n")
        osm_data <- readRDS(cache_file)
        osm_roads <- osm_data$roads
        osm_waterways <- osm_data$waterways
        osm_natural <- osm_data$natural
        osm_protected <- osm_data$protected
        osm_conservation <- osm_data$conservation
      } else {
        cat("  Downloading OSM data (will be cached)...\n")
      
      #	Major roads only (motorway, trunk, primary, secondary)
        osm_roads <- tryCatch({
          opq(bbox) %>%
            add_osm_feature(key = "highway", 
                            value = c("motorway", "trunk", "primary", "secondary")) %>%
            osmdata_sf()
        }, error = function(e) NULL)
      
      #	Major waterways only (river, canal)
        osm_waterways <- tryCatch({
          opq(bbox) %>%
            add_osm_feature(key = "waterway", value = c("river", "canal")) %>%
            osmdata_sf()
        }, error = function(e) NULL)
      
      #	Natural features (water, wood, forest)
        osm_natural <- tryCatch({
          opq(bbox) %>%
            add_osm_feature(key = "natural", value = c("water", "wood", "forest")) %>%
            osmdata_sf()
        }, error = function(e) NULL)
      
      #	National parks and protected area boundaries
        osm_protected <- tryCatch({
          opq(bbox) %>%
            add_osm_feature(key = "boundary", value = c("national_park", "protected_area")) %>%
            osmdata_sf()
        }, error = function(e) NULL)
      
      #	Conservation and forest landuse
        osm_conservation <- tryCatch({
          opq(bbox) %>%
            add_osm_feature(key = "landuse", value = c("conservation", "forest")) %>%
            osmdata_sf()
        }, error = function(e) NULL)
      
      #	Save to cache
        osm_data <- list(
          roads = osm_roads, 
          waterways = osm_waterways, 
          natural = osm_natural,
          protected = osm_protected,
          conservation = osm_conservation
        )
      saveRDS(osm_data, cache_file)
      cat("  Data cached to:", cache_file, "\n")
    }
  
    #	Build plot with base layer
      p <- ggplot() +
        geom_sf(data = world, fill = "antiquewhite1", color = "grey50", size = 0.3)
    
    #	Add protected areas (darker green)
      if (!is.null(osm_protected)) {
      if (!is.null(osm_protected$osm_multipolygons)) {
        p <- p + geom_sf(data = osm_protected$osm_multipolygons,
                         fill = "darkgreen", color = "forestgreen", alpha = 0.3, size = 0.4)
      }
      if (!is.null(osm_protected$osm_polygons)) {
        p <- p + geom_sf(data = osm_protected$osm_polygons,
                         fill = "darkgreen", color = "forestgreen", alpha = 0.3, size = 0.4)
      }
    }
    
    #	Add conservation landuse
      if (!is.null(osm_conservation)) {
      if (!is.null(osm_conservation$osm_multipolygons)) {
        p <- p + geom_sf(data = osm_conservation$osm_multipolygons,
                         fill = "forestgreen", color = NA, alpha = 0.3)
      }
      if (!is.null(osm_conservation$osm_polygons)) {
        p <- p + geom_sf(data = osm_conservation$osm_polygons,
                         fill = "forestgreen", color = NA, alpha = 0.3)
      }
    }
    
    #	Add OSM natural features (light green)
      if (!is.null(osm_natural)) {
      if (!is.null(osm_natural$osm_polygons)) {
        p <- p + geom_sf(data = osm_natural$osm_polygons,
                         fill = "lightgreen", color = NA, alpha = 0.3)
      }
      if (!is.null(osm_natural$osm_multipolygons)) {
        p <- p + geom_sf(data = osm_natural$osm_multipolygons,
                         fill = "lightgreen", color = NA, alpha = 0.3)
      }
    }
    
    #	Add OSM waterways
      if (!is.null(osm_waterways)) {
      if (!is.null(osm_waterways$osm_lines)) {
        p <- p + geom_sf(data = osm_waterways$osm_lines,
                         color = "steelblue", size = 0.5, alpha = 0.7)
      }
      if (!is.null(osm_waterways$osm_polygons)) {
        p <- p + geom_sf(data = osm_waterways$osm_polygons,
                         fill = "lightblue", color = "steelblue", size = 0.2)
      }
    }
    
    #	Add OSM roads
      if (!is.null(osm_roads)) {
      if (!is.null(osm_roads$osm_lines)) {
        p <- p + geom_sf(data = osm_roads$osm_lines,
                         color = "darkgray", size = 0.4, alpha = 0.5)
      }
    }
    
    #	Add community points
      p <- p +
        geom_point(data = project_data,
                   aes(x = lon, y = lat, fill = project),
                   size = 2,
                   alpha = 0.9,
                   shape = 21,
                   color = "black",
                   stroke = 0.6) +
        coord_sf(xlim = lon_lim, ylim = lat_lim, expand = FALSE) +
        common_theme +
        labs(
          title = paste0(proj, " Communities (n=", nrow(project_data), ")"),
          x = "Longitude",
          y = "Latitude",
          fill = "Project"
        ) +
        scale_fill_viridis_d(option = "turbo")
    
    #	Add labels for smaller projects
      if (nrow(project_data) <= 15) {
      p <- p + geom_text(data = project_data,
                         aes(x = lon, y = lat, label = community),
                         size = 2, hjust = -0.2, vjust = 0.5,
                         check_overlap = TRUE, fontface = "bold")
    }
    
    #	Display map
      print(p)
    
    #	Save map
      ggsave(paste0(fp, "/maps/", proj, "_map.png"), plot = p, 
            width = 10, height = 8, dpi = 300, bg = "white")
    
    cat("\n")
  }
  
cat("=== All maps created ===\n")

# Project Maps With Buffer Zone Circles
  # Buffer radii in kilometers
    buffer_radii_km <- c(1, 2, 3, 4, 5)

  cat("\n=== Creating maps with buffer zones for", length(projects), "projects ===\n\n")

  # Create map for each project with buffer circles
    for (proj in projects) {
      cat("Creating buffer map for:", proj, "\n")
      
      #	Filter data for this project
        project_data <- community_coords %>% filter(project == proj)
      
      #	Skip if no data
        if (nrow(project_data) == 0) {
          cat("  No data for project:", proj, "\n\n")
          next
        }
      
      #	Convert to sf object with geographic CRS
        project_sf <- st_as_sf(project_data, coords = c("lon", "lat"), crs = 4326)
      
      #	Transform to projected CRS for accurate buffering (meters)
        project_utm <- st_transform(project_sf, crs = 3857)
      
      #	Create buffer zones for each radius
        buffers_list <- list()
        for (radius_km in buffer_radii_km) {
          radius_m <- radius_km * 1000
          buffer_zone <- st_buffer(project_utm, dist = radius_m)
          buffer_zone <- st_transform(buffer_zone, crs = 4326)
          buffer_zone$radius_km <- radius_km
          buffers_list[[as.character(radius_km)]] <- buffer_zone
        }
        all_buffers <- do.call(rbind, buffers_list)
      
      #	Calculate bounding box with buffer
        lon_range <- range(project_data$lon)
        lat_range <- range(project_data$lat)
      
      #	Use 10km buffer for single communities (~0.09 degrees), larger for multi-community
        min_buffer <- if (nrow(project_data) == 1) 0.09 else 0.5
        
        lon_buffer <- max((lon_range[2] - lon_range[1]) * 0.3, min_buffer)
        lat_buffer <- max((lat_range[2] - lat_range[1]) * 0.3, min_buffer)
        lon_lim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
        lat_lim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
        bbox <- c(lon_lim[1], lat_lim[1], lon_lim[2], lat_lim[2])
      
      #	Check cache or download OSM data (all features in one query)
        cache_file <- file.path(cache_dir, paste0(proj, "_osm_v2.rds"))
      
        if (file.exists(cache_file)) {
          cat("  Loading cached OSM data\n")
          osm_data <- readRDS(cache_file)
          osm_roads <- osm_data$roads
          osm_waterways <- osm_data$waterways
          osm_natural <- osm_data$natural
          osm_protected <- osm_data$protected
          osm_conservation <- osm_data$conservation
        } else {
        cat("  Downloading OSM data (will be cached)...\n")
        
        #	Major roads only (motorway, trunk, primary, secondary)
          osm_roads <- tryCatch({
            opq(bbox) %>%
              add_osm_feature(key = "highway", 
                              value = c("motorway", "trunk", "primary", "secondary")) %>%
              osmdata_sf()
          }, error = function(e) NULL)
        
        #	Major waterways only (river, canal)
          osm_waterways <- tryCatch({
            opq(bbox) %>%
              add_osm_feature(key = "waterway", value = c("river", "canal")) %>%
              osmdata_sf()
          }, error = function(e) NULL)
        
        #	Natural features (water, wood, forest)
          osm_natural <- tryCatch({
            opq(bbox) %>%
              add_osm_feature(key = "natural", value = c("water", "wood", "forest")) %>%
              osmdata_sf()
          }, error = function(e) NULL)
        
        #	National parks and protected area boundaries
          osm_protected <- tryCatch({
            opq(bbox) %>%
              add_osm_feature(key = "boundary", value = c("national_park", "protected_area")) %>%
              osmdata_sf()
          }, error = function(e) NULL)
        
        #	Conservation and forest landuse
          osm_conservation <- tryCatch({
            opq(bbox) %>%
              add_osm_feature(key = "landuse", value = c("conservation", "forest")) %>%
              osmdata_sf()
          }, error = function(e) NULL)
        
        #	Save to cache
          osm_data <- list(
            roads = osm_roads, 
            waterways = osm_waterways, 
            natural = osm_natural,
            protected = osm_protected,
            conservation = osm_conservation
          )
          saveRDS(osm_data, cache_file)
          cat("  Data cached to:", cache_file, "\n")
      }
      
      #	Build plot with base layer
        p_buffer <- ggplot() +
          geom_sf(data = world, fill = "antiquewhite1", color = "grey50", size = 0.3)
      
      #	Add protected areas (darker green)
        if (!is.null(osm_protected)) {
          if (!is.null(osm_protected$osm_multipolygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_protected$osm_multipolygons,
                                           fill = "darkgreen", color = "forestgreen", alpha = 0.3, size = 0.4)
          }
          if (!is.null(osm_protected$osm_polygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_protected$osm_polygons,
                                           fill = "darkgreen", color = "forestgreen", alpha = 0.3, size = 0.4)
          }
        }
      
      #	Add conservation landuse
        if (!is.null(osm_conservation)) {
          if (!is.null(osm_conservation$osm_multipolygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_conservation$osm_multipolygons,
                                           fill = "forestgreen", color = NA, alpha = 0.3)
          }
          if (!is.null(osm_conservation$osm_polygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_conservation$osm_polygons,
                                           fill = "forestgreen", color = NA, alpha = 0.3)
          }
        }
      
      #	Add OSM natural features (light green)
        if (!is.null(osm_natural)) {
          if (!is.null(osm_natural$osm_polygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_natural$osm_polygons,
                                           fill = "lightgreen", color = NA, alpha = 0.3)
          }
          if (!is.null(osm_natural$osm_multipolygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_natural$osm_multipolygons,
                                           fill = "lightgreen", color = NA, alpha = 0.3)
          }
        }
      
      #	Add OSM waterways
        if (!is.null(osm_waterways)) {
          if (!is.null(osm_waterways$osm_lines)) {
            p_buffer <- p_buffer + geom_sf(data = osm_waterways$osm_lines,
                                           color = "steelblue", size = 0.5, alpha = 0.7)
          }
          if (!is.null(osm_waterways$osm_polygons)) {
            p_buffer <- p_buffer + geom_sf(data = osm_waterways$osm_polygons,
                                           fill = "lightblue", color = "steelblue", size = 0.2)
          }
      }
      
      #	Add OSM roads
        if (!is.null(osm_roads)) {
          if (!is.null(osm_roads$osm_lines)) {
            p_buffer <- p_buffer + geom_sf(data = osm_roads$osm_lines,
                                           color = "darkgray", size = 0.4, alpha = 0.5)
          }
      }
      
      #	Add buffer circles
        p_buffer <- p_buffer +
          geom_sf(data = all_buffers,
                  aes(fill = factor(radius_km)),
                  alpha = 0.2,
                  color = "black",
                  size = 0.3)
      
      #	Add community points
        p_buffer <- p_buffer +
          geom_point(data = project_data,
                     aes(x = lon, y = lat),
                     size = 2,
                     alpha = 0.9,
                     shape = 21,
                     fill = "red",
                     color = "black",
                     stroke = 0.6) +
          coord_sf(xlim = lon_lim, ylim = lat_lim, expand = FALSE) +
          common_theme +
          labs(
            title = paste0(proj, " Communities with Buffer Zones (n=", nrow(project_data), ")"),
            x = "Longitude",
            y = "Latitude",
            fill = "Buffer (km)"
          ) +
          scale_fill_viridis_d(option = "plasma")
      
      #	Add labels for smaller projects
        if (nrow(project_data) <= 15) {
          p_buffer <- p_buffer + geom_text(data = project_data,
                                           aes(x = lon, y = lat, label = community),
                                           size = 2, hjust = -0.2, vjust = 0.5,
                                           check_overlap = TRUE, fontface = "bold")
        }
      
      #	Display map
        print(p_buffer)
      
      #	Save map
        ggsave(paste0(fp, "/maps/", proj, "_buffer_map.png"), plot = p_buffer, 
               width = 10, height = 8, dpi = 300, bg = "white")
      
      cat("\n")
    }

cat("=== All buffer maps created ===\n")
cat("\nNote: OSM data cached in:", cache_dir, "\n")
cat("To force re-download, delete the cache directory\n")


