# Urbanicity Index for Anthropological Fieldsites - Vignette
# Developed as Part of the Population Ecology, Aging, and Health Network (PEcAHN)

# To run the code below, you will need to install the following packages.
  # install.packages("osmdata")
  # install.packages("sf")
  # install.packages("terra")
  # install.packages("raster")
  # install.packages("gdistance")
  # install.packages("httr")
  # install.packages("ggplot2")
  # install.packages("dplyr")

# Set Your File Path
# This is the path to the folder where you downloaded the test data and source functions.
# For example, mine is: "C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn/test_data"
  fp = "your_file_path"

# Load Test Data
  # The test_data_boundaries data table provides northeast, southeast, southwest, and northwest coordinates for seven communities.
    test_data_boundaries <- read.csv(paste0(fp, "/test_data_boundaries.csv"))
    
  # The test_data_5km data table provides a single coordinate for 38 communities.
    test_data_5km <- read.csv(paste0(fp, "/test_data_5km.csv"))
  
# Load Source Functions
  source(paste0(fp, "/urbanicity_source_functions.r"))

# Create Bounding Box for Each Test Dataset.
  # The bounding box defines the area that we are pulling data for.
  # For the test_data_boundaries example, we pre-specify the full boundary area with exact coordinates.
  # For the test_data_5km example, we define the boundary as 5km with in a single point.
  # The function defaults to 5km, but you can specify any distance by changing the distance_km value.
    
    bbox_boundaries <- create_bounding_boxes(test_data_boundaries)
    bbox_5km <- create_bounding_boxes(test_data_5km, distance_km = 5)
  
# Compute Urbanicity Measures
  # This function takes the bounding box (bbox) we created above as the main argument.
  # You also need to specify which metrics to calculate. The default is to calculate all metrics,
  # but you can choose not to compute specific metrics by changing TRUE to FALSE.
    
    # roads (percent of paved roads and ratio of paved-to-unpaved roads)
    # shops (number of formal shops / markets)
    # healthcare (distance to nearest healthcare facility)
    # transport (number of public transit stops)
    # financial (number of financial facilities [e.g., atms, banks])
    # schools (distance to nearest school)
    # cell_towers (number of cell towers)
    # buildings (building density [i.e., percent of area covered by buildings])
    # nighttime light (intensity of nighttime light)
    # population (population density)
    
  # Three additional arguments specify the location of the files used to compute the distance and nighttime light measures.
  # By default, if you do not supply a file path for the distance measures,
  # euclidean distance will be used to compute travel time instead of using the friction surface method.
    
    # friction_surface_path (file path to friction surface raster for computing travel time)
    # population_raster_path (file path to population raster for computing population density)
    # nighttime_light_path (file path to image for computing nighttime light intensity)
    
      test_results_boundaries <- compute_urbanicity_iterative(bbox_boundaries,
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
                                                              friction_surface_path = paste0(fp, "/friction_surface_walking.geotiff"),
                                                              population_raster_path = paste0(fp, "/pop_raster.tif"),
                                                              nighttime_light_path = paste0(fp, "/nighttime_lights.tif"))
      test_results_5km <- compute_urbanicity_iterative(bbox_5km,
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
                                                       friction_surface_path = paste0(fp, "/friction_surface_walking.geotiff"),
                                                       population_raster_path = paste0(fp, "/pop_raster.tif"),
                                                       nighttime_light_path = paste0(fp, "/nighttime_lights.tif"))
      
# Generate Summary Plots of the Results
  summary_plots_boundaries <- create_summary_plots(test_results_boundaries)
  summary_plots_5km <- create_summary_plots(test_results_5km)   
      
  # View Summary Plots for Continuous Measures
    summary_plots_5km$paved_to_unpaved_ratio
    summary_plots_5km$pct_paved_roads
    summary_plots_5km$n_shops
    summary_plots_5km$n_transport_stops
    summary_plots_5km$n_financial
    summary_plots_5km$building_density_pct
    summary_plots_5km$nighttime_light
    
      