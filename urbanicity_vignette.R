# Urbanicity Index for Anthropological Fieldsites - Vignette
# Developed as Part of the Population Ecology, Aging, and Health Network (PEcAHN)

# This vignette provides step-by-step instructions to compute urbanicity measures for anthropological fieldsites.
# Several public data sources are required to compute the urbanicity measures. They include:
  
  # 1. Friction surface data for calculating the travel time measures. For this vignette, we are using walking time raster
  # files from the Malaria Atlas Project. The file is available for download here: https://malariaatlas.org/project-resources/accessibility-to-healthcare/
  # And the associated paper is available here: https://doi.org/10.1038/s41591-020-1059-1

  # 2. The data for estimating population density are from the Gridded Population of the World NASA dataset available here:
  # https://www.earthdata.nasa.gov/data/projects/gpw

  # 3. The data for calculating light intensity are from the Earth at Night (Black Marble) 2016 Color Maps NASA dataset, which can be accessed here:
  # https://www.visibleearth.nasa.gov/images/144898/earth-at-night-black-marble-2016-color-maps/144944l

# All PEcAHN members can access these data files in the following Google Drive folder: https://drive.google.com/drive/u/0/folders/1C6RWI7yhAM9UDjwMc7w9nxW30W5u3wee
# Please download the files to your folder with the source code before completing the vignette.

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
    
  # The test_data_5km data table provides a single coordinate for 16 communities.
    test_data_5km <- read.csv(paste0(fp, "/test_data_5km.csv"))
  
# Load Source Functions
  source(paste0(fp, "/urbanicity_source_functions.r"))

# Create Bounding Box for Each Test Dataset.
  # The bounding box defines the area that we are pulling data for.
  # For the test_data_boundaries example, we pre-specify the full boundary area with exact coordinates.
  # For the test_data_5km example, we define the boundary as a 5km radius around a single point.
  # The function defaults to 5km, but you can specify any distance by changing the distance_km value.
    
    bbox_boundaries <- create_bounding_boxes(test_data_boundaries)
    bbox_5km <- create_bounding_boxes(test_data_5km, distance_km = 5)
  
# Compute Urbanicity Measures
  # This function takes the bounding box (bbox) we created above as the main argument.
  # You also need to specify which metrics to calculate. The default is to calculate all metrics
  # by using metrics = c("all"), but you can also specify specific metrics as shown below.
    
    # roads (percent of paved roads and ratio of paved-to-unpaved roads)
    # shops (number of formal shops / markets)
    # healthcare (travel time (walking) to nearest healthcare facility)
    # transport (number of public transit stops)
    # financial (number of financial facilities [e.g., atms, banks])
    # schools (travel time (walking) to nearest school)
    # cell_towers (number of cell towers)
    # buildings (building density [i.e., percent of area covered by buildings])
    # nighttime light (intensity of nighttime light)
    # population (population density)
    
  # The search_buffer argument sets the area used for the travel distance measures for each community.
  # the default is 1 degree (100 km).
    
  # Three additional arguments specify the location of the files used to compute the distance and nighttime light measures.
    # friction_surface_path (file path to friction surface raster for computing travel time)
    # population_raster_path (file path to population raster for computing population density)
    # nighttime_light_path (file path to image for computing nighttime light intensity)
    
    # Note: I will eventually simplify the output in the console but for now it helps with debugging.
    
        test_results_boundaries <- compute_urbanicity_iterative(
          bbox_boundaries,
          metrics = c("roads", "shops", "healthcare", "transport", "financial", "schools",
                      "cell_towers", "buildings", "nighttime_light", "population"),
          search_buffer = 1,
          friction_surface_path = paste0(fp, "/friction_surface_walking.geotiff"),
          population_raster_path = paste0(fp, "/pop_raster.tif"),
          nighttime_light_path = paste0(fp, "/nighttime_lights.tif")
        )
        
        test_results_5km <- compute_urbanicity_iterative(
          bbox_5km,
          metrics = c("roads", "shops", "healthcare", "transport", "financial", "schools",
                      "cell_towers", "buildings", "nighttime_light", "population"),
          search_buffer = 1,
          friction_surface_path = paste0(fp, "/friction_surface_walking.geotiff"),
          population_raster_path = paste0(fp, "/pop_raster.tif"),
          nighttime_light_path = paste0(fp, "/nighttime_lights.tif")
        )
      
# Generate Summary Plots of the Results
  summary_plots_boundaries <- create_summary_plots(test_results_boundaries)
  summary_plots_5km <- create_summary_plots(test_results_5km)   
      
  # View Summary Plots
    summary_plots_5km$paved_to_unpaved_ratio
    summary_plots_5km$pct_paved_roads
    summary_plots_5km$travel_time_paved_road_min
    summary_plots_5km$n_shops
    summary_plots_5km$n_transport_stops
    summary_plots_5km$travel_time_school_min
    summary_plots_5km$travel_time_healthcare_min
    summary_plots_5km$n_cell_towers
    summary_plots_5km$n_financial
    summary_plots_5km$building_density_pct
    summary_plots_5km$pop_density
    summary_plots_5km$nighttime_light
    
      