# Urbanicity Index for Anthropological Fieldsites - Vignette
# Developed as Part of the Population Ecology, Aging, and Health Network (PEcAHN)

# To run the code below, you will need to install two packages (osmdata and sf) for pulling data from OpenStreetMap.
# Uncomment and run the following if you need to install the packages:
  # install.packages("osmdata")
  # install.packages("sf")

# You will also need the following two packages for the plotting function to work. Uncomment and run as needed.
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
  # but you can also replace "all" below with any combination of the following to calculate only those metrics:
    
    # roads (percent of paved roads and ratio of paved-to-unpaved roads)
    # shops (number of formal shops / markets)
    # healthcare (presence of healthcare facilities [yes/no])
    # transport (number of public transit stops)
    # financial (number of financial facilities [e.g., atms, banks])
    # schools (presence of schools [yes/no])
    # cell_towers (number of cell towers)
    # buildings (building density [i.e., percent of area covered by buildings])
    
      test_results_boundaries <- compute_urbanicity_iterative(bbox_boundaries, metrics = c("all"))
      test_results_5km <- compute_urbanicity_iterative(bbox_5km, metrics = c("all"))
      
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
      