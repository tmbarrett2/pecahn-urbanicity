# Urbanicity Index for Anthropological Fieldsites - Vignette
# Developed as Part of the Population Ecology, Aging, and Health Network (PEcAHN)

# Set Your File Path
# This is the path to the folder where you downloaded the test data and source functions.
# For example, mine is: "C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn/test_data"
  fp = "C:/Users/tmbar/OneDrive - Duke University/Active Projects/pecahn/test_data"
  #fp = "your_file_path"

# Load Test Data
  # The test_data_boundaries example provides northeast, southeast, southwest, and northwest coordinates for seven communities
    test_data_boundaries <- read.csv(paste0(fp, "/test_data_boundaries.csv"))
    
  # The test_data_5km example provides a single coordinate for 67 communities.
    test_data_5km <- read.csv(paste0(fp, "/test_data_5km.csv"))
  
# Load Source Functions
  source(paste0(fp, "/urbanicity_source_functions.r"))

# Create Bounding Box for Each Test Dataset.
# The bounding box defines the area that we are pulling data for.
  bbox_boundaries <- create_bounding_boxes(test_data_boundaries)
  bbox_5km <- create_bounding_boxes(test_data_5km)
  
# Compute urbanicity measures...
  test_results_boundaries <- compute_urbanicity_iterative(bbox_boundaries, metrics = c("all"))
  test_results_5km <- compute_urbanicity_iterative(bbox_5km, metrics = c("all"))