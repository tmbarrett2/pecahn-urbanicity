## Urbanicity Index for Anthropological Field Sites
### Developed as part of the [Population Ecology, Aging, and Health Network (PEcAHN)](https://sites.duke.edu/pecahn/)

This repository contains the **`pecahnurbanicity`** R package, which provides tools to compute measures of urbanicity from publicly available geospatial data.

### Key Components
- **`pecahnurbanicity` package** – Core R functions for computing and summarizing urbanicity metrics.
- **`urbanicity_vignette.html`** – A step-by-step tutorial demonstrating how to install and use the package with example data.

**View the vignette:**  
[Urbanicity Index Vignette (HTML)](https://tmbarrett2.github.io/pecahn-urbanicity/urbanicity_vignette.html)

---

### Installation
You can install the package directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("tmbarrett2/pecahn-urbanicity/pecahnurbanicity")
```

### Data sources & Google Earth Engine setup
Raster data (friction surface, WorldPop population, VIIRS nighttime lights, and GHSL settlement
classification) are pulled directly from **Google Earth Engine** at runtime via the
[`rgee`](https://r-spatial.github.io/rgee/) package.
OpenStreetMap data (roads, facilities, buildings) are queried live via the Overpass API and cached locally.

You will need a Google Earth Engine account and a Google Cloud Project with the Earth Engine API enabled.
One-time setup:

```r
install.packages("rgee")
library(rgee)
ee_install()                                    # installs the EE Python API (run once)
ee_Initialize(project = "your-gcp-project-id")  # authenticate
```

Then pass your project to the package via the `ee_project` argument and request data by year, e.g.:

```r
library(pecahnurbanicity)
results <- compute_urbanicity_iterative(
  local_bboxes          = bbox_1km,
  regional_bboxes       = bbox_5km,
  metrics               = "all",
  ee_project            = "your-gcp-project-id",
  population_years      = c(2015, 2020),
  nighttime_light_years = c(2015, 2020)
)
```

See the [vignette](https://tmbarrett2.github.io/pecahn-urbanicity/urbanicity_vignette.html) for a full walkthrough.
