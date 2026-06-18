# Pinned software environment

Captured 2026-06-16 on the build machine (Windows 11, build 26200).

## R side (buffers, urbanicity run, dataset assembly)

R **4.4.3** (2025-02-28 ucrt), x86_64-w64-mingw32. Full package versions and
`sessionInfo()` are in [`sessionInfo.txt`](sessionInfo.txt). Key packages:

| Package | Version |
|---|---|
| pecahnurbanicity | 0.3.0 |
| sf | 1.1.1 |
| osmdata | 0.3.0 |
| rgee | 1.1.8 |
| reticulate | 1.46.0 |
| raster | 3.6.32 |
| gdistance | 1.6.5 |
| digest | 0.6.39 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| readr | 2.2.0 |
| tidyverse | 2.0.0 |
| usethis | 3.2.1 |

The `rgee` Python backend runs through `reticulate` against the conda env
`C:/Users/tmbar/AppData/Local/r-miniconda/envs/rgee/python.exe`
(see `04_compute_urbanicity_reference.Rmd`).

## ArcGIS side (sampling — `01_generate_reference_distribution.ipynb`)

| Component | Version |
|---|---|
| ArcGIS Pro | 3.7 |
| ArcGIS Pro Python (`arcgispro-py3` conda env) | 3.13.13 |
| `arcpy` | 3.7 |
| Spatial Analyst extension | required (checked out in the notebook) |

## OSM backend

| Component | Version / digest |
|---|---|
| Docker Engine | 29.5.3 (Docker Desktop) |
| Overpass image | `wiktorn/overpass-api@sha256:04d4e189ba98aa42619414b8dab2599e852d3e7add9b84001b6de1f22319c01c` |
| osmium-tool | 1.15.0 (libosmium 2.19.0) |
| OSM planet snapshot | 2026-06-08, MD5 `eb419b0f06089ef97e4b49f937ebf4a8` |

See [`osm_backend/BUILD.md`](osm_backend/BUILD.md) for the full build sequence.
