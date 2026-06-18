# PECAHN Urbanicity Reference Dataset

Reproducibility materials for the **urbanicity reference dataset** — a population-weighted, GHS-SMOD-stratified sample of 1,400 global points measured on the same urbanicity metrics as PEcAHN study sites. This folder documents how the sample is drawn, how its urbanicity metrics are computed, and how the per-point dataset is assembled.

The reference sample is a population-weighted, GHS-SMOD-stratified set of global locations carrying a design weight per point. Study-site coordinates are *not* part of the reference set — sites are measured on the same metrics downstream and compared against this dataset.

---

## Pipeline overview

```
GHS-SMOD + GHS-POP                     OSM planet (2026-06-08)
        │                                       │
        ▼                                       ▼
1. 01_generate_reference_distribution.ipynb   3. self-hosted Overpass backend
   → ref_pts.csv (1,400 pts)                     (osmium extract + Docker build)
        │                                       │
        ▼                                       │
2. 02_generate_local_buffers.R                  │
   → local_buffers.geojson  ────────────────────┤ (clips the local OSM layers)
        │                                       │
        └───────────────┬───────────────────────┘
                        ▼
        4. 04_compute_urbanicity_reference.Rmd
           → per-class urbanicity results (.rds)
                        │
                        ▼
        5. 05_assemble_reference_dataset.R
           → urbanicity_reference_dataset.rds (per-point metrics + weights)
```

*(The scripts are numbered in the order they run. The OSM backend is the one
step with no numbered script — it's a Docker / osmium build rather than an R or
Python file — so the prefixes jump from `02_` to `04_`; its artifacts live in
[`osm_backend/`](osm_backend/).)*

---

## Repository contents

| File | Runs | Purpose |
|---|---|---|
| `01_generate_reference_distribution.ipynb` | 1st | ArcGIS Pro / arcpy script. Draws the equalized stratified random sample from GHS-SMOD, masks to inhabited land via GHS-POP, extracts population, computes design weights, exports `ref_pts.csv`. |
| `02_generate_local_buffers.R` | 2nd | Buffers the reference points (3 km, equal-area) and dissolves them into `local_buffers.geojson`, used to clip the local OSM layers when building the self-hosted Overpass extract. |
| `osm_backend/` | 3rd | Self-hosted Overpass build (no numbered script): `docker-compose.yml`, `run_phase3_rest.ps1`, and `BUILD.md` (osmium command list, planet date+MD5, pinned image digest). |
| `04_compute_urbanicity_reference.Rmd` | 4th | Runs `pecahnurbanicity::compute_urbanicity_iterative()` (all metrics) over the 1,400 points in per-class batches, querying the local Overpass instance + Google Earth Engine, saving intermediate `.rds` results. |
| `05_assemble_reference_dataset.R` | 5th | Binds the per-class `.rds` batches, joins each point's design weight and SMOD class, and writes the assembled per-point reference dataset (`urbanicity_reference_dataset.rds`). |
| `ref_pts.csv` | — | The actual drawn sample (`Classified, GrndTruth, pop, w, lon, lat`). Shipped directly so the dataset does not depend on re-running the stochastic ArcGIS draw. See [`ref_pts_data_dictionary.md`](ref_pts_data_dictionary.md). |
| `ref_pts_data_dictionary.md` | — | Column-by-column description of `ref_pts.csv`, the SMOD class code→label mapping, and the zero-weight-point convention. |
| `environment.md` / `sessionInfo.txt` | — | Pinned software environment (R, ArcGIS, Docker/osmium) and the full R `sessionInfo()`. |
| `LICENSE` | — | MIT license for the code in this folder. |
| `urbanicity_reference_dataset.rds` | produced | Output of step 5: one row per reference point with its design weight, SMOD class, and all computed urbanicity metrics. Written once all seven classes are present. |

The assembled dataset is a tidy table — one row per reference point, carrying the
design weight, SMOD class, and every urbanicity metric computed in step 4. How to
summarise or compare study sites against it is left to the downstream analysis.

---

## Data sources & versions

| Input | Product / version | Notes |
|---|---|---|
| Settlement classes (sampling) | **GHS-SMOD** E2025, R2023A, Mollweide (ESRI:54009), 1 km, V2.0 (`GHS_SMOD_E2025_GLOBE_R2023A_54009_1000_V2_0.tif`) | Seven settlement classes (11/12/13/21/22/23/30); water (10) excluded. |
| Population (sampling + weights) | **GHS-POP** E2025, R2023A, Mollweide, 1 km, V1.0 (`GHS_POP_E2025_GLOBE_R2023A_54009_1000_V1_0.tif`) | Used both to mask to inhabited land and to weight the draw. |
| Road / facility / building features | **OpenStreetMap planet snapshot 2026-06-08** | Served from a self-hosted Overpass instance built from a tag-filtered + buffer-clipped extract (see "OSM backend"). |
| Travel-time friction | **MAP Global Friction Surface 2019** (walking-only), via Google Earth Engine | Used for facility travel-time metrics. |
| Population & nighttime-light rasters | via Google Earth Engine, years 2015–2025 | See the GEE asset IDs below. |

### Google Earth Engine asset IDs (used by `04_compute_urbanicity_reference.Rmd`)

These are the `pecahnurbanicity` 0.3.0 defaults (`R/ee_extract_rasters.R`),
overridable via the `*_asset` arguments on `compute_urbanicity_iterative()`:

| Role | Asset ID | Band | Scale |
|---|---|---|---|
| Friction surface | `projects/malariaatlasproject/assets/accessibility/friction_surface/2019_v5_1_walking_only` | `friction` | ~927.67 m |
| Population (WorldPop next-gen) | `projects/sat-io/open-datasets/WORLDPOP/pop` | `population` | native (~92.77 m fallback) |
| Nighttime lights (VIIRS DNB) | `NOAA/VIIRS/DNB/ANNUAL_V22` (fallback `NOAA/VIIRS/DNB/ANNUAL_V21`) | `average_masked` | ~463.83 m |
| Urban-center search (GHS-SMOD) | `JRC/GHSL/P2023A/GHS_SMOD_V2-0` (per-epoch `.../{YEAR}`) | `smod_code` | 1000 m |

---

## Reproducibility-critical parameters

**Sampling — `01_generate_reference_distribution.ipynb`**
- 200 points per class × 7 classes = **1,400 points**, `EQUALIZED_STRATIFIED_RANDOM`
- Minimum spacing **5,000 m** (reduces 1 km-cell autocorrelation)
- Random seed **`42 MERSENNE_TWISTER`**; output CRS World Mollweide (54009); snap raster = SMOD
- Class cell counts `N_h` taken from the **land-masked** raster
- Design weight: `w = pop × (N_h / n_h)`, clamped to ≥ 0

**Buffers — `02_generate_local_buffers.R`**
- **3 km** buffer in World Mollweide (covers the 1 km analysis box + projection-distortion margin), dissolved, written as WGS84 GeoJSON

**Urbanicity run — `04_compute_urbanicity_reference.Rmd`**
- `metrics = "all"`, `search_buffer = 1`, `batch_size = 5`
- Local bbox **1 km**, regional bbox **5 km** (`create_bounding_boxes` distances 1 and 5)
- Population & nighttime-light years: **2015–2025**
- All 1,400 points run against the **single** OSM snapshot for consistency

**Assembly — `05_assemble_reference_dataset.R`**
- Binds the per-class batch results and joins each point's design weight `w` and
  SMOD class back on via the `"lat, lon"` `community` key
- Writes only when all 1,400 points / seven classes are present (self-skips otherwise)
- `w = 0` points (155 uninhabited class-11 cells) are retained — they carry zero
  design weight, so any weighted summary gives them no mass, but their metrics remain

---

## OSM backend (self-hosted Overpass)

Because the run issues thousands of progressive feature queries, it uses a **local Overpass instance** rather than the public API (which rate-limits) — this also pins every query to one reproducible OSM snapshot.

The extract uses a **split-layer** strategy: roads and key amenities are kept worldwide (so progressive searches never clip), while the heavy local 1 km layers (buildings, shops, transport, financial, cell) are clipped to the point buffers before filtering. Build steps:

1. `osmium tags-filter` the planet → global layers (highways + amenity hospital/school/bank/atm).
2. `osmium extract -s smart -p local_buffers.geojson` → clip; then `tags-filter` → local layers.
3. `osmium merge` global + local → `planet-pecahn.osm.pbf`.
4. Build a `wiktorn/overpass-api` Docker instance (init mode) from the merged extract.

Served at `http://localhost:12345/api/interpreter`. Pinned tooling:

- OSM planet snapshot **2026-06-08**, MD5 `eb419b0f06089ef97e4b49f937ebf4a8`
- Docker image `wiktorn/overpass-api@sha256:04d4e189ba98aa42619414b8dab2599e852d3e7add9b84001b6de1f22319c01c`
- `osmium-tool 1.15.0` (libosmium 2.19.0)

Full command list and the `docker-compose.yml` are in [`osm_backend/`](osm_backend/).

---

## How to reproduce

1. **Sample** — open `01_generate_reference_distribution.ipynb` in the ArcGIS Pro Python window with the two GHS GeoTIFFs in the project folder; run to produce `ref_pts.csv`. *(Or skip and use the shipped CSV.)*
2. **Buffers** — set `INPUT_CSV` in `02_generate_local_buffers.R`, run to produce `local_buffers.geojson`.
3. **OSM backend** — build the self-hosted Overpass extract + instance (see `osm_backend/BUILD.md`); start it and confirm a test query returns features.
4. **Urbanicity** — in `04_compute_urbanicity_reference.Rmd`, set `osmdata::set_overpass_url("http://localhost:12345/api/interpreter")`, authenticate Earth Engine, run all classes (saves per-class `.rds` batches).
5. **Assemble** — point `RESULTS_DIR` in `05_assemble_reference_dataset.R` at the `.rds` batches and run it to bind them, join the design weights, and write `urbanicity_reference_dataset.rds`. The script self-skips the write until all seven classes are present.

---

## Software environment

See [`environment.md`](environment.md) for the full pinned environment; in brief:

- **ArcGIS Pro 3.7**, Python **3.13.13** (`arcgispro-py3` conda env), Spatial Analyst, `arcpy` 3.7
- **R 4.4.3** with `pecahnurbanicity` **0.3.0** (repo commit `99096fd`), `sf`, `osmdata`, `rgee`, `reticulate`, `tidyverse` — full versions in [`sessionInfo.txt`](sessionInfo.txt)
- **Docker Desktop** (Engine 29.5.3), `wiktorn/overpass-api` (digest above), `osmium-tool 1.15.0`

---

## Attribution & licensing

- **OpenStreetMap** data © OpenStreetMap contributors, licensed under the **ODbL**. Attribution required for any derived products.
- **GHSL** products (GHS-SMOD, GHS-POP) © European Union, European Commission **Joint Research Centre** — free to use with attribution.
- **MAP Global Friction Surface 2019** — cite the Malaria Atlas Project.
- Code in this folder is released under the MIT license — see [`LICENSE`](LICENSE).
