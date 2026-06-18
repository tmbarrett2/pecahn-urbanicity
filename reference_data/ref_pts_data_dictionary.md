# Data dictionary — `ref_pts.csv`

The drawn urbanicity **reference sample**: an equalized stratified random sample of
the global GHS-SMOD settlement gradient, with GHS-POP-based design weights. One row
per reference point.

- **Rows:** 1,400 (200 per settlement class × 7 classes)
- **Encoding:** UTF-8 (a BOM may precede the first header; read with `readr::read_csv()`
  or `read.csv(..., fileEncoding = "UTF-8-BOM")` to strip it)
- **CRS of `lon`/`lat`:** WGS84 (EPSG:4326)
- **Produced by:** `01_generate_reference_distribution.ipynb` (ArcGIS Pro / arcpy)

> Study-site coordinates are **not** part of this file — sites are measured on the same
> metrics downstream and compared against the reference dataset built from these points.

---

## Columns

| Column | Type | Units | Description |
|---|---|---|---|
| `Classified` | integer | GHS-SMOD code | Settlement class of the point's 1 km cell, taken from the classified (land-masked) GHS-SMOD raster. One of the seven codes in the mapping below. This is the stratification variable (200 points per code). |
| `GrndTruth` | integer | — | Ground-truth class field emitted by ArcGIS `CreateAccuracyAssessmentPoints`. It is a fixed placeholder of **`-1`** for every row (the tool is repurposed here only to draw the equalized stratified sample; no manual ground-truthing is performed). Carried through for provenance; **not used** downstream. |
| `pop` | double | persons per 1 km cell | GHS-POP E2025 population **count** sampled at the point's cell (`ExtractValuesToPoints`, `VALUE_ONLY`). Range in this draw: `0` – `59240.53`. `0` occurs at uninhabited cells (see note on `w`). |
| `w` | double | — (relative weight) | **Design weight**, `w = max(0, pop) × (N_h / n_h)`, where `N_h` is the land-masked cell count of class `h` and `n_h = 200` is the realized sample size of class `h`. Re-inflates each sampled point to the population its class represents across the globe, correcting the equalized (over-sampled rare classes) draw back to a population-weighted distribution. Range in this draw: `0` – `193103405.10`. |
| `lon` | double | decimal degrees | Longitude, WGS84 (`POINT_X`). Range: `-160.10` – `175.54`. |
| `lat` | double | decimal degrees | Latitude, WGS84 (`POINT_Y`). Range: `-41.92` – `80.68`. |

---

## GHS-SMOD class code → label

GHS Settlement Model (Degree of Urbanisation, level 2). Class **10 = Water** is
excluded before sampling, so it never appears in `Classified`.

| Code | GHS-SMOD label | Label used in `04_compute_urbanicity_reference.Rmd` (`project`) |
|---|---|---|
| `30` | Urban Centre grid cell | Urban Center |
| `23` | Dense Urban Cluster grid cell | Dense Urban Cluster |
| `22` | Semi-dense Urban Cluster grid cell | Semi-Dense Urban Cluster |
| `21` | Suburban / peri-urban grid cell | Suburban Grid Cell |
| `13` | Rural cluster grid cell | Rural Cluster |
| `12` | Low Density Rural grid cell | Low Density Rural Grid Cell |
| `11` | Very Low Density Rural grid cell | Very Low Density Grid Cell |
| `10` | Water grid cell | *(excluded from the sample)* |

These labels are assigned in `04_compute_urbanicity_reference.Rmd`
(`Classified` → `project`) and become the per-class batching key.

---

## Note on zero-weight points (`w = 0`)

**`w = 0` rows are valid and intentionally retained.** In this draw, **155** points
have `w = 0`, and **all 155 are class `11`** (Very Low Density Rural) cells whose
GHS-POP value is `0` — i.e. inhabitable land cells that nonetheless carry no
population in the E2025 epoch. Because the design weight is proportional to `pop`,
these points carry **zero design weight**, so any weighted summary of the reference
dataset gives them no mass. They are kept rather than dropped so that:

- the per-class sample sizes stay exactly `n_h = 200` (the `N_h / n_h` scale factor
  is computed against `n_h = 200`), and
- the file documents the full stratified draw, including the uninhabited tail of the
  very-low-density class.

Their urbanicity metrics are still computed by `04_compute_urbanicity_reference.Rmd`;
they simply have no influence on any weighted summary.

---

## Downstream join key

`04_compute_urbanicity_reference.Rmd` builds a `community` string as
`paste0(lat, ", ", lon)` (note the order: **lat, lon**) from this file via
`read.csv()`. `05_assemble_reference_dataset.R` reconstructs the identical
string the same way to join each point's design weight `w` back onto its computed
urbanicity metrics.
Because both sides format the same IEEE-754 doubles with R's default 15
significant digits, the keys match exactly — do not re-round `lon`/`lat`.
