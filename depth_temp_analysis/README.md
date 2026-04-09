# Depth vs temperature (segmented regression)

In **OffspringSizeRmax** this folder lives under **`Rscripts/depth_temp_analysis/`**. The same analysis is mirrored at the **repository root** as **`depth_temp_analysis/`** on **EBarrowclift/batoid-rmax-offspring-scaling** (paths inside the Rmd are portable).

`04_DepthTempAnalysis.Rmd` fits segmented regressions of median depth on median temperature for rays vs skates, builds diagnostic and publication-style figures, and compares 1- vs 2-breakpoint models.

## Dependencies (R packages)

| Package    | Role |
|------------|------|
| **here**   | Project-relative paths; knitr `root.dir` for knitting from subfolder. |
| **knitr**  | (bundled with R Markdown) used via `opts_knit$set(root.dir = …)`. |
| **segmented** | Breakpoint regression on `lm()` fits. |
| **dplyr** / **tibble** | Data manipulation; `tibble()` for panel labels. |
| **ape**    | `read.tree()` for the Stein et al. phylogeny (loaded; not yet used in downstream chunks). |
| **broom**  | `augment()` for fitted values. |
| **ggplot2** | All plots; viridis plasma scale for `rmax`. |
| **Cairo**  | Cairo graphics stack (PDF output uses `grDevices::cairo_pdf`). |

Unused in the current script but harmless if omitted later: none required beyond the table above.

## Data inputs (first match wins)

| File | Typical locations |
|------|-------------------|
| **stein-et-al-single.tree** | `data/stein-et-al-single.tree` or `data/processed/stein-et-al-single.tree` |
| **Batoid offspring / model table** | `data/batoid_model_offspring.csv`, or `data/processed/batoid_model_offspring.csv`, or `data/processed/Batoid_model_offspring.csv` |

Column subsets are hard-coded by index in the Rmd; the CSV must match the expected **OffspringSizeRmax** / **batoid** model table layout.

## Outputs

- **`Depth_Temp_Skates_Rays_singlecol.png`** (600 dpi) and **`.pdf`**  
- Written under **`figures/`** (this repo) or **`depth_temp_analysis/figures/`** when that directory exists (standalone copy on **batoid-rmax-offspring-scaling**).
