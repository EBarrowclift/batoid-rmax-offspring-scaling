**Data and analysis for 'Offspring size resolves a population growth paradox in rays and skates'**

------------------------------------------------------

This repository contains the data (single phylogenetic tree from Stein et al., 2018 and batoid model data frame) and analysis for model fitting for the paper 'Offspring size resolves a population growth paradox in rays and skates'. The preprint is available at https://doi.org/10.1101/2024.01.02.573919, and the paper is accepted for publication in *Fish and Fisheries*. The study investigated the scaling of offspring body mass with the maximum intrinsic rate of population increase (rmax) as well as adult body mass, temperature, and depth for 85 ray and skate species (Superorder Batoidea).

The repository contains a single script ```batoid-rmax-offspring-scaling.R``` required to fit phylogenetic generalised linear models with data available in the `/data` directory.

R script edited from https://github.com/EBarrowclift/batoid-rmax-scaling created for the paper Barrowclift et al.,, 2023. 'Tropical rays are intrinsically more sensitive to overfishing than the temperate skates'. https://doi.org/10.1016/j.biocon.2023.110003.

### Species richness mapping

The **`mapping/`** directory holds an R Markdown pipeline (hex-based global richness maps for skates, rays, and sharks, plus publication figures). Open the repository root (the `.Rproj` file) before rendering. See **`mapping/README.md`** for run order, dependencies, and **`data/processed/`** files.

The spatial caches under **`data/processed/`** (`.rds` files) are **too large to host on GitHub**. You can **request a copy** by emailing **dulvy@sfu.ca**.

### Depth vs temperature (segmented regression)

The **`depth_temp_analysis/`** directory contains **`04_DepthTempAnalysis.Rmd`**, which uses **`data/stein-et-al-single.tree`** and **`data/batoid_model_offspring.csv`** and writes figures to **`depth_temp_analysis/figures/`**. See **`depth_temp_analysis/README.md`** for R package dependencies and input details.

