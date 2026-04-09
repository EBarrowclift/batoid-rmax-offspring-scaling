**Data and analysis for 'Offspring size resolves a population growth paradox in rays and skates'**

------------------------------------------------------

This repository contains the data (single phylogenetic tree from Stein et al., 2018 and batoid model data frame) and analysis for model fitting for the paper 'Offspring size resolves a population growth paradox in rays and skates', published in [insert publication info]. The preprint is available at https://doi.org/10.1101/2024.01.02.573919. The study investigated the scaling of offspring body mass with the maximum intrinsic rate of population increase (rmax) as well as adult body mass, temperature, and depth for 85 ray and skate species (Superorder Batoidea).

The repository contains a single script ```batoid-rmax-offspring-scaling.R``` required to fit phylogenetic generalised linear models with data available in the `/data` directory.

R script edited from https://github.com/EBarrowclift/batoid-rmax-scaling created for the paper Barrowclift et al.,, 2023. 'Tropical rays are intrinsically more sensitive to overfishing than the temperate skates'. https://doi.org/10.1016/j.biocon.2023.110003.

### Species richness mapping

The **`mapping/`** directory holds an R Markdown pipeline (hex-based global richness maps for skates, rays, and sharks, plus publication figures). Open the repository root (the `.Rproj` file) before rendering. See **`mapping/README.md`** for run order, dependencies, and required **`data/processed/`** files.

