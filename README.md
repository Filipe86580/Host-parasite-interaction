# Host–Parasite Interaction Networks between Nematodes and Non-Volant Small Mammals in Brazil

**A systematic review and ecological network analysis**

## Overview

This repository contains scripts and data for the network analysis of host–parasite interactions between nematodes and non-volant small mammals in Brazil. The scripts enable modularity, nestedness, specialization, compound topology, and species role analysis using bipartite network approaches, as described in the corresponding systematic review.

## Script adapted 
Script adapted from Mello et al., 2013; 2019 and Pinheiro et al., 2022

Mello, M.A.R., Bezerra, E.L.S., Machado, I.C., 2013. Functional Roles of Centridini Oil Bees and Malpighiaceae Oil Flowers in Biome-wide Pollination Networks. Biotropica 45, 45–53. https://doi.org/10.1111/j.1744-7429.2012.00899.x

Mello, M.A.R., Felix, G.M., Pinheiro, R.B.P., Muylaert, R.L., Geiselman, C., Santana, S.E., Tschapka, M., Lotfi, N., Rodrigues, F.A., Stevens, R.D., 2019. Insights into the assembly rules of a continent-wide multilayer network. Nat Ecol Evol 3, 1525–1532. https://doi.org/10.1038/s41559-019-1002-3

Pinheiro, R.B.P., Felix, G.M.F., Lewinsohn, T.M., 2022. Hierarchical compound topology uncovers complex structure of species interaction networks. Journal of Animal Ecology 91, 2248–2260. https://doi.org/10.1111/1365-2656.13806


## Directory Structure

```
.
├── Supplement 2.R          # Main analysis workflow: data import, preprocessing, topology analysis, plotting
├── Supplement 3.R          # Custom igraph vertex shape (diamond)
├── Supplement 4.R          # Posterior probability function for modular networks
├── Supplement 5.R          # Restricted (area/modular) null model function
├── data/
│   └── Pasta2.csv          # Main host-parasite interaction matrix
├── figures/                # Output figures (created by analysis scripts)
├── results/                # Output statistics and results (created by analysis scripts)
```

## Setup

1. **Requirements**
   - R (≥ 4.0 recommended)
   - R packages: `igraph`, `bipartite`, `Rmisc`, `vegan`, `gdata`, `ggplot2`, `gridExtra`, `grid`

2. **Preparation**
   - Place all `Supplement *.R` scripts in your working directory.
   - Ensure `data/Pasta2.csv` exists for the main analysis.
   - Create output folders: `figures/`, `results/` if not present.

3. **Execution**
   - Open `Supplement 2.R` in RStudio or your R environment.
   - Run all code blocks (it will source dependencies).
   - Outputs such as network figures, statistics, species roles, and compound topology plots are saved to `figures/` and `results/`.

## Script Descriptions

- **Supplement 2.R**:  
  - Loads data and dependencies.  
  - Plots networks and matrix visualizations.  
  - Computes modularity (DIRTLPAwb+), specialization (H2’), nestedness (NODF/WNODA), and compound topology.  
  - Runs null models and permutation tests.  
  - Assigns and plots species roles.

- **Supplement 3.R**:  
  - Defines a custom "diamond" shape for igraph node plotting.

- **Supplement 4.R**:  
  - Function `PosteriorProb()` for calculating conditional probabilities of interactions in modular/networks as part of the null model construction.

- **Supplement 5.R**:  
  - Function `RestNullModel()` for area/modular-restricted null model generation based on network structure.

## Results

- **Topological Metrics**
  - Modularity (DIRTLPAwb+): overall network modularity and module count.
  - Specialization (H2’): degree of host/parasite specialization.
  - Nestedness (NODF, WNODA): nestedness at full network, within/between modules.
  - Compound topology visualizations: sorted and colored matrices to highlight modular/nested structure.
  - Species roles: Identifies hubs, connectors, peripherals, and kinless species by z-c scores.

Results and plots are exported as `.csv`, `.txt`, and image files to their respective folders.

## Citation

If you use this code or data, please cite the original publication (when available), and acknowledge the Ecological Synthesis Lab (SintECO).

## Contact / Contributing

For questions, requests, or contributions, please contact the authors or open an issue in this repository.
