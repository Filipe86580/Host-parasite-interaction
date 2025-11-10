# Host–Parasite Interaction Networks between Nematodes and Non-Volant Small Mammals in Brazil

**A systematic review and ecological network analysis**

## Overview

This repository contains scripts and data for the network analysis of host–parasite interactions between nematodes and non-volant small mammals in Brazil. The scripts enable modularity, nestedness, specialization, compound topology, and species role analysis using bipartite network approaches, as described in the corresponding systematic review.

**Authors**: Gabriel M. Felix, Rafael B. P. Pinheiro, Marco A. R. Mello  
**Lab**: Ecological Synthesis Lab (SintECO)

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
