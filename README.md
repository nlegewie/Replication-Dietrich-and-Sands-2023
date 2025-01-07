# Replication Dietrich and Sands 2023

## Authors

Nicolas Legewie, University of Erfurt  
Doron Shiffer-Sebba, Northwestern University  
Jannes Jacobsen, DeZIM Institute  
Yoav Goldstein, Independent Researcher  
Jörg Dollmann, DeZIM Institute and Mannheim Centre for European Social Research (MZES)  


## Purpose

This repository contains the code for the replication of Dietrich and Sands (2023) using data from a new context and new data collection technique (3DSR, see [Goldstein et al., 2023](https://doi.org/10.1177/00491241221147)), using 3D videos and computer vision models to estimate the effect of minoritized bystanders on pedestrian distance in Berlin, Germany.


## Abstract

Dietrich and Sands (2023) used New York City traffic camera footage to experimentally examine the effect of a pair of racialized confederate bystanders on the distance pedestrians maintained from those bystanders as they passed them on the street. Across their block-randomized experimental conditions, the authors found that pedestrians deviated by around 4 inches on average, maintaining larger distances from Black individuals as opposed to White individuals. Their point estimate was significant at the 5% level. In this conceptual replication we use data from a new context and new data collection technique, using 3D videos and computer vision models to estimate the effect of minoritized bystanders on pedestrian distance in Berlin, Germany. We directionally reproduce Dietrich and Sands’ main claim, finding that, writ large, pedestrians maintain a larger distance from Muslim bystanders of Middle Eastern descent wearing jalabiyas (a religious muslim garment), as opposed to White bystanders wearing jeans and t-shirts, in the German context. Using bootstrapped sub-samples from our data, our point estimates range from 0.39 inches (1 cm) to 9.84 inches (25 cm), with an average difference in distance of about 5.9 inches (15 cm). This roughly mirrors Dietrich and Sands’ finding of a 4 inch difference. However, our replication finds greater heterogeneity across locations, where different types of areas and different streets show opposite patterns. Our pooled results are only statistically significant at the 5% level when using a wild block bootstrap but are not significant when using clustered standard errors.


## Requirements

- R version [4.4.2 (2024-10-31)]
- Required packages (installed automatically via `pacman`):
  - `here` [1.0.1] - Project path management
  - `glue` [1.7.0] - String interpolation
  - `tidyverse` [2.0.0] - Data manipulation and visualization
  - `readxl` [1.4.3] - Excel file reading
  - `sp` [2.1-2] - Spatial data handling
  - `conflicted` [1.2.0] - Conflict resolution
  - `patchwork` [1.2.0] - Plot composition
  - `ggthemr` [1.1.0] - ggplot themes
  - `flextable` [0.9.4] - Table formatting
  - `broom` [1.0.5] - Model tidying
  - `sandwich` [3.1-0] - Robust standard errors
  - `lmtest` [0.9-40] - Testing linear regression models
  - `multiwayvcov` [1.2.3] - Multi-way clustering

Note: Package versions listed are the ones used in development. The code should work with newer versions as well.

## Data

### Input Data Requirements
The `data/raw/` directory should contain:
1. CSV files with trajectory data (format: `YYYYMMDD_HHMMSS_*.csv`)
2. `metacsv.xlsx` containing experimental conditions and block information

## Data Processing
The data processing script (`3DSR_racial_avoidance_data processing.R`):
1. Reads and cleans trajectory data
2. Creates analysis corridors
3. Calculates distances between participants and confederates
4. Generates validation plots
5. Creates random subsamples for analysis

The processed data (`3DSR_racial_avoidance_subsamples.rds`) contains:
- Minimum distances between participants and confederates
- Location information
- Experimental conditions
- Block and timing data


## Analysis Pipeline

1. **Analysis** (`3DSR_racial_avoidance_analysis.R`):
   - Analyzes distance distributions
   - Performs regression analyses with multiple specifications
   - Creates specification curves
   - Generates summary tables

Results are saved in the `outputs/figures` and `outputs/tables` directories.


## Usage

1. Clone this repository
2. Run scripts in order:
   ```R
   source("scripts/3DSR_racial_avoidance_data processing.R")
   source("scripts/3DSR_racial_avoidance_analysis.R")
   ```
3. If you want to test the pipeline without the existing processed data and outputs, you can delete the contents of the `data/processed/`, `outputs/figures/`, and `outputs/tables/` directories.


## Citation

[Add paper citation when available]

## Contact

For comments or questions, please contact: Nicolas Legewie at: nicolas (dot) legewie (at) uni-erfurt (dot) de
