# Analysis codes for Ushio et al. "Detecting and validating influential organisms for rice growth: An ecological network approach"

This repository contains analysis codes to reproduce the results in Ushio et al. (2023) bioRxiv https://doi.org/10.1101/xxxx.xx.xx.xxxxxx.

# License
See LICENSE.


# Analysis workflow
## Step 1. Detection of causal relationships between ecological community members and rice growth in 2017
- `01_2017_DNAxRice/`: Codes for Unified Information-theoretic Causality (UIC) analysis of quantitative environmental DNA data and rice growth in 2017.<br>

## Step 2. Summary of rice growth data in 2019
- `02_2019_Rice/`: Codes for summarizing and visualizing rice growth in 2019<br>

## Step 3. Sequence data processing of environmental DNA data in 2019
- `03_2019_eDNA/`: Codes for sequence data processing of quantitative environmental DNA data in 2019<br>

## Step 4. Detection of differentially expressed genes in rice in 2019
- `04_2019_riceRNA/`: Codes for RNA sequence data processing and detection of differentially expressed genes in rice leaves in 2019<br>

## Step 5. Format figures
- `FigCode/`: Format figures<br>
- `FigCode/Fig_Combine.Rmd`: Combine figures<br>


# Package versions
see text files in `00_SessionInfo/`


# Data availablity
Sequence data are deposited in DDBJ Sequence Read Archives (DRA). The accession numbers are as follows: DRA009658, DRA009659, DRA009660, and DRA009661 for eDNA data of ecological communities in 2017 (Ushio 2022 _Proceedings of the Royal Society B_ https://doi.org/10.1098/rspb.2021.2690), DRA015682, DRA015683, DRA015685, and DRA015686 for eDNA data of ecological communities in 2019, and DRAXXXXX for rice RNA expression data.
