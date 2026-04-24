# UK-Biobank-WMH-and-Structural-Covariance-Project
This repository contains R code used in Guocheng Jiang's first-author manuscript entitled **"Common vascular white matter lesions contribute to specific grey matter network disruptions"** If you are interested in applying or adapting any of the following scripts for your research, please cite the published manuscript once it becomes available.

# Project 1: Structural Covariance Network (SCN) simulator
**Related script: 1_SCN_Simulator.R**

In the first simulation, 3 × 3 SCN was built and step-functions were plotted for visualization. This simulation provides conceptual illustration of the progressive thresholding workflow using two simplified SCNs. Example outputs:
1. Correlation matrix plots for two 3 × 3 SCNs.
2. A step function plot showing the percentage of surviving connections across threshold values.

In the second simulation, we assessed of minimum effect size required for SCN connectivity comparisons using log rank tests and Cox-proportional hazards model. Simulations show the range of standardized effect size (Cohen’s d) required to achieve X% statistical power for detecting global differences in the network connectivity at different size of SCNs. Example output:
1. A small/medium/large effect size (Cohen's d) is required to achieve X% statistical power to detect a significant group difference in SCN connectivity given a network size of N nodes.

In the third simulation, we examined the effect of data stratification on the estimation of age-related differences in SCN connectivity. Example output:
1. At N folds, Cohen’s d estimates remained relatively consistent in direction across folds.
2. At 2N folds, simulation results showed substantially increased dispersion of effect size estimate. 

# Project 2:  Structural covariance network analysis: Comparisons of global and local backbone network connectivity.
**Related script: 2_SCN_Connectivity_Analysis.R**

The script compares grey matter SCN organization between two groups groups, including:
1. Construction of SCN adjacency matrices from regional brain volume data.
2. Global SCN connectivity comparison using threshold-survival curves and Cox proportional hazards models.
3. Identification of stable backbone network connections.
4. Visualization of preserved and WMH-vulnerable backbone edges using network and circular chord diagrams.

# Project 3:  Longitudinal associations between brain white matter lesions and grey matter volumetic changes.
**Related script: 3_Longitudinal_BLCS_Analysis**

The script uses bivariate latent change score models (BLCS) to examine longitudinal relationships between white matter hyperintensity (WMH) accumulation and regional grey matter volume change.
The example in this script analyze and tests for the following associations:
1. Whether baseline WMH volume predicts faster regional grey matter volume change,
2. Whether baseline regional grey matter volume predicts faster WMH volume change,
3. Whether baseline WMH volume and regional grey matter volume show significant covariance across the cohort,
4. Whether WMH volume change and regional grey matter volume change show significant covariance across the cohort.

# Dependency requirements:
The scripts were built based on R (Version 4.5.0). Required R code libraries were summarized within each R script.
