  <!-- badges: start -->
  [![R-CMD-check](https://github.com/Soheil-A-Neshat/MicroEcoTools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Soheil-A-Neshat/MicroEcoTools/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# MicroEcoTools
Theoretical Microbial Ecology Computational Tools
[![DOI](https://zenodo.org/badge/DOI/10.1111/2041-210X.70047.svg)](https://doi.org/10.1111/2041-210X.70047)


MicroEcoTools is an R package developed for microbial ecologists to apply ecological frameworks to microbial community data. This package helps analyze the effects of disturbances on microbial communities by assessing microbial diversity, community assembly mechanisms, and life-history strategies.

# Key Features
Null Model Analysis (NMA): Assesses the relative impacts of stochastic and deterministic assembly mechanisms in microbial communities.
Trait-Based Analysis: Applies Grime's CSR (Competitor, Stress-Tolerant, Ruderal) life-history strategy framework to microbial taxa and functional traits.
Diversity Analysis: Provides functions to calculate Hill numbers for diversity analysis.
Simulation Algorithms: Uses a multinomial distribution-based simulation algorithm to generate and analyze simulated communities.

## Installation
You can install the MicroEcoTools package directly from GitHub using the devtools package in R:

Install the devtools package if you haven't already:

install.packages("devtools")

Install MicroEcoTools from GitHub:

devtools::install_github("Soheil-A-Neshat/MicroEcoTools")

library(MicroEcoTools)

## Main Functions
Null Model Analysis (NMA):

NMA(): Performs null model analysis comparing observed communities to simulated communities under the null hypothesis.
Trait-Based Analysis:

CSR_assign(): Assigns competitor, stress-tolerant, and ruderal categories to traits based on experimental data.
CSR_Simulation(): Uses simulation to predict the accuracy of CSR assignments.
Statistical Testing:

Welch_ANOVA(): Conducts Welch's ANOVA, useful for datasets with unequal variances.
pairwise_welch(): Performs pairwise comparisons assuming non-constancy of variance.
Diversity Analysis:

Hill_diversity(): Calculates Hill numbers for diversity metrics.

## Example Datasets
The package includes two example datasets:

NMA_data: Genus-level count data from a study on activated sludge under varying disturbance frequencies.
CSR_data: Functional trait data categorized under the CSR framework.
These datasets can be loaded into your R session for practice and to explore the package's functionalities:

## Implementation Details
MicroEcoTools utilizes several R packages, including data.table, dplyr, ggplot2, vegan, and others, for data manipulation, statistical analysis, and visualization. The package also features extensive error-checking and reporting functionalities to ensure data quality and facilitate troubleshooting.

## References and Documentation
For a detailed description of the package's functions and additional examples, please refer to the official documentation available through R help function.
More details are available here: https://doi.org/10.1111/2041-210X.70047

## Limitations
Requires at least two biological replicates for analysis.
Simulation results may vary with fewer than 1,000 simulations.

## Unit Test
Unit tests was performed using testthat package and coverage assessed through covr:
  MicroEcoTools Coverage: 78.00%

## Citation
Neshat, S. A., Santillan, E., & Wuertz, S. (2025). MicroEcoTools: An R package for comprehensive theoretical microbial ecology analysis. Methods in Ecology and Evolution, 00, 1–9.
bioRxiv 2024.08.19.608598; doi: https://doi.org/10.1111/2041-210X.70047

