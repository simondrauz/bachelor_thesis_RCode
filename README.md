
# Bachelor Thesis R Code Repository

## Overview
This repository contains R and Stan code for modeling Bayesian Penalized Structured Regression (BPSR) models, as part of the Bayesian analysis conducted in the bachelor thesis "Bayesian Penalized Structured Regression for Probabilistic Forecasting of Conflict Fatalities". The code primarily focuses on implementing and analyzing these models using `rstan`, a package for running Stan models in R.

### Key Files and Directories
- **Dokumentation.xlsx**: Documentation related to the project (altough not thorough).
- **assessing_convergence.R**: Script for assessing the convergence and precision of the Bayesian models.
- **country_list.csv**: CSV file containing a list of countries and their identifiers.
- **crps_extraction.R**: Scripts for extracting CRPS (Continuous Ranked Probability Score) metric.
- **distribution_plots.R**: Script for creating different distribution plots.
- **fit_inspection.R**: Script for inspecting the fit of the models.
- **pit.r**: Script for caluclating and plotting the Probability Integral Transform (PIT) for count data distributions.
- **stan_drivers**: Directory containing R scripts for driving Stan models.
- **stan_models**: Contains Stan model files.
- **verfication_of_MRF_reparametrization_matrices.R**: Script for verifying Markov Random Field (MRF) reparametrization matrices.

## Installation
To use the code in this repository, you will need R and the `rstan` package. Addional depencies should be available for inspection in the 'renv.lock' file.

## Usage
- The Stan models in `stan_models` and can be run using `rstan`.
- R scripts are used for loading input data, running models, analyzing model outputs, and performing statistical checks.

## Additional Information
- The repository is part of a larger project that includes Python code for data preprocessing and analysis as well as the data used as input in these Stan models, located in a separate [repository](https://github.com/simondrauz/bachelor_thesis).
