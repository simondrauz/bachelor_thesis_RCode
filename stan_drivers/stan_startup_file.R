#This file is used to set overall specifications for the run

# LOAD NEEDED PACKAGES
library("rstan")
library('arrow')
library('dplyr')
library('parallel')

# Set Running variables
distribution_assumption = "zinb"
stan_models_dir_path = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/stan_models/"
# Set to "_enableParallel" if parallel computing should be activated
enable_parallel = ""
model = "model15"

feature_set = "feature_set1"
data_dir_path = "C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/data/"

stan_fit_dir = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/stan_fits/"

eval_months <- c(463)
lag_indicators <- c(9)
month_indicators <- c(7)
no_chains = 1
cores_for_par = 0
