# Load packages
library("splines")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library('arrow')
library('dplyr')

# Load data
df_raw <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/data_cm_features_allyears_ext.parquet")
df_raw <- df_raw_expanded_feature_space_std

df_raw <- df_raw[df_raw$ged_sb < 1000,]
df_raw <- df_raw[df_raw$ged_sb > 0,]

df_raw <- df_raw[df_raw$country_id %in% developed_countries_ids, ]
df_raw <- df_raw[df_raw$country_id %in% countries_in_transition_ids, ]
df_raw <- df_raw[df_raw$country_id %in% developing_countries_extended_without_LDC_ids, ]
df_raw <- df_raw[df_raw$country_id %in% least_developed_countries_ids, ]

# Subset the dataframe for the specified columns
df <- df_raw[-c(1:3), c('month_id', 'ged_sb', 'ged_sb_tlag_3', 'decay_ged_sb_5_tlag_3', 'decay_ged_sb_100_tlag_3', 'decay_ged_sb_500_tlag_3', 'wdi_sp_pop_totl_tlag_3')]
df <- df_raw[-c(1:3),]

# Remove Nan values from lags of countries which don't have obeservations for the first months and therefore missing lags at a later point
# Remove NaN values from 'ged_sb_tlag_3' column
df <- df[!is.na(df$ged_sb_tlag_3), ]

# Define training_data and evaluation_data
# Take training data till Oct. 2019
df_train <- df[df$month_id <= 478, ]
#Take Jan. 2020 as evaluation month
df_eval <- df[df$month_id == 481, ]

dataX1 <- df_train$ged_sb_tlag_3
dataX2 <- df_train$decay_ged_sb_5_tlag_3
dataX3 <- df_train$decay_ged_sb_100_tlag_3
dataX4 <- df_train$decay_ged_sb_500_tlag_3
dataX5 <- df_train$wdi_sp_pop_totl_tlag_3
spatial_data <- df_train %>% select(starts_with("country_"))
temporal_data <- df_train %>% select(Winter, Spring, Summer, Fall)
dataY <- df_train$ged_sb
no_data <- nrow(df_train)

dataX1_eval <- df_eval$ged_sb_tlag_3
dataX2_eval <- df_eval$decay_ged_sb_5_tlag_3
dataX3_eval <- df_eval$decay_ged_sb_100_tlag_3
dataX4_eval <- df_eval$decay_ged_sb_500_tlag_3
dataX5_eval <- df_eval$wdi_sp_pop_totl_tlag_3
spatial_data_eval <- df_eval %>% select(starts_with("country_"))
temporal_data_eval <- df_eval %>% select(Winter, Spring, Summer, Fall)
dataY_eval <- df_eval$ged_sb
no_data_eval <- nrow(df_eval)

no_countries = ncol(spatial_data)
no_seasons = ncol(temporal_data)


# PARAMETERS FROM PYTHON CODE
# Parameters related to splines
no_knots <- 20
spline_degree <- 3
no_basis <- no_knots + spline_degree
random_walk_order <- 1

# Parameters related to priors
tau_ig_alpha <- 0.001
tau_ig_beta <- 0.001
tau_ig_alpha_spatial <- 0.001
tau_ig_beta_spatial <- 0.001
tau_ig_alpha_temporal <- 0.001
tau_ig_beta_temporal <- 0.001
gamma_alpha_nb_alpha <- 1
gamma_alpha_gamma_beta_nb_alpha <- 1
gamma_beta_gamma_beta_nb_alpha <- 0.005
intercept_gaussian_mu <- 0
intercept_gaussian_sigma <- 5

# Simulation parameters
tune_default <- 2000
draws_default <- 1000
target_accept_default <- 0.95
no_chains_default <- 4

# Data for the Stan model
data_stan <- list(
  no_data = no_data,
  no_data_eval = no_data_eval,
  spline_degree = spline_degree,
  no_basis = no_basis,
  no_interior_knots = no_knots,
  random_walk_order = random_walk_order,
  covariate_data_X1 = dataX1,
  covariate_data_X1_eval = dataX1_eval,
  covariate_data_X2 = dataX2,
  covariate_data_X2_eval = dataX2_eval,
  covariate_data_X3 = dataX3,
  covariate_data_X3_eval = dataX3_eval,
  covariate_data_X4 = dataX4,
  covariate_data_X4_eval = dataX4_eval,
  covariate_data_X5 = dataX5,
  covariate_data_X5_eval = dataX5_eval,
  no_countries=no_countries,
  spatial_data = spatial_data,
  spatial_data_eval = spatial_data_eval,
  no_seasons=no_seasons,
  temporal_data = temporal_data,
  temporal_data_eval=temporal_data_eval,
  Y = as.integer(dataY),
  Y_eval = as.integer(dataY_eval),
  a_tau_squared = tau_ig_alpha,
  b_tau_squared = tau_ig_beta,
  a_tau_squared_spatial = tau_ig_alpha_spatial,
  b_tau_squared_spatial = tau_ig_beta_spatial,
  a_tau_squared_temporal = tau_ig_alpha_temporal,
  b_tau_squared_temporal = tau_ig_beta_temporal,
  a_alpha = gamma_alpha_nb_alpha,
  alpha_b_alpha = gamma_alpha_gamma_beta_nb_alpha,
  beta_b_alpha = gamma_beta_gamma_beta_nb_alpha,
  intercept_mu = intercept_gaussian_mu,
  intercept_sigma = intercept_gaussian_sigma
)

sm = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/stan_models/bayesian_psplines_gaussian_prior_integrated_functionalities_reparameterized_corrected_baseline_features.stan")

fit <- sampling(sm, data=data_stan, iter=tune_default + draws_default, warmup=tune_default, chains=no_chains_default, control=list(adapt_delta = target_accept_default))
fit_extract = rstan::extract(fit)
save(fit, file = "fit_versionNew_gedsb_between_0_and_1000_excluding.RData")
