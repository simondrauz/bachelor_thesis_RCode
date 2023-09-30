# Load packages
library("splines")
library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library('arrow')

# Load data
df_raw <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/data_cm_features_allyears.parquet")
df_raw_expanded_feature_space_non_std <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/data_cm_features_allyears_expanded_feature_space_not_std.parquet")
df_raw_expanded_feature_space_std <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/data_cm_features_allyears_expanded_feature_space_std.parquet")
df_raw <- df_raw_expanded_feature_space_non_std
df_raw <- df_raw_expanded_feature_space_std

df_raw <- df_raw[df_raw$ged_sb < 1000,]
df_raw <- df_raw[df_raw$ged_sb > 0,]

df_raw <- df_raw[df_raw$country_id %in% developed_countries_ids, ]
df_raw <- df_raw[df_raw$country_id %in% countries_in_transition_ids, ]
df_raw <- df_raw[df_raw$country_id %in% developing_countries_extended_without_LDC_ids, ]
df_raw <- df_raw[df_raw$country_id %in% least_developed_countries_ids, ]

# Subset the dataframe for the specified columns
df <- df_raw[-c(1:3), c('month_id', 'ged_sb', 'ged_sb_tlag_3', 'ged_sb_tsum_24_tlag_3')]

# Remove Nan values from lags of countries which don't have obeservations for the first months and therefore missing lags at a later point
# Remove NaN values from 'ged_sb_tlag_3' column
df <- df[!is.na(df$ged_sb_tlag_3), ]

# Remove NaN values from 'ged_sb_tsum_24_tlag_3' column
df <- df[!is.na(df$ged_sb_tsum_24_tlag_3), ]

# Define training_data and evaluation_data
# Take training data till Oct. 2019
df_train <- df[df$month_id <= 478, ]
#Take Jan. 2020 as evaluation month
df_eval <- df[df$month_id == 481, ]

dataX1 <- df_train$ged_sb_tlag_3
dataX2 <- df_train$ged_sb_tsum_24_tlag_3
dataY <- df_train$ged_sb
no_data <- nrow(df_train)

dataX1_eval <- df_eval$ged_sb_tlag_3
dataX2_eval <- df_eval$ged_sb_tsum_24_tlag_3
dataY_eval <- df_eval$ged_sb
no_data_eval <- nrow(df_eval)


# PARAMETERS FROM PYTHON CODE
# Parameters related to splines
no_knots_X1 <- 20
no_knots_X2 <- 20
spline_degree <- 3
no_basis_X1 <- no_knots_X1 + spline_degree
no_basis_X2 <- no_knots_X2 + spline_degree
random_walk_order_X1 <- 1
random_walk_order_X2 <- 1

# Parameters related to priors
tauX1_ig_alpha <- 0.001
tauX2_ig_alpha <- 0.001
tauX1_ig_beta <- 0.001
tauX2_ig_beta <- 0.001
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
  covariate_data_X1 = dataX1,
  covariate_data_X1_eval = dataX1_eval,
  no_basis_X1 = no_basis_X1,
  no_interior_knots_X1 = no_knots_X1,
  random_walk_order_X1 = random_walk_order_X1,
  covariate_data_X2 = dataX2,
  covariate_data_X2_eval = dataX2_eval,
  no_basis_X2 = no_basis_X2,
  no_interior_knots_X2 = no_knots_X2,
  random_walk_order_X2 = random_walk_order_X2,
  Y = as.integer(dataY),
  Y_eval = as.integer(dataY_eval),
  a_tau_X1_squared = tauX1_ig_alpha,
  a_tau_X2_squared = tauX2_ig_alpha,
  b_tau_X1_squared = tauX1_ig_beta,
  b_tau_X2_squared = tauX2_ig_beta,
  a_alpha = gamma_alpha_nb_alpha,
  alpha_b_alpha = gamma_alpha_gamma_beta_nb_alpha,
  beta_b_alpha = gamma_beta_gamma_beta_nb_alpha,
  intercept_mu = intercept_gaussian_mu,
  intercept_sigma = intercept_gaussian_sigma
)

sm = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/stan_models/bayesian_psplines_gaussian_prior_integrated_functionalities_reparameterized_corrected.stan")

fit <- sampling(sm, data=data_stan, iter=tune_default + draws_default, warmup=tune_default, chains=no_chains_default, control=list(adapt_delta = target_accept_default))
fit_extract = rstan::extract(fit)
save(fit, file = "fit_versionNew_gedsb_between_0_and_1000_excluding.RData")
