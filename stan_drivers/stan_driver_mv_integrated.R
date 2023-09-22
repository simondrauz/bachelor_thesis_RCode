# Load packages
library("splines")
library("rstan")
library('arrow')

# Load data
country = 'Syria'
df <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/fatality_data_Syria.parquet")
# Preliminary definition of evaluation data until better integration of cross validation

# Take the last 24 observations of training data
df_eval <- tail(df, 24)

# Select 3rd till 14th observation (equals year 2019 here) as evaluation data
df_eval <- df_eval[3:14, ]

# Cut the last 24 observations from the training data
df <- df[-(nrow(df) - 23):(nrow(df)),]

dataX1 <- df$ged_sb_tlag_1
dataX2 <- df$ged_sb_tsum_24
dataY <- df$ged_sb
no_data <- nrow(df)

dataX1_eval <- df_eval$ged_sb_tlag_1
dataX2_eval <- df_eval$ged_sb_tsum_24
dataY_eval <- df_eval$ged_sb
no_data_eval <- nrow(df_eval)


# PARAMETERS FROM PYTHON CODE
# Parameters related to splines
no_knots_X1 <- 10
no_knots_X2 <- 10
spline_degree <- 3
no_basis_X1 <- no_knots_X1 + spline_degree
no_basis_X2 <- no_knots_X2 + spline_degree

# Parameters related to priors
tauX1_ig_alpha_gaussian_sigma <- 1
tauX2_ig_alpha_gaussian_sigma <- 1
tauX1_ig_beta_gaussian_sigma <- 0.005
tauX2_ig_beta_gaussian_sigma <- 0.005
gamma_alpha_nb_alpha <- 0.1
gamma_beta_nb_alpha <- 0.1
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
  covariate_data_X2 = dataX2,
  covariate_data_X2_eval = dataX2_eval,
  no_basis_X2 = no_basis_X2,
  no_interior_knots_X2 = no_knots_X2,
  Y = as.integer(dataY),
  Y_eval = as.integer(dataY_eval),
  a_tau_X1 = tauX1_ig_alpha_gaussian_sigma,
  a_tau_X2 = tauX2_ig_alpha_gaussian_sigma,
  b_tau_X1 = tauX1_ig_beta_gaussian_sigma,
  b_tau_X2 = tauX2_ig_beta_gaussian_sigma,
  a_alpha = gamma_alpha_nb_alpha,
  b_alpha = gamma_beta_nb_alpha,
  intercept_mu = intercept_gaussian_mu,
  intercept_sigma = intercept_gaussian_sigma
)

sm = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/bayesian_psplines_gaussian_prior_integrated_functionalities.stan")

fit <- sampling(sm, data=data_stan, iter=tune_default + draws_default, warmup=tune_default, chains=no_chains_default, control=list(adapt_delta = target_accept_default))
