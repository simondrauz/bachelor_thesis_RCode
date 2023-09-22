# Priod predicitve Sampling
set.seed(12345)
library(rstan)
library(arrow)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
sm_ppc = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/bayesian_psplines_gaussian_prior_integrated_functionalities_reparameterized_ppc.stan")

# VERSION WHICH WAS SIMULATED BY CHATGPT
# Set Up Prior Distributions for Data
prior_data <- list(
  no_data = 100,
  spline_degree = 3,
  no_basis_X1 = 10,
  no_interior_knots_X1 = 7,
  random_walk_order_X1 = 1,
  no_basis_X2 = 10,
  no_interior_knots_X2 = 7,
  random_walk_order_X2 = 1,
  a_tau_X1 = 2,
  b_tau_X1 = 1,
  a_tau_X2 = 2,
  b_tau_X2 = 1,
  a_alpha = 2,
  b_alpha = 1,
  intercept_mu = 0,
  intercept_sigma = 1,
  covariate_data_X1 = runif(100, 0, 1),  # Mock data for covariate X1
  covariate_data_X1_eval = runif(50, 0, 1),  # Mock evaluation data for covariate X1
  covariate_data_X2 = runif(100, 0, 1),  # Mock data for covariate X2
  covariate_data_X2_eval = runif(50, 0, 1)  # Mock evaluation data for covariate X2
)

# VERSION WHICH IS CLOSE TO OUT MODEL
# Load data
df_raw <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/data_cm_features_allyears.parquet")
# Subset the dataframe for the specified columns
df <- df_raw[-c(1:3), c('month_id', 'ged_sb', 'ged_sb_tlag_3', 'ged_sb_tsum_24_tlag_3')]

# Remove Nan values from lags of countries which don't have obeservations for the first months and therefore missing lags at a later point
# Remove NaN values from 'ged_sb_tlag_3' column
df <- df[!is.na(df$ged_sb_tlag_3), ]

# Remove NaN values from 'ged_sb_tsum_24_tlag_3' column
df <- df[!is.na(df$ged_sb_tsum_24_tlag_3), ]


# Define training_data and evaluation_data
# 120 months (10 years)
df_train <- df[df$month_id == 478, ]
df_eval <- df[df$month_id == 481, ]

dataX1 <- df_train$ged_sb_tlag_3
dataX2 <- df_train$ged_sb_tsum_24_tlag_3
dataY <- df_train$ged_sb
no_data <- nrow(df_train)


# PARAMETERS FROM PYTHON CODE
# Parameters related to splines
no_knots_X1 <- 10
no_knots_X2 <- 10
spline_degree <- 3
no_basis_X1 <- no_knots_X1 + spline_degree
no_basis_X2 <- no_knots_X2 + spline_degree
random_walk_order_X1 <- 1
random_walk_order_X2 <- 1
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
prior_data <- list(
  no_data = no_data,
  spline_degree = spline_degree,
  covariate_data_X1 = dataX1,
  no_basis_X1 = no_basis_X1,
  no_interior_knots_X1 = no_knots_X1,
  random_walk_order_X1 = random_walk_order_X1,
  covariate_data_X2 = dataX2,
  no_basis_X2 = no_basis_X2,
  no_interior_knots_X2 = no_knots_X2,
  random_walk_order_X2 = random_walk_order_X2,
  a_tau_X1 = tauX1_ig_alpha_gaussian_sigma,
  a_tau_X2 = tauX2_ig_alpha_gaussian_sigma,
  b_tau_X1 = tauX1_ig_beta_gaussian_sigma,
  b_tau_X2 = tauX2_ig_beta_gaussian_sigma,
  a_alpha = gamma_alpha_nb_alpha,
  b_alpha = gamma_beta_nb_alpha,
  intercept_mu = intercept_gaussian_mu,
  intercept_sigma = intercept_gaussian_sigma
)

# Run Prior Predictive Sampling
prior_samples <- sampling(sm_ppc, data = prior_data, algorithm = "Fixed_param", iter = 1000)

# Extract the samples for analysis.
prior_draws <- rstan::extract(prior_samples)

library(ggplot2)

# Assuming y_sim is what you want to plot
y_sim_draws <- prior_draws$y_sim

# Convert to data frame for ggplot2
# Convert the matrix to a long-format data frame
y_sim_df <- as.data.frame(as.table(as.matrix(y_sim_draws)))

# Here we compare the distribution of y values (which resulted from the priors and input data) 
for (column_index in 1:no_data) {
  # Extract one column (e.g., the first column) from y_sim_draws
  one_column_draws <- y_sim_draws[, column_index]
  
  # Create a data frame for ggplot2
  one_column_df <- data.frame(Value = one_column_draws)
  
  # Get the actual data observation for this column_index
  actual_data_value <- dataY[column_index]
  
  # Plot using ggplot2
  p <- ggplot(one_column_df, aes(x = Value)) +
    geom_histogram(aes(y = (..count..)/sum(..count..) * 100), binwidth = 1, fill = "blue", alpha = 0.7) +
    geom_vline(aes(xintercept = actual_data_value), color = "red", linetype = "dashed", size = 1) +
    xlab(paste("y_sim values for column: ", column_index)) +
    ylab("Frequency (%)") +
    ggtitle(paste("Prior Draws of y_sim for Column: ", column_index))
  
  print(p)  # Explicitly print the plot in a loop
}
