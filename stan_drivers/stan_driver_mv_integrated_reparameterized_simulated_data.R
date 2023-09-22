# Load packages

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(MASS)
library(dplyr)
library(tidyr)
library(zoo)

random_seed = 12345
set.seed(random_seed)
# Function to generate simulated negative binomial data with random alpha
generate_negbinom_data <- function(n, mu, alpha_shape, alpha_rate) {
  alpha_rv <- rgamma(n, shape = alpha_shape, rate = alpha_rate)
  print("Summary of Alpha value to generate neg. bin. data")
  print(summary(alpha_rv))
  y <- rnbinom(n, mu = mu, size = alpha_rv)
  print("Summary of generated neg. bin. data")
  print(summary(y))
  return(y)
}

# Function to extend data with new columns
extend_data <- function(df) {
  df <- df %>%
    mutate(ged_sb_tlag_3 = lag(ged_sb, 3)) %>%
    mutate(ged_sb_tsum_24 = rollapply(ged_sb, width = 24, FUN = mean, align = "right", fill = NA)) %>%
    mutate(ged_sb_tsum_24_tlag_3 = lag(ged_sb_tsum_24, 3)) %>%
    mutate(ged_sb_tlag_3 = as.vector(scale(ged_sb_tlag_3, center = TRUE, scale = TRUE)),
           ged_sb_tsum_24_tlag_3 = as.vector(scale(ged_sb_tsum_24_tlag_3, center = TRUE, scale = TRUE)))
  
  return(df)
}
# We generate taget data based on a negative binomial distribution with mean = 1 which should in a way reflect the 
# mu in the model with mu=exp(n*) and n* a linear combination of standerdized variables
# The n generated data points are a result of n neg bin distributions with different alphas from a gamma distribution
# Generate 1100 training data points
generated_data <- generate_negbinom_data(n = 1100, mu = 1, alpha_shape = 0.1, alpha_rate = 0.1)

# Create a data frame with index and original data
df_single <- data.frame(index = seq(1, length(generated_data)), ged_sb = generated_data)

# Extend the data frame with new columns
df_extended <- extend_data(df_single)
df_extended <- df_extended %>% rename(month_id = index)

# View the first few rows to check for NaN or NA values
head(df_extended, 30)
# alpha_rv <- rgamma(1100, shape = 1, rate = 0.1)
# y <- rnbinom(1100, mu = 1, size = alpha_rv)

# summary(alpha_rv)

df_raw <- df_extended

# Subset the dataframe for the specified columns
df <- df_raw[-c(1:3), c('month_id', 'ged_sb', 'ged_sb_tlag_3', 'ged_sb_tsum_24_tlag_3')]

# Remove Nan values from lags of countries which don't have obeservations for the first months and therefore missing lags at a later point
# Remove NaN values from 'ged_sb_tlag_3' column
df <- df[!is.na(df$ged_sb_tlag_3), ]

# Remove NaN values from 'ged_sb_tsum_24_tlag_3' column
df <- df[!is.na(df$ged_sb_tsum_24_tlag_3), ]

# Define training_data and evaluation_data
df_train <- df[df$month_id <= 1000, ]
df_eval <- df[df$month_id > 1000, ]

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
  a_tau_X1 = tauX1_ig_alpha_gaussian_sigma,
  a_tau_X2 = tauX2_ig_alpha_gaussian_sigma,
  b_tau_X1 = tauX1_ig_beta_gaussian_sigma,
  b_tau_X2 = tauX2_ig_beta_gaussian_sigma,
  a_alpha = gamma_alpha_nb_alpha,
  b_alpha = gamma_beta_nb_alpha,
  intercept_mu = intercept_gaussian_mu,
  intercept_sigma = intercept_gaussian_sigma
)

sm = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/bayesian_psplines_gaussian_prior_integrated_functionalities_reparameterized.stan")

fit <- sampling(sm, data=data_stan, iter=tune_default + draws_default, warmup=tune_default, chains=no_chains_default, control=list(adapt_delta = target_accept_default))

fit_extract = rstan::extract(fit)
y_pred_train = fit_extract$y_pred_train


save(fit, file = "fit_bayesian_psplines_gaussian_prior_integrated_functionalities_reparameterized_simulated_data.RData")

stan_hist(fit, pars="alpha")
stan_hist(fit, pars="tau_X1")
stan_hist(fit, pars="tau_X2")

stan_hist(fit, pars="intercept")
stan_hist(fit, pars="spline_coefficients_X1_penalized", binwidth = 0.2)
stan_hist(fit, pars="spline_coefficients_X2_penalized", binwidth = 0.2)
stan_hist(fit, pars="spline_coefficients_X1_non_penalized", binwidth = 0.2)
stan_hist(fit, pars="spline_coefficients_X2_non_penalized", binwidth = 0.2)

# Plot first twenty posterior predicitve samples for evaluation set
samples = rstan::extract(fit, pars="y_pred_eval")
subset_samples = samples$y_pred_eval[, 1:20]
df_samples = as.data.frame(subset_samples)
for(i in 1:20) {
  # Subset to get only the i-th element for each sample
  subset_samples <- samples$y_pred_eval[, i]
  
  # Create a data frame for ggplot2
  df <- data.frame(value = subset_samples)
  
  # Use ggplot2 to plot the histogram
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 1) +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(paste("Posterior Predictive Distribution for y_pred_eval[", i, "]", sep = "")) +
    ylab("Percentage") +
    xlab("Value")
  
  print(p)
}

rstan::summary(fit, pars = "alpha")
rstan::summary(fit, pars = "tau_X1")
rstan::summary(fit, pars = "tau_X2")
rstan::summary(fit, pars = "spline_coefficients_X1_penalized")
rstan::summary(fit, pars = "mu_eval")
stan_diag(fit)

for (column_index in 1:20) {
  # Extract one column (e.g., the first column) from y_pred_train
  one_column_draws <- y_pred_train[, column_index]
  
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
