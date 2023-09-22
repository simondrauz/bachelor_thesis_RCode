# Load packages
library("splines")
library("rstan")
library('arrow')

# Define functions
generate_penalty_matrix <- function(no_coefficients) {
  # Create a no_coefficients x no_coefficients matrix filled with zeros
  mat <- matrix(0, nrow=no_coefficients, ncol=no_coefficients)
  
  # Fill the main diagonal with 2s
  diag(mat) <- 2
  
  # Modify the [1,1] and [no_coefficients,no_coefficients] elements to be 1
  mat[1,1] <- 1
  mat[no_coefficients,no_coefficients] <- 1
  
  # Add -1s to the sub-diagonals above and below the main diagonal
  if (no_coefficients > 1) {
    mat[cbind(2:no_coefficients, 1:(no_coefficients-1))] <- -1   # Sub-diagonal below
    mat[cbind(1:(no_coefficients-1), 2:no_coefficients)] <- -1   # Sub-diagonal above
  }
  
  return(mat)
}
standardize <- function(x) {
  (x - mean(x)) / sd(x)
}
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
dataY <- standardize(dataY)
N <- nrow(df)

dataX1_eval <- df_eval$ged_sb_tlag_1
dataX2_eval <- df_eval$ged_sb_tsum_24
dataY_eval <- df_eval$ged_sb
dataY_eval <- standardize(dataY_eval)
N_eval <- nrow(df_eval)


# PARAMETERS FROM PYTHON CODE
# Parameters related to splines
num_knots_X1 <- 10
num_knots_X2 <- 10
spline_degree <- 3

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

# Determine the knots for equally spaced knots
knot_list_X1 <- seq(min(dataX1), max(dataX1), length.out = num_knots_X1 + 2)[-c(1, num_knots_X1 + 2)]
knot_list_X2 <- seq(min(dataX2), max(dataX2), length.out = num_knots_X2 + 2)[-c(1, num_knots_X2 + 2)]

# Create basis matrices
basis_X1 <- bs(dataX1, knots = knot_list_X1, degree = spline_degree, intercept = FALSE)
basis_X2 <- bs(dataX2, knots = knot_list_X2, degree = spline_degree, intercept = FALSE)

# Determine the knots for equally spaced knots for df_eval
knot_list_X1_eval <- seq(min(dataX1_eval), max(dataX1_eval), length.out = num_knots_X1 + 2)[-c(1, num_knots_X1 + 2)]
knot_list_X2_eval <- seq(min(dataX2_eval), max(dataX2_eval), length.out = num_knots_X2 + 2)[-c(1, num_knots_X2 + 2)]

# Create basis matrices for df_eval
basis_X1_eval <- bs(dataX1_eval, knots = knot_list_X1_eval, degree = spline_degree, intercept = FALSE)
basis_X2_eval <- bs(dataX2_eval, knots = knot_list_X2_eval, degree = spline_degree, intercept = FALSE)


# Assuming the generate_penalty_matrix function exists
K_X1 <- generate_penalty_matrix(ncol(basis_X1))
K_X2 <- generate_penalty_matrix(ncol(basis_X2))

# Data for the Stan model
data_mv <- list(
  N = N,
  num_basis_X1 = ncol(basis_X1),
  num_basis_X2 = ncol(basis_X2),
  basis_X1 = basis_X1,
  basis_X2 = basis_X2,
  Y = as.integer(dataY),
  K_X1 = K_X1,
  K_X2 = K_X2,
  a_tau_X1 = tauX1_ig_alpha_gaussian_sigma,
  a_tau_X2 = tauX2_ig_alpha_gaussian_sigma,
  b_tau_X1 = tauX1_ig_beta_gaussian_sigma,
  b_tau_X2 = tauX2_ig_beta_gaussian_sigma,
  a_alpha = gamma_alpha_nb_alpha,
  b_alpha = gamma_beta_nb_alpha,
  intercept_mu = intercept_gaussian_mu,
  intercept_sigma = intercept_gaussian_sigma,
  N_eval = N_eval,
  basis_X1_eval = basis_X1_eval,
  basis_X2_eval = basis_X2_eval,
  Y_eval = as.integer(dataY_eval)
)
sm = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/spline_test_model_mv_std.stan")

fit <- sampling(sm, data=data_mv, iter=tune_default + draws_default, warmup=tune_default, chains=no_chains_default, control=list(adapt_delta = target_accept_default))
