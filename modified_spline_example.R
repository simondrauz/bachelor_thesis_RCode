library("splines")
library("rstan")

num_knots <- 11
num_interior_knots <- 9
spline_degree <- 3
num_basis <- num_interior_knots + spline_degree
num_basis_base_model <- num_knots + spline_degree - 1

X <- seq(from=-5, to=5, by=.1)
num_data <- length(X)
knots <- seq(-5, 5, 1)

a0 <- 0.2
a <- rnorm(num_basis_base_model, 0, 1)
B_true <- t(bs(X, knots=knots[-c(1, length(knots))], degree=spline_degree, intercept = TRUE))
B_true_wi <- t(bs(X, knots=knots[-c(1, length(knots))], degree=spline_degree, intercept = FALSE))
Y_true <- as.vector(a0*X + a%*%B_true)
Y <- Y_true + rnorm(length(X), 0, 0.2)

modified_stan_data <- list(
  num_data = num_data,
  num_interior_knots = num_interior_knots,
  num_basis = num_basis,
  spline_degree = spline_degree,
  X = X,
  Y = Y
)
basis_modified_stan_data <- list(
  num_data = num_data,
  num_knots = num_knots,
  num_basis = num_basis_base_model,
  knots = knots,
  spline_degree = spline_degree,
  X = X,
  Y = Y
)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile the modified Stan model (assuming its name is "b_spline_retrieval.stan")
spline_model_for_B <- stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/modified_spline_example_code.stan")
basis_spline_model_for_B <- stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/basis_modified_spline_example_code.stan")

# Run for one iteration only, no actual sampling
fit_spline_for_B <- sampling(spline_model_for_B, data = modified_stan_data, iter = 1, chains = 1, algorithm = "Fixed_param", seed = 123)
basis_fit_spline_for_B <- sampling(basis_spline_model_for_B, data = basis_modified_stan_data, iter = 1, chains = 1, algorithm = "Fixed_param", seed = 123)

# Extract the generated B matrix from sophisticated model
extracted_data = extract(fit_spline_for_B, permuted = TRUE)
B_stan_without_icpt_modified = extracted_data$B_generated_without_icpt[1,,]
penalty_matrix_generated = extracted_data$penalty_matrix_generated[1,,]

# Extract the generated B matrix from basic model
extracted_data = extract(basis_fit_spline_for_B, permuted = TRUE)
basis_B_stan_without_icpt_modified = extracted_data$B_generated_without_icpt[1,,]
