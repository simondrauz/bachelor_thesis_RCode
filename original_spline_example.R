library("splines")
library("rstan")

num_knots <- 10
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1

X <- seq(from=-5, to=5, by=.1)
num_data <- length(X)
knots <- unname(quantile(X,probs=seq(from=0, to=1, length.out = num_knots)))

a0 <- 0.2
a <- rnorm(num_basis, 0, 1)
B_true <- t(bs(X, knots = knots[2:9], degree=spline_degree, intercept = TRUE))
B_true_wi <- t(bs(X, knots=knots[2:9], degree=spline_degree, intercept = FALSE))
Y_true <- as.vector(a0*X + a%*%B_true)
Y <- Y_true + rnorm(length(X), 0, 0.2)

stan_data <- list(
  num_data = num_data,
  num_knots = num_knots,
  knots = knots,
  spline_degree = spline_degree,
  X = X,
  Y = Y
)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile the modified Stan model (assuming its name is "b_spline_retrieval.stan")
spline_model_for_B_original <- stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/RModels/original_spline_example_code.stan")

# Run for one iteration only, no actual sampling
fit_spline_for_B_original <- sampling(spline_model_for_B, data = stan_data, iter = 1, chains = 1, algorithm = "Fixed_param", seed = 123)

# Extract the generated B matrix
extracted_data = extract(fit_spline_for_B_original, permuted = TRUE)
B_stan_with_icpt = extracted_data$B_generated_with_icpt[1,,]
B_stan_without_icpt = extracted_data$B_generated_without_icpt[1,,]
