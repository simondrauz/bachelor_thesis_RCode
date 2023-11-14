# This file is meant to enable runnung different models simultainously
# AT the moment the specified model has
#   - baseline_features (set 1), non-linear
#   - spatial unstructured
#   - temporal coefficients season
#   - negative binomial distribution

# DEFINE FUNCTIONS
fit_model_for_eval_month <- function(stan_model, df, eval_month, lag_indicator, 
                                     no_knots, spline_degree, random_walk_order,  
                                     tune_iterations, draw_iterations, no_chains, 
                                     target_accept_rate, cores_for_par) {
  stan_data_model = create_stan_data_model(
    df = df,
    eval_month = eval_month,
    lag_indicator = lag_indicator,
    no_knots = no_knots,
    spline_degree = spline_degree,
    random_walk_order = random_walk_order
  )
  
  fit <- fit_model(
    stan_model = stan_model,
    stan_data = stan_data_model,
    tune_iterations = tune_iterations,
    draw_iterations = draw_iterations,
    no_chains = no_chains,
    target_accept_rate = target_accept_rate,
    cores_for_par = cores_for_par
  )
  return(fit)
}

create_stan_data_model <- function(df, eval_month, lag_indicator, no_knots, 
                                   spline_degree, random_walk_order ){
  # DOCSTRING:
  #   This function creates the stan data object for the model:
  #   - baseline_features (set 1), non-linear
  #   - spatial unstructured
  #   - temporal coefficients season
  #   - negative binomial distribution
  
  # Select the training data
  last_train_month = eval_month-lag_indicator
  df_train <- df[df$month_id <= last_train_month, ]
  # Determine length of trainin data set
  no_data <- nrow(df_train)
  # Select covariate data from df for training
  
  dataX1 <- df_train[[paste0('ged_sb_tlag_', lag_indicator)]]
  dataX2 <- df_train[[paste0('ged_sb_tlag_', (lag_indicator + 1))]]
  dataX3 <- df_train[[paste0('ged_sb_tlag_', (lag_indicator + 2))]]
  dataX4 <- df_train[[paste0('decay_ged_sb_5_tlag_', lag_indicator)]]
  dataX5 <- df_train[[paste0('decay_ged_sb_100_tlag_', lag_indicator)]]
  dataX6 <- df_train[[paste0('decay_ged_sb_500_tlag_', lag_indicator)]]
  dataX7 <- df_train[[paste0('ged_sb_tsum_24_tlag_', lag_indicator)]]
  # Select spatial data from df for training
  spatial_data <- df_train %>% select(starts_with("country_"))
  # Drop country_id column from spatial data
  spatial_data <- subset(spatial_data, select = -c(country_id))
  # Select target data from df for training
  dataY <- df_train$ged_sb
  
  # Select the evaluation data
  eval_months = c(eval_month, eval_month+12, eval_month+24, eval_month+36)
  df_eval <- df[df$month_id %in% eval_months, ]
  df_eval <- df_eval[order(df_eval$month_id, df_eval$country_id), ]
  # Determine length of trainin data set
  no_data_eval <- nrow(df_eval)
  # Select covariate data from df for evaluation
  dataX1_eval <- df_eval[[paste0('ged_sb_tlag_', lag_indicator)]]
  dataX2_eval <- df_eval[[paste0('ged_sb_tlag_', (lag_indicator + 1))]]
  dataX3_eval <- df_eval[[paste0('ged_sb_tlag_', (lag_indicator + 2))]]
  dataX4_eval <- df_eval[[paste0('decay_ged_sb_5_tlag_', lag_indicator)]]
  dataX5_eval <- df_eval[[paste0('decay_ged_sb_100_tlag_', lag_indicator)]]
  dataX6_eval <- df_eval[[paste0('decay_ged_sb_500_tlag_', lag_indicator)]]
  dataX7_eval <- df_eval[[paste0('ged_sb_tsum_24_tlag_', lag_indicator)]]
  # Select spatial data from df for training
  spatial_data_eval <- df_eval %>% select(starts_with("country_"))
  # Drop country_id column from spatial data
  spatial_data_eval <- subset(spatial_data_eval, select = -c(country_id))
  # Select target data from df for training
  dataY_eval <- df_eval$ged_sb
  
  # Determine factor relevevant for defining regression coefficients dimension
  no_countries = ncol(spatial_data)
  
  # Define parameters related to splines
  no_knots <- no_knots
  spline_degree <- spline_degree
  no_basis <- no_knots + spline_degree
  random_walk_order <- random_walk_order
  
  # Define parameters for Priors
  # Parameters of tau ~ IG(a. b) for regularization
  tau_ig_alpha <- 0.001
  tau_ig_beta <- 0.001
  tau_ig_alpha_spatial <- 0.001
  tau_ig_beta_spatial <- 0.001
  
  # Parameters for alpha ~ G(a, b) with y ~ NB(mu, alpha)
  gamma_alpha_nb_alpha <- 1
  gamma_alpha_gamma_beta_nb_alpha <- 1
  gamma_beta_gamma_beta_nb_alpha <- 0.005
  
  
  # Put information into stan_data object
  stan_data <- list(
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
    covariate_data_X6 = dataX6,
    covariate_data_X6_eval = dataX6_eval,
    covariate_data_X7 = dataX7,
    covariate_data_X7_eval = dataX7_eval,
    
    no_countries=no_countries,
    spatial_data = spatial_data,
    spatial_data_eval = spatial_data_eval,
    Y = as.integer(dataY),
    Y_eval = as.integer(dataY_eval),
    a_tau_squared = tau_ig_alpha,
    b_tau_squared = tau_ig_beta,
    a_tau_squared_spatial = tau_ig_alpha_spatial,
    b_tau_squared_spatial = tau_ig_beta_spatial,
    a_alpha = gamma_alpha_nb_alpha,
    alpha_b_alpha = gamma_alpha_gamma_beta_nb_alpha,
    beta_b_alpha = gamma_beta_gamma_beta_nb_alpha
  )
  return(stan_data)
}
fit_model <- function(stan_model, stan_data, tune_iterations, draw_iterations, no_chains, 
                      target_accept_rate, cores_for_par){
  # DOCSTRING
  #   This function is meant to initiate the fit of the passed model, data, and sampling specifications
  options(mc.cores = (no_chains + cores_for_par))
  fit <- sampling(stan_model, data=stan_data, iter=tune_iterations + draw_iterations, warmup=tune_iterations, chains=no_chains, 
                  control=list(adapt_delta = target_accept_rate), refresh=300)
  return(fit)
}
# SET RSTAN SAMPLING OPTIONS
rstan_options(auto_write = TRUE)


# READ DATA
data_path <- paste0(data_dir_path, "cm_features_allyears_", feature_set, ".parquet")
df_raw <- read_parquet(data_path)

# Compile Stan Model
model_path = paste0(stan_models_dir_path, model, "_", distribution_assumption, enable_parallel, ".stan")
stan_model = stan_model(model_path)
# List of evaluation months you want to process (for now we test with Jan2018-Apr2018)
# Start fitting of stan_model for each of the specified evaluation months
samplingfunction <- function(x){
  if(x %in% 1:length(eval_months)){
    res <- fit_model_for_eval_month(stan_model = stan_model, 
                                    df = df_raw, 
                                    eval_month = eval_months[x],
                                    lag_indicator = lag_indicators[x], 
                                    no_knots = 20, 
                                    spline_degree = 3, 
                                    random_walk_order = 1,  
                                    tune_iterations = 1000, 
                                    draw_iterations = 500, 
                                    no_chains = no_chains, 
                                    target_accept_rate = 0.8,
                                    cores_for_par = cores_for_par)
    return(res)
  } else {
    stop("Unexpected value of x!")
  }
}

# Retrieve results
cl <- makeCluster(((no_chains + cores_for_par) * length(eval_months)) + 1)
# Load the dplyr library (or magrittr for just the pipe) in each worker
clusterEvalQ(cl, {
  library(dplyr)
  library(rstan)
})
clusterExport(cl,c('stan_model','df_raw', 'fit_model_for_eval_month', 'create_stan_data_model', 'fit_model', 'eval_months', 'lag_indicators', 'no_chains', 'cores_for_par'))
out <- parLapply(cl, c(1:length(eval_months)),samplingfunction)

for (i in 1:length(eval_months)) {
  month_indicator = month_indicators[i]
  assign(paste0("fit_", month.abb[month_indicator]), out[[i]])
  save(list = paste0("fit_", month.abb[month_indicator]), 
       file = paste0(stan_fit_dir, "fit_", model, "_", distribution_assumption, "_", feature_set, "_",  month.abb[month_indicator], "2018_composed.RData"), 
       compress = TRUE)
  rm(list = paste0("fit_", month.abb[month_indicator]))  # remove from memory after saving
}

# Stop the cluster
stopCluster(cl)

fit_test <- fit_model_for_eval_month(stan_model = stan_model, 
                                     df = df_raw, 
                                     eval_month = 457,
                                     lag_indicator = 3, 
                                     no_knots=20, 
                                     spline_degree=3, 
                                     random_walk_order=1,  
                                     tune_iterations=1000, 
                                     draw_iterations=500, 
                                     no_chains=1, 
                                     target_accept_rate=0.8,
                                     cores_for_par=3)


