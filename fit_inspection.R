# Load necessary libraries
library(StanHeaders)
library(rstan)
library(ggplot2)
library(reshape2)
library(gridExtra)

library(tidyr)
library(dplyr)

library(bayesplot)

library(scoringRules)
library(SpecsVerification)


# Functions

# Plots samples histograms for all countries at the actuals observation together with the observation
ppc_hist_plus_obs <- function(fit_extract, df_eval, country_mapping){
  
  # Loop through the observations (number of columns in fit_extract$y_pred_eval)
  for(i in 1:ncol(fit_extract$y_pred_eval)){
    
    # Extract the 2000 iterations for the i-th observation
    posterior_samples <- fit_extract$y_pred_eval[,i]
    
    # Get the observed value for the i-th observation
    observed_value <- df_eval$ged_sb[i]
    
    # Get the country name using the country_id
    country_name <- country_mapping$name[country_mapping$country_id == df_eval$country_id[i]]
    
    # Plot
    hist(posterior_samples, probability = TRUE, 
         main = paste("Posterior Predictive Samples vs Observation for", country_name),
         xlab = "Value", ylab = "Density", col = "lightblue")
    
    # Add observed value as a vertical line
    abline(v = observed_value, col = "red", lwd = 2)
    
    legend("topright", legend = c("Observation", "Posterior Samples"),
           fill = c("red", "lightblue"))
    
    # Pause for the user to view each plot
    readline(prompt="Press [enter] to see the next plot.")
  }
}

# Plots KDE and CDF for all countries at the actuals observation together with the observation
ppc_kde_and_cdf_plus_obs <- function(fit_extract, df_eval, country_mapping) {
  
  # Loop through the observations (number of columns in fit_extract$y_pred_eval)
  for(i in 1:ncol(fit_extract$y_pred_eval)) {
    
    # Extract the samples for the i-th observation
    posterior_samples <- fit_extract$y_pred_eval[, i]
    
    # Create a data frame for ggplot2
    posterior_samples_df <- data.frame(Value = posterior_samples)
    
    # Get the observed value for the i-th observation
    observed_value <- df_eval$ged_sb[i] # Changed from get_sp to ged_sb
    
    # Check if observed_value is numeric and not NA
    if(is.numeric(observed_value) & !is.na(observed_value)) {
      
      # Get the country name using the country_id
      country_name <- country_mapping$name[country_mapping$country_id == df_eval$country_id[i]]
      
      # KDE plot
      max_density <- max(density(posterior_samples)$y)
      p1 <- ggplot(posterior_samples_df, aes(x = Value)) + 
        geom_density(aes(y = after_stat(density)), fill = "darkblue", alpha = 0.5) + 
        geom_vline(aes(xintercept = observed_value), color = "red", linetype = "dashed") +
        geom_text(aes(x = max(posterior_samples), y = max_density, label = paste("Obs:", round(observed_value, 2))), hjust = 1, color = "red") +
        ggtitle(paste("KDE for", country_name)) +
        xlab("Value") +
        ylab("Density")
      
      # CDF plot
      p2 <- ggplot(posterior_samples_df, aes(x = Value)) + 
        stat_ecdf(geom = "step", color = "darkblue") + 
        geom_vline(aes(xintercept = observed_value), color = "red", linetype = "dashed") +
        geom_text(aes(x = max(posterior_samples), y = 1, label = paste("Obs:", round(observed_value, 2))), hjust = 1, color = "red") +
        ggtitle(paste("CDF for", country_name)) +
        xlab("Value") +
        ylab("Cumulative Probability")
      
      # Display plots side-by-side
      grid.arrange(p1, p2, ncol = 2)
      
      # Pause for the user to view each plot
      readline(prompt="Press [enter] to see the next plot.")
    } else {
      print(paste("Skipping plot for index", i, "because observed_value is not numeric or is NA"))
    }
  }
}

# Compute the CRPS score for each country (not working)
compute_CRPS_with_package <- function(fit_extract, df_eval, country_mapping) {
  
  # Initialize an empty vector to store the CRPS scores
  crps_scores <- numeric(0)
  country_names <- character(0)
  
  # Loop through the observations (number of columns in fit_extract$y_pred_eval)
  for(i in 1:ncol(fit_extract$y_pred_eval)) {
    
    # Extract the samples for the i-th observation
    posterior_samples <- fit_extract$y_pred_eval[, i]
    
    # Get the observed value for the i-th observation
    observed_value <- df_eval$ged_sb[i]
    
    # Check if observed_value is numeric and not NA
    if(is.numeric(observed_value) & !is.na(observed_value)) {
      
      # Calculate the CRPS for the i-th observation
      crps_score <- SpecsVerification::EnsCrps(posterior_samples, observed_value)
      
      # Append the CRPS score to the vector
      crps_scores <- c(crps_scores, crps_score)
      
      # Get the country name using the country_id
      country_name <- country_mapping$name[country_mapping$country_id == df_eval$country_id[i]]
      
      # Append the country name to the vector
      country_names <- c(country_names, country_name)
      
    } else {
      print(paste("Skipping CRPS calculation for index", i, "because observed_value is not numeric or is NA"))
    }
  }
  
  # Create a data frame to store CRPS scores and corresponding country names
  crps_df <- data.frame(Country = country_names, CRPS = crps_scores)
  
  return(crps_df)
}


# Load Fit
load("stan_fits/fit_fatalities_bigger_0.RData")
load("stan_fits/fit_zinb_all_fatalities.RData")
fit <- fit_parallel
# Extract sample from fit
posterior_incl_warmup_combined <- rstan::extract(fit_zinb_parallel, inc_warmup = TRUE, permuted = TRUE)
# Extract sample from fit without warmup
posterior_excl_warmup_combined <- rstan::extract(fit_zinb_parallel, inc_warmup = FALSE, permuted = TRUE)
pp_eval <-posterior_excl_warmup_combined$y_pred_eval
# Inclusive warm_up
posterior_incl_warmup <- rstan::extract(fit_parallel, inc_warmup = TRUE, permuted = FALSE)

# Seems to only take the samples, not the warmup
posterior <- as.matrix(fit_parallel)

# Get fit information including convergence statistics as a table
fit_summary <- summary(fit_parallel)
fit_summary_table <- summary(fit_parallel)$summary

# Load the country mapping from the CSV file
country_mapping <- read.csv("country_list.csv", stringsAsFactors = FALSE)


# Plots
ppc_hist_plus_obs(fit_extract, df_eval, country_mapping)
ppc_kde_and_cdf_plus_obs(posterior_incl_warmup_combined, df_eval, country_mapping)

#CRPS
crps_df <- compute_CRPS_with_package(posterior_incl_warmup_combined, df_eval, country_mapping)


# Trace plots of warm-up +  sampling for selected parameters using rstan
rstan::traceplot(fit_zinb_parallel, pars= "alpha", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= c("tau_squared", "tau_squared_spatial", "tau_squared_temporal"), inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "intercept", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "b_alpha", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= c("intercept","alpha", "b_alpha", "tau_squared", "tau_squared_spatial", "tau_squared_temporal"), inc_warmup=TRUE)

rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X1_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X1_non_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X2_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X2_non_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X3_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X3_non_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X4_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X4_non_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X5_penalized", inc_warmup=TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "spline_coefficients_X5_non_penalized", inc_warmup=TRUE)

rstan::traceplot(fit_zinb_parallel, pars= "spatial_coefficients", inc_warmup = TRUE)
rstan::traceplot(fit_zinb_parallel, pars= "temporal_coefficients", inc_warmup = TRUE)


# Trace plot for sampling iterations using bayesplot
bayesplot::mcmc_trace(posterior, pars= c("tau_squared", "tau_squared_spatial", "tau_squared_temporal"))


# Trace plot for warm_up + sampling, bt some issue with the passed data containing null values
bayesplot::mcmc_trace(posterior_incl_warmup, pars= c("tau_squared", "tau_squared_spatial", "tau_squared_temporal"), n_warmup=1250)

# Density plots for selected parameters using rstan
rstan::stan_dens(fit, pars = c("mu_eval", "alpha"))

# Density plots for selected parameters using bayesplot with medians and 80% intervals
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
bayesplot::mcmc_areas(posterior,
           pars = c("tau_squared", "tau_squared_spatial", "tau_squared_temporal"),
           prob = 0.8) + plot_title
bayesplot::mcmc_areas(posterior,
                      pars = c("alpha", "b_alpha"),
                      prob = 0.8) + plot_title
bayesplot::mcmc_areas(posterior,
                      pars = c("intercept"),
                      prob = 0.8) + plot_title

# For spline_coefficients_X1_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X1_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X1_penalized")

# For spline_coefficients_X1_non_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X1_non_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X1_non_penalized")

# For spline_coefficients_X2_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X2_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X2_penalized")

# For spline_coefficients_X2_non_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X2_non_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X2_non_penalized")

# For spline_coefficients_X3_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X3_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X3_penalized")

# For spline_coefficients_X3_non_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X3_non_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X3_non_penalized")

# For spline_coefficients_X4_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X4_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X4_penalized")

# For spline_coefficients_X4_non_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X4_non_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X4_non_penalized")

# For spline_coefficients_X5_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X5_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X5_penalized")

# For spline_coefficients_X5_non_penalized
bayesplot::mcmc_areas(posterior, 
                      pars = c("spline_coefficients_X5_non_penalized"), 
                      prob = 0.8) + 
  ggtitle("Traceplot for spline_coefficients_X5_non_penalized")




# scatter plot also showing divergences (don't really know how valuabel)
color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(fit_parallel),
  pars = c("alpha", "tau_squared"), 
  np = nuts_params(fit_parallel), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

# Nuts energy plot
color_scheme_set("red")
np <- nuts_params(fit_zinb_parallel)
mcmc_nuts_energy(np) + ggtitle("NUTS Energy Diagnostic")

# Obtain number of divergences
divergences <- rstan::get_sampler_params(fit_parallel, inc_warmup = FALSE)$divergent__
num_divergences <- sum(divergences)
print(num_divergences)
rstan::check_hmc_diagnostics(fit_parallel)


# EXPLORATIVE
# Draft to plot all posterior samples in training data of one country against actuals data in training data of one country
# Specify the country_id
country_id <- "220"

# Subset the relevant observations
relevant_obs <- spatial_data %>%
  filter(!!sym(paste0("country_", country_id)) == 1)

# Get indices of the relevant observations
indices <- which(!!sym(paste0("country_", country_id))(spatial_data) == 1)



# Subset the posterior samples
posterior_country <- posterior_samples[, relevant_obs$row_number]

# Gather the posterior samples
posterior_country_long <- as.data.frame(posterior_country) %>%
  gather(key = "sample", value = "predicted")

# Add the observed data (ged_sb) to the long format dataframe
posterior_country_long$observed <- rep(relevant_obs$ged_sb, nrow(posterior_country))

# Plot
ggplot(posterior_country_long, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  theme_minimal() +
  labs(title = paste("Posterior Predictive vs Observed for Country", country_id),
       x = "Observed Data (ged_sb)", y = "Posterior Predictive Sample")

# Extract samples
posterior_samples <- rstan::extractact(fit_parallel, permuted = TRUE)
y_rep_train <- fit_extract$y_pred_train
dim(y_rep_train)
y_rep_eval <- fit_extract$y_pred_eval
dim(y_rep_eval)
color_scheme_set("brightblue")
ppc_dens_overlay(df_train$ged_sb, y_rep_train_truncated[1:500, ]) + xlim(500, 2500)

# Step 1: Remove zeros from each column and store in a list
y_rep_train_no_zeros <- apply(y_rep_train, 2, function(x) x[x != 0])

# Step 2: Determine the minimum number of non-zero elements across all columns
min_length <- min(sapply(y_rep_train_no_zeros, length))

# Step 3: Truncate each column vector to that length and reconstruct the matrix
y_rep_train_truncated <- do.call(cbind, lapply(y_rep_train_no_zeros, function(x) x[1:min_length]))


posterior_samples_zinb <- rstan::extract(fit_zinb_parallel, inc_warmup = FALSE, permuted = TRUE)$y_pred_eval
posterior_samples <- posterior_incl_warmup_combined$y_pred_eval

write.csv(posterior_samples_zinb, "C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/posterior_samples_zinb.csv")
write.csv(posterior_samples, "C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/posterior_samples.csv")

write.csv(df_eval, "C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/temp/df_eval.csv")
