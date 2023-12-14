library(rstan) # or the appropriate package for your fits
library(bayesplot) # for NUTS energy plots
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(arrow)

# Define the models and their associated fit files
# models <- c('model13_nb_feature_set1', 'model3_zinb_feature_set1', 'model1_zinb_feature_set1', 
#            'model4_zinb_feature_set1', 'model15_zinb_feature_set1', 'model19_zinb_feature_set3', 
#            'model23_zinb_feature_set4')

models <- c('model15_zinb_feature_set1')

# models <- c('model3_zinb_feature_set1', 'model1_zinb_feature_set1', 
#            'model4_zinb_feature_set1', 'model15_zinb_feature_set1')

# Define the list of scientific identifiers for the models
# scientific_model_identifiers <- c('$M_{1}$', '$M_{2}$', '$M_{3}$', 
#                                  '$M_{4}$', '$M_{5}$', '$M_{6}$', '$M_{7}$')
scientific_model_identifiers <- c('$M_{5}$')
# scientific_model_identifiers <- c('$M_{2}$', '$M_{3}$', 
#                                 '$M_{4}$', '$M_{5}$')
# Initialize a dataframe to store results
results <- data.frame(
  Model = character(),
  Month = character(),
  MaxRhat = numeric(),
  MinESS = numeric(),
  TopRhatParams = character(),            # To store top Rhat parameters as a string
  BottomESSParams = character(),          # To store bottom ESS parameters as a string
  NumberRhatAbove1.1All = integer(),        # Number of Rhat > 1.1 for all parameters
  NumNumberESSBelowThresholdAll = integer(),         # Number of ESS < 100 for all parameters
  NumberRhatAbove1.1Filtered = integer(),   # Number of Rhat > 1.1 excluding specific parameters
  NumNumberESSBelowThresholdFiltered = integer(),    # Number of ESS < 100 excluding specific parameters
  RhatStatsAll = character(),             # String with quantiles and mean for all Rhat values
  RhatStatsFiltered = character(),        # String with quantiles and mean for filtered Rhat values
  ESSStatsAll = character(),              # String with quantiles and mean for all ESS values
  ESSStatsFiltered = character(),         # String with quantiles and mean for filtered ESS values
  AvgTotalTimePerChain = character(),
  NumberOfChains = integer(),
  ConvergenceIssue = logical()            # Indicator of convergence issue
)

process_fit <- function(fit, model_name, model_identifier, month, generate_trace_plots = FALSE, countries = NULL, indices = NULL) {
  print('Analyse chains and sampling times...')

  num_chains <- fit@sim$chains
  
  # Get the elapsed time for each chain
  elapsed_time <- get_elapsed_time(fit)
  # Calculate the average total time per chain (including warm-up and sampling)
  avg_total_time_per_chain <- mean(elapsed_time[, 1]) + mean(elapsed_time[,2])
  # Convert average total time per chain to hours, minutes, and seconds
  avg_total_time_hms <- convert_to_hms(avg_total_time_per_chain)
  
  print(paste0('Checking Convergence Statistics...'))
  fit_summary <- summary(fit)$summary
  param_names <- rownames(fit_summary)
  
  rhat_values <- fit_summary[,"Rhat"]
  ess_values <- fit_summary[,"n_eff"]
  
  # Handle NaN values
  rhat_values <- na.omit(rhat_values)
  ess_values <- na.omit(ess_values)
  
  # Sort and get top 10 Rhat values and names
  sorted_rhat <- sort(rhat_values, decreasing = TRUE)
  top_rhat_params <- names(sorted_rhat)[1:10]
  top_rhat_values <- sorted_rhat[1:10]
  
  # Combine parameter names with their Rhat values
  top_rhat_combined <- paste(top_rhat_params, top_rhat_values, sep = ", ", collapse = "; ")
  
  # Sort and get bottom 10 ESS values and names
  sorted_ess <- sort(ess_values)
  bottom_ess_params <- names(sorted_ess)[1:10]
  bottom_ess_values <- sorted_ess[1:10]
  
  # Combine parameter names with their ESS values
  bottom_ess_combined <- paste(bottom_ess_params, bottom_ess_values, sep = ", ", collapse = "; ")
  
  ess_threshold = num_chains * 10
  # Numbers including all parameters
  num_rhat_above_1_1_all <- sum(rhat_values > 1.1)
  num_ess_below_threshold_all <- sum(ess_values < ess_threshold)
  
  # Exclude specific parameters for Number calculation
  # Determine which parameters to exclude based on model_name
  if (model_name == 'model13_nb_feature_set1') {
    exclude_pattern <- "mu|mu_eval"
  } else if (model_name == 'model4_zinb_feature_set1') {
    exclude_pattern <- "mu|mu_eval|pi|pi_eval"
  } else {
    exclude_pattern <- "mu|mu_eval|pi|pi_eval"
  }
  
  # Exclude specific parameters for Number calculation
  excluded_params <- grep(exclude_pattern, param_names, value = TRUE)
  filtered_rhat_values <- rhat_values[!names(rhat_values) %in% excluded_params]
  filtered_ess_values <- ess_values[!names(ess_values) %in% excluded_params]
  
  num_rhat_above_1_1_filtered <- sum(filtered_rhat_values > 1.1)
  num_ess_below_threshold_filtered <- sum(filtered_ess_values < ess_threshold)
  
  # Function to calculate quantiles and mean, and format as string
  summarize_stats <- function(values) {
    quantiles <- quantile(values, probs = c(0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95))
    mean_val <- mean(values)
    stats_str <- paste("Q(0.05):", quantiles[1], 
                       "Q(0.15):", quantiles[2], 
                       "Q(0.25):", quantiles[3], 
                       "Q(0.5):", quantiles[4], 
                       "Q(0.75):", quantiles[5], 
                       "Q(0.85):", quantiles[6], 
                       "Q(0.95):", quantiles[7], 
                       "Mean:", mean_val, 
                       sep = ", ", collapse = "; ")
    return(stats_str)
  }
  
  
  # Compute and format statistics for Rhat
  rhat_stats_all <- summarize_stats(rhat_values)
  rhat_stats_filtered <- summarize_stats(filtered_rhat_values)
  
  # Compute and format statistics for ESS
  ess_stats_all <- summarize_stats(ess_values)
  ess_stats_filtered <- summarize_stats(filtered_ess_values)
  
  # Check for convergence issues
  convergence_issue <- num_rhat_above_1_1_all > 0 | num_ess_below_threshold_all > 0
  
  # Add results to the dataframe
  results <<- rbind(results, data.frame(
    Model = model_name, 
    Month = month,
    MaxRhat = max(rhat_values, na.rm = TRUE), 
    MinESS = min(ess_values, na.rm = TRUE),
    TopRhatParams = top_rhat_combined,
    BottomESSParams = bottom_ess_combined,
    NumberRhatAbove1.1All = num_rhat_above_1_1_all,
    NumNumberESSBelowThresholdAll = num_ess_below_threshold_all,
    NumberRhatAbove1.1Filtered = num_rhat_above_1_1_filtered,
    NumNumberESSBelowThresholdFiltered = num_ess_below_threshold_filtered,
    RhatStatsAll = rhat_stats_all,
    RhatStatsFiltered = rhat_stats_filtered,
    ESSStatsAll = ess_stats_all,
    ESSStatsFiltered = ess_stats_filtered,
    AvgTotalTimePerChain = avg_total_time_hms,
    NumberOfChains = num_chains,
    ConvergenceIssue = convergence_issue
  ))
  print('Creating NUTS Energy plot...')
  # Plot NUTS energy
  np <- nuts_params(fit)
  # Modify the plot titles to include the scientific identifier
  # Assuming model_identifier is something like 'M[1]'
  # Assuming model_identifier is a LaTeX-like string
  # Use latex2exp to render the LaTeX-style title
  energy_plot <- mcmc_nuts_energy(np) + 
    ggtitle(latex2exp::TeX(paste("NUTS Energy Plot - ", model_identifier, " - ", month)))
  
  plot_file <- sprintf("stan_fits/NUTS_energy_plots/energy_plot_%s_%s.png", model_name, month)
  ggsave(plot_file, energy_plot, width = 8, height = 6)  
  # Call plot_trace_plots if requested
  if (generate_trace_plots && !is.null(countries) && !is.null(indices)) {
    plot_trace_plots(fit, model_name, model_identifier, month, countries, indices)
  }
}




plot_trace_plots <- function(fit, model_name, model_identifier, month, countries, indices) {
  for (i in seq_along(countries)) {
    country <- countries[i]
    index <- indices[i]
    
    # Determine the year based on index value
    year <- ifelse(index <= 94, 2018,
                   ifelse(index <= 188, 2019,
                          ifelse(index <= 282, 2020, 2021)))
    
    if (model_name=='model13_nb_feature_set1'){
      pars_to_plot <- c(paste0("mu_eval[", index, "]"), 'alpha')
      pars_title_name <- c('\\mu', '\\alpha')
    }
    else if(model_name=='model4_zinb_feature_set1'){
      pars_to_plot <- c(paste0("mu_eval[", index, "]"), 'pi_eval', 'alpha')
      pars_title_name <- c('\\mu', '\\pi', '\\alpha')
    }
    else{
      pars_to_plot <- c(paste0("mu_eval[", index, "]"), paste0("pi_eval[", index, "]"), 'alpha')
      pars_title_name <- c('\\mu', '\\pi', '\\alpha')
    }

    
    for (j in seq_along(pars_to_plot)) {
      par <- pars_to_plot[j]
      par_title <- pars_title_name[j]
      
      if (par == 'alpha') {
        title_separate <- latex2exp::TeX(sprintf("Trace Plot for %s - %s - %s - %s", par_title, month, year, model_identifier))
      } else {
        title_separate <- latex2exp::TeX(sprintf("Trace Plot for %s - %s - %s - %s - %s", par_title, country, month, year, model_identifier))
      }
      
      p_separate <- rstan::traceplot(fit, pars = par, inc_warmup = TRUE) + ggtitle(title_separate)
      
      
      ggsave(sprintf("stan_fits/trace_plots/separate_trace_plot_%s_%s_%s_%s_%s.png", model_name, month, year, country, par), 
             p_separate, 
             width = 8, 
             height = 4)
    }
  }
}

convert_to_hms <- function(seconds) {
  hours <- floor(seconds / 3600)
  minutes <- floor((seconds %% 3600) / 60)
  seconds <- round(seconds %% 60)
  return(sprintf("%02d:%02d:%02d", hours, minutes, seconds))
}


countries <- c('Saudi Arabia', 'Syria', 'Ethiopia', 'Cameroon', 'Spain')
# Indices corresponding to SA 2018, Syr 2019, Eth 2020, Cam 2021, Spa 2018
country_indices <- c(53, 179, 211, 310, 14)
# Process each model and fit
for (i in seq_along(models)) {
  model <- models[i]
  model_identifier <- scientific_model_identifiers[i]
  for (month in c("Jan")) {
    print(paste0('Reading Fit of ', model, ' and ', month, '...'))
    file_path <- sprintf("E:/stan_fits/%s/fit_%s_%s2018_composed.RData", model, model, month)
    load(file_path)
    process_fit(get(paste0('fit_', month)), model, model_identifier, month, generate_trace_plots = TRUE, countries = countries, indices = country_indices )
    print('Removing Fit Object...')
    # Optionally, remove the fit from memory
    # Correct way to remove the fit object from memory
    rm(list = paste0('fit_', month))
  }
}


write_parquet(results, 'stan_fits/results_convergence_assessment.parquet')
c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

#mcmc_trace(fit_Jan, pars = "mu_eval[53]", 
#           main = expression(paste("Trace plot for ", alpha, " - ", country, " - ", month, " - ", year)))
#mu_eval_data <- as.data.frame(fit_Jan)$mu_eval[53]
rm(fit_Jan)
load('E:/stan_fits/model15_zinb_feature_set1/fit_model15_zinb_feature_set1_Jan2018_composed.RData')
rstan::traceplot(fit_Jan, pars= "alpha", inc_warmup=TRUE) + ggtitle("Your Desired Title Here")
for (i in seq_along(models)){
  print(i)
}
results_prev <- results
# Alternatively, you can use the format function when printing or outputting the dataframe
results_formatted <- data.frame(lapply(results, function(x) {
  if(is.numeric(x)) format(x, scientific = FALSE) else x
}))
summary_example <- summary(fit_Jan)$summary
elapsed_time <- get_elapsed_time(fit_Jan)
mean <- mean(elapsed_time[, 1]) + mean(elapsed_time[,2])

total_time_per_chain <- rowSums(elapsed_time)
avg_total_time_per_chain <- mean(total_time_per_chain)
