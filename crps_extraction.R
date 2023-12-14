library(scoringRules)
library(arrow)
library(rstan)

stan_fit_dir = "E:/stan_fits/"
data_dir_path = "C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/data/"
feature_set = 'feature_set4'
data_path <- paste0(data_dir_path, "cm_features_allyears_", feature_set, ".parquet")
data <- read_parquet(data_path)
eval_months <- c(457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468)
# eval_months <- c(457, 459, 460, 461, 462, 463, 464, 465, 466, 468)

# eval_months <- c(457, 459, 460, 461)

eval_months_names <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
# eval_months_names <- c('Jan', 'Mar', 'Apr', 'May')
# eval_months_names <- c('Jan', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Dec')

model_name = 'model15_zinb_feature_set1'
crps_dfs <- list()
pps_list <- list()


get_eval_id <- function(eval_month, df){
  # Select the evaluation data
  eval_months = c(eval_month, eval_month + 12, eval_month + 24, eval_month + 36)
  df_eval <- df[df$month_id %in% eval_months, ]
  df_eval <- df_eval[order(df_eval$month_id, df_eval$country_id), ]
  
  # Select only month_id and country_id columns
  # Using base R:
  df_id <- df_eval[, c("month_id", "country_id", "ged_sb")]
  
  # Using dplyr (if you prefer and have dplyr loaded):
  # df_id <- dplyr::select(df_eval, month_id, country_id)
  
  return(df_id)
}

create_pps_dataframe <- function(model_fit, eval_id) {
  # Extract using the extract function
  y_eval_samples <- extract(model_fit, pars = "y_pred_eval")$y_pred_eval
  
  # Identify rows with NaN values and print their indices
  nan_rows <- which(apply(y_eval_samples, 1, function(x) any(is.nan(x))))
  
  # Check if there are any NaN values
  if (length(nan_rows) > 0) {
    print("Indices of rows with NaN values:")
    print(nan_rows)
    
    # Drop rows with NaN values from y_eval_samples
    y_eval_samples <- y_eval_samples[-nan_rows, ]
  } else {
    print("No NaN values found.")
  }
  
  
  # Adjust draw_indexes based on the new number of draws
  draw_indexes <- seq(1, nrow(y_eval_samples), by = 2)
  
  # Initialize an empty dataframe for the result
  result_df <- data.frame(matrix(ncol = length(draw_indexes) + 3, nrow = nrow(eval_id)))
  names(result_df) <- c("month_id", "country_id", "ged_sb", paste0("draw_", 1:length(draw_indexes)))
  
  # Populate month_id and country_id
  result_df$month_id <- eval_id$month_id
  result_df$country_id <- eval_id$country_id
  result_df$ged_sb <- eval_id$ged_sb
  
  # Iterate over each column in the y_eval_samples matrix
  for (j in 1:ncol(y_eval_samples)) {
    # Extract every second draw for this month-country combination
    result_df[j, 4:(length(draw_indexes) + 3)] <- y_eval_samples[draw_indexes, j]
  }
  
  return(result_df)
}

calculate_crps <- function(result_df) {
  # Initialize an empty dataframe for the result
  crps_df <- data.frame(matrix(ncol = 3, nrow = nrow(result_df)))
  names(crps_df) <- c("month_id", "country_id", "crps")
  
  # Loop through each row in the dataframe
  for (i in 1:nrow(result_df)) {
    crps_df$month_id[i] <- result_df$month_id[i]
    crps_df$country_id[i] <- result_df$country_id[i]
    
    # Extract the observed value
    observed <- result_df$ged_sb[i]
    
    # Extract the vector of predictive samples for the row
    predictive_samples <- as.numeric(result_df[i, 4:ncol(result_df)])
    
    # Calculate CRPS for the row and store it
    # The verification package's crps function can handle ensemble forecasts directly
    crps_df$crps[i] <- scoringRules::crps_sample(y = observed, dat = predictive_samples)
  }
  
  # Return the dataframe with CRPS values
  return(crps_df)
}
compute_crps_averages <- function(crps_dfs) {
  # Define the starting month_id for January of the first year and the number of months in a year
  start_month_id <- 457
  year_length <- 12
  
  # Create a template for the months and years
  months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Year Average')
  years <- 2018:2021
  
  # Initialize an empty dataframe
  results_df <- data.frame(matrix(ncol = length(months) + 1, nrow = length(years)))
  colnames(results_df) <- c('Year', months)
  results_df$Year <- years
  
  # Calculate the mean CRPS scores for each month of each year
  for (year in years) {
    yearly_crps <- numeric(year_length)
    for (month_index in 1:year_length) {
      month_id <- start_month_id + (year - years[1]) * year_length + (month_index - 1)
      # Extract the relevant CRPS data frame from crps_dfs list
      crps_df <- crps_dfs[[month_index]]
      # Calculate mean CRPS for the month and store in the vector
      yearly_crps[month_index] <- mean(crps_df[crps_df$month_id == month_id, 'crps'], na.rm = TRUE)
    }
    # Add the monthly CRPS to the dataframe
    results_df[year - years[1] + 1, 2:(year_length+1)] <- yearly_crps
    # Calculate and add the Year Average (annual average) CRPS to the dataframe
    results_df[year - years[1] + 1, 'Year Average'] <- mean(yearly_crps, na.rm = TRUE)
  }
  
  return(results_df)
}

compute_crps_averages_by_country <- function(crps_dfs) {
  # Define the starting month_id for January of the first year and the number of months in a year
  start_month_id <- 457
  year_length <- 12
  
  # Initialize an empty list to hold dataframes for each country
  results_list <- list()
  
  # Unique country IDs in the data (assuming all dfs have the same set of countries)
  country_ids <- unique(crps_dfs[[1]]$country_id)
  
  # Loop over each country ID
  for (country_id in country_ids) {
    # Template for the months and years
    months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Year Average')
    years <- 2018:2021
    
    # Initialize an empty dataframe for this country
    country_df <- data.frame(matrix(ncol = length(months) + 2, nrow = length(years)))
    colnames(country_df) <- c('Year', 'Country_ID', months)
    country_df$Year <- years
    country_df$Country_ID <- country_id
    
    # Calculate the mean CRPS scores for each month of each year
    for (year in years) {
      yearly_crps <- numeric(year_length)
      for (month_index in 1:year_length) {
        month_id <- start_month_id + (year - years[1]) * year_length + (month_index - 1)
        # Extract the relevant CRPS data frame from crps_dfs list
        crps_df <- crps_dfs[[month_index]]
        # Filter by country_id and calculate mean CRPS for the month
        yearly_crps[month_index] <- mean(crps_df[crps_df$month_id == month_id & crps_df$country_id == country_id, 'crps'], na.rm = TRUE)
      }
      # Add the monthly CRPS to the dataframe
      country_df[year - years[1] + 1, 3:(year_length+2)] <- yearly_crps
      # Calculate and add the Year Average (annual average) CRPS to the dataframe
      country_df[year - years[1] + 1, 'Year Average'] <- mean(yearly_crps, na.rm = TRUE)
    }
    
    # Add the country's dataframe to the results list
    results_list[[as.character(country_id)]] <- country_df
  }
  
  # Combine all country dataframes into one
  results_df <- do.call(rbind, results_list)
  return(results_df)
}

concat_pps_list <- function(pps_list) {
  # Find the maximum number of columns across all data frames
  max_cols <- max(sapply(pps_list, ncol))
  
  # Standardize and concatenate the data frames
  standardized_pps_list <- lapply(pps_list, function(df) {
    # Get the number of missing columns
    missing_cols <- max_cols - ncol(df)
    cols_df <- ncol(df)
    
    # Add missing columns with NA if needed
    if (missing_cols > 0) {
      for (i in 1:missing_cols) {
        df[paste0("draw_", cols_df -3 + i)] <- NA
      }
    }
    
    return(df)
  })
  
  # Concatenate all data frames in the list
  concatenated_pps <- do.call(rbind, standardized_pps_list)
  
  return(concatenated_pps)
}







for (i in 1:length(eval_months)){
  eval_month = eval_months[i]
  eval_month_name = eval_months_names[i]
  print(paste0('Loading fit of ', eval_month_name, '...'))
  load(paste0(stan_fit_dir, model_name, '/', 'fit_', model_name, '_', eval_month_name, '2018_composed.RData'))
  print('Extracting CRPS scores...')
  eval_id <- get_eval_id(eval_month, data)
  pps <- create_pps_dataframe(get(paste0('fit_', eval_month_name)), eval_id)
  pps_list[[i]] <- pps
  crps_dfs[[i]] <- calculate_crps(pps)
  print(paste0('CRPS for ', eval_month_name, ' has been ectracted. Remove fit object...'))
  rm(list = paste0('fit_', eval_month_name))
}

crps_averages_df <- compute_crps_averages(crps_dfs)

crps_averages_by_country_df <- compute_crps_averages_by_country(crps_dfs)

# Example usage:
concatenated_pps <- concat_pps_list(pps_list)

arrow::write_parquet(crps_averages_df, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_df_averages.parquet"))
arrow::write_parquet(crps_averages_by_country_df, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_averages_by_country_df.parquet"))
arrow::write_parquet(concatenated_pps, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_posterior_predicitve_samples.parquet"))

# load(r"(C:\Users\Uwe Drauz\RProjects\bachelor_thesis\stan_fits\model4_zinb_feature_set1\fit_model4_zinb_feature_set1_Jan2018_composed.RData)")
# rstan::traceplot(fit_Jan, pars= "mu", inc_warmup=TRUE)

model15_results <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/model15_zinb_feature_set1_crps_df_averages.parquet")
model_15_scores <- read_parquet(paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_df_averages.parquet"))

load(paste0(stan_fit_dir, model_name, '/', 'fit_', model_name, '_', 'Jan', '2018_composed.RData'))

fit_summary <- summary(fit_Jan)

model15_results <- read_parquet("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/model15_zinb_feature_set1_crps_df_averages.parquet")
model_15_scores <- read_parquet(paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_df_averages.parquet"))

load(paste0(stan_fit_dir, model_name, '/', 'fit_', model_name, '_', 'Aug', '2018_composed.RData'))

fit_summary <- summary(fit_Jan)

drop <-c('Year', 'Feb', 'Nov', 'Year Average')
results_23 <- subset(crps_averages_df, select = -c('Feb', 'Nov'))
df = crps_averages_df[,!(names(crps_averages_df) %in% drop)]
year_avergae <- rowSums(df)/10
mean(year_avergae)
load(paste0(stan_fit_dir, model_name,'_run2', '/', 'fit_', model_name, '_', 'Aug', '2018_composed.RData'))
