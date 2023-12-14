library(scoringRules)
library(arrow)
library(rstan)

data_dir_path = "C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/data/"
feature_set = 'feature_set1'
data_path <- paste0(data_dir_path, "cm_features_allyears_", feature_set, ".parquet")
data <- read_parquet(data_path)
eval_months <- c(457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468)
eval_months_names <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
model_name = 'model15_zinb_feature_set1'
crps_train_dfs <- list()
pps_train_list <- list()
crps_eval_dfs <- list()
pps_eval_list <- list()
last_training_month = 454


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


create_pps_dataframe <- function(model_fit, id, pps_name) {
  # Extract using the extract function
  if (pps_name == 'y_pred_train') {
    y_samples <- extract(model_fit, pars = 'y_pred_train')$y_pred_train
  } else if (pps_name == 'y_pred_eval') {
    y_samples <- extract(model_fit, pars = 'y_pred_eval')$y_pred_eval
  }
  # Identify rows with NaN values and print their indices
  nan_rows <- which(apply(y_samples, 1, function(x) any(is.nan(x))))
  
  # Check if there are any NaN values
  if (length(nan_rows) > 0) {
    print("Indices of rows with NaN values:")
    print(nan_rows)
    
    # Drop rows with NaN values from y_samples
    y_samples <- y_samples[-nan_rows, ]
  } else {
    print("No NaN values found.")
  }
  
  
  # Adjust draw_indexes based on the new number of draws
  draw_indexes <- seq(1, nrow(y_samples), by = 2)
  
  # Initialize an empty dataframe for the result
  result_df <- data.frame(matrix(ncol = length(draw_indexes) + 3, nrow = nrow(id)))
  names(result_df) <- c("month_id", "country_id", "ged_sb", paste0("draw_", 1:length(draw_indexes)))
  
  # Populate month_id and country_id
  result_df$month_id <- id$month_id
  result_df$country_id <- id$country_id
  result_df$ged_sb <- id$ged_sb
  
  # Iterate over each column in the y_samples matrix
  for (j in 1:ncol(y_samples)) {
    # Extract every second draw for this month-country combination
    result_df[j, 4:(length(draw_indexes) + 3)] <- y_samples[draw_indexes, j]
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

compute_crps_averages <- function(crps_dfs, start_year, end_year, start_month) {
  # Define the starting month_id for January of the first year and the number of months in a year
  start_month_id <- start_month
  year_length <- 12
  
  # Create a template for the months and years
  months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Year Average')
  years <- start_year:end_year
  
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

compute_crps_averages_by_country <- function(crps_dfs, start_year, end_year, start_month) {
  # Define the starting month_id for January of the first year and the number of months in a year
  start_month_id <- start_month
  year_length <- 12
  
  # Initialize an empty list to hold dataframes for each country
  results_list <- list()
  
  # Unique country IDs in the data (assuming all dfs have the same set of countries)
  country_ids <- unique(crps_dfs[[1]]$country_id)
  
  # Loop over each country ID
  for (country_id in country_ids) {
    # Template for the months and years
    months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Year Average')
    years <- start_year:end_year
    
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


# PERFORM EXTRACTION FOR TRAINING DATA
for (i in 1:length(eval_months)){
  eval_month = eval_months[i]
  eval_month_name = eval_months_names[i]
  train_id <- subset(data, month_id <= last_training_month, select= c(month_id, country_id, ged_sb))
  print(paste0('Loading fit of ', eval_month_name, '...'))
  load(paste0("stan_fits/", model_name, '/', 'fit_', model_name, '_', eval_month_name, '2018_composed.RData'))
  print('Extracting CRPS scores...')
  pps_train <- create_pps_dataframe(get(paste0('fit_', eval_month_name)), train_id, "y_pred_train")
  pps_train_list[[i]] <- pps_train
  crps_train_dfs[[i]] <- calculate_crps(pps_train)
  print(paste0('CRPS for ', eval_month_name, ' has been ectracted. Remove fit object...'))
  rm(list = paste0('fit_', eval_month_name))
}

crps_train_averages_df <- compute_crps_averages(crps_train_dfs, 1990, 2017, 121)
crps_train_averages_by_country_df <- compute_crps_averages_by_country(crps_train_dfs, 1990, 2017, 121)
concatenated_train_pps <- concat_pps_list(pps_train_list)

arrow::write_parquet(crps_train_averages_df, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_train_df_averages.parquet"))
arrow::write_parquet(crps_train_averages_by_country_df, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_train_averages_by_country_df.parquet"))
arrow::write_parquet(concatenated_train_pps, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_train_posterior_predicitve_samples.parquet"))


### PERFORM EXTRACTION FOR EVALUATION DATA
for (i in 1:length(eval_months)){
  eval_month = eval_months[i]
  eval_month_name = eval_months_names[i]
  print(paste0('Loading fit of ', eval_month_name, '...'))
  load(paste0("stan_fits/", model_name, '/', 'fit_', model_name, '_', eval_month_name, '2018_composed.RData'))
  print('Extracting CRPS scores...')
  eval_id <- get_eval_id(eval_month, data)
  pps_eval <- create_pps_dataframe(get(paste0('fit_', eval_month_name)), eval_id, 'y_pred_eval')
  pps_eval_list[[i]] <- pps_eval
  crps_eval_dfs[[i]] <- calculate_crps(pps_eval)
  print(paste0('CRPS for ', eval_month_name, ' has been ectracted. Remove fit object...'))
  rm(list = paste0('fit_', eval_month_name))
}

crps_eval_averages_df <- compute_crps_averages(crps_eval_dfs, 2018, 2021, 457)
crps_eval_averages_by_country_df <- compute_crps_averages_by_country(crps_eval_dfs, 2018, 2021, 457)
concatenated_eval_pps <- concat_pps_list(pps_eval_list)

arrow::write_parquet(crps_averages_df, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_eval_df_averages.parquet"))
arrow::write_parquet(crps_averages_by_country_df, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_eval_averages_by_country_df.parquet"))
arrow::write_parquet(concatenated_pps, paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_eval_posterior_predicitve_samples.parquet"))

# crps_eval_averages_df <- read_parquet(paste0("C:/Users/Uwe Drauz/Documents/bachelor_thesis_local/personal_competition_data/Results/", model_name, "_crps_df_averages.parquet"))
