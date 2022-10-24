######################################################################################
# GAUSSIAN PROCESS                                                                   #
#                                                                                    #
# Task: Fit a GP to interpolate whole parameter space using OM simulations.          #
# Input: Summary table of parameter values and indicator of interest                 #
# Output: GP fitted to each constrained arms                                         #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch                                               #                                      
######################################################################################

# Load function
require(hetGP)

# ---------------------------------------------------------
# Parent function to perform all GP steps.
# ---------------------------------------------------------
run_gp <- function(pm) {
  
  # Create a list elements
  indices <- list()
  
  # Message in the consoles
  message("* Fitting GP models")
  
  # Load the full summary data
  param_table_name <- param_set_name(all_samples = TRUE)
  summary_file <- paste0(param_table_name, "_", ".txt")
  summary_path <- paste0(pm$pth$summarised, summary_file)
  indicator_table <- read.table(summary_path, sep = "\t", header = TRUE, as.is = TRUE)
  
  # Loop through all the different arms
  for (this_IC50_SP_R in pm$settings$IC50_SP_R) {
    for (this_IC50_SP in pm$settings$IC50_SP) {
      for (this_Age in pm$settings$Age) {
        for (this_round in pm$settings$Number_round) {
          for (this_reduction in pm$settings$Coverage_reduction) {
            for (this_season in names(pm$settings$seasonality)) {
              
              # ---- select the data specific to each arm -----
              
              # Construct file name based on arm values when constrained analysis
              setting_name <- paste(pm$opts$om_outcome,
                                    this_season,
                                    "IC50_SP_R", this_IC50_SP_R,
                                    "IC50_SP", this_IC50_SP,
                                    "age", this_Age,
                                    "Coverage_reduction", this_reduction,
                                    "Number_round", this_round,
                                    sep = "_") 
              
              sample_name <-paste(setting_name, "sample", pm$sample_num, sep = "_")
              
              # Reduce the table to only consider this arm
              setting_idx <- indicator_table$season == this_season &
                indicator_table$IC50_SP_R == this_IC50_SP_R &
                indicator_table$IC50_SP == this_IC50_SP &
                indicator_table$Age == this_Age &
                indicator_table$Coverage_reduction == this_reduction &
                indicator_table$Number_round == this_round
              
              # Message to console to know which arm is done
              message("  - ", setting_name)
              
              # Select the scenario of the specific arms
              setting_table <- indicator_table[setting_idx,]
              
              # Remove annoying row numbers
              rownames(setting_table) <- NULL
              
              # Save all the scenario run in this arms
              setting_file <-paste0(pm$pth$gp_samples, sample_name, ".txt")
              write.table(setting_table, setting_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) # save the table
              
              # ---- Train the model ----
              
              # Select input variables of interest
              all_progs <- pm$prog$prog_names ## select only parameter and indicators
              input_data <- setting_table
              
              # Train a GP model on this data using a cross validation approach
              gp_result <- split_train(pm, input_data, 0.2) # see function bellow
              
              # Save the fitted GP model to file - adaptive sample iteration
              gp_file <- paste0(pm$pth$gp_samples, sample_name, "_gp.RData")
              saveRDS(gp_result, file = gp_file)
              
              # Save the fitted GP model to file - latest fitted model
              gp_file <- paste0(pm$pth$gp_models, setting_name, "_gp.RData")
              saveRDS(gp_result, file = gp_file)
              
              # Save the fitted GP in a list to return to the user
              indices$sample_name <- gp_result
              names(indices)[names(indices) == "sample_name"] <- setting_name
            }
          }
        }
      }
    }
  }
  
  # Return indices
  return(indices)
}


# ---------------------------------------------------------
# Simple split of train-test data to train the GP.
# ---------------------------------------------------------
split_train <- function(pm, input_data, p) {
  
  # Define the parameter name
  all_progs <- pm$prog$prog_names
  
  # Select all scenario with 1 seed
  input_data_1 <- input_data[input_data$seed == 1,]
  
  # Calculate the total number of scenario
  n_samples <- nrow(input_data_1)
  
  # Randomly sample the scenario to be part of the test and training data set
  n_test <- round(n_samples * p)
  test_idx <- sort(sample(1:n_samples, n_test))
  train_idx <- setdiff(1:n_samples, test_idx)
  
  # Use the indices to separate the data
  train_data <- input_data_1[train_idx,]
  test_data <- input_data_1[test_idx,]
  
  # Select all the other seed to be part of the good data set
  train_data <- input_data[input_data$Access %in% train_data$Access & input_data$eir %in% train_data$eir,]
  TEST_DATA <- input_data[input_data$Access %in% test_data$Access & input_data$eir %in% test_data$eir,]
  
  # Select the  parameter columns (except dosage)+ indicators columns only
  train_data <- train_data[, c(all_progs[all_progs != "Dosage_long"], "Indicator")]
  TEST_DATA <- TEST_DATA[, c(all_progs[all_progs != "Dosage_long"], "Indicator")]
  
  # Remove NA in the trained data set
  train_data <- train_data[rowSums(is.na(train_data)) < 1,]
  
  # For the test data set calculate the mean selection coefficient (indicator) across seed
  test_data_2 <- test_data[FALSE, c(all_progs[all_progs != "Dosage_long"], "Indicator")]
  for (i in 1:length(test_data$Access)) {
    test_data_2[i,] <- colMeans(TEST_DATA[test_data$Access[i] == TEST_DATA$Access & test_data$half_life_long[i] == TEST_DATA$half_life_long,], na.rm = TRUE)
  }
  
  # Remove NA in the trained data set
  test_data_2 <- test_data_2[rowSums(is.na(test_data_2)) < 1,]
  
  # Train model using these training data
  trained_model <- gp_train(train_data) # see function bellow
  
  # Assess the model with both testing and training data
  gp_output <- gp_test(trained_model, train_data, test_data_2) # see function bellow
  
  # Append  to this output
  gp_output$gp_model <- trained_model
  
  return(gp_output)
}


# -----------------------------------------------------------
# Train a GP regression with hetGP given a training data set.
# -----------------------------------------------------------
gp_train <- function(train_data) {
  
  # Define indices of parameters and response
  n_params <- ncol(train_data) - 1
  param_idx <- 1:n_params
  response_idx <- n_params + 1
  
  # Find and remove duplicates
  prep_data <- find_reps(
    X = as.matrix(train_data[, param_idx]),
    Z = as.matrix(train_data[, response_idx]),
    rescale = FALSE,
    normalize = FALSE)
  
  # Prepare input for mleHetGP model
  X <- list(
    X0 = as.matrix(prep_data$X0),
    Z0 = as.matrix(prep_data$Z0),
    mult = prep_data$mult)
  
  # Run GP model
  gp_model <- mleHetGP(X = X, Z = prep_data$Z, covtype = "Matern5_2")
  
  return(gp_model)
}

# ---------------------------------------------------------
# Run GP with test data to determine model error.
# ---------------------------------------------------------
gp_test <- function(gp_model, train_data, test_data_2) {
  
  # Define indices of parameters and response
  n_params <- ncol(train_data) - 1
  param_idx <- 1:n_params
  response_idx <- n_params + 1
  
  # Apply the model on the test data
  predict_test <- predict(x = as.matrix(test_data_2[, param_idx]), object = gp_model)
  predict_test <- predict_test$mean
  
  # The actual response associated with the test data
  actual_test <- test_data_2[, response_idx]
  
  # Repeat this for the training data
  predict_train <- predict(x = as.matrix(train_data[, param_idx]), object = gp_model)
  predict_train <- predict_train$mean
  
  # The actual response associated with the test data
  actual_train <- train_data[, response_idx]
  
  # Concatenate output into list
  gp_output <- list(
    actual_test = actual_test,
    predict_test = predict_test,
    actual_train = actual_train,
    predict_train = predict_train)
  
  # return the prediction and true value
  return(gp_output)
}
