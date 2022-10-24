#######################################################################################
# Post process results of the global sensitivity analysis                             #
#                                                                                     #
# Task : Organised the result of the sensitivity analysis from a list to a table      #
# Input: Output of the sensitivity analysis in a list format                          #
# output: Output of the sensitivity analysis in a table format                        #
#                                                                                     #
# Author:  thiery.masserey@swisstph.ch                                                #                                      
#######################################################################################

#-----------------------------------------------------------------
# function to organize the Sobol indices   from a list to a table.
#-----------------------------------------------------------------
Post_process_sensitivity <- function(Results_SA) {
  
  # Select the parameter names
  parameters <- pm$prog$Parameter[pm$prog$Parameter != "Dosage_long"]
  
  # ---- Define the different variable of our table and their length ----
  
  # Setting names
  Setting_names <- rep(0,
                      length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                      length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # Parameter names
  factors <- rep(parameters,
                 length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                 length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # If first/ total order indices
  First <- rep(0,
              length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
              length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # Value of the Sobol indices
  Effect <- rep(0,
               length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
               length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # Maximum values of the Sobol indices
  MAX <- rep(0,
            length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
            length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # Minimum values of the Sobol indices
  MIN <- rep(0,
             length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
             length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # Define the value for the factors that were constrained
  IC50_SP_R <- rep(0,
                   length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                   length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  IC50_SP <- rep(0,
                 length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                 length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  AGE <- rep(0,
             length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
             length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  Coverage_reduction <- rep(0,
                            length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                            length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  Number_round <- rep(0,
                      length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                      length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  Seasonality <- rep(0,
                     length(pm$settings$IC50_SP_R) * length(pm$settings$IC50_SP),
                     length(pm$settings$Age) * length(pm$settings$Number_round) * length(parameters) * 2 * length(pm$settings$Coverage_reduction) * length(pm$settings$seasonality))
  
  # Merge all the information into a table
  data <- cbind(Setting_names,
                as.character(factors),
                First,
                Effect,
                MAX,
                MIN,
                IC50_SP_R,
                IC50_SP,
                AGE,
                Coverage_reduction,
                Number_round,
                Seasonality)
  
  colnames(data) <-c("Setting",
                     "Factor",
                     "First",
                     "Effect",
                     "MAX",
                     "MIN",
                     "IC50_SP_R",
                     "IC50_SP",
                     "AGE",
                     "Coverage_reduction",
                     "Number_round",
                      "Seasonality")
  
  #---- Restructure the table ----
  
  # Define the number of setting post processed for the loop
  n_settings <- 1
  
  # Do a loop across each arm
  for (this_IC50_SP_R in pm$settings$IC50_SP_R) {
    for (this_IC50_SP in pm$settings$IC50_SP) {
      for (this_Age in pm$settings$Age) {
        for (this_round in pm$settings$Number_round) {
          for (this_reduction in pm$settings$Coverage_reduction) {
            for (this_season in names(pm$settings$seasonality)) {
              
              # Define the arm name
              setting_name <- paste(pm$opts$om_outcome,
                                    this_season,
                                    "IC50_SP_R", this_IC50_SP_R,
                                    "IC50_SP", this_IC50_SP,
                                    "age", this_Age,
                                    "Coverage_reduction", this_reduction,
                                    "Number_round", this_round,
                                    sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
              # Load the data
              SA_file <-paste0(pm$pth$results, setting_name, "_gp.RData")
              SA <- readRDS(SA_file)
              
              # Other options:
              # SA<-Results_SA[setting_name]
              # SA<-list.ungroup(SA, level = 1L)
              
              # Transfer the  data from the sensitivity analysis output into the data frame
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 1] <- setting_name
              data[n_settings:(n_settings + length(parameters) - 1), 3] <- SA$S$original[1:length(parameters)]
              data[(n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1), 3] <- SA$T$original[1:length(parameters)]
              data[c(n_settings:(n_settings + length(parameters) - 1)), 4] <- "First"
              data[c((n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1)), 4] <- "Total"
              data[n_settings:(n_settings + length(parameters) - 1), 5] <- SA$S$`max. c.i.`[1:length(parameters)]
              data[(n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1), 5] <- SA$T$`max. c.i.`[1:length(parameters)]
              data[n_settings:(n_settings + length(parameters) - 1), 6] <- SA$S$`min. c.i.`[1:length(parameters)]
              data[(n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1), 6] <- SA$T$`min. c.i.`[1:length(parameters)]
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 7] <- this_IC50_SP_R
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 8] <- this_IC50_SP
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 9] <- this_Age
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 10] <- this_reduction
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 11] <- this_round
              data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 12] <- this_season
              
              # Update the number of setting we went across with the post processing
              n_settings <- n_settings + length(parameters) * 2
            }
          }
        }
      }
    }
  }
  
  # Finalize the data frame
  data <- as.data.frame(data)
  
  # Update the type of variable
  data$First <- as.numeric(levels(data$First))[data$First]
  data$MAX <- as.numeric(levels(data$MAX))[data$MAX]
  data$MIN <- as.numeric(levels(data$MIN))[data$MIN]
  data$First[data$First <= 0] <- 0
  
  # Save the data frame of Sobol indices for each arm and parameter
  name_file <- paste0(pm$pth$results, "Sobol_indices", ".txt")
  write.table(data, file = name_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # Return the  dataframe of sobol indices for each arm and parameter
  return(data)

}

#------------------------------------------- ---------------------------
# function to organise the effect of each factor from a list to a table.
#-----------------------------------------------------------------------
Post_process_sensitivity_2 <- function(Results_SA) {
  
  # Select the parameter name
  parameters <- pm$prog$Parameter[pm$prog$Parameter != "Dosage_long"]
  
  # Create a parameter table that will contains the estimate selection coefficient over the parameter range for each arm
  Quantil_final_final <- matrix(NA, nrow = 0, ncol = 5)
  colnames(Quantil_final_final) <- c("L", "M", "U", "x", "G")
  
  # Loop across each arm
  for (this_IC50_SP_R in pm$settings$IC50_SP_R) {
    for (this_IC50_SP in pm$settings$IC50_SP) {
      for (this_Age in pm$settings$Age) {
        for (this_round in pm$settings$Number_round) {
          for (this_reduction in pm$settings$Coverage_reduction) {
            for (this_season in names(pm$settings$seasonality)) {
              
              # Define the name of the arm
              setting_name <- paste(pm$opts$om_outcome,
                                    this_season,
                                    "IC50_SP_R", this_IC50_SP_R,
                                    "IC50_SP", this_IC50_SP,
                                    "age", this_Age,
                                    "Coverage_reduction", this_reduction,
                                    "Number_round", this_round,
                                    sep = "_")
              
              # Load the results of the sensitivity analysis for this arm
              SA_file <-paste0(pm$pth$results, setting_name, "_gp.RData")
              SA <- readRDS(SA_file)
              
              # other options:
              # SA<-Results_SA[setting_name]
              # SA<-list.ungroup(SA, level = 1L)
              
              # Select the predicted selection coefficient into a vector
              Y <- SA$y
              
              # Select the input parameter into a table
              X <- matrix(NA, ncol = length(SA$X), nrow = length(Y))
              colnames(X) <- colnames(SA$X)
              for (i in 1:length(SA$X)) {
                # i<-6
                X[, i] <- SA$X[, i]
              }
              
              # For each parameter, divide the parameter into N fraction, and estimate the mean selection coefficient in this fraction
              # Define the number of fraction
              N <- 40
              NN <- 40
              
              # Create a vector that will contain the boundary of fraction of the parameter space + median selection coefficient
              Quantil_final <- matrix(NA, nrow = 0, ncol = 5)
              
              # Loop across each parameters
              for (i in 1:length(SA$X)) {
                
                # Create a data frame of parameter value and  selection coefficient
                DF <- as.data.frame(cbind(X[, i], Y))
                colnames(DF)[1:2] <- c("X", "Y")
                
                # Cut the parameter space into N piece
                splitFitness <- seq(min(X[, i]), max(X[, i]), length = N)
                
                # Create a table for this specific parameter
                Quantiles_1 <- matrix(NA, nrow = N, ncol = 3)
                
                # Estimate the median selection coefficient in each fraction of the parameter space
                for (k1 in 1:(NN - 1)) {
                  thisXY <- DF[X[, i] > splitFitness[k1] & X[, i] < splitFitness[k1 + 1], ]
                  Quantiles_1[k1, ] <- quantile(thisXY$Y, c(0.25, 0.5, 0.75))
                }
                
                # Transform data into data frame
                Quantiles_1 <- as.data.frame(cbind(Quantiles_1, 1:40))
                
                # Select the parameter name
                Quantiles_1$G <- colnames(X)[i]
                
                # Add columns name
                colnames(Quantiles_1) <- c("L", "M", "U", "x", "G")
                
                # Update the table to contain the data of all parameters
                Quantil_final <- rbind(Quantil_final, Quantiles_1)
              }
              
              # Add information about the arm
              Quantil_final$IC50_SP_R <- this_IC50_SP_R
              Quantil_final$IC50_SP <- this_IC50_SP
              Quantil_final$Age <- this_Age
              Quantil_final$Number_round <- this_round
              Quantil_final$Coverage_reduction <- this_reduction
              Quantil_final$Seasonality <- this_season
              
              # Update the table to contain the data of all arm
              Quantil_final_final <- rbind(Quantil_final_final, Quantil_final)
            }
          }
        }
      }
    }
  }
  
  # Save the table of mean selection coefficient over the parameter range of each parameter, for each arm
  name_file <- paste0(pm$pth$results, "Factors_effect", ".txt")
  write.table(Quantil_final_final, file = name_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # Return the table of mean selection coefficient over the parameter range of each parameter, for each arm
  return(Quantil_final_final)
}
