###################################################################################
# Perform the sensitivity analysis                                                #
#                                                                                 #
# Task : Perform the global sensitivity analysis of the GP                        #
# Input: fitted GP                                                                #
# Output: Sobol indices + direction of effect                                     #
#                                                                                 #
# authors: thiery.masserey@swisstph.ch                                            #                   
###################################################################################


# Load the package
library("sensitivity")
library("multisensi")
library("lhs")
library("hetGP")
library("tgp")
library("dplyr")
library("ggplot2")
library("reshape2")
library("gridExtra")
library("hrbrthemes")
library("grid")
library("gridExtra")
library("rlist")

#--------------------------------------------------------------------
# function to prepare data for sensitivity analysis and save results.
#--------------------------------------------------------------------
sensitivity_analysis <- function(pm) {
  
  # Create a list to store the results of the SA
  indices <- list()
  
  # Define the parameter range
  param_ranges <- cbind(pm$prog$min, pm$prog$max)
  
  # Define the name of the parameter
  row.names(param_ranges) <- pm$prog$prog_names
  
  # Define the column names
  colnames(param_ranges) <- c("min", "max")
  param_ranges <- param_ranges[rownames(param_ranges) != "Dosage_long",]
  
  # Defined the number of sampled points
  sa_n <- pm$opts$sa_n
  
  # Loop across each arms
  for (this_IC50_SP_R in pm$settings$IC50_SP_R) {
    for (this_IC50_SP in pm$settings$IC50_SP) {
      for (this_Age in pm$settings$Age) {
        for (this_round in pm$settings$Number_round) {
          for (this_reduction in pm$settings$Coverage_reduction) {
            for (this_season in names(pm$settings$seasonality)) {
              
              # Define the name of the arm
              setting_name <- paste( pm$opts$om_outcome,
                                     this_season,
                                     "IC50_SP_R", this_IC50_SP_R,
                                     "IC50_SP", this_IC50_SP,
                                     "age", this_Age,
                                     "Coverage_reduction",this_reduction, 
                                     "Number_round", this_round,
                                     sep = "_") 
              
              # Message to the console
              message("  - ", setting_name)
              
              # Load the GP of the define arm
              gp_file <- paste0(pm$pth$gp_models, setting_name, "_gp.RData")
              trained_model <- readRDS(gp_file)
              trained_model <- trained_model$gp_model
              
              # Define parameter ranges
              param_ranges_2 <- param_ranges
              
              # Do the sensitivity analysis
              sobol_indices <- calc_sobol_idx(trained_model, param_ranges_2, sa_n) # see function bellow
              
              # Save the results in a data frame
              SA_file <- paste0(pm$pth$results, setting_name, "_gp.RData")
              saveRDS(sobol_indices, file = SA_file)
              
              # Save the results as an output of this function
              indices$setting_name <- sobol_indices
              names(indices)[names(indices) == "setting_name"] <- setting_name
            }
          }
        }
      }
    }
  }
  
  # Return the Sobol indices
  return(indices)
  
}

#----------------------------------------------
# function to perform the sensitivity analysis.
#----------------------------------------------
calc_sobol_idx <- function(trained_model, param_ranges_2, sa_n) {
  
  # Wrapper function for the GP_model prediction
  GP_f <- function(X) {
    out <- predict(x = as.matrix(X), object = trained_model)
    return(out$mean)
  }
  
  # Construct the two random lhs samples
  X1 <- as.data.frame(lhs(sa_n, as.matrix(param_ranges_2)))
  X2 <- as.data.frame(lhs(sa_n, as.matrix(param_ranges_2)))
  
  # Define the columns name
  colnames(X1) <- pm$prog$Parameter[pm$prog$Parameter != "Dosage_long"]
  colnames(X2) <- pm$prog$Parameter[pm$prog$Parameter != "Dosage_long"]
  
  # Estimate the Sobol indices
  SA <- soboljansen(model = GP_f,
                as.data.frame(X1),
                as.data.frame(X2),
                nboot = 1500)
  
  # Return the sobole indices + direction of effect
  return(SA)
}
