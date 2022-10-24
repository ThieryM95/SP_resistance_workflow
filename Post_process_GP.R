#######################################################################
# Post process of the GP prediction for a better visualization        #
#                                                                     #
# Task : Organised the gp accuracy result from a list to a data frame #
# Input: GP prediction and true value from  the train data set        #
# Output: Table of the GP test and prediction for each arm            #
#                                                                     #
# Author:  thiery.masserey@swisstph.ch                                #
#######################################################################


# Load package
library(SummarizedExperiment)

#--------------------------------------------------------------------------
# function to Reorganise the gp accuracy result from a list to a dataframe.
#--------------------------------------------------------------------------
Post_process_GP <- function(Results_gp) {
  
  # Create the data frame
  precision_final = precision_3 = precision_2 <- matrix(NA, nrow = 0, ncol = 9)
  
  # Add the columns name
  colnames(precision_2)[1:9] <- c("Test_True",
                                  "Test_predicted",
                                  "iteration",
                                  "seasonality",
                                  "IC50_SP_R",
                                  "IC50_SP",
                                  "Age",
                                  "Number_round",
                                  "Coverage_reduction")
  
  # Loop across the different iteration and arms
  for (iter in 1:length(Results_gp)) {
   z <- Results_gp[[iter]] # select data of the ith iteration
    for (this_IC50_SP_R in pm$settings$IC50_SP_R) {
      for (this_IC50_SP in pm$settings$IC50_SP) {
        for (this_Age in pm$settings$Age) {
          for (this_round in pm$settings$Number_round) {
            for (this_reduction in pm$settings$Coverage_reduction) {
              for (this_season in names(pm$settings$seasonality)) {
                
                # Define the arm name
                setting_name = paste(pm$opts$om_outcome,
                                     this_season,
                                     "IC50_SP_R", this_IC50_SP_R,
                                     "IC50_SP", this_IC50_SP,
                                     "age", this_Age,
                                     "Coverage_reduction", this_reduction,
                                     "Number_round", this_round,
                                      sep = "_") 
                
                # Select the data from the arm
                precision <- z[setting_name]
                
                # Change format
                precision <- list.ungroup(precision, level = 1L)
                
                # Select the variable true value and predicted value
                Test_True <- precision$actual_test
                Test_predicted <- precision$predict_test
                
                # Create the data set with all variables
                precision_2 <- data.frame(Test_True,
                                          Test_predicted,
                                          iter,
                                          this_season,
                                          this_IC50_SP_R,
                                          this_IC50_SP,
                                          this_Age,
                                          this_reduction,
                                          this_round)
                
                # Create a data set that contain all the data from all the arms
                precision_3 <- rbind(precision_3, precision_2)
                
              }
            }
          }
        }
      }
    }
    
  }
  
  # Define columns name
  colnames(precision_3)[1:9] <- c("Test_True",
                                  "Test_predicted",
                                  "iteration",
                                  "seasonality",
                                  "IC50_SP_R",
                                  "IC50_SP",
                                  "Age",
                                  "Coverage_reduction",
                                  "Number_round")
  
  # Save the data set
  name_file = paste0(pm$pth$gp_models, "Precision", ".txt")
  write.table(
    precision_3,
    file = name_file,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE)
  
  # Return the data set
  return(precision_final)
}
