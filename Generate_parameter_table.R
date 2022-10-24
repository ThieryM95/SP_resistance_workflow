######################################################################################
#      Generate the parameter tables for the simulation that will be run             #
#                                                                                    #
#                                                                                    #
# Task: This function create the table of input parameter for each simulation        #
#                                                                                    #  
# Input: Number of simulations, number of seeds, constrained parameter values,       #
#        parameter range                                                             #
#                                                                                    #
# Output: A table of the parameter value of each parameter for each simulation       #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch                                               #
######################################################################################

# --------------------------------------------
# Generate the input table for the simulation.
# --------------------------------------------
generate_param_table <- function(pm, param_table = NULL) {
  
  # Message in console
  message("* Generating parameter table")
  
  # File name for this parameter table (see adaptive_sampling.R)
  param_table_name <- param_set_name(sample_num = pm$sample_num)
  param_file <- paste0(pm$pth$param_table, param_table_name, ".txt")
  
  # If sample number is == 0 the parameter table is generate via Latin hyper cube sampling
  if (pm$sample_num == 0) {
    
    # Continue if file doesn't exist OR we want to overwrite
    if (!file.exists(param_file) | pm$opts$do_resample == TRUE) {
      
      # ---- Generate Latin Hyper cube sampling for intervention coverage ----
      
      # Load the define parameter  ranges into a matrix
      cov_ranges <- matrix(c(pm$prog$min, pm$prog$max), ncol = 2)
      
      # Generate Latin hyper cube sample between the coverage limits
      param_table <- as.data.frame(lhs(pm$opts$lhc_samples, cov_ranges))
      
      # Variable names of parameters must match @parameter@ fields in base xml
      colnames(param_table) <- paste(pm$prog$prog_names)
      
      # ---- Define the parameter that were constrained (the different arm) ----
      
      # Parameter range when we aimed to estimate the prophylactic period
      if (pm$opts$Type_analysis == "Propylaxis") {
        
        # Define the level and value for each parameter
        eir <- pm$settings$eir
        season <- pm$settings$seasonality
        Access <- pm$settings$Access
        Diagnostic <- pm$settings$Diagnostic
        
        # Merge into full factorial table
        setting_table <-Reduce(merge,
                 list(
                   names(season),
                   as.data.frame(eir),
                   as.data.frame(Access),
                   as.data.frame(Diagnostic)))
        
        # Give name to columns
        colnames(setting_table) <- c("seasonality", "eir", "Access", "Diagnostic") 
      }
      
      # Parameter range when we aimed to estimate the rate of spread of resistant genotype
      if (pm$opts$Type_analysis == "SMC") {
        
        # Define the level and value for each parameter
        IC50_SP_R <- pm$settings$IC50_SP_R
        IC50_SP <- pm$settings$IC50_SP
        Age <- pm$settings$Age
        Number_round <- pm$settings$Number_round
        Coverage_reduction <- pm$settings$Coverage_reduction
        season <- pm$settings$seasonality
        
        # Merge into full factorial table
        setting_table <-Reduce(merge,
            list(
              names(season),
              as.data.frame(IC50_SP_R),
              as.data.frame(IC50_SP),
              as.data.frame(Age),
              as.data.frame(Number_round),
              as.data.frame(Coverage_reduction)))
        
        # Give name to columns
        colnames(setting_table) <- c("seasonality", "IC50_SP_R", "IC50_SP", "Age", "Number_round", "Coverage_reduction")
      }
      
      # Parameter range when we aimed to estimate the effectiveness of SMC resistant parasite
      if (pm$opts$Type_analysis == "Efficacy") {
        
        # Define the level and value for each parameter
        IC50_SP <- pm$settings$IC50_SP
        Age <- pm$settings$Age
        Number_round <- pm$settings$Number_round
        Coverage_reduction <- pm$settings$Coverage_reduction
        season <- pm$settings$seasonality
        
        # Merge into full factorial table
        setting_table <- Reduce(merge,
            list(
              names(season),
              as.data.frame(IC50_SP),
              as.data.frame(Age),
              as.data.frame(Number_round),
              as.data.frame(Coverage_reduction)))
        
        # Give name to columns
        colnames(setting_table) <- c("seasonality", "IC50_SP", "Age", "Number_round", "Coverage_reduction")
      }
      
      # Parameter range when we aimed to replicate Zongo trail
      if (pm$opts$Type_analysis == "Trial") {
        
        # Define the level and value for each parameter
        IC50_SP_R <- pm$settings$IC50_SP_R
        IC50_SP <- pm$settings$IC50_SP
        Coverage <- pm$settings$Coverage
        
        # Merge into full factorial table
        setting_table <- Reduce(merge,
                 list(
                   as.data.frame(IC50_SP_R),
                   as.data.frame(IC50_SP),
                   as.data.frame(Coverage)))
        
        # Give name to columns
        colnames(setting_table) <- c("IC50_SP_R", "IC50_SP", "Coverage")
      
      }
      
      # Merge the setting and parameter table
      param_table <- merge(setting_table, param_table)
      
      # Incorporate a scenario name column
      scenario_name <- paste("scenario", 1:nrow(param_table), sep = "_")
      param_table <- cbind(scenario_name, param_table)
      
    }
    
    # Repeat the scenarios for a predefined number of seeds
    seed <- 1:pm$opts$n_seeds
    param_table <- merge(param_table, as.data.frame(seed))
    
    # Write table to file
    write.table(
      param_table,
      file = param_file,
      sep = "\t",
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE)
  
  # Case where sample number is greater than 0 the parameter table is generate via adaptive sampling    
  } else {
    
    # Write extended parameter table to a new file
    write.table(
      param_table,
      file = param_file,
      sep = "\t",
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE)
    
  }
  
}
