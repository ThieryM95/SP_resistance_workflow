######################################################################################
# Option.R                                                                           #
#                                                                                    #
# Task: Set key simulation options: analysis type, parameter spaces                  #
# number of simulations and seeds, directories                                       #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch                                               #                                      
######################################################################################

# Download function in the environments
source("myRfunctions.R")
require(plyr)
require(dplyr)
require(xlsx)

# ------------------------------------------
# Set options for a locally called analysis.
# ------------------------------------------
set_options <-function(do_step = NA, sample_num = 0, quiet = FALSE) {
    
    #----general information----
    pm <- list()
    
    # Define analyse name
    pm$analysis_name <- "SMC"
    
    # Define the user name
    pm$user <- "masthi00"
    
    # Iteration
    pm$sample_num <- sample_num
    
    # ---- Sampling options ----
    
    opts <- list()
    
    # Number of latin hypercube samples to generate (prior to GP)
    opts$lhc_samples <- 1
    
    # Number of different seeds to simulate
    opts$n_seeds <- 50
    
    # Flag for sampling EIR values
    opts$do_resample <- TRUE
    
    # ---- Post processing  ----
    
    # Define the type of analysis
    opts$Type_analysis <- "Trial" # option are: Prophylaxis or SMC or Efficacy or Trial
    
    # Define if same or different children recieve SMC at different round
    opts$Type_coverage <- "Fixed" # option are : Random/Fixed (NB: Random is only for SMC or efficacy or Trial)
    
    # Define the outcome
    opts$om_outcome <- "Spread"
    
    # ---- Gaussian process  ----
    
    # Select GP kernel function
    opts$gp_kernel <- "Gaussian" # "Gaussian", "Matern5_2", or "Matern3_2"
    
    # Proportion of data withheld from training set
    opts$gp_test_split <- 0.2 
    
    # Maximum number of iteration of optimization algorithm
    opts$gp_max_iter <- 10000
    
    # ---- Adaptive sampling  ----
    
    # Maximum number of adaptive sampling attempts
    opts$sampling_max_iter <- 10
    
    # Number of points to re-simulate per arm
    opts$n_adaptive_samples <- 100
    
    # Quantitative threshold for accepting GP performance
    opts$stop_criteria <- 0.99 
    
    # ---- Option for sensitivity analysis
    
    # Number of sampling points
    opts$sa_n <- 10000
    
    # Create output directories
    pm$opts <- opts
    
    # Specify parameter range and which one are constrained
    pm <- constrained_parameter(pm, opts$Type_analysis) # see function bellow
    pm <- variable_parameter(pm, opts$Type_analysis) # see function bellow
    
    # Define the directory
    pm <- set_dirs(pm, pm$opts$Type_analysis) # check directories.R
    
    # ---- Display key options and return ----
    
    # Number of scenarios defined - just to inform the user
    n_param <- 0
    
    for (i in 1:length(pm$settings)) {
      
      n_param[i] <- length(pm$settings[[i]])
      
    }
    
    # Define the number of setting  
    n_settings <- prod(n_param)
    
    # Estimate the number of scenario  
    pm$opts$n_total_jobs <- pm$opts$lhc_samples * n_settings * pm$opts$n_seeds
    
    # Only display if flag is on
    if (quiet == FALSE) {
        
      message(" - Total number of arms: ",
            format(n_settings, big.mark = ","))
      
      message(" - Total number of simulations: ",
            format(pm$opts$n_total_jobs, big.mark = ","))
    }
   
   # Return the output   
   return(pm)
}

# ------------------------------------------------------
# Define the parameter values of constrained parameters.
# ------------------------------------------------------
constrained_parameter <- function(pm, type_analysis) {
  
  # Create a list with all the arm to simulate
  settings <- list()
  
  # Define the different seasonality profile via fourrier coefficient
  sesonality1 <- c(0, 0, 0, 0)
  sesonality2 <-c(-1.1688319842020671, -1.9682456652323406, -0.23417717218399048, 0.3833462397257487)
  sesonality3 <-c(-1.0296205679575603, -1.7833550771077473, -0.7692119280497233, 1.332314173380534) # 0.07+10 test 10
  
  # Setting when aimed to estimate the prophylactic period
  if (type_analysis == "Propylaxis") {
    
    # Define level of EIR
    EIR <- c(5,  50, 100, 150, 500)
    
    # Seasonality pattern
    Seasonality <- c("sesonality1")
    
    # Level of access to treatment
    Access <- c(0)
    
    # Diagnostic detection limit
    Diagnostic <- c(20)
    
    # Merge all the information
    setting_data <-Reduce(merge,
               list(
               as.data.frame(EIR),
               as.data.frame(Access),
               as.data.frame(Diagnostic),
               as.data.frame(Seasonality)))
    
    # Give name to the columns
    colnames(setting_data) <- c("eir", "Access", "Diagnostic", "seasonality")
    
    # Save the information into variable setting
    settings$eir <- unique(setting_data$eir)
    settings$Access <- unique(setting_data$Access)
    settings$Diagnostic <- unique(setting_data$Diagnostic)
    settings$seasonality <- data.frame(sesonality1)
  }
  
  # Setting when aimed to estimate the selection coefficient
  if (type_analysis == "SMC") {
    
    # Seasonality profile
    Seasonality <- c("sesonality2", "sesonality3")
    
    # Reduction of coverage at each round of adaptive sampling
    Coverage_reduction <- c(0, 0.1)
    
    # Number of round of adaptive sampling (NB: 5.4 = 4 + 1 before, 4.5= 4+ 1 after)
    Number_round <- c(4, 5.4, 4.5)
    
    # Maximum age targeted by SMC
    Age <- c(5, 10)
    
    # EC50 of the resistant genotype
    IC50_SP_R <- c(24.20)
    
    # EC50 of the sensitive genotype
    IC50_SP <- c(2.39)
    
    # Merge all the information
    setting_data <- Reduce( merge,
        list(
          as.data.frame(IC50_SP_R),
          as.data.frame(IC50_SP),
          as.data.frame(Age),
          as.data.frame(Number_round),
          as.data.frame(Coverage_reduction),
          as.data.frame(Seasonality)))
    
    # Name the columns
    colnames(setting_data) <-c(
        "IC50_SP_R",
        "IC50_SP",
        "Age",
        "Number_round",
        "Coverage_reduction",
        "seasonality")
    
    # Save the information into variable setting
    settings$IC50_SP_R <- unique(setting_data$IC50_SP_R)
    settings$IC50_SP <- unique(setting_data$IC50_SP)
    settings$Age <- unique(Age)
    settings$Number_round <- unique(setting_data$Number_round)
    settings$Coverage_reduction <- unique(setting_data$Coverage_reduction)
    settings$seasonality <- data.frame(sesonality2, sesonality3)
    
  }
  
  # Setting when aimed to assessed the efficacy of SMC
  if (type_analysis == "Efficacy") {
    
    # Seasonality Profile
    Seasonality <- c("sesonality2")
    
    # Reduction of coverage at each round of adaptive sampling
    Coverage_reduction <- c(0.1)
    
    # Number of round of adaptive sampling (NB: 5.4 = 4 + 1 before, 4.5= 4+ 1 after)
    Number_round <- c(4)
    
    # Maximum age targeted by SMC
    Age <- c(5)
    
    # EC50 of the sensitive genotype
    IC50_SP <- c(60.12)
    
    # Merge all the information
    setting_data <-
      Reduce(merge,
        list(
          as.data.frame(IC50_SP),
          as.data.frame(Age),
          as.data.frame(Number_round),
          as.data.frame(Coverage_reduction),
          as.data.frame(Seasonality)))
    
    # Name the columns
    colnames(setting_data) <- c(
        "IC50_SP",
        "Age",
        "Number_round",
        "Coverage_reduction",
        "seasonality")
    
    # Save the information into variable setting
    settings$IC50_SP_R <- unique(setting_data$IC50_SP_R)
    settings$IC50_SP <- unique(setting_data$IC50_SP)
    settings$Age <- unique(Age)
    settings$Number_round <- unique(setting_data$Number_round)
    settings$Coverage_reduction <-unique(setting_data$Coverage_reduction)
    settings$seasonality <- data.frame(sesonality2)
  
  }
  
  # Setting when aimed to replicate the trial of Zongo et al (2015)
  if (type_analysis == "Trial") {
    
    # EC50 of the resistant genotype
    IC50_SP_R <- c(0.5)
    
    # EC50 of the sensitive genotype
    IC50_SP <- c(2.39, 0.5)
    
    # SMC coverage
    Coverage <- c(0, 1) #0 control group, 1 trial group
    
    # Merge all the information
    setting_data <- Reduce(merge,
             list(
               as.data.frame(IC50_SP_R),
               as.data.frame(IC50_SP),
               as.data.frame(Coverage)))
    
    # Name the columns
    colnames(setting_data) <- c("IC50_SP_R", "IC50_SP", "Coverage")
    
    # Save the information into variable setting
    settings$IC50_SP_R <- unique(setting_data$IC50_SP_R)
    settings$IC50_SP <- unique(setting_data$IC50_SP)
    settings$Coverage <- unique(setting_data$Coverage)
  }
  
  # Append settings to pm list
  pm$settings <- settings
  
  # Return pm
  return(pm)
}

# -------------------------------------------------------------------
# Define the parameter space for parameters that are not constrained.
# -------------------------------------------------------------------
variable_parameter <- function(pm, type_analysis) {
  
  # Parameter space if aimed to estimate the prophylactic period
  if (type_analysis == "Propylaxis") {
    
    # Parameter name
    Parameter <- c("IC50_SP")
    
    # Maximum values
    max <- c(0.02) # 0.02, 0.3,100
    
    # Minimum values
    min <- c(0.0005) # 0.0005,0.002,0.01
  }
  
  # Parameter space if aimed to identify the key driver of SMC-resistance
  if (type_analysis == "SMC") {
    
    # Parameter names
    Parameter <- c("Coverage", "Access", "eir", "half_life_long", "Dosage_long") # NB: dosae long do not vary at the end
    
    # Maximum values
    max <- c(1, 0.5, 500, 21, 30)
    
    # Minimum values
    min <- c(0.7, 0.04, 5, 7, 30)
  }
  
  if (type_analysis == "Efficacy") {
  
    # Parameter names
    Parameter <- c("Coverage", "Access", "eir", "half_life_long", "Dosage_long") # NB: dosae long do not vary at the end
    
    # Maximum values
    max <- c(1, 0.04, 5, 15, 30)
    
    # Minimum values
    min <- c(0.9, 0.5, 500, 15, 30)
  }
  
  if (type_analysis == "Trial") {
    
    # Parameter names
    Parameter <- c("eir") # NB: dosage long do not vary at the end
    
    # Maximum values
    max <- c(350)
    
    # Minimum values
    min <- c(350)
  }
  
  # Merge the information into a dataframe
  program_data <- data.frame(Parameter, max, min)
  
  # Convert dataframe to list
  prog <- as.list(program_data)
  
  # Names
  prog$prog_names <- Parameter
  
  # Easy access number of programs
  prog$n_progs <- length(prog$prog_names)
  
  # Append program details to pm list
  pm$prog <- prog
  
  # Return pm
  return(pm)
}
