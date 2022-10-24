######################################################################################
# SET DIRECTORIES                                                                    #
#                                                                                    #
# Set and get directories in one place in the name of consistency and ease.          #
# Creates any directories that do not currently exist.                               #
#                                                                                    #
# Outputs: A list of relevant directories (within pm$pth) which can be referenced    #
# elsewhere                                                                          #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch                                               #
######################################################################################

# Load the functions
source("myRfunctions.R")

# --------------------------------------------
# Define paths for project inputs and outputs.
# --------------------------------------------
set_dirs <- function(pm, type_analysis){
  
  # Initiate file path lists
  pth <- out <- list()
  
  # Define the sample number (for adaptive sampling)
  sample <- paste0("sample_", pm$sample_num)
  
  # ---- Code and resource locations ----
  
  # Base path to code repositories
  base_stm <- file.path("/scicore/home/penny/masthi00")
  
  # Parent paths to all input files relating to this project
  pth$code <- file.path(base_stm, "smc_resistance") # working directory
  input <- file.path(pth$code, "SIM_FOLDER") # file with input in the working folder
  
  # ---- Select the good base XML files based on analysis----
  
  # Base file if aimed to estimate the prophylactic period
  if (type_analysis == "Propylaxis") {
    pth$xml_base <- file.path(input, "Propylaxis_new_dosage_SP_only.xml")
  }
  
  # Base file if aimed to estimate the rate of spread of SP resistant parasites
  if (type_analysis == "SMC") {
    pth$xml_base <- file.path(input, "SMC_SP_fourier_new_dosage.xml")
  }
  
  # Base file if aimed to estimate the rate of spread of SP resistant parasites with random coverage
  if (type_analysis == "SMC" & pm$opts$Type_coverage=="Random") {
    pth$xml_base <- file.path(input, "SMC_SP_fourier_new_dosage_random.xml")
  }
  
  # Base file if aimed to estimate the efficacy of SMC
  if (type_analysis == "Efficacy") {
    pth$xml_base <- file.path(input, "SMC_efficacy.xml")
  }
  
  # Base file if aimed to estimate the efficacy of SMC with random coverage
  if (type_analysis == "Efficacy" & pm$opts$Type_coverage=="Random") {
    pth$xml_base <- file.path(input, "SMC_efficacy_random.xml")
  }
  
  # Base file if aimed to replicate the trials of Zongo
  if (type_analysis == "Trial") {
    pth$xml_base <- file.path(input, "Zongo_4.xml")
  }
  
  # Path to output form scicore (if simulation run or not)
  pth$log_files <- file.path(pth$code, "JOB_OUT")
  
  # Path to OpenMalaria resource files
  pth$om_files <- file.path(base_stm, "OM_schema43_0")
  
  # ---- Output directories and files ----
  
  # Parent path to all output files relating to this project
  file_stem <- file.path(base_stm, "OUT_SMC", "ZONGO_FINAL_50_seeds")
  
  # Path to parameter table
  out$param_table <- file.path(file_stem, "0_parameters")
  
  # Paths to  simulation folders
  out$sim_files <- file.path(file_stem, "1_sim_files")
  
  # Paths to seed files
  out$sim_sed <- file.path(out$sim_files, "sed_files", sample)
  
  # Paths to xml files
  out$sim_xml <- file.path(out$sim_files, "xml_files", sample)
  
  # Paths to Outputs folders
  sim_output <- file.path(file_stem, "2_sim_outputs")
  
  # Path to OpenMalaria raw outputs
  out$sim_out <- file.path(sim_output, "raw_output", sample)
  
  # Path to Post processed Outputs
  out$processed <- file.path(sim_output, "processed", sample)
  
  # Path to summary outputs
  out$summarised <- file.path(sim_output, "summarised")
  
  # Paths to  GP fitted during every round of adaptative sampling
  out$gp_samples <- file.path(file_stem, "3_gp_models")
  
  # Path to the GP fitted during the last round of adaptative sampling (use for the sensitivity analysis)
  out$gp_models <- file.path(file_stem, "3_gp_models", "models")
  
  # Paths to sensitivity analysis
  out$results <- file.path(file_stem, "4_results")
  out$results_gp <- file.path(file_stem, "4_results")
  out$results_optim <- file.path(file_stem, "4_results")
  
  # Make all output directories
  make_out_dirs(out)  # see bellow
  
  # Append paths to pm list
  pm <- append_dirs(pm, pth, out)  # see bellow
  
  # Return pm
  return(pm)
}

# ---------------------------------------------------------
# Make all output directories if they do not already exist.
# ---------------------------------------------------------
make_out_dirs <- function(out) {
  
  # Extract all path names in list
  pth_names <- names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth <- out[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth)) {
      dir.create(this_pth, recursive = TRUE)
    }
    
  # Close directory loop  
  } 
}

# ---------------------------------------------------------
# Concatenate separators and append directories to pm list.
# ---------------------------------------------------------
append_dirs <- function(pm, pth, out) {
  
  # Extract all path names in list
  pth_names <- names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth <- out[[pth_name]]
    
    # Add a file separator to end of output paths
    out[[pth_name]] <- paste0(this_pth, file_sep())
  }
  
  # Concatenate lists
  pm$pth <- c(pth, out)
  
  return(pm)
}
