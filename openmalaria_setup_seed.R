########################################################################################
# SIMULATION SETUP                                                                     #
#                                                                                      #
# TASK: For each simulations, create the seed pattern and xml file                     #
# INPUTS: Parameter table                                                              #
# OUTPUTS: Seed files and then XML file                                                #
# Written by thiery.masserey@swisstph.ch                                               #
########################################################################################

# Load functions required
source("myRfunctions.R")
require(pracma)
require(stringr)
require(tgp)
require(qdap)
require(xlsx)

# -------------------------------------------------------------------
# Generate XML seed replacement patterns for full country simulations.
# -------------------------------------------------------------------
simulation_sed_patterns <- function(pm) {
  
  # Message in the console
  message("  - Generating scenarios")
  
  # Load the table that contain the time of deployment of SMC for the different seasonality pattern and EIR level
  TIMING <- read.table("/scicore/home/penny/masthi00/smc_resistance/Data_time_2.txt", sep = "\t", header = T)
  
  # File path to parameter table
  param_table_name <- param_set_name(sample_num = pm$sample_num)
  param_file <- paste0(pm$pth$param_table, param_table_name, ".txt")
  
  # Read the parameter table generated by generate_param_table
  if (file.exists(param_file)) {
    param_table <-
      read.table(param_file,
                 sep = "\t",
                 header = TRUE,
                 as.is = TRUE)
  } else {
  
  # Throw an error if file cannot be found
  stop("Parameter table does not seem to exist")
  
  }
  
  # Load parameter name used for the analysis
  name_para <- colnames(param_table)
  
  # Name of parameter used as input parameter in OpenMalaria
  param_names <-  c(name_para[3:length(name_para)], "a_1", "b_1", "a_2", "b_2") # change
  
  # Name of parameter if aimed to estimate the spread of SP resistant genotype
    if (pm$opts$Type_analysis == "SMC") {
    param_names <- c(name_para[3:5], name_para[8:13], paste0("coverage_", 0:5), "a_1", "b_1", "a_2", "b_2", "time_0_1", "time_1_1", "time_2_1", "time_3_1", "time_4_1", "time_5_1", "DAY")
  }
  
  # Name of parameter if aimed to estimate the spread of SP resistant genotype with constant coverage
  if (pm$opts$Type_analysis == "SMC" & pm$opts$Type_coverage == "Random") {
    param_names <- c(name_para[3:5], name_para[9:13], paste0("coverage_", 0:5), "a_1", "b_1", "a_2", "b_2", "time_0_1", "time_1_1", "time_2_1", "time_3_1", "time_4_1", "time_5_1", "DAY")
  }
  
  # Name of parameter if aimed to estimate the efficacy of SMC
  if (pm$opts$Type_analysis == "Efficacy") {
    param_names <- c(name_para[3:4], name_para[7:12], paste0("coverage_", 0:5), "a_1", "b_1", "a_2", "b_2", "time_0_1", "time_1_1", "time_2_1", "time_3_1", "time_4_1", "time_5_1", "DAY")
  }
  
  # Name of parameter if aimed to estimate the efficacy of SMC with constant coverage
  if (pm$opts$Type_analysis == "Efficacy" &
      pm$opts$Type_coverage == "Random") {
    param_names <- c(name_para[3:4], name_para[8:12], paste0("coverage_", 0:5), "a_1", "b_1", "a_2", "b_2", "time_0_1", "time_1_1", "time_2_1", "time_3_1", "time_4_1", "time_5_1", "DAY")
  }
  
  # Name of parameter if aimed to replicate trial
  if (pm$opts$Type_analysis == "Trial") {
    param_names <- c(name_para[2:6])
  }
  
  # Name of parameter used as Output of Openmalaria for the Global sensitivity analysis
  param_names_2 <- c(name_para[3:length(name_para)], "seasonality") # change
  
  if (pm$opts$Type_analysis == "Trial") {
    param_names_2 <- param_names
  }
  
  # Create two lists of scenario
  # List of scenario without adjustment of their values (table used for the sensitivity analysis later)
  Liste_scenario <- data.frame(matrix(ncol = length(param_names_2) + 1, nrow = nrow(param_table)))
  colnames(Liste_scenario) <- c(param_names_2, "scenario_name")
  
  # Lists of scenario with parameter adjusted for the input of OpenMalaria
  Liste_scenario_adjusted <- data.frame(matrix(ncol = length(param_names) + 1,nrow = nrow(param_table)))
  colnames(Liste_scenario_adjusted) <- c(param_names, "scenario_name")
  
  # Initiate progress bar
  pb <- txtProgressBar(min = 0, max = nrow(param_table), initial = 0, width = 100, style = 3)
  
  # Loop through parameter table rows
  for (i in 1:nrow(param_table)) {
    
    # Select the ith scenario
    this_scen <- param_table[i,]
    
    # Define the seasonality pattern
    this_season <- this_scen$seasonality
    
    #---- Create the table used for the sensitivity analysis (without adjustment)----#
    param_values_2 <- c(this_scen[3:length(name_para)])
    
    
    # ---- Adjust the parameter for the input table ---- #
    
    # Adjustment needed when we estimate the selection coefficient of SMC-resistant parasite
    if (pm$opts$Type_analysis == "SMC" | pm$opts$Type_analysis == "Efficacy") {
      
      # Adjust half_life of long acting drug of ACT
      beta <- log(2) / this_scen$half_life_long
      a12 <- 8.46242774566474
      a21 <- 3.3058064516129035
      m <- 0.25
      W <- 60
      KK <- (beta ^ 2 - a12 * beta - a21 * beta) / (beta - a21)
      this_scen$half_life_long <- KK * (W ^ m)
      
      # Adjust coverage of each round of SMC
      # for setting with low seasonality pattern
      if (this_scen$seasonality == "sesonality2") {
      
        # If four rounds of SMC
        if (this_scen$Number_round == 4) {
              this_scen$coverage_0 <- 0
              this_scen$coverage_1 <- 1
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_4 <- 1 * (1 - this_scen$Coverage_reduction) ^ 3
              this_scen$coverage_5 <- 0
        }
        
        # If four rounds of SMC + 1 before
        if (this_scen$Number_round == 5.4) {
              this_scen$coverage_0 <- 1
              this_scen$coverage_1 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 3
              this_scen$coverage_4 <- 1 * (1 - this_scen$Coverage_reduction) ^ 4
              this_scen$coverage_5 <- 0
        }
        
        # If four rounds of SMC + 1 later
        if (this_scen$Number_round == 4.5) {
              this_scen$coverage_0 <- 0
              this_scen$coverage_1 <- 1
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_4 <- 1 * (1 - this_scen$Coverage_reduction) ^ 3
              this_scen$coverage_5 <- 1 * (1 - this_scen$Coverage_reduction) ^ 4
        }
      }
      
      # For setting with high seasonality pattern
      if (this_scen$seasonality == "sesonality3") {
        
        # If three rounds of SMC
        if (this_scen$Number_round == 4) {
              this_scen$coverage_0 <- 0
              this_scen$coverage_1 <- 1
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_4 <- 0
              this_scen$coverage_5 <- 0
        }
        
        # If three rounds of SMC + 1 before
        if (this_scen$Number_round == 5.4) {
              this_scen$coverage_0 <- 1
              this_scen$coverage_1 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 3
              this_scen$coverage_4 <- 0
              this_scen$coverage_5 <- 0
        }
        
        # If three rounds of SMC + 1 later
        if (this_scen$Number_round == 4.5) {
              this_scen$coverage_0 <- 0
              this_scen$coverage_1 <- 1
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_4 <- 1 * (1 - this_scen$Coverage_reduction) ^ 3
              this_scen$coverage_5 <- 0
        }
      }
      
      # For setting with high seasonality pattern
      if (this_scen$seasonality == "sesonality4") {
              this_scen$coverage_0 <- 0
              this_scen$coverage_1 <- 1
              this_scen$coverage_2 <- 1 * (1 - this_scen$Coverage_reduction)
              this_scen$coverage_3 <- 1 * (1 - this_scen$Coverage_reduction) ^ 2
              this_scen$coverage_4 <- 0
              this_scen$coverage_5 <- 0
      }
      
      # In case of random coverage modify the coverage at each round
      if (pm$opts$Type_coverage == "Random") {
              this_scen$coverage_0 <- this_scen$coverage_0 * this_scen$Coverage
              this_scen$coverage_1 <- this_scen$coverage_1 * this_scen$Coverage
              this_scen$coverage_2 <- this_scen$coverage_2 * this_scen$Coverage
              this_scen$coverage_3 <- this_scen$coverage_3 * this_scen$Coverage
              this_scen$coverage_4 <- this_scen$coverage_4 * this_scen$Coverage
              this_scen$coverage_5 <- this_scen$coverage_5 * this_scen$Coverage
      }
      
      # Adjust timing of SMC deployment deepening on seasonality profile and EIR
      
      # Adjust timing for setting with low seasonality profile
      if (this_scen$seasonality == "sesonality2") {
        Timming_deploy <- TIMING[TIMING$Sesonality == "sesonality2",][c(9:15)]
      }
      
      # Adjust timing for setting with high seasonality profile
      if (this_scen$seasonality == "sesonality3") {
        
        # Adjust timing for setting with low seasonality profile and low EIR
        if (this_scen$eir <= 50) {
          Timming_deploy <- TIMING[TIMING$Sesonality == "sesonality3" & TIMING$EIR == round(this_scen$eir, digits = 0),][c(9:15)]
        
          } else {
          
          # Adjust timing for setting with low seasonality profile and high EIR
          Timming_deploy <- TIMING[TIMING$Sesonality == "sesonality3" & TIMING$EIR == 50,][c(9:15)]
        }
      }
      
      # Adjust timing for setting with high seasonality profile
      if (this_scen$seasonality == "sesonality4") {
        
        # Adjust timing for setting with low seasonality profile and low EIR
        Timming_deploy <- TIMING[TIMING$Sesonality == "sesonality3" & TIMING$EIR == 50,][c(9:15)]
        Timming_deploy$time_0_1 <- 7
        Timming_deploy$time_1_1 <- 8
        Timming_deploy$time_2_1 <- 9
        Timming_deploy$time_3_1 <- 10
        Timming_deploy$time_4_1 <- 11
        Timming_deploy$time_5_1 <- 12
        Timming_deploy$DAY <- 15
      }
      
      # Adjust EIR
      
      # Adjust input EIR for setting with low level of seasonality
      if (this_scen$seasonality == "sesonality2") {
        this_scen$eir <- adjust_EIR(this_scen$eir, this_scen$Access)
      }
      
      # Adjust input EIR for setting with high level of seasonality
      if (this_scen$seasonality == "sesonality3") {
        this_scen$eir <- adjust_EIR_2(this_scen$eir, this_scen$Access)
      }
    }
    
    # Adjustment needed when we estimate the length of the prophylactic period
    if (pm$opts$Type_analysis == "Propylaxis") {
      
      # Adjust input EIR for setting with low level of seasonality
      if (this_scen$seasonality == "sesonality2") {
        this_scen$eir <- adjust_EIR(this_scen$eir, 0.04)
      }
      
      # Adjust input EIR for setting with high level of seasonality
      if (this_scen$seasonality == "sesonality3") {
        this_scen$eir <- adjust_EIR_2(this_scen$eir, 0.04)
      }
    }
    
    # Load the fourrier coefficient for the seasonality pattern
    seasonality_eir <- pm$settings$seasonality[[this_season]]
    
    # Select all the input parameter that will be used in Openmalaria if the analysis aimed to estimate the selection coefficient of SMC resistant parasite
    if (pm$opts$Type_analysis == "SMC") {
          param_values <- c(this_scen[, c(3:5, 8:19)])
          param_values <- c(param_values, seasonality_eir, Timming_deploy) ##
    }
    
    if (pm$opts$Type_analysis == "SMC" &
        pm$opts$Type_coverage == "Random") {
          param_values <- c(this_scen[, c(3:5, 9:19)])
          param_values <- c(param_values, seasonality_eir, Timming_deploy) ##
    }
    
    # Select all the input parameter that will be used in Openmalaria if the analysis aimed to estimate the prophylactic period
    if (pm$opts$Type_analysis == "Propylaxis") {
          param_values <- c(this_scen[, c(3:7)])
          param_values <- c(param_values, seasonality_eir)
    }
    
    # Select all the input parameter that will be used in Openmalaria if the analysis aimed to estimate the efficacy
    if (pm$opts$Type_analysis == "Efficacy") {
          param_values <- c(this_scen[, c(3:4, 7:18)])
          param_values <- c(param_values, seasonality_eir, Timming_deploy) ##
    }
    
    # Select all the input parameter that will be used in Openmalaria if the analysis aimed to estimate the efficacy if random coverage
    if (pm$opts$Type_analysis == "Efficacy"  & pm$opts$Type_coverage == "Random") {
          param_values <- c(this_scen[, c(3:4, 8:18)])
          param_values <- c(param_values, seasonality_eir, Timming_deploy) ##
    }
    
    # Select all the input parameter that will be used in Openmalaria if the analysis aimed to replicate Zongo trial
        if (pm$opts$Type_analysis == "Trial") {
          param_values_2 <- c(this_scen[2:length(name_para)])
          param_values <- c(this_scen[, c(2:6)])
        }
    
    # ---- Save the seed file and parameter tables ---- #
    
    # Format parameter values to not use scientific notation
    param_format <- format(param_values, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
    
    # Create sed command replacement pattern for each parameter
    sed_patterns <- paste0("s$@", param_names, "@$", param_format, "$g;")
    
    # Save file
    file_name <- paste0(this_scen$scenario_name, "_", this_scen$seed, ".txt")
    file_path <- paste0(pm$pth$sim_sed, file_name)
    write.table(sed_patterns,
      file = file_path,
      quote = FALSE,
      col.names = FALSE,
      row.names = FALSE)
    
    # Save all info on the scenario (save without and with adjusted parameters)
    if (pm$opts$Type_analysis == "SMC") {
      Liste_scenario[i,] <- c(param_values_2[1:(length(param_values_2) - 4)], this_season, file_name)
      Liste_scenario_adjusted[i,] <- c(param_values[1:(length(param_values))], file_name)
    }
    
    if (pm$opts$Type_analysis == "Propylaxis") {
      Liste_scenario[i,] <- c(param_values_2, this_season, file_name)
      Liste_scenario_adjusted[i,] <- c(param_values, file_name)
    }
    
    if (pm$opts$Type_analysis == "Efficacy") {
      Liste_scenario[i,] <- c(param_values_2[1:(length(param_values_2))], this_season, file_name)
      Liste_scenario_adjusted[i,] <- c(param_values[1:(length(param_values))], file_name)
    }
    
    if (pm$opts$Type_analysis == "Trial") {
      Liste_scenario[i,] <- c(param_values_2, file_name)
      Liste_scenario_adjusted[i,] <- c(param_values, file_name)
    }
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
  
  # Close seed loop    
  }
  
  # Close progress bar
  close(pb)
  
  # Save the list of scenario with adjustement and without adjustement
  file_name <- paste0(pm$pth$sim_sed, "scenarios_", ".txt")
  write.table(Liste_scenario, file = file_name, row.names = FALSE)
  file_name <- paste0(pm$pth$sim_sed, "scenarios_adjusted_", ".txt")
  write.table(Liste_scenario_adjusted,
              file = file_name,
              row.names = FALSE)
  
  # Return the number of scenario
  return(nrow(param_table))
  
}

# -------------------------------------------------------------
# Generate scenario xml files from base file using sed patterns.
# -------------------------------------------------------------
generate_xml_files <- function(pm,
                               n_jobs,
                               file_path = NA,
                               sed_path = NA,
                               xml_path = NA) {
  # Message to the console
  message("  - Generating xml files")
  
  # Create a new log file for the cluster jobs
  log_file <- file.path(pm$pth$log_files, "scicore_log.txt")
  create_bash_log(log_file)
  
  # Add country sub directory to file paths
  sed_country <- paste0(sed_path)
  xml_country <- paste0(xml_path)
  
  # Create sequence of xml to be run at the same time
  n_seq <- 100
  n_submit <- ceil(n_jobs / n_seq)
  
  # Message to the console
  message("   ~ ", thou_sep(n_submit), " sets of ", thou_sep(n_seq))
  
  # Construct slurm array command for running in parallel
  slurm_array <- paste0("--array=1-", n_submit, "%", 300)
  
  # Concatenate system command
  sys_command <- paste("sbatch", slurm_array, "bash_make_xmls.sh", n_seq, file_path, sed_country,  xml_country, log_file)
  
  # Invoke this command
  system(sys_command)
  
  # Wait for all cluster jobs to complete
  wait_for_jobs(pm, log_file, n_jobs)
  
}
