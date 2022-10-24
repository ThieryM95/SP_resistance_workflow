##################################################################
# Estimate the measrument of interest in each simulation         #
#                                                                #
# Input: Post-processed data  of OpenMalaria                     #
# Output: Table of parameter input of each simulation and        #
#         the estimated indicator of interest                    #
#                                                                #
# authors: thiery.masserey@swisstph.ch                           #
##################################################################

# Load function
library(grDevices)
library(pkr)

# ------------------------------------------------------------------------------
# Estimate the selection coefficient or prophylactic period for each simulation.
# ------------------------------------------------------------------------------
SummaryResults <- function(pm) {
  
  # Message in the console
  message("  - Calulate the spread or esthablihsment ")
  
  # Load the list of scenario
  Scenario_liste <- read.table(file.path(pm$pth$sim_sed, "scenarios_.txt"), header = T)
  
  # Load the list of timing of deployment to make sure it occure as expected
  TIMING <- read.table("/scicore/home/penny/masthi00/smc_resistance/Data_time_2.txt", sep = "\t", header = T)
  
  # Load the list of Output files
  Output <- list.files(pm$pth$processed)
  
  # Define the time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step
  Number_survey_years <- 365 / time_step
  
  # List of new variable that are estimated during the analysis that aimed to estimate the spread of drug resistance
  if (pm$opts$Type_analysis == "SMC") {
    Scenario_liste$Indicator <- NA
    Scenario_liste$Indicator_1 <- NA
    Scenario_liste$Input_EIR <- NA
    Scenario_liste$Simulated_EIR <- NA
    Scenario_liste$time_eir <- NA
    Scenario_liste$time_eir_2 <- NA
  }
  
  # List of new variable that are estimated during the analysis that aimed to estimate the prophylactic periode
  if (pm$opts$Type_analysis == "Propylaxis") {
    Scenario_liste$Indicator_40_b <- NA
    Scenario_liste$Indicator_45_b <- NA
    Scenario_liste$Indicator_120 <- NA
    Scenario_liste$Indicator_120_infect <- NA
    Scenario_liste$Simulated_EIR <- NA
  }
  
  # List of new variable that are estimated during the analysis that aimed to estimate the efficacy
  if (pm$opts$Type_analysis == "Efficacy") {
    Scenario_liste$Uncomplicated_0  <- NA
    Scenario_liste$Uncomplicated_1  <- NA
    Scenario_liste$Uncomplicated_0_year <- NA
    Scenario_liste$Uncomplicated_1_year <- NA
    Scenario_liste$Uncomplicated_0_round  <- NA
    Scenario_liste$Uncomplicated_1_round  <- NA
    Scenario_liste$Person_0 <- NA
    Scenario_liste$Person_1 <- NA
    Scenario_liste$Person_1_round <- NA
    Scenario_liste$Input_EIR <- NA
    Scenario_liste$Simulated_EIR <- NA
    Scenario_liste$time_eir <- NA
    Scenario_liste$time_eir_2 <- NA
  }
  
  # List of new variable that are estimated during the analysis that aimed to replicate Zongo trial
  if (pm$opts$Type_analysis == "Trial") {
    Scenario_liste$prevalance_before <- NA
    Scenario_liste$prevalance_after <- NA
    Scenario_liste$Input_EIR <- NA
    Scenario_liste$Simulated_EIR <- NA
    Scenario_liste$Simulated_EIR_2 <- NA
    UNCOMPLICATE <- NULL
  }
  
  # Initiate progress bar
  pb <- txtProgressBar(min = 0, max = length(Output), initial = 0, width = 100, style = 3)
  
  # Loop across each simulation
  for (i in 1:length(Output)) {
    
    # Download the data
    Output_data_name <- Output[i]
    Output_data_file_path <- file.path(pm$pth$processed, Output_data_name)
    Output_data <- read.table(Output_data_file_path, sep = ";")
    name_scenario <- gsub("PostProcess_", "", Output_data_name)
    name_scenario <- gsub("_out", "", name_scenario)
    
    # Variable estimate when we estimate the selection coefficient of SMC resistant parasite
    if (pm$opts$Type_analysis == "SMC") {
      
      # Estimate the selection coefficient
      Scenario_liste$Indicator[Scenario_liste$scenario_name == name_scenario] <- spread_R2(Output_data) # see function bellow
      
      # Estimate the selection coefficient via the moving average
      Scenario_liste$Indicator_1[Scenario_liste$scenario_name == name_scenario] <- spread_R5(Output_data) # see function bellow
      
      # Estimate the input EIR
      Scenario_liste$Input_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$inputEIR[Output_data$Survey >= (20 * Number_survey_years) & Output_data$Survey <= (30 * Number_survey_years)] * 365 / 5)
      
      # Estimate the simulated EIR
      Scenario_liste$Simulated_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$simulatedEIR[Output_data$Survey >= (20 * Number_survey_years) & Output_data$Survey <= (30 * Number_survey_years)] * 365 / 5)
      
      # Estimate when the EIR reach it maximum
      Scenario_liste$time_eir[Scenario_liste$scenario_name == name_scenario] <- Time_eir_max(Output_data) # see function bellow
      
      # Estimate the time at which SMC was deployed for the high seasonality setting
      if (Scenario_liste$eir[Scenario_liste$scenario_name == name_scenario] <= 50 & Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality3") {
        Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario] <- TIMING$time_2[TIMING$EIR == round(Scenario_liste$eir[Scenario_liste$scenario_name == name_scenario])] - 0.7
      }
      
      # Estimate the time at which SMC was deployed for the high seasonality setting
      if (Scenario_liste$eir[Scenario_liste$scenario_name == name_scenario] >= 50 & Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality3") {
        Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario] <- TIMING$time_2[TIMING$EIR == 50] - 0.7
      }
      
      # Estimate the time at which SMC was deployed for the moderate seasonality setting
      if (Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality2") {
        Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario] <- TIMING$time_2[TIMING$EIR == 1] - 0.7
      }
      
      # If the peak of EIR is more than one month away from the time during which the second round of SMC is deploy delete the simulation
      if (abs(Scenario_liste$time_eir[Scenario_liste$scenario_name == name_scenario] - Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario]) >= 1) {
        Scenario_liste$Indicator[Scenario_liste$scenario_name == name_scenario] <- NA
      }
    }
    
    # Variable estimated when we estimate the prophylactic period
    if (pm$opts$Type_analysis == "Propylaxis") { 
      
      # Estimate the prophylactic period over 40 day based on the number of patent individual
      Scenario_liste$Indicator_40_b[Scenario_liste$scenario_name == name_scenario] <- time_prophylaxies_40_b(Output_data)
      Scenario_liste$Indicator_45_b[Scenario_liste$scenario_name == name_scenario] <- time_prophylaxies_45_b(Output_data)
      
      # Estimate the prophylactic period over 120 days based on the number of patent individual
      Scenario_liste$Indicator_120[Scenario_liste$scenario_name == name_scenario] <- time_prophylaxies_120(Output_data)
      
      # Estimate the prophylactic period over 120 days based on the number of infected individuals
      Scenario_liste$Indicator_120_infect[Scenario_liste$scenario_name == name_scenario] <- time_prophylaxies_120_infect(Output_data)
      
      # Estimate the simulated EIR
      Scenario_liste$Simulated_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$simulatedEIR[Output_data$Survey >= (30 * Number_survey_years) & Output_data$Survey <= (40 * Number_survey_years)] * 365 / 5)
      
    }
    
    # Variable estimate when we estimate efficacy of SMC
    if (pm$opts$Type_analysis == "Efficacy") {
      
      # Estimate the incidence of uncomplicated malaria during the five months of high transmission when SMC is not deployed
      Scenario_liste$Uncomplicated_0[Scenario_liste$scenario_name == name_scenario] <- Uncomplicated_0(Output_data) # see function bellow
      
      # Estimate  the incidence of of uncomplicated malaria during the five month of high transmission when SMC is deployed
      Scenario_liste$Uncomplicated_1[Scenario_liste$scenario_name == name_scenario] <- Uncomplicated_1(Output_data) # see function bellow
      
      # Estimate the incidence of uncomplicated malaria during one year when SMC is not deployed
      Scenario_liste$Uncomplicated_0_year[Scenario_liste$scenario_name == name_scenario] <- Uncomplicated_0_year(Output_data) # see function bellow
      
      # Estimate the incidence of uncomplicated malaria during one year when SMC is deployed
      Scenario_liste$Uncomplicated_1_year[Scenario_liste$scenario_name == name_scenario] <- Uncomplicated_1_year(Output_data) # see function bellow
      
      # Estimate the incidence of uncomplicated malaria during the five month of high transmission when SMC is not deployed
      Scenario_liste$Uncomplicated_0_round[Scenario_liste$scenario_name == name_scenario] <- Uncomplicated_0_round(Output_data) # see function bellow
      
      # Estimate  the incidence of of uncomplicated malaria during the five month of high transmission when SMC is deployed
      Scenario_liste$Uncomplicated_1_round[Scenario_liste$scenario_name == name_scenario] <- Uncomplicated_1_round(Output_data) # see function bellow
      
      # Estimate the number of person at risk when no SMC is deployed.
      Scenario_liste$Person_0[Scenario_liste$scenario_name == name_scenario] <- Person_0(Output_data) # see function bellow
      
      # Estimate the number of person at risk when SMC is deployed.
      Scenario_liste$Person_1[Scenario_liste$scenario_name == name_scenario] <- Person_1(Output_data) # see function bellow
      
      # Estimate the number of person at risk when SMC is deployed over one years.
      Scenario_liste$Person_1_round[Scenario_liste$scenario_name == name_scenario] <- Person_1_round(Output_data) # see function bellow
      
      # Estimate the input EIR
      Scenario_liste$Input_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$inputEIR[Output_data$Survey >= (20 * Number_survey_years) & Output_data$Survey <= (30 * Number_survey_years)] * 365 / 5)
      
      # Estimate the simulated EIR
      Scenario_liste$Simulated_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$simulatedEIR[Output_data$Survey >= (20 * Number_survey_years) & Output_data$Survey <= (30 * Number_survey_years)] * 365 / 5)
      
      # Estimate when the EIR reach it maximum
      Scenario_liste$time_eir[Scenario_liste$scenario_name == name_scenario] <- Time_eir_max(Output_data) # see function bellow
      
      # Estimate the time at which SMC was deployed in the high seasonality setting
      if (Scenario_liste$eir[Scenario_liste$scenario_name == name_scenario] <= 50 & Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality3") {
        Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario] <- TIMING$time_2[TIMING$EIR == round(Scenario_liste$eir[Scenario_liste$scenario_name == name_scenario])] - 0.7
      }
      
      # Estimate the time at which SMC was deployed in the high seasonality setting
      if (Scenario_liste$eir[Scenario_liste$scenario_name == name_scenario] >= 50 & Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality3") {
        Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario] <- TIMING$time_2[TIMING$EIR == 50] - 0.7
      }
      
      # Estimate the time at which SMC was deployed in the moderate seasonality setting
      if (Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality2") {
        Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario] <- TIMING$time_2[TIMING$EIR == 1] - 0.7
      }
      
      # If the peak of EIR is more than one month away from the time during which the second round of SMC is deploy delete the simulation
      if (abs(Scenario_liste$time_eir[Scenario_liste$scenario_name == name_scenario] - Scenario_liste$time_eir_2[Scenario_liste$scenario_name == name_scenario]) >= 1) {
        Scenario_liste$Uncomplicated_1[Scenario_liste$scenario_name == name_scenario] <- NA
        Scenario_liste$Uncomplicated_1_year[Scenario_liste$scenario_name == name_scenario] <- NA
      }
    }
    
    # Variable estimate when we aim to replicate the Zongo trial
    if (pm$opts$Type_analysis == "Trial") {
      
      # Assess the number of uncomplicated malaria
      UNCOMPLICATE[[name_scenario]] <- Output_data$nUncomp_2[Output_data$Survey >= 2906 & Output_data$Survey <= 2936] + Output_data$nUncomp_3[Output_data$Survey >= 2906 & Output_data$Survey <= 2936]
      
      # Estimate the input EIR
      Scenario_liste$Input_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$inputEIR[Output_data$Survey >= (20 * Number_survey_years) & Output_data$Survey <= (30 * Number_survey_years)] * 365 / 5)
      
      # Estimate the simulated EIR
      Scenario_liste$Simulated_EIR[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$simulatedEIR[Output_data$Survey >= (20 * Number_survey_years) & Output_data$Survey <= (30 * Number_survey_years)] * 365 / 5)
      Scenario_liste$Simulated_EIR_2[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$simulatedEIR[Output_data$Survey >= (39 * Number_survey_years) & Output_data$Survey <= (40 * Number_survey_years)] * 365 / 5)
      
      # Define time of SMC
      Time_MDA <- c(2894, 2900, 2906)
      
      # For the test group assess the prevalence in August, and November
      if (Scenario_liste$Coverage[Scenario_liste$scenario_name == name_scenario] == 1) {
        Scenario_liste$prevalance_before[Scenario_liste$scenario_name == name_scenario] <- (Output_data$nPatent_2[Time_MDA[1] - 1] + Output_data$nPatent_3[Time_MDA[1] - 1]) / (Output_data$nHost_2[Time_MDA[1] - 1] + Output_data$nHost_3[Time_MDA[1] - 1])
        Scenario_liste$prevalance_after[Scenario_liste$scenario_name == name_scenario] <- (Output_data$nPatent_2[Time_MDA[3] + 6] + Output_data$nPatent_3[Time_MDA[3] + 6]) / (Output_data$nHost_2[Time_MDA[3] + 6] + Output_data$nHost_3[Time_MDA[3] + 6])
      
      # For the control group assess the prevalence in September, and December  
      } else{
        Scenario_liste$prevalance_before[Scenario_liste$scenario_name == name_scenario] <- (Output_data$nPatent_2[Time_MDA[2] - 1] + Output_data$nPatent_3[Time_MDA[2] - 1]) / (Output_data$nHost_2[Time_MDA[2] - 1] + Output_data$nHost_3[Time_MDA[2] - 1])
        Scenario_liste$prevalance_after[Scenario_liste$scenario_name == name_scenario] <- (Output_data$nPatent_2[Time_MDA[3] + 12] + Output_data$nPatent_3[Time_MDA[3] + 12]) / (Output_data$nHost_2[Time_MDA[3] + 12] + Output_data$nHost_3[Time_MDA[3] + 12])
      }
      
    }
    
    Scenario_liste$Liste[Scenario_liste$scenario_name == name_scenario] <- i
    
    # Update the progress bar
    setTxtProgressBar(pb, i)
    
    # end the loop
  }
  
  # Close the progress bar
  close(pb)
  
  # Save the dataset in the summarized folder
  name_sumarised <- param_set_name(sample_num = pm$sample_num)
  summary_file <- paste0(name_sumarised, ".txt")
  summary_path <- paste0(pm$pth$summarised, summary_file)
  write.table( Scenario_liste, file = summary_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  # If zongo trial, save the UNCOMPLCATED dataset too
  if (pm$opts$Type_analysis == "Trial") {
    saveRDS(UNCOMPLICATE, file = paste0(pm$pth$summarised, "Data_delay.RData"))
  }
  
  # Merge the summarized results of each iterations
  full_table_name <- param_set_name(all_samples = TRUE)
  
  # Construct file name and path
  full_summary_file <- paste0(full_table_name, "_", ".txt")
  full_summary_path <- paste0(pm$pth$summarised, full_summary_file)
  
  # Initiate 'all sample' table
  full_param_table <- NULL
  
  # Process only necessary for adaptive sampling summary
  if (pm$sample_num > 0) {
    
    # Loop through previously generated samples
    for (i in (1:pm$sample_num - 1)) {
      
      # Construct sample results summary file name and path
      sample_file <- paste0(param_set_name(sample_num = i), "_", ".txt")
      
      # Load up the param table associated with this sample number
      sample_param_table <- read.table(paste0(pm$pth$summarised, sample_file), sep = "\t", header = TRUE, as.is = TRUE)
      full_param_table <- rbind(full_param_table, sample_param_table)
      
    }
  }
  
  # Concatenate the param table we have just summarized
  full_param_table <- rbind(full_param_table, Scenario_liste)
  
  # Overwrite scenario names in this 'all samples' file
  full_param_table$scenario_name <- paste0("scenario_", 1:nrow(full_param_table))
  
  # (Over)-write full result tables across all samples
  write.table(full_param_table, full_summary_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
} # end the function


#--------------------------------------------------------------
# function to estimate the time at which the EIR reach it peak.
#--------------------------------------------------------------

Time_eir_max <- function(Output_data_2) {
  
  # Define the timestep
  time_step <- 5
  Number_survey_generation <- 60 / time_step
  Number_survey_years <- 365 / time_step
  
  # Select the data over a year
  TIME <- Output_data_2[Output_data_2$Survey >= (29 * Number_survey_years) & Output_data_2$Survey < (30 * Number_survey_years),]
  
  # Select the time at which the EIR is at it maximum
  TIME_MAX <- TIME$Survey[TIME$simulatedEIR == max(TIME$simulatedEIR)]
  
  # Convert the data into month
  TIME_MAX <- (TIME_MAX / Number_survey_years - 29)
  TIME_MAX_2 <- TIME_MAX * 365 / 30.4
  
  # Return the time
  return(TIME_MAX_2)
  
}

#--------------------------------------------------------------------------------------------------------------
# Estimate the mean time at which individual become patent after the deployment of SMC over a period of 40 day.
#--------------------------------------------------------------------------------------------------------------

time_prophylaxies_40_b <- function(Output_data_2) {
  
  # Define the time step
  time_step <- 5
  Number_survey_years <- 365 / time_step
  
  # Select the time at which we deploy SMC
  Time_MDA <- Output_data_2$Survey[Output_data_2$nMDAs_2 >= 1]
  time_mean <- 0
  mean_rinfection <- 0
  
  # Loop across each SMC deployment and estimate the median time at which  patient become patent
  for (i in 1:length(Time_MDA)) {
    
    # Estimate how many individual are patent at each time step
    time_before <- Output_data_2$nPatent_2[Time_MDA[i] - 5 / 5]
    time_0 <- Output_data_2$nPatent_2[Time_MDA[i]]
    time_5 <- Output_data_2$nPatent_2[Time_MDA[i] + 5 / 5]
    time_10 <- Output_data_2$nPatent_2[Time_MDA[i] + 10 / 5]
    time_15 <- Output_data_2$nPatent_2[Time_MDA[i] + 15 / 5]
    time_20 <- Output_data_2$nPatent_2[Time_MDA[i] + 20 / 5]
    time_25 <- Output_data_2$nPatent_2[Time_MDA[i] + 25 / 5]
    time_30 <- Output_data_2$nPatent_2[Time_MDA[i] + 30 / 5]
    time_35 <- Output_data_2$nPatent_2[Time_MDA[i] + 35 / 5]
    time_40 <- Output_data_2$nPatent_2[Time_MDA[i] + 40 / 5]
    
    time_5 <- ifelse(time_5 < time_before, time_5, time_before)
    time_10 <- ifelse(time_10 < time_before, time_10, time_before)
    time_15 <- ifelse(time_15 < time_before, time_15, time_before)
    time_20 <- ifelse(time_20 < time_before, time_20, time_before)
    time_25 <- ifelse(time_25 < time_before, time_25, time_before)
    time_30 <- ifelse(time_30 < time_before, time_30, time_before)
    time_35 <- ifelse(time_35 < time_before, time_35, time_before)
    time_40 <- ifelse(time_40 < time_before, time_40, time_before)
    
    
    # Estimate the number of newly patent individual at each time step (can be only positive)
    Number_5 <- time_5 #ifelse(time_5 < 0, 0, time_5)
    Number_10 <- ifelse(time_10 - time_5 < 0, 0, time_10 - time_5)
    Number_15 <- ifelse(time_15 - time_10 < 0, 0, time_15 - time_10)
    Number_20 <- ifelse(time_20 - time_15 < 0, 0, time_20 - time_15)
    Number_25 <- ifelse(time_25 - time_20 < 0, 0, time_25 - time_20)
    Number_30 <- ifelse(time_30 - time_25 < 0, 0, time_30 - time_25)
    Number_35 <- ifelse(time_35 - time_30 < 0, 0, time_35 - time_30)
    Number_40 <- ifelse(time_40 - time_35 < 0, 0, time_40 - time_35)
    
    # Estimate the median time at which patient become patent
    time_mean[i] <-
      median(c(rep(0, Number_5),
        rep(5, Number_10),
        rep(10, Number_15),
        rep(15, Number_20),
        rep(20, Number_25),
        rep(25, Number_30),
        rep(30, Number_35),
        rep(35, Number_40)))
    
    # If the number of patent patient at time 40 is still low, then the prophylactic periode is defined to be 40
    if ((time_40 - 10) <= mean(c(time_5, time_35)) &
        time_40 <= (time_before / 2)) {
      time_mean[i] <- 40
    }
  }
  
  # Estimate the mean of the median time at which patient become patent again
  mean_rinfection <- mean(time_mean)
  
  # Return the mean
  return(mean_rinfection)
}

#--------------------------------------------------------------------------------------------------------------
# Estimate the mean time at which individual become patent after the deployment of SMC over a period of 45 day.
#--------------------------------------------------------------------------------------------------------------

time_prophylaxies_45_b <- function(Output_data_2) {
  
  # Define the time step
  time_step <- 5
  Number_survey_years <- 365 / time_step
  
  # Select the time at which we deploy SMC
  Time_MDA <- Output_data_2$Survey[Output_data_2$nMDAs_2 >= 1]

  time_mean <- 0
  mean_rinfection <- 0
  
  # Loop across each SMC deployment and estimate the median time at which  patient become patent
  for (i in 1:length(Time_MDA)) {
    
    # Estimate how many individual are patent at each time step
    time_before <- Output_data_2$nPatent_2[Time_MDA[i] - 5 / 5]
    time_0 <- Output_data_2$nPatent_2[Time_MDA[i]]
    time_5 <- Output_data_2$nPatent_2[Time_MDA[i] + 5 / 5]
    time_10 <- Output_data_2$nPatent_2[Time_MDA[i] + 10 / 5]
    time_15 <- Output_data_2$nPatent_2[Time_MDA[i] + 15 / 5]
    time_20 <- Output_data_2$nPatent_2[Time_MDA[i] + 20 / 5]
    time_25 <- Output_data_2$nPatent_2[Time_MDA[i] + 25 / 5]
    time_30 <- Output_data_2$nPatent_2[Time_MDA[i] + 30 / 5]
    time_35 <- Output_data_2$nPatent_2[Time_MDA[i] + 35 / 5]
    time_40 <- Output_data_2$nPatent_2[Time_MDA[i] + 40 / 5]
    time_45 <- Output_data_2$nPatent_2[Time_MDA[i] + 45 / 5]
    
    # Estimate the number of newly patent individual at each time step (can be only positive)
    Number_5 <- ifelse(time_5 < 0, 0, time_5)
    Number_10 <- ifelse(time_10 - time_5 < 0, 0, time_10 - time_5)
    Number_15 <- ifelse(time_15 - time_10 < 0, 0, time_15 - time_10)
    Number_20 <- ifelse(time_20 - time_15 < 0, 0, time_20 - time_15)
    Number_25 <- ifelse(time_25 - time_20 < 0, 0, time_25 - time_20)
    Number_30 <- ifelse(time_30 - time_25 < 0, 0, time_30 - time_25)
    Number_35 <- ifelse(time_35 - time_30 < 0, 0, time_35 - time_30)
    Number_40 <- ifelse(time_40 - time_35 < 0, 0, time_40 - time_35)
    Number_45 <- ifelse(time_45 - time_40 < 0, 0, time_45 - time_40)
    
    # Estimate the median time at which patient become patent
    time_mean[i] <- median(c(
        rep(0, Number_5),
        rep(5, Number_10),
        rep(10, Number_15),
        rep(15, Number_20),
        rep(20, Number_25),
        rep(25, Number_30),
        rep(30, Number_35),
        rep(35, Number_40),
        rep(40, Number_45)))
    
    # If the number of patent patient at time 40 is still low, then the prophylactic period is defined to be 45
    if ((time_45 - 10) <= mean(c(time_5, time_40)) &
        time_45 <= (time_before / 2)) {
      time_mean[i] <- 45
    }
  }
  
  # Estimate the mean of the median time at which patient become patent again
  mean_rinfection <- mean(time_mean)
  
  # Return the mean
  return(mean_rinfection)
}

#---------------------------------------------------------------------------------------------------------------
# Estimate the mean time at which individual become patent after the deployment of SMC over a period of 120 day.
#---------------------------------------------------------------------------------------------------------------

time_prophylaxies_120 <- function(Output_data_2) {
  
  # Define the time step
  time_step <- 5
  Number_survey_years <- 365 / time_step
  
  # Select the time at which we deploy SMC
  Time_MDA <- Output_data_2$Survey[Output_data_2$nMDAs_2 >= 1]
  
  time_mean <- 0
  mean_rinfection <- 0
  
  # Loop across each SMC deployment and estimate the median time at which  patient become patent
  for (i in 1:length(Time_MDA)) {
    
    # Estime how many individual are patent at each time step
    time_0 <- Output_data_2$nPatent_2[Time_MDA[i] - 5 / 5]
    time_5 <- Output_data_2$nPatent_2[Time_MDA[i] + 5 / 5]
    time_10 <- Output_data_2$nPatent_2[Time_MDA[i] + 10 / 5]
    time_15 <- Output_data_2$nPatent_2[Time_MDA[i] + 15 / 5]
    time_20 <- Output_data_2$nPatent_2[Time_MDA[i] + 20 / 5]
    time_25 <- Output_data_2$nPatent_2[Time_MDA[i] + 25 / 5]
    time_30 <- Output_data_2$nPatent_2[Time_MDA[i] + 30 / 5]
    time_35 <- Output_data_2$nPatent_2[Time_MDA[i] + 35 / 5]
    time_40 <- Output_data_2$nPatent_2[Time_MDA[i] + 40 / 5]
    time_45 <- Output_data_2$nPatent_2[Time_MDA[i] + 45 / 5]
    time_50 <- Output_data_2$nPatent_2[Time_MDA[i] + 50 / 5]
    time_55 <- Output_data_2$nPatent_2[Time_MDA[i] + 55 / 5]
    time_60 <- Output_data_2$nPatent_2[Time_MDA[i] + 60 / 5]
    time_65 <- Output_data_2$nPatent_2[Time_MDA[i] + 65 / 5]
    time_70 <- Output_data_2$nPatent_2[Time_MDA[i] + 70 / 5]
    time_75 <- Output_data_2$nPatent_2[Time_MDA[i] + 75 / 5]
    time_80 <- Output_data_2$nPatent_2[Time_MDA[i] + 80 / 5]
    time_85 <- Output_data_2$nPatent_2[Time_MDA[i] + 85 / 5]
    time_90 <- Output_data_2$nPatent_2[Time_MDA[i] + 90 / 5]
    time_95 <- Output_data_2$nPatent_2[Time_MDA[i] + 95 / 5]
    time_100 <- Output_data_2$nPatent_2[Time_MDA[i] + 100 / 5]
    time_105 <- Output_data_2$nPatent_2[Time_MDA[i] + 105 / 5]
    time_110 <- Output_data_2$nPatent_2[Time_MDA[i] + 110 / 5]
    time_115 <- Output_data_2$nPatent_2[Time_MDA[i] + 115 / 5]
    time_120 <- Output_data_2$nPatent_2[Time_MDA[i] + 120 / 5]
    
    # Adjust that this value cannot be higher than the number of people patent before SMC deployment
    time_5 <- ifelse(time_5 < time_0, time_5, time_0)
    time_10 <- ifelse(time_10 < time_0, time_10, time_0)
    time_15 <- ifelse(time_15 < time_0, time_15, time_0)
    time_20 <- ifelse(time_20 < time_0, time_20, time_0)
    time_25 <- ifelse(time_25 < time_0, time_25, time_0)
    time_30 <- ifelse(time_30 < time_0, time_30, time_0)
    time_35 <- ifelse(time_35 < time_0, time_35, time_0)
    time_40 <- ifelse(time_40 < time_0, time_40, time_0)
    time_45 <- ifelse(time_45 < time_0, time_45, time_0)
    time_50 <- ifelse(time_50 < time_0, time_50, time_0)
    time_55 <- ifelse(time_55 < time_0, time_55, time_0)
    time_60 <- ifelse(time_60 < time_0, time_60, time_0)
    time_65 <- ifelse(time_65 < time_0, time_65, time_0)
    time_70 <- ifelse(time_70 < time_0, time_70, time_0)
    time_75 <- ifelse(time_75 < time_0, time_75, time_0)
    time_80 <- ifelse(time_80 < time_0, time_80, time_0)
    time_85 <- ifelse(time_85 < time_0, time_85, time_0)
    time_90 <- ifelse(time_90 < time_0, time_90, time_0)
    time_95 <- ifelse(time_95 < time_0, time_95, time_0)
    time_100 <- ifelse(time_100 < time_0, time_100, time_0)
    time_105 <- ifelse(time_105 < time_0, time_105, time_0)
    time_110 <- ifelse(time_110 < time_0, time_110, time_0)
    time_115 <- ifelse(time_115 < time_0, time_115, time_0)
    time_120 <- ifelse(time_120 < time_0, time_120, time_0)
    
    # Estime the number of nelwy patent individual at each timestep (can be only positive)
    Number_5 <- ifelse(time_5 < 0, 0, time_5)
    Number_10 <- ifelse(time_10 - time_5 < 0, 0, time_10 - time_5)
    Number_15 <- ifelse(time_15 - time_10 < 0, 0, time_15 - time_10)
    Number_20 <- ifelse(time_20 - time_15 < 0, 0, time_20 - time_15)
    Number_25 <- ifelse(time_25 - time_20 < 0, 0, time_25 - time_20)
    Number_30 <- ifelse(time_30 - time_25 < 0, 0, time_30 - time_25)
    Number_35 <- ifelse(time_35 - time_30 < 0, 0, time_35 - time_30)
    Number_40 <- ifelse(time_40 - time_35 < 0, 0, time_40 - time_35)
    Number_45 <- ifelse(time_45 - time_40 < 0, 0, time_45 - time_40)
    Number_50 <- ifelse(time_50 - time_45 < 0, 0, time_50 - time_45)
    Number_55 <- ifelse(time_55 - time_50 < 0, 0, time_55 - time_50)
    Number_60 <- ifelse(time_60 - time_55 < 0, 0, time_60 - time_55)
    Number_65 <- ifelse(time_65 - time_60 < 0, 0, time_65 - time_60)
    Number_70 <- ifelse(time_70 - time_65 < 0, 0, time_70 - time_65)
    Number_75 <- ifelse(time_75 - time_70 < 0, 0, time_75 - time_70)
    Number_80 <- ifelse(time_80 - time_75 < 0, 0, time_80 - time_75)
    Number_85 <- ifelse(time_85 - time_80 < 0, 0, time_85 - time_80)
    Number_90 <- ifelse(time_90 - time_85 < 0, 0, time_90 - time_85)
    Number_95 <- ifelse(time_95 - time_90 < 0, 0, time_95 - time_90)
    Number_100 <- ifelse(time_100 - time_95 < 0, 0, time_100 - time_95)
    Number_105 <- ifelse(time_105 - time_100 < 0, 0, time_105 - time_100)
    Number_110 <- ifelse(time_110 - time_105 < 0, 0, time_110 - time_105)
    Number_115 <- ifelse(time_115 - time_110 < 0, 0, time_115 - time_110)
    Number_120 <- ifelse(time_120 - time_115 < 0, 0, time_120 - time_115)
    
    # estimate the median time at which people become patent again
    time_mean[i] <- median(
        c(rep(0, Number_5),
          rep(5, Number_10),
          rep(10, Number_15),
          rep(15, Number_20),
          rep(20, Number_25),
          rep(25, Number_30),
          rep(30, Number_35),
          rep(35, Number_40),
          rep(40, Number_45),
          rep(45, Number_50),
          rep(50, Number_55),
          rep(55, Number_60),
          rep(60, Number_65),
          rep(65, Number_70),
          rep(70, Number_75),
          rep(75, Number_80),
          rep(80, Number_85),
          rep(85, Number_90),
          rep(90, Number_95),
          rep(95, Number_100),
          rep(100, Number_105),
          rep(105, Number_110),
          rep(110, Number_115),
          rep(115, Number_120)
        )
      )
    
  }
  
  # Estimate the mean across each SMC deployment
  mean_rinfection <- mean(time_mean)
  
  # Return the mean
  return(mean_rinfection)
}

#-----------------------------------------------------------------------------------------------------------------
# Estimate the mean time at which individual become infected after the deployment of SMC over a period of 120 day.
#-----------------------------------------------------------------------------------------------------------------

time_prophylaxies_120_infect <- function(Output_data_2) {
  
  # Define the time step
  time_step <- 5
  Number_survey_years <- 365 / time_step
  
  # Select the time at which we deploy SMC
  Time_MDA <- Output_data_2$Survey[Output_data_2$nMDAs_2 >= 1]
  
  time_mean <- 0
  mean_rinfection <- 0
  
  # Loop across each SMC deployment and estimate the median time at which  patient become infected
  for (i in 1:length(Time_MDA)) {
    
    # Estimate how many individual are infected at each time step
    time_0 <- Output_data_2$nInfect_2[Time_MDA[i] - 5 / 5]
    time_5 <- Output_data_2$nInfect_2[Time_MDA[i] + 5 / 5]
    time_10 <- Output_data_2$nInfect_2[Time_MDA[i] + 10 / 5]
    time_15 <- Output_data_2$nInfect_2[Time_MDA[i] + 15 / 5]
    time_20 <- Output_data_2$nInfect_2[Time_MDA[i] + 20 / 5]
    time_25 <- Output_data_2$nInfect_2[Time_MDA[i] + 25 / 5]
    time_30 <- Output_data_2$nInfect_2[Time_MDA[i] + 30 / 5]
    time_35 <- Output_data_2$nInfect_2[Time_MDA[i] + 35 / 5]
    time_40 <- Output_data_2$nInfect_2[Time_MDA[i] + 40 / 5]
    time_45 <- Output_data_2$nInfect_2[Time_MDA[i] + 45 / 5]
    time_50 <- Output_data_2$nInfect_2[Time_MDA[i] + 50 / 5]
    time_55 <- Output_data_2$nInfect_2[Time_MDA[i] + 55 / 5]
    time_60 <- Output_data_2$nInfect_2[Time_MDA[i] + 60 / 5]
    time_65 <- Output_data_2$nInfect_2[Time_MDA[i] + 65 / 5]
    time_70 <- Output_data_2$nInfect_2[Time_MDA[i] + 70 / 5]
    time_75 <- Output_data_2$nInfect_2[Time_MDA[i] + 75 / 5]
    time_80 <- Output_data_2$nInfect_2[Time_MDA[i] + 80 / 5]
    time_85 <- Output_data_2$nInfect_2[Time_MDA[i] + 85 / 5]
    time_90 <- Output_data_2$nInfect_2[Time_MDA[i] + 90 / 5]
    time_95 <- Output_data_2$nInfect_2[Time_MDA[i] + 95 / 5]
    time_100 <- Output_data_2$nInfect_2[Time_MDA[i] + 100 / 5]
    time_105 <- Output_data_2$nInfect_2[Time_MDA[i] + 105 / 5]
    time_110 <- Output_data_2$nInfect_2[Time_MDA[i] + 110 / 5]
    time_115 <- Output_data_2$nInfect_2[Time_MDA[i] + 115 / 5]
    time_120 <- Output_data_2$nInfect_2[Time_MDA[i] + 120 / 5]
    
    # Adjust that this value cannot be higher than the number of people patent before SMC deployment
    time_5 <- ifelse(time_5 < time_0, time_5, time_0)
    time_10 <- ifelse(time_10 < time_0, time_10, time_0)
    time_15 <- ifelse(time_15 < time_0, time_15, time_0)
    time_20 <- ifelse(time_20 < time_0, time_20, time_0)
    time_25 <- ifelse(time_25 < time_0, time_25, time_0)
    time_30 <- ifelse(time_30 < time_0, time_30, time_0)
    time_35 <- ifelse(time_35 < time_0, time_35, time_0)
    time_40 <- ifelse(time_40 < time_0, time_40, time_0)
    time_45 <- ifelse(time_45 < time_0, time_45, time_0)
    time_50 <- ifelse(time_50 < time_0, time_50, time_0)
    time_55 <- ifelse(time_55 < time_0, time_55, time_0)
    time_60 <- ifelse(time_60 < time_0, time_60, time_0)
    time_65 <- ifelse(time_65 < time_0, time_65, time_0)
    time_70 <- ifelse(time_70 < time_0, time_70, time_0)
    time_75 <- ifelse(time_75 < time_0, time_75, time_0)
    time_80 <- ifelse(time_80 < time_0, time_80, time_0)
    time_85 <- ifelse(time_85 < time_0, time_85, time_0)
    time_90 <- ifelse(time_90 < time_0, time_90, time_0)
    time_95 <- ifelse(time_95 < time_0, time_95, time_0)
    time_100 <- ifelse(time_100 < time_0, time_100, time_0)
    time_105 <- ifelse(time_105 < time_0, time_105, time_0)
    time_110 <- ifelse(time_110 < time_0, time_110, time_0)
    time_115 <- ifelse(time_115 < time_0, time_115, time_0)
    time_120 <- ifelse(time_120 < time_0, time_120, time_0)
    
    # Estime the number of nelwy infected individual at each timestep (can be only positive)
    Number_5 <- ifelse(time_5 < 0, 0, time_5)
    Number_10 <- ifelse(time_10 - time_5 < 0, 0, time_10 - time_5)
    Number_15 <- ifelse(time_15 - time_10 < 0, 0, time_15 - time_10)
    Number_20 <- ifelse(time_20 - time_15 < 0, 0, time_20 - time_15)
    Number_25 <- ifelse(time_25 - time_20 < 0, 0, time_25 - time_20)
    Number_30 <- ifelse(time_30 - time_25 < 0, 0, time_30 - time_25)
    Number_35 <- ifelse(time_35 - time_30 < 0, 0, time_35 - time_30)
    Number_40 <- ifelse(time_40 - time_35 < 0, 0, time_40 - time_35)
    Number_45 <- ifelse(time_45 - time_40 < 0, 0, time_45 - time_40)
    Number_50 <- ifelse(time_50 - time_45 < 0, 0, time_50 - time_45)
    Number_55 <- ifelse(time_55 - time_50 < 0, 0, time_55 - time_50)
    Number_60 <- ifelse(time_60 - time_55 < 0, 0, time_60 - time_55)
    Number_65 <- ifelse(time_65 - time_60 < 0, 0, time_65 - time_60)
    Number_70 <- ifelse(time_70 - time_65 < 0, 0, time_70 - time_65)
    Number_75 <- ifelse(time_75 - time_70 < 0, 0, time_75 - time_70)
    Number_80 <- ifelse(time_80 - time_75 < 0, 0, time_80 - time_75)
    Number_85 <- ifelse(time_85 - time_80 < 0, 0, time_85 - time_80)
    Number_90 <- ifelse(time_90 - time_85 < 0, 0, time_90 - time_85)
    Number_95 <- ifelse(time_95 - time_90 < 0, 0, time_95 - time_90)
    Number_100 <- ifelse(time_100 - time_95 < 0, 0, time_100 - time_95)
    Number_105 <- ifelse(time_105 - time_100 < 0, 0, time_105 - time_100)
    Number_110 <- ifelse(time_110 - time_105 < 0, 0, time_110 - time_105)
    Number_115 <- ifelse(time_115 - time_110 < 0, 0, time_115 - time_110)
    Number_120 <- ifelse(time_120 - time_115 < 0, 0, time_120 - time_115)
    
    # Estimate the median time of reinfection
    time_mean[i] <- median(
        c(rep(0, Number_5),
          rep(5, Number_10),
          rep(10, Number_15),
          rep(15, Number_20),
          rep(20, Number_25),
          rep(25, Number_30),
          rep(30, Number_35),
          rep(35, Number_40),
          rep(40, Number_45),
          rep(45, Number_50),
          rep(50, Number_55),
          rep(55, Number_60),
          rep(60, Number_65),
          rep(65, Number_70),
          rep(70, Number_75),
          rep(75, Number_80),
          rep(80, Number_85),
          rep(85, Number_90),
          rep(90, Number_95),
          rep(95, Number_100),
          rep(100, Number_105),
          rep(105, Number_110),
          rep(110, Number_115),
          rep(115, Number_120)
        )
      )
    
  }
  
  # Estimate the mean across each SMC deployment
  mean_rinfection <- mean(time_mean)
  
  # Return the mean
  return(mean_rinfection)
}

#--------------------------------------------------------------
# Estimate the selection coefficient of the resistant genotype.
#--------------------------------------------------------------

spread_R2 <- function(Output_data) {
  
  # Define the time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Select the time at which SMC is deployed
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  time_start <- Time_first_SMC
  
  # Select the end of the regression (9 year later)
  time_end <- Time_first_SMC + 9 * Number_survey_years # define time to stop the regression
  
  # Select the data within this boundary
  time_spread <- Output_data$Survey[Output_data$Survey > time_start & Output_data$Survey < time_end] / Number_survey_years # time spread is converted in years
  R_num <- Output_data$innoculationsPerAgeGroup_R1[Output_data$Survey > time_start & 
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_R2[Output_data$Survey > time_start &     
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_R3[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_R4[Output_data$Survey > time_start &
              Output_data$Survey < time_end]
  
  S_num <- Output_data$innoculationsPerAgeGroup_S1[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_S2[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_S3[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_S4[Output_data$Survey > time_start &
              Output_data$Survey < time_end]
  
  # Estimate the frequency of resistant genotype
  Inoculation_R <- R_num / (S_num + R_num)
  
  # Delete measurement that have to low frequency
  Measurment_R <- Inoculation_R[Inoculation_R <= 0.9 & Inoculation_R >= 0.3]
  time_spread <- time_spread[1:length(Measurment_R)]
  
  # Estimate the selection coefficient if we have at least two measurement
  if (length(Measurment_R) >= 2 & is.na(sum(Measurment_R)) == F) {
    
    # Estimate the slope fo the regression
    slope <- lm(log(Measurment_R / (1 - Measurment_R)) ~ time_spread)$coefficients[2]
    
    # Adjust the slope to be in number of parasite generation
    slope <- slope / 6
    
  } else {
    slope <- NA
  }
  
  # Return the slope
  return(slope)
}

#----------------------------------------------------------------------------------
# Estimate the selection coefficient of the resistant genotype via a moving average.
#----------------------------------------------------------------------------------

spread_R5 <- function(Output_data) {
  
  # Define the time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  time_start <- Time_first_SMC
  
  # Define the time at which SMC end (stop the regression)
  time_end <- Time_first_SMC + 9 * Number_survey_years #
  
  # Select the data within this boundary
  time_spread <- Output_data$Survey[Output_data$Survey > time_start & Output_data$Survey < time_end] / Number_survey_years # time spread is converted in years
  
  R_num <- Output_data$innoculationsPerAgeGroup_R1[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_R2[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_R3[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_R4[Output_data$Survey > time_start &
              Output_data$Survey < time_end]
  
  S_num <- Output_data$innoculationsPerAgeGroup_S1[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_S2[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_S3[Output_data$Survey > time_start &
              Output_data$Survey < time_end] + Output_data$innoculationsPerAgeGroup_S4[Output_data$Survey > time_start &
                                                                                                                                                                                                                                                                               Output_data$Survey < time_end]
  
  # Estimate the frequency of resistant genotype
  Inoculation_R <- R_num / (S_num + R_num)
  
  # Estimate the moving average of the frequency over 1 year (36*5 day before and 36*5 day after)
  Measurment_MA <- 0
  time_spread_MA <- 0
  for (i in ((1 + 36):(length(Inoculation_R) - 36))) {
    Measurment_MA[i - 36] <- mean(Inoculation_R[(i - 36):(i + 36)])
    time_spread_MA[i - 36] <- time_spread[i]
  }
  
  # Remove estimation in which the frequency of one genotype is to low
  Measurment_MA <- Measurment_MA[Measurment_MA >= 0.3 & Measurment_MA <= 0.9]
  time_spread_MA <- time_spread_MA[1:length(Measurment_MA)]
  
  # If at least two measurement estimate the selection coefficient
  if (length(Measurment_MA) >= 2 & is.na(sum(Measurment_MA)) == F) {
    
    # Estimate the slope of the regression
    slope <- lm(log(Measurment_MA / (1 - Measurment_MA)) ~ time_spread_MA)$coefficients[2]
    
    # Adjust the slope to be in generation time
    slope <- slope / 6
    
    } else {
      slope <- NA
    }
  
  # Return the slope
  return(slope)
}

#-----------------------------------------------
# Estimate the Efficacy when no SMC is deployed.
#-----------------------------------------------

Uncomplicated_0 <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1:10) {
    
    # Define the start and end
    time_start[i] <- Time_first_SMC - i * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
  }
  
  # Do the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  
  # Estimate number of individual at risk
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the protective effectiveness
  return(Mean_Uncomplicated / Mean_Person)

}

#-----------------------------------------------
# Estimate the Efficacy when  SMC is deployed.
#-----------------------------------------------

Uncomplicated_1 <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1:10) {
    
    # Define the start and end
    time_start[i] <- Time_first_SMC + (i - 1) * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the protective efficacy
  return(Mean_Uncomplicated / Mean_Person)
}

#--------------------------------------------------------------
# Estimate the Efficacy over one years when no SMC is deployed.
#---------------------------------------------------------------

Uncomplicated_0_year <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # define the time at which SMC begin
  Time_first_SMC <- 20 * Number_survey_years
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1:10) {
    
    # Define the start and end
    time_start[i] <- Time_first_SMC + (i - 1) * Number_survey_years
    time_end[i] <- time_start[i] + 73
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  Mean_Person <- mean(Person) * 12 # follow up for 12 months
  
  # Return the protective efficacy
  return(Mean_Uncomplicated / Mean_Person)
}

#-----------------------------------------------
# Estimate the Efficacy over one years when  SMC is deployed.
#-----------------------------------------------

Uncomplicated_1_year <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- 30 * Number_survey_years
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1:10) {
    
    # Define the start and end
    time_start[i] <- Time_first_SMC + (i - 1) * Number_survey_years
    time_end[i] <- time_start[i] + 73
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  Mean_Person <- mean(Person) * 12 #follow up for 12 months
  
  # Return the protective efficacy
  return(Mean_Uncomplicated / Mean_Person)
}

#-----------------------------------------------
# Estimate the Efficacy when no SMC is deployed over one years
#-----------------------------------------------

Uncomplicated_0_round <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1) {
    
    # Define the start and end
    time_start[i] <- Time_first_SMC - i * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  
  # Estimate number of individual at risk
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the protective efficacy
  return(Mean_Uncomplicated / Mean_Person)
}

#-----------------------------------------------
# Estimate the Efficacy when  SMC is deployed over one year
#-----------------------------------------------

Uncomplicated_1_round <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1) {
    
    # Define the start and end
    time_start[i] <- Time_first_SMC + (i - 1) * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the protective efficacy
  return(Mean_Uncomplicated / Mean_Person)
}

#-----------------------------------------------
# Estimate the number of person when no SMC is deployed.
#-----------------------------------------------

Person_0 <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1:10) {
    
    # Define the begining and the end
    time_start[i] <- Time_first_SMC - i * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  
  # Estimate the number of individual at risk
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the number of individual at risk
  return(Mean_Person)
}

#-----------------------------------------------
# Estimate the number of person when  SMC is deployed.
#-----------------------------------------------

Person_1 <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 year sooner
  for (i in 1:10) {
    
    # Defin the begining and the end
    time_start[i] <- Time_first_SMC + (i - 1) * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    
  }
  
  # Estimate the mean 
  Mean_Uncomplicated <- mean(Uncomplicated)
  
  # Estimate the number of individual at risk
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the number of individual at risk
  return(Mean_Person)
}

#------------------------------------------------------------------
# Estimate the number of person when SMC is deployed over one year
#------------------------------------------------------------------

Person_1_round <- function(Output_data) {
  
  # Define time step
  time_step <- 5
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step
  
  # Define the time at which SMC begin
  Time_first_SMC <- Output_data$Survey[Output_data$nMDAs_2 >= 1][1] # time start after one generation
  
  time_start <- 0
  time_end <- 0
  Uncomplicated <- 0
  Person <- 0
  
  # Select the time 10 years sooner
  for (i in 1) {
    
    # Define the begining and the end
    time_start[i] <- Time_first_SMC + (i - 1) * Number_survey_years
    time_end[i] <- time_start[i] + 30 * 4 / time_step
    
    # Estimate the total number of uncomplicated case during the four month of SMC deployment
    Uncomplicated[i] <- sum(Output_data$nUncomp_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    Person[i] <- mean(Output_data$nHost_2[Output_data$Survey >= time_start[i] & Output_data$Survey <= time_end[i]])
    
  }
  
  # Estimate the mean
  Mean_Uncomplicated <- mean(Uncomplicated)
  
  # Estimate the number of individual at risk
  Mean_Person <- mean(Person) * 4 #follow up of for month
  
  # Return the number of individual at risk
  return(Mean_Person)
}
