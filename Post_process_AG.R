######################################################################################
#     Post process of the data                                                       #
#                                                                                    #
# Task : Reorganize the output of OpenMalaria in a table                             #
#        that display all the measurement at each time step                          #
#                                                                                    #
# Input: Outpout file form Openmalaria                                               #
#                                                                                    #
# Output: Table of each measurement values at each time step                         #
#                                                                                    #
# author: thiery.masserey@swisstph.ch                                                #
######################################################################################

#--------------------------------------------------------
# function to post-process the outputdata of OpenMalaria.
#--------------------------------------------------------
Postprocess <- function(pm) {
  
  # Message to the console
  message("  - Post Processing ")
  
  # Define list of output files
  Output <- list.files(pm$pth$sim_out) # Make a list of all the output file
  
  # Initiate progress bar
  pb <- txtProgressBar(min = 0, max = length(Output), initial = 0, width = 100, style = 3)
  
  # Start the loop that will do the post process for each output file
  for (i in 1:length(Output)) {
    
    # Define file name, file path, and download it i<-1
    Output_data_name <- Output[i]
    Output_data_file_path <- file.path(pm$pth$sim_out, Output_data_name)
    Output_data <- read.table(Output_data_file_path)
    
    # Extract the survey Time in a new dataframe
    a <- as.data.frame(table(Output_data$V1))
    Survey <- as.numeric(a$Var1)
    Output_data_2 <- as.data.frame(Survey)
    
    # Extract the survey measure number
    Indicator <- Output_data$V3[Output_data$V1 == 1]
    
    # Extract information about age group and resistant/sensitive parasite
    Age_groupe <- Output_data$V2[Output_data$V1 == 1]
    
    # Define the meaning of each measurement
    Indicators_meaning <- c(
      
      "nHost_1", # 0 The number of host in each age group
      "nHost_2",
      "nHost_3",
      "nHost_4",
      
      "nInfect_1", # 1 The number of human hosts with an infection
      "nInfect_2",
      "nInfect_3",
      "nInfect_4",
      
      "nExpectd_1", # 2 The expected number of infected hosts
      "nExpectd_2",
      "nExpectd_3",
      "nExpectd_4",
      
      "nPatent_1", # 3 The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
      "nPatent_2",
      "nPatent_3",
      "nPatent_4",
      
      "nTransmit", # 7 Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes.
      
      "nTreatments1_1", # 11 The number of treatments (1st line)
      "nTreatments1_2",
      "nTreatments1_3",
      "nTreatments1_4",
      
      "nUncomp_1", # 14 The number of episodes uncomplicated
      "nUncomp_2",
      "nUncomp_3",
      "nUncomp_4",
      
      "innoculationsPerAgeGroup_S1", # 30 The number of inoculation per genotype
      "innoculationsPerAgeGroup_R1",
      "innoculationsPerAgeGroup_S2",
      "innoculationsPerAgeGroup_R2",
      "innoculationsPerAgeGroup_S3",
      "innoculationsPerAgeGroup_R3",
      "innoculationsPerAgeGroup_S4",
      "innoculationsPerAgeGroup_R4",
      
      "Vector_Nv", # 32 Host seeking mosquito population size at this time step, species 1
      
      "inputEIR", # 35 The input EIR
      
      "simulatedEIR", # 36 The simulated EIR
      
      "nMDAs_1", # 52 The number of drug doses given via mass deployment (MDA or screen&treat)
      "nMDAs_2",
      "nMDAs_3",
      "nMDAs_4",
      
      "nTreatDiagnostics_1", # 64 The number of diagnostic tests performed (if in the health system description, useDiagnosticUC="true").
      "nTreatDiagnostics_2",
      "nTreatDiagnostics_3",
      "nTreatDiagnostics_4",
      
      "nInfectByGenotype_S1", # 69 The number of infection by genotype
      "nInfectByGenotype_R1",
      "nInfectByGenotype_S2",
      "nInfectByGenotype_R2",
      "nInfectByGenotype_S3",
      "nInfectByGenotype_R3",
      "nInfectByGenotype_S4",
      "nInfectByGenotype_R4",
      
      "nPatentByGenotype_S1", # 70 The number of patent infection by genotype
      "nPatentByGenotype_R1",
      "nPatentByGenotype_S2",
      "nPatentByGenotype_R2",
      "nPatentByGenotype_S3",
      "nPatentByGenotype_R3",
      "nPatentByGenotype_S4",
      "nPatentByGenotype_R4",
      
      "nHostDrugConcNonZero_A1", # 72 For each drug type in the pharmacology section of the XML, report the number of humans with non-zero concentration
      "nHostDrugConcNonZero_B1",
      "nHostDrugConcNonZero_C1",
      "nHostDrugConcNonZero_D1",
      "nHostDrugConcNonZero_E1",
      "nHostDrugConcNonZero_A2", 
      "nHostDrugConcNonZero_B2",
      "nHostDrugConcNonZero_C2",
      "nHostDrugConcNonZero_D2",
      "nHostDrugConcNonZero_E2",
      "nHostDrugConcNonZero_A3",
      "nHostDrugConcNonZero_B3",
      "nHostDrugConcNonZero_C3",
      "nHostDrugConcNonZero_D3",
      "nHostDrugConcNonZero_E3",
      "nHostDrugConcNonZero_A4", 
      "nHostDrugConcNonZero_B4",
      "nHostDrugConcNonZero_C4",
      "nHostDrugConcNonZero_D4",
      "nHostDrugConcNonZero_E4",
      
      "sumLogDrugConcNonZero_A1", # 73 For each drug type in the pharmacology section of the XML, report the sum of the natural logarithm of the drug concentration in hosts with non-zero concentration.
      "sumLogDrugConcNonZero_B1",
      "sumLogDrugConcNonZero_C1",
      "sumLogDrugConcNonZero_D1",
      "sumLogDrugConcNonZero_E1",
      "sumLogDrugConcNonZero_A2", 
      "sumLogDrugConcNonZero_B2",
      "sumLogDrugConcNonZero_C2",
      "sumLogDrugConcNonZero_D2",
      "sumLogDrugConcNonZero_E2",
      "sumLogDrugConcNonZero_A3", 
      "sumLogDrugConcNonZero_B3",
      "sumLogDrugConcNonZero_C3",
      "sumLogDrugConcNonZero_D3",
      "sumLogDrugConcNonZero_E3",
      "sumLogDrugConcNonZero_A4", 
      "sumLogDrugConcNonZero_B4",
      "sumLogDrugConcNonZero_C4",
      "sumLogDrugConcNonZero_D4",
      "sumLogDrugConcNonZero_E4")
    
    # Define the meaning of each measurement if aimed to estimate the length of the prophylactic period
    if (pm$opts$Type_analysis == "Propylaxis") {
      Indicators_meaning <- c(
        
        "nHost_1",# 0 The number of host in each age group
        "nHost_2",
        "nHost_3",
        "nHost_4",
        
        "nInfect_1", # 1 The number of human hosts with an infection
        "nInfect_2",
        "nInfect_3",
        "nInfect_4",
        
        "nExpectd_1", # 2 The expected number of infected hosts
        "nExpectd_2",
        "nExpectd_3",
        "nExpectd_4",
        
        "nPatent_1", #  3 The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
        "nPatent_2",
        "nPatent_3",
        "nPatent_4",
        
        "nTransmit", # 7 Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes.
        
        "nTreatments1_1", # 11 The number of treatments (1st line)
        "nTreatments1_2",
        "nTreatments1_3",
        "nTreatments1_4",
        
        "nUncomp_1", # 14 The number of episodes uncomplicated1
        "nUncomp_2",
        "nUncomp_3",
        "nUncomp_4",
        
        "Vector_Nv", # 32 Host seeking mosquito population size at this time step, species 1
        
        "inputEIR", # 35 The input EIR
        
        "simulatedEIR", # 36 The simulated EIR
        
        "nMDAs_1", # 52 Number of drug doses given via mass deployment (MDA or screen&treat)
        "nMDAs_2",
        "nMDAs_3",
        "nMDAs_4",
        
        "nTreatDiagnostics_1", # 64 The number of diagnostic tests performed (if in the health system description, use DiagnosticUC="true").
        "nTreatDiagnostics_2",
        "nTreatDiagnostics_3",
        "nTreatDiagnostics_4",
        
        "nHostDrugConcNonZero_A1", # 72 For each drug type in the pharmacology section of the XML, report the number of humans with non-zero concentration
        "nHostDrugConcNonZero_B1",
        "nHostDrugConcNonZero_C1",
        "nHostDrugConcNonZero_D1",
        "nHostDrugConcNonZero_A2", 
        "nHostDrugConcNonZero_B2",
        "nHostDrugConcNonZero_C2",
        "nHostDrugConcNonZero_D2",
        "nHostDrugConcNonZero_A3", 
        "nHostDrugConcNonZero_B3",
        "nHostDrugConcNonZero_C3",
        "nHostDrugConcNonZero_D3",
        "nHostDrugConcNonZero_A4", 
        "nHostDrugConcNonZero_B4",
        "nHostDrugConcNonZero_C4",
        "nHostDrugConcNonZero_D4",
        
        "sumLogDrugConcNonZero_A1", # 73 For each drug type in the pharmacology section of the XML, report the sum of the natural logarithm of the drug concentration in hosts with non-zero concentration.
        "sumLogDrugConcNonZero_B1",
        "sumLogDrugConcNonZero_C1",
        "sumLogDrugConcNonZero_D1",
        "sumLogDrugConcNonZero_A2",
        "sumLogDrugConcNonZero_B2",
        "sumLogDrugConcNonZero_C2",
        "sumLogDrugConcNonZero_D2",
        "sumLogDrugConcNonZero_A3",
        "sumLogDrugConcNonZero_B3",
        "sumLogDrugConcNonZero_C3",
        "sumLogDrugConcNonZero_D3",
        "sumLogDrugConcNonZero_A4",
        "sumLogDrugConcNonZero_B4",
        "sumLogDrugConcNonZero_C4",
        "sumLogDrugConcNonZero_D4"
      )
    }
    
    # Define the meaning of each measurement if aimed to replicate the trial of Zongo et al
    if (pm$opts$Type_analysis == "Trial") {
      Indicators_meaning <- c(
        
        "nHost_1", # 0 The number of host in each age group
        "nHost_2",
        "nHost_3",
        "nHost_4",
        "nHost_5",
        "nHost_6",
        "nHost_7",
        
        "nInfect_1", # 1 The number of human hosts with an infection
        "nInfect_2",
        "nInfect_3",
        "nInfect_4",
        "nInfect_5",
        "nInfect_6",
        "nInfect_7",
        
        "nExpectd_1", # 2 The expected number of infected hosts
        "nExpectd_2",
        "nExpectd_3",
        "nExpectd_4",
        "nExpectd_5",
        "nExpectd_6",
        "nExpectd_7",
        
        "nPatent_1", # 3 The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
        "nPatent_2",
        "nPatent_3",
        "nPatent_4",
        "nPatent_5",
        "nPatent_6",
        "nPatent_7",
        
        "nTransmit", # 7 Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes.
        
        "nTreatments1_1", #  11 The number of treatments (1st line)
        "nTreatments1_2",
        "nTreatments1_3",
        "nTreatments1_4",
        "nTreatments1_5",
        "nTreatments1_6",
        "nTreatments1_7",
        
        "nUncomp_1", # 14 The number of episodes uncomplicated1
        "nUncomp_2",
        "nUncomp_3",
        "nUncomp_4",
        "nUncomp_5",
        "nUncomp_6",
        "nUncomp_7",
        
        "nSevere_1", # 15 The number of sever episodes
        "nSevere_2",
        "nSevere_3",
        "nSevere_4",
        "nSevere_5",
        "nSevere_6",
        "nSevere_7",
        
        "nDirDeaths_1", # 15 The number of death
        "nDirDeaths_2",
        "nDirDeaths_3",
        "nDirDeaths_4",
        "nDirDeaths_5",
        "nDirDeaths_6",
        "nDirDeaths_7",
        
        "innoculationsPerAgeGroup_S1", # 30 The number of inoculation by Genotype
        "innoculationsPerAgeGroup_R1",
        "innoculationsPerAgeGroup_S2",
        "innoculationsPerAgeGroup_R2",
        "innoculationsPerAgeGroup_S3",
        "innoculationsPerAgeGroup_R3",
        "innoculationsPerAgeGroup_S4",
        "innoculationsPerAgeGroup_R4",
        "innoculationsPerAgeGroup_S5",
        "innoculationsPerAgeGroup_R5",
        "innoculationsPerAgeGroup_S6",
        "innoculationsPerAgeGroup_R6",
        "innoculationsPerAgeGroup_S7",
        "innoculationsPerAgeGroup_R7",
        
        "Vector_Nv_1", # 32 Host seeking mosquito population size at this time step, species 1
        "Vector_Nv_2",
        "Vector_Nv_3",
        "Vector_Nv_4",
        
        "inputEIR", # 35 The input EIR
        "simulatedEIR", # 36 The simulated EIR
        
        "nMassITNs_1", # 44 The number of ITN given via mass deployment
        "nMassITNs_2",
        "nMassITNs_3",
        "nMassITNs_4",
        "nMassITNs_5",
        "nMassITNs_6",
        "nMassITNs_7",
        
        "nMDAs_1", # 52 The number of drug doses given via mass deployment (MDA or screen&treat)
        "nMDAs_2",
        "nMDAs_3",
        "nMDAs_4",
        "nMDAs_5",
        "nMDAs_6",
        "nMDAs_7",
        
        "nCtsMDA_1", # 59 The number of drug doses given via age based deployment (MDA or screen&treat)
        "nCtsMDA_2",
        "nCtsMDA_3",
        "nCtsMDA_4",
        "nCtsMDA_5",
        "nCtsMDA_6",
        "nCtsMDA_7",
        
        "nInfectByGenotype_S1", # 69 The number of infection by Genotype
        "nInfectByGenotype_R1",
        "nInfectByGenotype_S2",
        "nInfectByGenotype_R2",
        "nInfectByGenotype_S3",
        "nInfectByGenotype_R3",
        "nInfectByGenotype_S4",
        "nInfectByGenotype_R4",
        "nInfectByGenotype_S5",
        "nInfectByGenotype_R5",
        "nInfectByGenotype_S6",
        "nInfectByGenotype_R6",
        "nInfectByGenotype_S7",
        "nInfectByGenotype_R7",
        
        "nPatentByGenotype_S1", # 70 The number of patent infection by Genotype
        "nPatentByGenotype_R1",
        "nPatentByGenotype_S2",
        "nPatentByGenotype_R2",
        "nPatentByGenotype_S3",
        "nPatentByGenotype_R3",
        "nPatentByGenotype_S4",
        "nPatentByGenotype_R4",
        "nPatentByGenotype_S5",
        "nPatentByGenotype_R5",
        "nPatentByGenotype_S6",
        "nPatentByGenotype_R6",
        "nPatentByGenotype_S7",
        "nPatentByGenotype_R7",
        
        "nHostDrugConcNonZero_A1", # 72 For each drug type in the pharmacology section of the XML, report the number of humans with non-zero concentration
        "nHostDrugConcNonZero_B1",
        "nHostDrugConcNonZero_C1",
        "nHostDrugConcNonZero_D1",
        "nHostDrugConcNonZero_E1",
        "nHostDrugConcNonZero_A2", 
        "nHostDrugConcNonZero_B2",
        "nHostDrugConcNonZero_C2",
        "nHostDrugConcNonZero_D2",
        "nHostDrugConcNonZero_E2",
        "nHostDrugConcNonZero_A3",
        "nHostDrugConcNonZero_B3",
        "nHostDrugConcNonZero_C3",
        "nHostDrugConcNonZero_D3",
        "nHostDrugConcNonZero_E3",
        "nHostDrugConcNonZero_A4", 
        "nHostDrugConcNonZero_B4",
        "nHostDrugConcNonZero_C4",
        "nHostDrugConcNonZero_D4",
        "nHostDrugConcNonZero_E4",
        "nHostDrugConcNonZero_A5", 
        "nHostDrugConcNonZero_B5",
        "nHostDrugConcNonZero_C5",
        "nHostDrugConcNonZero_D5",
        "nHostDrugConcNonZero_E5",
        "nHostDrugConcNonZero_A6", 
        "nHostDrugConcNonZero_B6",
        "nHostDrugConcNonZero_C6",
        "nHostDrugConcNonZero_D6",
        "nHostDrugConcNonZero_E6",
        "nHostDrugConcNonZero_A7",
        "nHostDrugConcNonZero_B7",
        "nHostDrugConcNonZero_C7",
        "nHostDrugConcNonZero_D7",
        "nHostDrugConcNonZero_E7",
        
        "sumLogDrugConcNonZero_A1", # 73 For each drug type in the pharmacology section of the XML, report the sum of the natural logarithm of the drug concentration in hosts with non-zero concentration.
        "sumLogDrugConcNonZero_B1",
        "sumLogDrugConcNonZero_C1",
        "sumLogDrugConcNonZero_D1",
        "sumLogDrugConcNonZero_E1",
        "sumLogDrugConcNonZero_A2", 
        "sumLogDrugConcNonZero_B2",
        "sumLogDrugConcNonZero_C2",
        "sumLogDrugConcNonZero_D2",
        "sumLogDrugConcNonZero_E2",
        "sumLogDrugConcNonZero_A3", 
        "sumLogDrugConcNonZero_B3",
        "sumLogDrugConcNonZero_C3",
        "sumLogDrugConcNonZero_D3",
        "sumLogDrugConcNonZero_E3",
        "sumLogDrugConcNonZero_A4", 
        "sumLogDrugConcNonZero_B4",
        "sumLogDrugConcNonZero_C4",
        "sumLogDrugConcNonZero_D4",
        "sumLogDrugConcNonZero_E4",
        "sumLogDrugConcNonZero_A5",
        "sumLogDrugConcNonZero_B5",
        "sumLogDrugConcNonZero_C5",
        "sumLogDrugConcNonZero_D5",
        "sumLogDrugConcNonZero_E5",
        "sumLogDrugConcNonZero_A6",
        "sumLogDrugConcNonZero_B6",
        "sumLogDrugConcNonZero_C6",
        "sumLogDrugConcNonZero_D6",
        "sumLogDrugConcNonZero_E6",
        "sumLogDrugConcNonZero_A7", 
        "sumLogDrugConcNonZero_B7",
        "sumLogDrugConcNonZero_C7",
        "sumLogDrugConcNonZero_D7",
        "sumLogDrugConcNonZero_E7"
      )
    }
    
    # Extract the information for each measurement at each timestep, and save it in the new dataframe (Output_data_2)
    for (j in 1:(length(Indicator))) {
      Output_data_2[j + 1] <- Output_data$V4[Output_data$V3 == Indicator[j] & Output_data$V2 == Age_groupe[j]] # put in the j +1 colone (due that we have one colone of the time survey), the value of the jth indicator
      colnames(Output_data_2)[j + 1] <- Indicators_meaning[j] # give the name of the indicators
    }
    
    # Save the new dataset in the postprocess folder
    Postprocess_data_name <- paste0("PostProcess_", Output_data_name)
    Postprocess_data_path <- file.path(pm$pth$processed, Postprocess_data_name)
    write.table(Output_data_2, file = Postprocess_data_path, sep = ";")
    
    # Update progress bar
    setTxtProgressBar(pb, i)
    
  }
  
  # close the progress bar
  close(pb)
}
