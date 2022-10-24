#################################################################################
# Workflow to run four different analysis (need to be defin in Option.R):       #
#                                                                               #
# 1) Assess the rate of spread of parasite resistant to SP when used for SMC    #
#    (use all the steps of the workflow)                                        #
#                                                                               #
# 2) Assess the prophylactic period confer by SP for different EC50 values of SP#
#    (use step 1 to 7 of the workflow)                                          #
#                                                                               #
# 3) Assess the protective effectiveness of SMC against clinical cases          #
#    (use step 1 to 7 of the workflow)                                          #
#                                                                               #
# 4) Replicate the trial of Zongo et al. (2015).                                #
#    (use step 1 to 7 of the workflow)                                          #
#                                                                               #
# Run on OpenMalaria V43.0                                                      #
#                                                                               #
# Author: thiery Masserey (thiery.masserey@swisstph.ch)                         #
#################################################################################

# Clear global environment
rm(list = ls())

# Set the working directory
setwd("/scicore/home/penny/masthi00/smc_resistance")

# Download functions in the environment
library("tgp")
library("styler")
source("myRfunctions.R")
source("directories.R")
source("Option.R")
source("Generate_parameter_table.R")
source("openmalaria_setup_seed.R")
source("openmalaria_run.R")
source("Post_process_AG.R")
source("summary_results.R") #  CHECK WHICH function can be deleted ##############
source("GP.R")
source("adaptative_sample.R")
source("sensitivity_analysis.R")
source("Post_process_GP.R")
source("Post_process_sensitivity_analysis.R")

# Tidy up
# Clear console
if (interactive())
  clc() # see myRfunctions.R

# Close figures
if (interactive())
  clf() # see myRfunctions.R


#########################
# Start of the workflow #
#########################

# 1) Create directory, and options :
pm <- set_options() # See Option.R

# 2) Generate Latin Hypercube Sample and list of scenario
generate_param_table(pm, param_table = as_param_table) # see generate_parameter_table_2.R

# 3) Generate seed patterns
n_jobs_unique <-  simulation_sed_patterns(pm) # see Open_malaria_setup.R

message("   ~ Total number of scenarios: ", 
        thou_sep(n_jobs_unique),
        " (with seed)")

# 4) Generate scenario xml files from the base xml
generate_xml_files(
  pm,
  n_jobs = n_jobs_unique,
  file_path = pm$pth$xml_base,
  sed_path = pm$pth$sim_sed,
  xml_path = pm$pth$sim_xml) # see Open_malaria_setup.R

# 5) Run simulation
run_model(
  pm,
  n_jobs_unique,
  xml_path = pm$pth$sim_xml,
  sim_path = pm$pth$sim_out) # see Openmalaria_run.R

# 6) Post processing
Postprocess(pm) # see Post_process_AG.R

# 7) Summarized the results
SummaryResults(pm) # see Summarized_results.R

# 8) fit the Gaussian process
Results_gp_1 <- run_gp(pm) # see run GP.R

# 9) Run adaptive sampling
Results_gp <- run_adaptive_sampling(pm) # see adaptive_sample.R

# 10) Post process to visualize the results of the GP in a table format
Precision_final <- Post_process_GP(Results_gp) # see Post_process_GP.R

# 11) Sensitivity analysis
Results_SA <- sensitivity_analysis(pm) # see sensitivity_analysis.R

# 12) Post process the result of the sensitivity analysis in a table format
# Function to transform the sobol indies into a table format
data <- Post_process_sensitivity(Results_SA) # see Post_process_sensitivity_analysis.R

# Function to transform the direction of the effect of parameter on the selection coeffcient in a table format
Quantil_final_final <- Post_process_sensitivity_2(Results_SA) # see Post_process_sensitivity_analysis.R

# see the visualize folder to visualize the results