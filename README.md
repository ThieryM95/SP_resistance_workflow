This is the first release of the analysis code for:
Seasonal malaria chemoprevention potentially remains effective with the spread of Plasmodium falciparum parasites resistant to sulfadoxine-pyrimethamine: a modelling study   

Thiery Masserey1,2, Tamsin Lee1,2, Sherrie L Kelly1,2, Ian M Hastings3, Melissa   A Penny1,2 
1 Swiss Tropical and Public Health Institute, Basel, Switzerland 
2 University of Basel, Basel, Switzerland 	
3 Liverpool School of Tropical Medicine, Liverpool, UK
Correspondence to: Prof Melissa A. Penny, melissa.penny@unibas.ch


In the study mentioned above, we have developed a disease modelling approach with emulator-based global sensitivity analysis to quantify the influence of multiple factors on the spread of Plasmodium falciparum parasites resistant to sulfadoxine-pyrimethamine (SP) (quintuple mutant) and predicted the time needed for the quintuple mutant to spread from 1% to 50% of inoculations for several SMC deployment strategies. We also estimated the impact of this spread on SMC effectiveness against clinical malaria.
This analysis was performed using an individual-based malaria model (https://github.com/SwissTPH/openmalaria/wiki). 

The folder "WF_SMC" contains the code used to:
1) Assess the rate of spread of parasites resistant to SP when used for seasonal malaria chemoprevention (SMC) and quantify the influence of multiple factors on the rate of spread via global sensitivity analysis. 
   The analysis involved the following steps:
	(i) randomly sampling combinations of parameters of the factors of interest;
	(ii) simulating and estimating the rate of spread of the SP-resistant mutant for each parameter combination in OpenMalaria;
	(iii) training an HGP to learn the relationship between the input and output with iterative improvements to fitting through adaptive sampling;
	(iv) performing a global sensitivity analysis based on the Sobol variance decomposition. 
   The script "launch.R" calls all the necessary functions to perform all the analysis steps specified in the script "Option.R". The code asks to specify the analysis details in the script "Option.R". 
                                                                              
2) Assess the prophylactic period conferred by SP for different EC50 values of SP. The script "launch.R" calls all the necessary functions to perform all the analysis steps specified in the script "Option.R". The code asks to specify the analysis details in the script "Option.R". 

3) Assess the protective effectiveness of SMC against clinical cases. The script "launch.R" calls all the necessary functions to perform all the analysis steps specified in the script "Option.R". The code asks to specify the analysis details in the script "Option.R". 
         
4) Replicate the trial of Zongo et al. (2015). The script "launch.R" calls all the necessary functions to perform all the analysis steps specified in the script "Option.R". The code asks to specify the analysis details in the script "Option.R". 
            
Note that the code uses the folder structure and working directories used by the researchers. The file paths will have to be adjusted.