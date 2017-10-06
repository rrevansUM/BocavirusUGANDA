# HH-transmission-model
Code for transmission hazard model used in Petrie et al. to examine household and community transmission of influenza using the HIVE study. Code is written in R and includes three main files:
* THmodel_premcmc_commented.R - setup file to organize the data to feed into the MCMC/model
* THmodel_mcmc_commented.R - the actual file that runs the MCMC
* THmodel_output_commented.R - post-processing to make nice output (MCMC diagnostics & outcomes)
* THmodel_infectionSimulator_commented.R - simulates community and household infections using parameters estimated from MCMC, weekly community infection proxy data, and household structure and covariate data (without observed infection outcome data)
* THmodel_exampleData.csv - example data set with individual ID (ID), household ID (hhID), age categories (AGECAT2 == 9-17 years; AGECAT3 >= 18 years), high risk health status (HR), vaccination status (vax), vaccination date (season day of follow-up; 0 == vaccinated prior to follow up; 999 == unvaccinated), and duration of influenza season follow up (endfudate). This data can be used to simulate infections using THmodel_infectionSimulator_commented, the resulting data formatted using TH_premcmc_commented, and model parameters estimated using THmodel_mcmc_commented. 
*  THmodel_exampleComData.csv - example data set with weekly community infection proxy data standardized to the peak week of infection. This data set is used both parameter estimation as well as simulation.
