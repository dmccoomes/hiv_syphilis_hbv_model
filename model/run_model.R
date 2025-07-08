
#### Title: Run HIV-syphilis-HBV maternal testing models
#### Purpose: Run testing models for HIV, syphilis, and HBV among pregnant populations and collect output
#### Author: David Coomes, Akash Malhotra
#### Date: June 26, 2025


## This script will read in a parameter worksheet with population parameters by country and population
## then calculate maternal transmission of HIV, syphilis, and HBV depending on parameters in the worksheet,
## adjustable parameters, and testing algorithm. 
## It will output a table of deaths, infections, and costs, and calculate ICERs


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table",
  "readxl"
)

rm(list = ls())

source(here("model", "supplemental", "life_table.R"))
source(here("model", "syphilis_model", "syphilis_model.R"))
source(here("model", "hiv_model", "hiv_model_final.R"))
source(here("model", "child_model", "infant_model.R"))
source(here("model", "model_output", "hiv_syphilis_hbv_model_output.R"))



## Set up testing algorithm
  # this one compares a single early test to repeat testing using dual HIV-syphilis tests and RDT plus lab
testing_algorithm <- list(
  two_dual_test = list(c(1,1), c("dual", "dual")),
  two_rdt_test = list(c(1,1), c("rdt", "rdt")),
  one_dual_test = list(c(1,0), c("dual", "dual")),
  one_rdt_test = list(c(1,0), c("rdt", "rdt"))
)


out <- cost_effective_func(testing_algorithm, syph_test = TRUE, baseline = "one_rdt_test",
                           here("parameters", "inp_par_dmc.xlsx"), "HIV_Syphilis_Dual_Preg",
                           "Kenya", "Pregnant population")


write.csv(out, here("output", "dual_model_output.csv"))




### NEEDS
# Get CE output
# Wrap all HIV functions into a single function
# Wrap all syphilis functions into a single function
# Wrap HIV and syphilis functions into a single function
# set up sensitivity analysis
# set up stochastic process - wrap this into a function



