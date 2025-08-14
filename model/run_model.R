
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
source(here("model", "population_model", "infant_model.R"))
source(here("model", "population_model", "maternal_model.R"))
source(here("model", "model_output", "hiv_syphilis_hbv_model_output.R"))
source(here("model", "model_output", "triplex_output.R"))



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




#############
# TRIPLEX MODEL OUTPUT
#############


rm(list = ls())

source(here("model", "supplemental", "life_table.R"))
source(here("model", "syphilis_model", "syphilis_model.R"))
source(here("model", "hiv_model", "hiv_model_final.R"))
source(here("model", "population_model", "infant_model.R"))
source(here("model", "population_model", "maternal_model.R"))
source(here("model", "model_output", "hiv_syphilis_hbv_model_output.R"))
source(here("model", "model_output", "triplex_output.R"))


# First, get population parameters
pop_param <- pop_param_func(here("parameters", "inp_par_dmc.xlsx"), "HIV_Syphilis_Dual_Preg",
                            "Kenya", "Pregnant population")

list2env(pop_param, envir = parent.frame())

# Load mortality file
mort_df <- read_csv(here("parameters", "kenya_mort_ihme.csv"))




triplex_out <- triplex_out_func(late_ANC_test = 0,
                             syph_test = TRUE,
                             test_type_1 = "rdt",
                             test_type_2 = "rdt")

child_hiv_syph <- child_db$child_outcomes %>%
  select(-c(discounted_ylls, discounted_ylds,
            discounted_dalys))
# 
# write.csv(child_hiv_syph, here("output", "dual_child_model_base.csv"))


# Get different outputs
  # Base case with 1 antenatal test using rdt and lab
triplex_out_base <- triplex_out_func(first_ANC_test = 0,
                                     late_ANC_test = 0,
                                     syph_test = TRUE,
                                     test_type_1 = "rdt",
                                     test_type_2 = "rdt")

triplex_out_test_hivrdt <- triplex_out_func(first_ANC_test = 1, 
                                            late_ANC_test = 0,
                                            syph_test = TRUE,
                                            test_type_1 = "rdt",
                                            test_type_2 = "rdt")

triplex_out_test_dualrdt <- triplex_out_func(first_ANC_test = 1, 
                                            late_ANC_test = 0,
                                            syph_test = TRUE,
                                            test_type_1 = "dual",
                                            test_type_2 = "dual")

# Put them together and write to CSV
child_model <- triplex_out_base$child_outcomes %>%
  mutate(test = "no test") %>%
  bind_rows(triplex_out_test_hivrdt$child_outcomes %>%
              mutate(test = "HIV RDT")) %>%
  bind_rows(triplex_out_test_dualrdt$child_outcomes %>%
              mutate(test = "HIV-syph dual RDT"))

# Write to csv
write_csv(child_model, here("output", "hiv_syphilis_child_outcomes.csv"))


# Put together maternal outcomes and write to csv
maternal_model <- triplex_out_base$maternal_outcomes %>%
  mutate(test = "no test") %>%
  bind_rows(triplex_out_test_hivrdt$maternal_outcomes %>%
              mutate(test = "HIV RDT")) %>%
  bind_rows(triplex_out_test_dualrdt$maternal_outcomes %>%
              mutate(test = "HIV-syph dual RDT"))

# Write to csv
write_csv(maternal_model, here("output", "hiv_syphilis_maternal_outcomes.csv"))






