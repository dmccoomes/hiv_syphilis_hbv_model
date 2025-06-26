
#### Title: HIV-syphilis-HBV maternal testing model output
#### Purpose: Run testig models for HIV, syphilis, and HBV among pregnant populations and collect output
#### Author: David Coomes, Akash Malhotra
#### Date: June 3, 2025


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table",
  "readxl"
)

rm(list = ls())

source(here("code", "supplemental", "life_table.R"))
source(here("code", "syphilis_model", "syphilis_model.R"))
source(here("code", "hiv_model", "hiv_model_final.R"))
source(here("code", "child_model", "infant_model.R"))


## This script will read in a parameter worksheet with population parameters by country and population
  ## then calculate maternal transmission of HIV, syphilis, and HBV depending on parameters in the worksheet,
  ## adjustable parameters, and testing algorithm. 
  ## It will output a table of deaths, infections, and costs, and calculate ICERs





##############
# Parameters #
##############

# Reading in parameters to set up model - PUT IN FUNCTION LATER
# Read in population parameters
pop_parameters <- readxl::read_xlsx(here("parameters", "inp_par_dmc.xlsx"), 
                                    sheet = "pop_characteristics") %>%
  mutate(across(value:value_ub,  ~as.double(.)))

# Read in test characteristics
test_parameters <- readxl::read_xlsx(here("parameters", "inp_par_dmc.xlsx"), 
                                     sheet = "test_characteristics")


# only for the selected combinations
selected_params <- pop_parameters %>%
  filter(model == "HIV_Syphilis_Dual_Preg", #will make a loop later
         country == "Kenya", #will make a loop later
         target_pop == "Pregnant population") %>% #will make a loop later
  bind_rows(test_parameters)

# Set up all parameters for the selected population
for (i in 1:nrow(selected_params)) {
  param_name <- as.character(selected_params$var_name[i])
  assign(param_name, selected_params$value[i])
  # assign(paste0(param_name, "_lb"), selected_params$value_lb[i])  # Assign lower bound
  # assign(paste0(param_name, "_ub"), selected_params$value_ub[i])  # Assign upper bound
}





##############
# Run models #
##############

# # Maternal syphilis model
# maternal_syph_df <- syph_maternal_func(first_ANC_test = 1, late_ANC_test = 1,
#                                        test_type_1 = "lab", test_type_2 = "lab")
# 
# adverse_df <- syph_adverse_outcomes_func(maternal_syph_df, test_type_1 = "lab", test_type_2 = "lab") 
# syph_infant_deaths <- syph_infant_deaths_func(adverse_df)
# 
# syph_infant_outcomes_1yr <- syph_infant_outcomes_1yr_func(adverse_df)
# 
# syph_cost_df <- syph_cost_func(test_type_1="lab", test_type_2 = "lab", maternal_syph_df = maternal_syph_df, adverse_df = adverse_df)
# 
# 
# 
# ## Run HIV model
# # maternal data frame
# maternal_hiv_df <- hiv_maternal_func(first_ANC_test = 1, late_ANC_test = 1,
#                                      test_type_1 = "rdt", test_type_2 = "rdt",
#                                      syph_test = TRUE, 
#                                      syph_infant_deaths = syph_infant_deaths)
# # person weeks data frame
# person_weeks <- pw_func(maternal_hiv_df)
# # Maternal costs - make sure these match the same input as the maternal DF ************
# maternal_costs <- mat_cost_func(person_weeks, maternal_hiv_df, first_ANC_test = 1, late_ANC_test = 1,
#                                 test_type_1 = "rdt", test_type_2 = "rdt")
# # Infant infections model
# infant_infections <- inf_infect_func(person_weeks)
# # lifetable
# lifetable <- lifetable_output()
# 
# # Get total infant outcomes
# infant_outcomes <- infant_outcomes_func(infant_infections, maternal_hiv_df)
# total_hiv_infected_infants <- infant_outcomes$total_infected_infants
# total_infant_treatment_costs <- infant_outcomes$total_infant_treatment_costs
# infant_deaths <- infant_outcomes$infant_deaths
# 
# # Infant table for 20 years
# child_outcomes <- infant_hiv_syp_model(maternal_hiv_df, infant_outcomes$total_infected_infants,
#                                        infant_outcomes$total_infant_treatment_costs,
#                                        syph_infant_outcomes_1yr,
#                                        infant_outcomes$infant_deaths,
#                                        syph_cost_df,
#                                        lifetable)
# 
# 
# total_infant_dalys <- total_infant_dalys_func(child_outcomes)
# 
# 
# # Get cost-effectiveness output
# total_costs <- tot_mat_costs_func(maternal_costs) + total_infant_costs_func(child_outcomes) +
#   sum(unlist(syph_cost_df), na.rm=TRUE)
# 
# total_hiv_infected_infants <- infant_outcomes$total_infected_infants
# 
# syphilis_outcomes <- rowSums(adverse_df)
# 
# 








# Creating function for HIV-syphilis dual model

hiv_syph_dual_func <- function(first_ANC_test=1, late_ANC_test,      # 1=test, 0=no test
                               syph_test=TRUE,                    # TRUE = syphilis test
                               first_ANC_test_type,               # test types: "dual", "rdt"
                               late_ANC_test_type) {                      
  
  ### Purpose: this function takes in parameters for a population, 
    ### runs HIV and syphilis functions to get population outputs
    ### and calculates costs, outcomes, and ICERs
  ### Input:
  ### Output:
  
  # Setting up 
  if (first_ANC_test_type == "dual") {
    first_syphilis_test <- "dual"
    first_hiv_test <- "dual"
  } else if (first_ANC_test_type == "rdt") {
    first_syphilis_test <- "lab"
    first_hiv_test <- "rdt"
  } else {
    print("First ANC test type not selected")
  }
  
  if (late_ANC_test_type == "dual") {
    late_syphilis_test <- "dual"
    late_hiv_test <- "dual"
  } else if (late_ANC_test_type == "rdt") {
    late_syphilis_test <- "lab"
    late_hiv_test <- "rdt"
  } else {
    print("Late ANC test type not selected")
  }
  
  
  ## Run syphilis model - choices are "lab" or "rdt"
  maternal_syph_df <- syph_maternal_func(first_ANC_test = first_ANC_test, 
                                         late_ANC_test = late_ANC_test,
                                         test_type_1 = first_syphilis_test, 
                                         test_type_2 = late_syphilis_test)
  
  adverse_df <- syph_adverse_outcomes_func(maternal_syph_df, 
                                           test_type_1 = first_syphilis_test,
                                           test_type_2 = late_syphilis_test)
  
  syph_infant_deaths <- syph_infant_deaths_func(adverse_df)
  
  syph_infant_outcomes_1yr <- syph_infant_outcomes_1yr_func(adverse_df)
  
  syph_cost_df <- syph_cost_func(test_type_1 = first_syphilis_test,
                                 test_type_2 = late_syphilis_test,
                                 maternal_syph_df = maternal_syph_df, 
                                 adverse_df = adverse_df)
  
  ## Run HIV model
  maternal_hiv_df <- hiv_maternal_func(first_ANC_test = first_ANC_test, 
                                       late_ANC_test = late_ANC_test, 
                                       test_type_1 = first_hiv_test, test_type_2 = late_hiv_test, 
                                       syph_test = syph_test,
                                       syph_infant_deaths = syph_infant_deaths)
  
  # person weeks data frame
  person_weeks <- pw_func(maternal_hiv_df)
  # Maternal costs - make sure these match the same input as the maternal DF ************
  maternal_costs <- mat_cost_func(person_weeks, maternal_hiv_df,
                                  first_ANC_test = first_ANC_test,
                                  late_ANC_test = late_ANC_test, 
                                  test_type_1 = first_hiv_test, 
                                  test_type_2 = late_hiv_test)
  
  # Infant infections model
  infant_infections <- inf_infect_func(person_weeks_df = person_weeks)
  # lifetable
  lifetable <- lifetable_output()
  
  # Get total infant outcomes
  infant_outcomes <- infant_outcomes_func(infant_infections, maternal_hiv_df)
  total_infected_infants <- infant_outcomes$total_infected_infants
  total_infant_treatment_costs <- infant_outcomes$total_infant_treatment_costs
  infant_deaths <- infant_outcomes$infant_deaths
  
  ######### OUTPUT
  
  # Infant table for 20 years
  child_outcomes <- infant_hiv_syp_model(maternal_hiv_df, infant_outcomes$total_infected_infants,
                                         infant_outcomes$total_infant_treatment_costs,
                                         syph_infant_outcomes_1yr,
                                         infant_outcomes$infant_deaths,
                                         syph_cost_df,
                                         lifetable)
  
  total_infant_dalys <- total_infant_dalys_func(child_outcomes)
  
  # Get cost-effectiveness output 
  total_costs <- tot_mat_costs_func(maternal_costs) + total_infant_costs_func(child_outcomes) + 
    sum(unlist(syph_cost_df), na.rm=TRUE)
  
  total_hiv_infected_infants <- infant_outcomes$total_infected_infants
  
  syphilis_outcomes <- rowSums(adverse_df)

  ### Getting data frame of output
  output <- NULL
  output$total_costs <- total_costs
  output$hiv_infected_infants <- total_hiv_infected_infants
  output$syph_preterm <- syphilis_outcomes["premature"]
  output$syph_asymptomatic <- syphilis_outcomes["asympt"]
  output$syph_clinical_congenital_syph <- syphilis_outcomes["congenital"]
  output$syph_stillbirths <- syphilis_outcomes["stillbirths"]
  output$syph_neonatal_death <- syphilis_outcomes["neonatal"]
  output$infant_deaths <- infant_deaths
  output$dalys <- total_infant_dalys_func(child_outcomes)

  return(unlist(output))

}



two_dual_test <- hiv_syph_dual_func(1,1,syph_test = TRUE, "dual", "dual")
two_rdt_test <- hiv_syph_dual_func(1,1,syph_test = TRUE, "rdt", "rdt")
single_rdt <- hiv_syph_dual_func(1,0,syph_test = TRUE, "rdt", "rdt")
single_dual <- hiv_syph_dual_func(1,0,syph_test = TRUE, "dual", "dual")

comparison <- rbind(single_rdt, single_dual, two_dual_test, two_rdt_test) %>%
  as.data.frame()



## Getting ICERs

comparison <- rbind(comparison[rownames(comparison) == "single_rdt", ],
                    comparison[rownames(comparison) != "single_rdt", ] %>%
                      arrange(total_costs)) %>%
  mutate(incremental_cost = total_costs - lag(total_costs),
         incremental_dalys = lag(dalys) - dalys,
         icer = incremental_cost / incremental_dalys,
         syph_total = syph_asymptomatic.asympt +
                          syph_clinical_congenital_syph.congenital +
                          syph_stillbirths.stillbirths +
                          syph_neonatal_death.neonatal,
         icer2 = ifelse(icer > lead(icer), "Dom", round(icer, 2)),
         icer2 = ifelse(is.na(icer2), 
                        round((total_costs - lag(total_costs, 2)) / (lag(dalys, 2) - dalys) ,2), icer2))

comparison_out <- comparison %>%
  select(total_costs, hiv_infected_infants, syph_total, dalys, incremental_cost, incremental_dalys, icer2)




cost_effective_func <- function()





#### Function here to read in parameters, split population up by HIV and syphilis status, and push output to df









### NEEDS
  # Get CE output
  # Wrap all HIV functions into a single function
  # Wrap all syphilis functions into a single function
  # Wrap HIV and syphilis functions into a single function
  # set up sensitivity analysis
  # set up stochastic process - wrap this into a function
















