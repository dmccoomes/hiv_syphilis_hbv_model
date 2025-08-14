
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


#### This script includes functions for downloading population parameters, 
  #### running all models for the HIV-syphilis dual testing,
  #### and collecting the output into a cost-effectiveness table



##############
# Parameters #
##############

# Creating function for population parameters
pop_param_func <- function(input_workbook, model, country, target_pop) {
  
  ### Function description
  
  # Read in population characteristics
  pop_parameters <- readxl::read_xlsx(input_workbook, 
                                      sheet = "pop_characteristics") %>%
    mutate(across(value:value_ub,  ~as.double(.)))
  
  # Read in test characteristics
  test_parameters <- readxl::read_xlsx(input_workbook, 
                                       sheet = "test_characteristics")
  
  # Read in adult model characteristics
  adult_parameters <- readxl::read_xlsx(input_workbook, 
                                        sheet = "adult_model_characteristics")
  
  # only for the selected combinations
  selected_params <- pop_parameters %>%
    bind_rows(adult_parameters) %>%
    filter(model == model, #will make a loop later
           country == country, #will make a loop later
           target_pop == target_pop) %>% #will make a loop later
    bind_rows(test_parameters) 
    
  # Set up all parameters for the selected population
  for (i in 1:nrow(selected_params)) {
    param_name <- as.character(selected_params$var_name[i])
    assign(param_name, selected_params$value[i])
    # assign(paste0(param_name, "_lb"), selected_params$value_lb[i])  # Assign lower bound
    # assign(paste0(param_name, "_ub"), selected_params$value_ub[i])  # Assign upper bound
  }
  
  env_list <- mget(ls())
  return(env_list)

}


##############
# Run models #
##############

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
    sum(syph_cost_df$test_cost, syph_cost_df$treatment_cost_fp, syph_cost_df$treatment_cost_tp,
        syph_cost_df$cong_syph_costs, na.rm=TRUE)
  
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







#### Function here to read in parameters, split population up by HIV and syphilis status, and push output to df
cost_effective_func <- function(testing_algorithm, 
                                syph_test = TRUE,
                                baseline = NULL,
                                # Define the population
                                param_workbook, model, country, target_pop) {
  
  ### Function description:
  
  ## Set up empty list for output
  output <- list(NULL)
  
  pop_param <- pop_param_func(param_workbook, model, country, target_pop)
  # Get specified population parameters
  list2env(pop_param, envir = parent.frame())
  

  ## Run models for each testing algorithm
  for (x in 1:length(testing_algorithm)) {

      out <- hiv_syph_dual_func(first_ANC_test = testing_algorithm[[x]][[1]][[1]],
                       late_ANC_test = testing_algorithm[[x]][[1]][[2]],
                       first_ANC_test_type = testing_algorithm[[x]][[2]][[1]],
                       late_ANC_test_type = testing_algorithm[[x]][[2]][[2]],
                       syph_test = TRUE)
    
    output[[names(testing_algorithm)[x]]] <- out

  }

  ## Combine output for all algorithms 
  comparison <- do.call(rbind, output) %>%
    as.data.frame()

  ## Getting ICERs
  if (!is.null(baseline)) {
  comparison <- rbind(comparison[rownames(comparison) == baseline, ],
                      comparison[rownames(comparison) != baseline, ] %>%
                        arrange(total_costs))
  } 
  
  comparison <- comparison %>%
    mutate(incremental_cost = total_costs - lag(total_costs),
           incremental_dalys = lag(dalys) - dalys,
           icer = incremental_cost / incremental_dalys,
           syph_total = syph_asymptomatic.asympt +
             syph_clinical_congenital_syph.congenital +
             syph_stillbirths.stillbirths +
             syph_neonatal_death.neonatal)
  
  ## Keeping only output that we need
  comparison_out <- comparison %>%
    select(total_costs, hiv_infected_infants, syph_total, dalys, incremental_cost, incremental_dalys, icer)
  
  return(comparison_out)
  
}



























