
#### Title: HIV-syphilis-HBV maternal testing model output
#### Purpose: Run testig models for HIV, syphilis, and HBV among pregnant populations and collect output
#### Author: David Coomes, Akash Malhotra
#### Date: June 3, 2025


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table",
  "readxl",
  "withr"
)


#### This script includes functions for downloading population parameters, 
  #### running all models for the HIV-syphilis dual testing,
  #### and collecting the output into a cost-effectiveness table

# Function to assign variable to global environment inside function and then revert after function is complete
local_assign <- function(name, value, env = parent.frame()) {
  old <- get(name, envir = env, inherits = TRUE)  # save old value
  do.call("assign", list(name, value, envir = env))  # assign new value
  withr::defer(assign(name, old, envir = env), envir = parent.frame())  # restore
}


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
  
  # Setting up test type
  if (first_ANC_test_type == "dual") {
    first_syphilis_test_type <- "dual"
    first_hiv_test_type <- "dual"
  } else if (first_ANC_test_type == "rdt") {
    first_syphilis_test_type <- "lab"
    first_hiv_test_type <- "rdt"
  } else {
    print("First ANC test type not selected")
  }
  
  if (late_ANC_test_type == "dual") {
    late_syphilis_test_type <- "dual"
    late_hiv_test_type <- "dual"
  } else if (late_ANC_test_type == "rdt") {
    late_syphilis_test_type <- "lab"
    late_hiv_test_type <- "rdt"
  } else {
    print("Late ANC test type not selected")
  }
  
  ## Run syphilis model - choices are "lab" or "rdt"
  maternal_syph_list <- syph_maternal_func(first_ANC_test = first_ANC_test, 
                                         late_ANC_test = late_ANC_test,
                                         test_type_1 = first_syphilis_test_type, 
                                         test_type_2 = late_syphilis_test_type)
  
  adverse_df <- syph_adverse_outcomes_func(maternal_syph_list, 
                                           test_type_1 = first_syphilis_test_type,
                                           test_type_2 = late_syphilis_test_type)
  
  syph_infant_deaths <- syph_infant_deaths_func(adverse_df)
  
  syph_infant_outcomes_1yr <- syph_infant_outcomes_1yr_func(adverse_df)
  
  syph_cost_df <- syph_cost_func(test_type_1 = first_syphilis_test_type,
                                 test_type_2 = late_syphilis_test_type,
                                 maternal_syph_list = maternal_syph_list, 
                                 adverse_df = adverse_df)
  
  ## Run HIV model
  maternal_hiv_df <- hiv_maternal_func(first_ANC_test = first_ANC_test, 
                                       late_ANC_test = late_ANC_test, 
                                       test_type_1 = first_hiv_test_type, test_type_2 = late_hiv_test_type, 
                                       syph_test = syph_test,
                                       syph_infant_deaths = syph_infant_deaths)
  
  # person weeks data frame
  person_weeks <- pw_func(maternal_hiv_df)
  # Maternal costs - make sure these match the same input as the maternal HIV DF ************
  maternal_costs <- mat_cost_func(person_weeks, maternal_hiv_df,
                                  first_ANC_test = first_ANC_test,
                                  late_ANC_test = late_ANC_test, 
                                  test_type_1 = first_hiv_test_type, 
                                  test_type_2 = late_hiv_test_type)
  
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









# first_ANC_test_type <- "dual"
# late_ANC_test_type <- "dual"
# first_ANC_test <- 1
# late_ANC_test <- 0
# syph_test <- TRUE
# first_syphilis_test_type <- "dual"
# late_syphilis_test_type <- "dual"
# first_syph_test <- 1
# late_syph_test <- 0
# prev_syph <- FALSE

# Creating wrapper function to split population based on known HIV status and previous syphilis infection status
  # those that are HIV positive and known status (on treatment) do not get an HIV test - syphilis lab test only
  # Those that have tested positive for syphilis do not get a treponemal syphilis test (including dual RDT) - HIV RDT only + syph lab
  # Those that are both HIV positive and have tested positive for syphilis do not get a trep test and no HIV test
hiv_syph_dual_func2 <- function(first_ANC_test=1, late_ANC_test,      # 1=test, 0=no test
                                syph_test=TRUE,                    # TRUE = syphilis test
                                first_ANC_test_type,               # test types: "dual", "rdt"
                                late_ANC_test_type) {          # splitting based on HIV and syphilis status
  
  # split population up into bins
  hiv_syph_pos_pct <- (hiv_prev * pct_ART_female) * (previous_syph_prev - syp_prev) * prob_syph_hiv_coinf
  hiv_pos_pct <- hiv_prev * pct_ART_female - hiv_syph_pos_pct                    # Those with HIV who know their status and are on ART
  syph_pos_pct <- (previous_syph_prev - syp_prev) - hiv_syph_pos_pct             # Those with previous syphilis infection
  neg_pct <- 1 - (hiv_pos_pct + syph_pos_pct + hiv_syph_pos_pct)                 # Those with no known infections
  
  full_pop <- population_size
  base_hiv_prev <- hiv_prev
  
  # Running all populations through testing regimen 
  # No previous infections - get whatever tests are part of the intervention
  local_assign("population_size", full_pop * neg_pct, env = globalenv())    # temporarily assigning 
  local_assign("hiv_prev", 
               ((base_hiv_prev - hiv_pos_pct - hiv_syph_pos_pct) / (neg_pct)) * (neg_pct / (neg_pct + syph_pos_pct)),
               env = globalenv())
  neg_pop <- hiv_syph_dual_func(first_ANC_test = first_ANC_test,
                                late_ANC_test = late_ANC_test, 
                                syph_test = syph_test, 
                                first_ANC_test_type = first_ANC_test_type,
                                late_ANC_test_type = late_ANC_test_type)
  
  # Previous HIV infection - no HIV test, just syphilis labs
  local_assign("population_size", full_pop * hiv_pos_pct, env = globalenv())    # temporarily assigning 
  local_assign("hiv_prev", 1, env = globalenv())
  hiv_pop <- hiv_syph_dual_func(first_ANC_test = first_ANC_test,
                                late_ANC_test = late_ANC_test, 
                                syph_test = syph_test, 
                                first_ANC_test_type = "rdt",
                                late_ANC_test_type = "rdt",
                                previous_HIV_infection = TRUE)
  
  # Previous syphilis infection - no dual test, HIV RDT plus syphilis labs
  local_assign("population_size", full_pop * syph_pos_pct, env = globalenv())    # temporarily assigning 
  local_assign("hiv_prev", ((base_hiv_prev - hiv_pos_pct - hiv_syph_pos_pct) / (syph_pos_pct)) * (syph_pos_pct / (neg_pct + syph_pos_pct)),
               env = globalenv())
  syph_pop <- hiv_syph_dual_func(first_ANC_test = first_ANC_test,
                                late_ANC_test = late_ANC_test, 
                                syph_test = syph_test, 
                                first_ANC_test_type = "rdt",
                                late_ANC_test_type = "rdt",
                                prev_syph = TRUE)
  
  # HIV infection with previous syphilis infection - no HIV test and syphilis labs
  local_assign("population_size", full_pop * hiv_syph_pos_pct, env = globalenv())    # temporarily assigning 
  local_assign("hiv_prev", 1, env = globalenv())
  hiv_syph_pop <- hiv_syph_dual_func(first_ANC_test = first_ANC_test,
                                     late_ANC_test = late_ANC_test, 
                                     syph_test = syph_test, 
                                     first_ANC_test_type = "rdt",
                                     late_ANC_test_type = "rdt",
                                     prev_syph = TRUE,
                                     previous_HIV_infection = TRUE)
  
  out <- neg_pop + hiv_pop + syph_pop + hiv_syph_pos_pct
  return(out)
  #return(list("neg_pop" = neg_pop, "hiv_pop" = hiv_pop, "syph_pop" = syph_pop, "hiv_syph_pop" = hiv_syph_pop))
  
}

# Current issues
  # We are not testing the HIV pos population for HIV, but from my model's standpoint this means we are not treating either
    # Ideas? 



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




### Testing new function that splits populations
#basic_test <- hiv_syph_dual_func(1,0, TRUE, "dual", "dual")
#update <- hiv_syph_dual_func2(1, 0, TRUE, "dual", "dual")


# one_dual <- hiv_syph_dual_func(1,0, TRUE, "dual", "dual")
# two_dual <- hiv_syph_dual_func(1,1, TRUE, "dual", "dual")
# 
# one_dual2 <- hiv_syph_dual_func2(1, 0, TRUE, "dual", "dual")
# two_dual2 <- hiv_syph_dual_func2(1, 1, TRUE, "dual", "dual")


## Notes:
  # This doesn't work (splitting population up and running through the HIV model because
    # when I decide not to test someone )












