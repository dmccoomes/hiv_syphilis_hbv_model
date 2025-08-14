
#### Title: Triplex output data
#### Purpose: Output HIV and syphilis population counts for triplex model
#### Author: David Coomes
#### Date: July 7, 2025


#### Purpose ####
### This script runs HIV-syphilis dual models and outputs dataframes of child and adult populations 
  ### by HIV and syphilis status, and treatment
  

# DMC 7/7/25: Testing function
# first_ANC_test <- 1
# late_ANC_test <- 0
# first_syphilis_test <- "lab"
# late_syphilis_test <- "dual"
# first_hiv_test <- "rdt"
# late_hiv_test <- "dual"
# syph_test <- 1


# Get maternal costs split out by testing costs and treatment costs
mat_cost_split_func <- function(maternal_costs_df) {
  
  ### this function takes the maternal costs dataframe and splits it up between 
  ### testing costs and treatment costs
  ### Input: maternal costs df that comes from the HIV model
  ### Output: list with test costs and treatment costs separated out
  
  test_cost <- sum(maternal_costs_df$hiv_neg, maternal_costs_df$ahiv_pos_nart,
                   maternal_costs_df$chiv_pos_uk, na.rm=TRUE)
  treat_cost <- sum(maternal_costs_df$ahiv_pos_art, maternal_costs_df$chiv_pos_art,
                    maternal_costs_df$chiv_pos_nart, na.rm=TRUE)
  
  out <- list("test_cost" = test_cost, "treat_cost" = treat_cost)
  
  return(out)
  
}


## DMC 8/12/25 - This function is pulled into function below
  # They both use a lot of the same input so we may be able to streamline the two functions
  # this function is similar to one in the infant model, but we have a few extra outputs
infant_hiv_syp_model_tri <- function(maternal_hiv_df,
                                     total_hiv_infected_infants,               # this is output from HIV infant outcomes
                                     total_infant_treatment_costs,
                                     syph_infant_outcomes_1yr,
                                     infant_deaths,
                                     syph_cost_df,
                                     lifetable) {
  
  ### this function uses input from HIV and syphilis models to create child outcomes
  ## input: maternal HIV df, infant HIV df, infant syphilis df
  ## output: 
  
  

  # Setting up empty matrix
  infant_matrix <- matrix(nrow = 21, ncol = 13,
                          dimnames = list(NULL,
                                          c("age", "hiv_infected_infants", "syphilis_infected_infants",
                                            "coinfected_infants", "lbw_infants", "healthy_infants",
                                            "syph_dead", "hiv_dead", "bg_dead", "discounted_ylls",
                                            "discounted_ylds", "discounted_dalys",
                                            "discounted_costs")))

  infant_df <- as.data.frame(infant_matrix)

  infants_1yr <- population_size - sum(maternal_hiv_df$deaths[1:deliveryweek], na.rm=TRUE) - infant_deaths
  pct_hiv_infected_infants <- total_hiv_infected_infants / infants_1yr
  pct_syph_infected_infants <- (syph_infant_outcomes_1yr$cong_syphilis + syph_infant_outcomes_1yr$asympt) / infants_1yr

  hiv_infected_infants <- total_hiv_infected_infants * (1-pct_syph_infected_infants)

  # Getting syphilis infected infants
  syphilis_infected_infants <- (syph_infant_outcomes_1yr$cong_syphilis + syph_infant_outcomes_1yr$asympt) *
    (1-pct_hiv_infected_infants)
  # This is different than excel sheet - I think they are not counting correctly
  coinfected_infants <- total_hiv_infected_infants * pct_syph_infected_infants + syphilis_infected_infants * pct_hiv_infected_infants
  lbw_infants <- syph_infant_outcomes_1yr$premature
  syph_dead <- syph_infant_outcomes_1yr$stillbirths + syph_infant_outcomes_1yr$neonatal
  hiv_dead <- 0
  bg_dead <- 0

  # Getting non-infected infants
  healthy_infants <- population_size - sum(maternal_hiv_df$deaths[1:(deliveryweek-1)]) -
    sum(total_hiv_infected_infants, syphilis_infected_infants, coinfected_infants, lbw_infants, syph_dead, hiv_dead, bg_dead)
  # Setting up first week parameters
  infant_df[1, ] <- c(NA, hiv_infected_infants, syphilis_infected_infants, coinfected_infants, lbw_infants, healthy_infants,
                      syph_dead, hiv_dead, bg_dead, rep(NA, 4))
  # Setting up age to 19
  infant_df[ ,1] <- c(NA, seq(0, 19, 1))
  # Adding in initial treatment costs at age 0
  infant_df[2, "discounted_costs"] <- (total_infant_treatment_costs / (1+disc_rate)^(infant_df[3,1]+1)) + syph_cost_df$cong_syph_costs


  # Filling out infant model
  # Get vectors of probability of death
  prob_death_hiv <- rep(lifetable$annual_cond_prob_hiv_death, lifetable$interval_length)
  prob_death_nhiv <- rep(lifetable$annual_cond_prob_nonhiv_death, lifetable$interval_length)

  for (x in 2:21) {
    infant_df$hiv_infected_infants[x] <- infant_df$hiv_infected_infants[x-1] -
      infant_df$hiv_infected_infants[x-1] * prob_death_hiv[x-1]

    infant_df$syphilis_infected_infants[x] <- infant_df$syphilis_infected_infants[x-1] -
      infant_df$syphilis_infected_infants[x-1] * prob_death_nhiv[x-1]

    infant_df$coinfected_infants[x] <- infant_df$coinfected_infants[x-1] -
      infant_df$coinfected_infants[x-1] * prob_death_hiv[x-1]

    infant_df$lbw_infants[x] <- infant_df$lbw_infants[x-1] -
      infant_df$lbw_infants[x-1] * prob_death_nhiv[x-1]

    # Updating for full model
    infant_df$syph_dead[x] <- infant_df$syph_dead[x-1] + 0

    infant_df$hiv_dead[x] <- infant_df$hiv_dead[x-1] +
      sum(infant_df$hiv_infected_infants[x-1], coinfected_infants, na.rm=TRUE) * prob_death_hiv[x-1]

    infant_df$bg_dead[x] <- infant_df$bg_dead[x-1] +
      sum(infant_df$syphilis_infected_infants[x-1], infant_df$lbw_infants[x-1], infant_df$healthy_infants[x-1], na.rm=TRUE) *
      prob_death_nhiv[x-1]

    infant_df$healthy_infants[x] <- infant_df$healthy_infants[x-1] -
      infant_df$healthy_infants[x-1] * prob_death_nhiv[x-1]

    infant_df$discounted_ylls[x] <- (infant_df$syph_dead[x] + infant_df$hiv_dead[x] + infant_df$bg_dead[x])/(1+disc_rate)^(infant_df$age[x]+1)

    infant_df$discounted_ylds[x] <-
      # hiv infected infants
      (infant_df$hiv_infected_infants[x] * (inf_ART_coverage*avg_ART_adherence*dw_HIVtx) +
         infant_df$hiv_infected_infants[x] * (1-inf_ART_coverage*avg_ART_adherence) *
         (dw_HIVnotx*propLE_HIV_notx+dw_AIDSnotx*(1-propLE_HIV_notx)) +
         # syphilis infected infants
         ifelse(infant_df$age[x] < 3, infant_df$syphilis_infected_infants[x] * dw_syp, 0) +
         # coinfected infants
         infant_df$coinfected_infants[x]*(inf_ART_coverage*avg_ART_adherence*ifelse(infant_df$age[x]<3, `dw_syp+HIVtx`, dw_HIVtx) +
                                            (1-inf_ART_coverage*avg_ART_adherence)*
                                            (ifelse(infant_df$age[x]<3, `dw_syp+HIVnotx` , dw_HIVnotx)*propLE_HIV_notx) +
                                            ifelse(infant_df$age[x]<3, `dw_syp+AIDSnotx`, dw_AIDSnotx) * (1-propLE_HIV_notx)) +
         # low birthweight infants
         ifelse(infant_df$age[x] < 1, infant_df$lbw_infants[x]*dw_lbw, 0) ) /
      (1+disc_rate)^(infant_df$age[x]+1)


    infant_df$discounted_dalys[x] <- infant_df$discounted_ylls[x] + infant_df$discounted_ylds[x]

    infant_df$discounted_costs[x] <- ifelse(infant_df$age[x]==0, total_infant_treatment_costs,
                                            ((infant_df$hiv_infected_infants[x] + infant_df$coinfected_infants[x]) *
                                               inf_ART_coverage*cost_infant_ART_after_2weeks*
                                               52) / (1+disc_rate)^(infant_df$age[x]+1))

  }

  return(infant_df)

}




# Creating function for HIV-syphilis dual model for triplex output

triplex_out_func <- function(first_ANC_test=1, late_ANC_test,      # 1=test, 0=no test
                               syph_test=TRUE,                    # TRUE = syphilis test
                               test_type_1,               # test types: "dual", "rdt"
                               test_type_2) {                      
  
  ### Purpose: this function takes in parameters for a population, 
  ### runs HIV and syphilis functions to get population outputs
  ### and calculates costs, outcomes, and ICERs
  ### Input: first and second ANC test, plus test type, inclusion of syphilis test
  ### Output: child outcomes, maternal outcomes, maternal costs
  
  # Setting up 
  if (test_type_1 == "dual") {
    first_syphilis_test <- "dual"
    first_hiv_test <- "dual"
  } else if (test_type_1 == "rdt") {
    first_syphilis_test <- "lab"
    first_hiv_test <- "rdt"
  } else {
    print("First ANC test type not selected")
  }
  
  if (test_type_2 == "dual") {
    late_syphilis_test <- "dual"
    late_hiv_test <- "dual"
  } else if (test_type_2 == "rdt") {
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
  child_outcomes <- infant_hiv_syp_model_tri(maternal_hiv_df, infant_outcomes$total_infected_infants,
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
  
  child_outcomes$age <- seq(0, 20, 1)
  
  maternal_costs_split <- mat_cost_split_func(maternal_costs)
  
  # Getting maternal df
  mat_df <- maternal_output_func(first_ANC_test = first_ANC_test, late_ANC_test = late_ANC_test, 
                                 test_type_1 = test_type_1, test_type_2 = test_type_2,
                                 syph_test = FALSE,
                                 syph_infant_deaths = syph_infant_deaths)
  
  out <- list("child_outcomes" = child_outcomes, 
              "maternal_costs" = maternal_costs_split, 
              "maternal_outcomes" = mat_df)

  return(out)
  
}










# child_db <- triplex_out_func(late_ANC_test = 0, 
#                              syph_test = TRUE,
#                              first_ANC_test_type = "rdt", 
#                              late_ANC_test_type = "rdt")
# 
# child_hiv_syph <- child_db$child_outcomes %>%
#   select(-c(discounted_ylls, discounted_ylds,
#             discounted_dalys)) 
# 
# write.csv(child_hiv_syph, here("output", "dual_child_model_base.csv"))













