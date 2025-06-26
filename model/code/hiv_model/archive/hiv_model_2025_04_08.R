
#### Title: HIV retesting model
#### Purpose: Create HIV retesting model based on Excel models
#### Author: David Coomes, Akash Malhotra
#### Date: April 8, 2025


# Load libraries
pacman::p_load(
  "tidyverse", 
  "here",
  "data.table"
)


hiv_model <- function() {
  
  ### This function takes in parameters and creates output of HIV and treatment states
  ### using a Markov model
  ### Inputs include paramters
  ### Outputs are a dataframe of HIV and treatment states, and mortality
  
  # Precompute common expressions
  survival_rate <- ifelse(start_week %in% c(delivery:(delivery+early_pp_visit-1)), 
                          (1-mort_maternal_perinatal), (1-mort_adult_female))
  
  
  
}




hiv_neg <- NULL
hiv_neg[1] <- population_size * neg
incidence <- NULL
week <- NULL


hiv_neg_func <- function(data, PrEP = 0) {
  
  # First, set up some constant terms across stages
  ## 
  
  prep_use <- PrEP
  
  hiv_neg <- data
  
  # Update HIV negative population
  for (x in 2:(delivery+end_of_model)) {
    
    week <- x

    # Set up mortality - use maternal mortality from delivery to 14 weeks pp, general mortality for all other times
    mortality <- ifelse(week >= delivery & week < (delivery + early_pp_visit), mort_maternal_perinatal, mort_adult_females)
    
    # Set up incidence depending on time period
    if (week > 0 & week <= second_ANC) {
      incidence == inc_preg_early 
    } else if (week > second_ANC & week <= delivery) {
      incidence == inc_preg_late
    } else if (week > delivery & week <= (delivery+pp_14weeks_visit)) {
      incidence == inc_pp_early
    } else if (week > (delivery+pp_14weeks_visit) & week <= (delivery+pp_mid_visit)) {
      incidence == inc_pp_mid 
    } else if (week > (delivery+pp_mid_visit) & week <= end_of_model) {
      incidence == inc_pp_late
    } else {
      print("out of bounds")
    }
    
    # incidence <- case_when(week > 0 & week <= second_ANC ~ inc_preg_early,
    #                        week > second_ANC & week <= delivery ~ inc_preg_late,
    #                        week > delivery & week <= (delivery+pp_14weeks_visit) ~ inc_pp_early,
    #                        week > (delivery+pp_14weeks_visit) & week <= (delivery+pp_mid_visit) ~ inc_pp_mid,
    #                        week > (delivery+pp_mid_visit) & week <= end_of_model ~ inc_pp_late)

    # Set up force of infection
    infection_force <- (1-incidence * (1-prep_use*red_prep))
    
    # Fill out vector of values
    #hiv_neg[x] <- hiv_neg[x-1] * (1-mortality) * infection_force
    week[x] <- week 
    incidence[x] <- incidence

  }
  
  df <- bind_cols(week, incidence)
  # 
  # #return(hiv_neg)
  # return(df)
  return(df) 
}

hiv_neg <- hiv_neg_func(data = hiv_neg)





hiv_neg_func <- function(PrEP = 0) {
  
  prep_use <- PrEP
  
  # Update HIV negative population
  for (x in 2:(delivery+end_of_model)) {
    
    print(x)
    
  }
    
    week <- x
    
    # Set up mortality - use maternal mortality from delivery to 14 weeks pp, general mortality for all other times
    mortality <- ifelse(week >= delivery & week < (delivery + early_pp_visit), mort_maternal_perinatal, mort_adult_females)
    
    # Set up incidence depending on time period
    if (week > 0 & week <= second_ANC) {
      incidence == inc_preg_early 
    } else if (week > second_ANC & week <= delivery) {
      incidence == inc_preg_late
    } else if (week > delivery & week <= (delivery+pp_14weeks_visit)) {
      incidence == inc_pp_early
    } else if (week > (delivery+pp_14weeks_visit) & week <= (delivery+pp_mid_visit)) {
      incidence == inc_pp_mid 
    } else if (week > (delivery+pp_mid_visit) & week <= end_of_model) {
      incidence == inc_pp_late
    } else {
      print("out of bounds")
    }
    
    # Set up force of infection
    infection_force <- (1-incidence * (1-prep_use*red_prep))
    
    # Fill out vector of values
    #hiv_neg[x] <- hiv_neg[x-1] * (1-mortality) * infection_force
    week[x] <- x 
    incidence[x] <- incidence
    
  }
  
  df <- bind_cols(week, incidence)
  # 
  # #return(hiv_neg)
  # return(df)
  return(df) 
}

hiv_neg <- hiv_neg_func()











