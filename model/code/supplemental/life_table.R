

#### Title: Life table for HIV retesting model 
#### Purpose: Create HIV retesting infant model based on Excel models
#### Author: David Coomes
#### Date: January 21, 2025


## Notes: 
  # I don't think the cumulative and conditional mortality rates need to be in the life table
  # We want to use this script to output the lifetable, and then we will use that to calculate 
    # probability of survival to 1 year for HIV+ infants not on ART (surv_hiv)
  # Then we use the cumulative and conditional mortality parameters to calculate the 
    # probability of survival to 1 year for HIV- infants or HIV+ infants not on ART (surv_nhiv)


library(tidyverse)
library(here)


### Input needs

  # cum_mort_year1 : cumulative mortality for first year (0-1) 
  # cum_mort_year2 : cumulative mortality for second year (0-1)
  # inf_ART_coverage : Infant ART coverage (0-1)
  # avg_ART_adherence : average ART adherence (0-1)
  # conditional mortality by age


### Output
  
  # Lifetable dataframe used for infant_model input


### Function

lifetable_output <- function() {

  ### INPUTS: cumulative mortality for HIV+ children for the first 2 years
  ### conditional mortality for HIV- children by age group 
  ### (<1 year, 1-4 years, 5-9 years, 10-14 years, 15-19 years). 
    ### This is a conc string that includes 5 numbers - one for each age group
  ### Also requires infant ART coverage and average ART adherence
  
  ### OUTPUTS: lifetable that includes annual conditional probability of death for 
  ### HIV- and HIV+ children by age group from 0-19 years
  
  required_vars <- c("cum_mort_year1", "cum_mort_year2", "inf_ART_coverage", "avg_ART_adherence", 
                     "cond_mort_hiv_neg_1", "cond_mort_hiv_neg_1_4", "cond_mort_hiv_neg_5_9", "cond_mort_hiv_neg_10_14", "cond_mort_hiv_neg_15_19")
  missing_vars <- required_vars[!sapply(required_vars, exists, envir = parent.frame())]
  
  if (length(missing_vars) > 0) {
    message("Missing variables: ", paste(missing_vars, collapse = ", "))
    return(invisible(NULL))  # Exit early
  }
  
  # Starting with cumulative mortality parameters for HIV positive infants by year
  cond_mort_year1 <- cum_mort_year1
  cond_mort_year2 <- 1 - (cum_mort_year2-1)/(-1*(1-cond_mort_year1))
  
  #put conditional mortality into a concatenated string
  cond_mort_hiv_neg <- c(cond_mort_hiv_neg_1, cond_mort_hiv_neg_1_4, cond_mort_hiv_neg_5_9, cond_mort_hiv_neg_10_14, cond_mort_hiv_neg_15_19)
  
  # Creating HIV negative life table
  cond_mort_hiv_neg <- cond_mort_hiv_neg
  cond_surv_hiv_neg <- 1 - cond_mort_hiv_neg
  cum_surv_hiv_neg <- c(cond_surv_hiv_neg[1])
  
  for (i in 2:5) {
    cum_surv_hiv_neg [i] = cond_surv_hiv_neg[i] * cum_surv_hiv_neg[i-1]
  }
  
  hiv_neg_lt <- data.frame(cond_mort = cond_mort_hiv_neg, cond_surv = cond_surv_hiv_neg, cum_surv = cum_surv_hiv_neg)
  rownames(hiv_neg_lt) <- c("<1 year", "1-4 years", "5-9 years", "10-14 years", "15-19 years")
  hiv_neg_lt <- hiv_neg_lt %>%
    mutate(cum_mort = 1 - cum_surv)
  
  
  # Creating conditional probability of mortality table
  cond_mort_hiv_pos <- c(cum_mort_year1*(1-inf_ART_coverage*avg_ART_adherence), 
                         cond_mort_year2*(1-inf_ART_coverage*avg_ART_adherence),
                         0, 0, 0)
  cond_prob_nonhiv_mort <- cond_mort_hiv_neg
  interval_length <- c(1, 4, 5, 5, 5)
  year <- c(0, 1, 5, 10, 15)
  
  cond_prob_mort_lt <- data.frame(cond_prob_hiv_mort = cond_mort_hiv_pos,
                                  cond_prob_nhiv_mort = cond_prob_nonhiv_mort,
                                  interval_length = interval_length, 
                                  year = year)
  cond_prob_mort_lt <- cond_prob_mort_lt %>%
    mutate(cond_prob_tot_mort = 1-(1-cond_prob_hiv_mort)*(1-cond_prob_nhiv_mort))
  rownames(cond_prob_mort_lt) <- c("<1 year", "1-4 years", "5-9 years", "10-14 years", "15-19 years")
  
  
  # Creating HIV positive life table 
  cond_mort_hiv_pos <- cond_prob_mort_lt$cond_prob_tot_mort
  cond_surv_hiv_pos <- 1 - cond_mort_hiv_pos
  
  cum_surv_hiv_pos <- c(cond_surv_hiv_pos[1])
  for (i in 2:5) {
    cum_surv_hiv_pos [i] = cond_surv_hiv_pos[i] * cum_surv_hiv_pos[i-1]
  }
  
  hiv_pos_lt <- data.frame(cond_mort = cond_mort_hiv_pos, cond_surv = cond_surv_hiv_pos, cum_surv = cum_surv_hiv_pos)
  rownames(hiv_pos_lt) <- c("<1 year", "1-4 years", "5-9 years", "10-14 years", "15-19 years")
  hiv_pos_lt <- hiv_pos_lt %>%
    mutate(cum_mort = 1 - cum_surv)
  
  
  # Getting final table
  life_table <- data.frame(cond_prob_nonhiv_mort = cond_prob_mort_lt$cond_prob_nhiv_mort,
                           cond_prob_tot_mort = cond_prob_mort_lt$cond_prob_tot_mort, 
                           interval_length = cond_prob_mort_lt$interval_length)
  
  life_table <- life_table %>%
    mutate(annual_cond_prob_hiv_death = 1 - exp((log(1-cond_prob_tot_mort)) / interval_length),
           annual_cond_prob_nonhiv_death = 1 - exp((log(1-cond_prob_nonhiv_mort)) / interval_length))
  life_table$year <- c(0, 1, 5, 10, 15)
  
  
  return(life_table)
  
  
}











