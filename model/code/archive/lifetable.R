

#### Title: Life table for HIV retesting model 
#### Purpose: Create HIV retesting infant model based on Excel models
#### Author: David Coomes
#### Date: January 21, 2025


library(tidyverse)
library(here)

source(here("code", "input_parameters.R"))


### Lifetable calculations


if (country == "Kenya") {

  # HIV mortality rate calculations
    # these might be country specific and should go into the parameter section?
  cum_mort_year1 <- 0.352
  cond_mort_year1 <- cum_mort_year1
  cum_mort_year2 <- 0.525
  cond_mort_year2 <- 1 - (cum_mort_year2-1)/(-1*(1-cond_mort_year1))
  
  
  # Creating HIV negative life table
  cond_mort_hiv_neg <- c(0.0365, 0.0145, 0.007, 0.0055, 0.0095)
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
  cond_mort_hiv_pos <- c(cum_mort_year1*(1-inf_ART_coverage_base*avg_ART_adherence), 
                         cond_mort_year2*(1-inf_ART_coverage_base*avg_ART_adherence),
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
  
  # Describe
  surv_nhiv <- 1 - life_table$annual_cond_prob_nonhiv_death[1]

  
} else if (country == "South Africa") {
  
  # HIV mortality rate calculations
  # these might be country specific and should go into the parameter section?
  cum_mort_year1 <- 0.352
  cond_mort_year1 <- cum_mort_year1
  cum_mort_year2 <- 0.525
  cond_mort_year2 <- 1 - (cum_mort_year2-1)/(-1*(1-cond_mort_year1))
  
  
  # Creating HIV negative life table
  cond_mort_hiv_neg <- c(0.0335, 0.0095, 0.004, 0.0065, 0.0075)
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
  cond_mort_hiv_pos <- c(cum_mort_year1*(1-inf_ART_coverage_base*avg_ART_adherence), 
                         cond_mort_year2*(1-inf_ART_coverage_base*avg_ART_adherence),
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
  
  # Describe
  surv_nhiv <- 1 - life_table$annual_cond_prob_nonhiv_death[1]
  
  
  
} else if (country == "Ukraine") {
  
  # HIV mortality rate calculations
  # these might be country specific and should go into the parameter section?
  cum_mort_year1 <- 0.352
  cond_mort_year1 <- cum_mort_year1
  cum_mort_year2 <- 0.525
  cond_mort_year2 <- 1 - (cum_mort_year2-1)/(-1*(1-cond_mort_year1))
  
  
  # Creating HIV negative life table
  cond_mort_hiv_neg <- c(0.0085, 0.001, 0.001, 0.001, 0.003)
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
  cond_mort_hiv_pos <- c(cum_mort_year1*(1-inf_ART_coverage_base*avg_ART_adherence), 
                         cond_mort_year2*(1-inf_ART_coverage_base*avg_ART_adherence),
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
  
  # Describe
  surv_nhiv <- 1 - life_table$annual_cond_prob_nonhiv_death[1]
  
} else {
  
  print("No population selected")
  
}





