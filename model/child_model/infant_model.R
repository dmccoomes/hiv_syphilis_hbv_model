
#### Title: HIV retesting infant model 
#### Purpose: Create HIV retesting infant model based on Excel models
#### Author: David Coomes
#### Date: January 21, 2025


pacman::p_load(
  "EpiModel", 
  "tidyverse", 
  "here",
  "data.table"
)


infant_hiv_syp_model <- function(maternal_hiv_df, 
                                 total_hiv_infected_infants,               # this is output from HIV infant outcomes
                                 total_infant_treatment_costs,
                                 syph_infant_outcomes_1yr,
                                 infant_deaths,
                                 syph_cost_df, 
                                 lifetable) {
  
  # Setting up empty matrix
  infant_matrix <- matrix(nrow = 21, ncol = 11, 
                            dimnames = list(NULL, 
                                            c("age", "hiv_infected_infants", "syphilis_infected_infants",
                                              "coinfected_infants", "lbw_infants", "healthy_infants", 
                                              "dead", "discounted_ylls", 
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
  dead <- syph_infant_outcomes_1yr$stillbirths + syph_infant_outcomes_1yr$neonatal

  # Getting non-infected infants
  healthy_infants <- population_size - sum(maternal_hiv_df$deaths[1:(deliveryweek-1)]) - 
    sum(total_hiv_infected_infants, syphilis_infected_infants, coinfected_infants, lbw_infants, dead)
  # Setting up first week parameters
  infant_df[1, ] <- c(NA, hiv_infected_infants, syphilis_infected_infants, coinfected_infants, lbw_infants, healthy_infants, dead, rep(NA, 4))
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
    
    infant_df$dead[x] <- infant_df$dead[x-1] +
      sum(infant_df$hiv_infected_infants[x-1], coinfected_infants, na.rm=TRUE) * prob_death_hiv[x-1] +
      sum(infant_df$syphilis_infected_infants[x-1], infant_df$lbw_infants[x-1], infant_df$healthy_infants[x-1], na.rm=TRUE) * prob_death_nhiv[x-1]
    
    infant_df$healthy_infants[x] <- infant_df$healthy_infants[x-1] -
      infant_df$healthy_infants[x-1] * prob_death_nhiv[x-1]
    
    infant_df$discounted_ylls[x] <- infant_df$dead[x]/(1+disc_rate)^(infant_df$age[x]+1)
    
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






#################
# Cost and DALY #
#################

# Function for getting total infant costs
total_infant_costs_func <- function(infant_df) {
  
  ### Function description:
  
  total_infant_costs <- sum(infant_df$discounted_costs, na.rm=TRUE)
  return(total_infant_costs)
}


# Function for getting total infant DALYs
total_infant_dalys_func <- function(infant_df) {
  
  ### Function description:
  
  total_infant_dalys <- sum(infant_df$discounted_dalys, na.rm=TRUE)
  return(total_infant_dalys)
}


