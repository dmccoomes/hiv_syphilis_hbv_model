






# Define a function to update maternal_df efficiently
update_stage <- function(start_week, end_week, inc_rate, test_prob, att_prob, test_sens, stage_number) {
  for (x in start_week:end_week) {
    maternal_df$weeks[x] <<- x
    
    # Print iteration details
    print(paste("Week:", x))
    
    # Precompute common expressions
    survival_rate <- (1 - mort_adult_female)
    
    if (stage_number == 0) {
      
      infection_factor <- (1 - (inc_rate * (1 - (PrEP * red_prep))))
      
      # HIV negative progression
      maternal_df$hiv_neg[x] <<- maternal_df$hiv_neg[x-1] * survival_rate * infection_factor
      
      # Acute HIV infection progression
      maternal_df$ahiv_pos_nart[x] <<- maternal_df$hiv_neg[x-1] * survival_rate * inc_rate * (1 - (PrEP * red_prep)) +
        maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 - 1/duration_of_recent_infection) +
        maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 - (1/duration_of_recent_infection)) * (1 - ART_drop_factor)
      
      # ART progression
      maternal_df$ahiv_pos_art[x] <<- maternal_df$ahiv_pos_art[x-1] * survival_rate *(1-(1/duration_of_recent_infection)) * ART_drop_factor 
      
      # Chronic HIV progression
      maternal_df$chiv_pos_uk[x] <<- maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 / duration_of_recent_infection)  +
        maternal_df$chiv_pos_uk[x-1] * survival_rate 
      
      maternal_df$chiv_pos_nart[x] <<- maternal_df$chiv_pos_nart[x-1] * survival_rate +
        maternal_df$chiv_pos_art[x-1] * survival_rate * (1 - ART_drop_factor)
      
      maternal_df$chiv_pos_art[x] <<- maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 / duration_of_recent_infection) +
        maternal_df$chiv_pos_art[x-1] * survival_rate * ART_drop_factor
      
    } else if (stage_number %in% c(1,2,3,4,5,6,7,8) )  {
      
      y <- ifelse(x == 19 | x == 36 | x == 39 | x == 40 | x == 45 | x == 53 | x == 65 | x == 78, 1, 0)
      
      infection_factor <- (1 - (inc_rate * (1 - (PrEP * red_prep))))
      
      # HIV negative progression
      maternal_df$hiv_neg[x] <<- maternal_df$hiv_neg[x-1] * survival_rate * infection_factor
      
      # Acute HIV infection progression
      maternal_df$ahiv_pos_nart[x] <<- maternal_df$hiv_neg[x-1] * survival_rate * (1-infection_factor) +
        (maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 - 1/duration_of_recent_infection) *
           (1 - test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results * p_ART) *y) + 
        (maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 - 1/duration_of_recent_infection) * (1-y)) +
        maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 - (1/duration_of_recent_infection)) * (1 - ART_drop_factor)
      
      # ART progression
      maternal_df$ahiv_pos_art[x] <<- maternal_df$ahiv_pos_art[x-1] * survival_rate *(1-(1/duration_of_recent_infection)) * ART_drop_factor +
        ((maternal_df$ahiv_pos_nart[x-1] * survival_rate * ( 1-1/duration_of_recent_infection) * 
            test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results * p_ART) * y) 
      
      
      # Chronic HIV progression
      maternal_df$chiv_pos_uk[x] <<- (maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1/duration_of_recent_infection) *
                                        (1 - test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results * p_ART) *y) + 
        (maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1/duration_of_recent_infection) * (1-y)) +
        (maternal_df$chiv_pos_uk[x-1] * survival_rate *
           (1 - (test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results))*y) +
        (maternal_df$chiv_pos_uk[x-1] * survival_rate * (1-y)) 
      
      maternal_df$chiv_pos_nart[x] <<- (maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 / duration_of_recent_infection) *
                                          test_prob * att_prob * test_accept * (1 - stockout) * p_results * test_sens * (1 - p_ART) * y) +
        (maternal_df$chiv_pos_uk[x-1] * survival_rate * 
           test_prob * att_prob * test_accept * (1 - stockout) * p_results * test_sens * (1 - p_ART) * y ) +
        maternal_df$chiv_pos_nart[x-1] * survival_rate +
        maternal_df$chiv_pos_art[x-1] * survival_rate * (1 - ART_drop_factor)
      
      maternal_df$chiv_pos_art[x] <<- ((maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 / duration_of_recent_infection) *
                                          test_prob * att_prob * test_sens * test_accept * (1 - stockout) * p_results * p_ART) * y) +
        maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 / duration_of_recent_infection) +
        ((maternal_df$chiv_pos_uk[x-1] * survival_rate * 
            test_prob * att_prob * test_sens * test_accept *  (1 - stockout) * p_results * p_ART) *y ) +
        maternal_df$chiv_pos_art[x-1] * survival_rate * ART_drop_factor
    }
    
    # Mortality calculations
    maternal_df$deaths[x] <<- mort_adult_female * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <<- stage_number
  }
}






