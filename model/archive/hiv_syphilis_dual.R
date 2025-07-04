
#### Title: HIV syphilis retesting model
#### Purpose: Create HIV retesting model based on Excel models
#### Author: Akash Malhotra, David Coomes
#### Date: January 27, 2025


#### Title: HIV retesting model
#### Purpose: Create HIV retesting model based on Excel models
#### Author: David Coomes
#### Date: January 9, 2025


# Load libraries
pacman::p_load(
  "EpiModel", 
  "tidyverse", 
  "here",
  "data.table",
  "dplyr"
)

parameters <- read.csv("inp_par.csv")

# only for the selected combinations
selected_params <- parameters %>%
  filter(model == "HIV_Syphilis_Dual_Maternal", #will make a loop later
         country == "South Africa", #will make a loop later
         target_pop == "Maternal") #will make a loop later


for (i in 1:nrow(selected_params)) {
  param_name <- as.character(selected_params$var_name[i])
  
  assign(param_name, selected_params$value[i])
  assign(paste0(param_name, "_lb"), selected_params$value_lb[i])  # Assign lower bound
  assign(paste0(param_name, "_ub"), selected_params$value_ub[i])  # Assign upper bound
}

# Starting states
neg <- (1-hiv_prev) - (inc_preg_early*duration_infection)
recent <- duration_infection*inc_preg_early
established_uk_nart <- hiv_prev*(1-status_known)
established_k_nart <- hiv_prev*status_known*(1-p_ART)
established_art <- hiv_prev*status_known*p_ART

# set the test sensitivity by selecting e3 or e4
test_sens_e <- ifelse(test_gen==3, sens_e3, sens_e4)
test_sens_c <- ifelse(test_gen==3, sens_c3, sens_c4)
test_spec <- ifelse(test_gen==3, spec3, spec4)

# Initialize maternal matrix with column names
maternal_df <- data.frame(
  weeks = integer(92),
  hiv_neg = numeric(92),
  ahiv_pos_nart = numeric(92),
  ahiv_pos_art = numeric(92),
  chiv_pos_uk = numeric(92),
  chiv_pos_nart = numeric(92),
  chiv_pos_art = numeric(92),
  deaths = numeric(92),
  stage = character(92)
)

# Set up first week parameters
maternal_df[1, ] <- c(1, population * neg, population * recent, 0, 
                      population * established_uk_nart, population * established_k_nart, 
                      population * established_art, 0, "Stage 0: Onset of pregnancy to first ANC")

maternal_df$hiv_neg <- as.numeric(maternal_df$hiv_neg)
maternal_df$ahiv_pos_nart <- as.numeric(maternal_df$ahiv_pos_nart)
maternal_df$ahiv_pos_art <- as.numeric(maternal_df$ahiv_pos_art)
maternal_df$chiv_pos_uk <- as.numeric(maternal_df$chiv_pos_uk)
maternal_df$chiv_pos_nart <- as.numeric(maternal_df$chiv_pos_nart)
maternal_df$chiv_pos_art <- as.numeric(maternal_df$chiv_pos_art)
maternal_df$deaths <- as.numeric(maternal_df$deaths)
# Define a function to update maternal_df efficiently
update_stage <- function(start_week, end_week, inc_rate, test_prob, att_prob, test_sens, stage_name) {
  for (x in start_week:end_week) {
    maternal_df$weeks[x] <<- x
    
    # Print iteration details
    print(paste("Week:", x))
    
    # Precompute common expressions
    survival_rate <- (1 - mort_adult_female)
    infection_factor <- (1 - (inc_rate * (1 - (PrEP * red_prep))))
    
    # HIV negative progression
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * survival_rate * infection_factor
    
    # Acute HIV infection progression
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * survival_rate * inc_rate * (1 - (PrEP * red_prep)) +
      maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 - 1/duration_infection) * 
      (1 - test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results * p_ART) +
      maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 - (1/duration_infection)) * (1 - ART_drop_factor)
    
    # ART progression
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 - (1/duration_infection)) * ART_drop_factor +
      maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 - 1/duration_infection) * 
      test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results * p_ART
    
    # Chronic HIV progression
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 / duration_infection) * 
      (1 - (test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results)) +
      maternal_df$chiv_pos_uk[x-1] * survival_rate * 
      (1 - (test_prob * att_prob * (1 - stockout) * test_accept * test_sens * p_results))
    
    maternal_df$chiv_pos_nart[x] <- maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 / duration_infection) * 
      test_prob * att_prob * test_accept * (1 - stockout) * p_results * test_sens * (1 - p_ART) +
      maternal_df$chiv_pos_uk[x-1] * survival_rate * test_prob * att_prob * test_accept * 
      (1 - stockout) * p_results * test_sens * (1 - p_ART) +
      maternal_df$chiv_pos_nart[x-1] * survival_rate +
      maternal_df$chiv_pos_art[x-1] * survival_rate * (1 - ART_drop_factor)
    
    maternal_df$chiv_pos_art[x] <- maternal_df$ahiv_pos_nart[x-1] * survival_rate * (1 / duration_infection) * 
      test_prob * att_prob * test_sens * test_accept * (1 - stockout) * p_results * p_ART +
      maternal_df$ahiv_pos_art[x-1] * survival_rate * (1 / duration_infection) +
      maternal_df$chiv_pos_uk[x-1] * survival_rate * test_prob * att_prob * test_sens * test_accept * 
      (1 - stockout) * p_results * p_ART +
      maternal_df$chiv_pos_art[x-1] * survival_rate * ART_drop_factor
    
    # Mortality calculations
    maternal_df$deaths[x] <- mort_adult_female * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- stage_name
  }
}



# Update different pregnancy stages
update_stage(2, first_ANC - 1, inc_preg_early, first_ANC_test, att_firstANC, test_sens_e, "Stage 0: Onset of pregnancy to first ANC")
update_stage(first_ANC, second_ANC - 1, inc_preg_late, first_ANC_test, att_firstANC, test_sens_c, "Stage 1: First ANC to second ANC")
update_stage(second_ANC, delivery - 1, inc_preg_late, late_ANC_test, att_lategest, test_sens_e, "Stage 2: Second ANC to Delivery")
update_stage(delivery, delivery, inc_preg_late, delivery_test_prob, att_delivery, test_sens_e, "Stage 3: Delivery")
update_stage(delivery + 1, delivery + early_pp_visit - 1, inc_pp_early, pp_early_test, att_6wk, test_sens_e, "Stage 4: 1-6 weeks post-partum")
update_stage(delivery + early_pp_visit, delivery + pp_14weeks_visit - 1, inc_pp_mid, pp_14weeks_test, att_14wk, test_sens_e, "Stage 5: 6-14 weeks post-partum")
update_stage(delivery + pp_14weeks_visit, delivery + pp_mid_visit - 1, inc_pp_mid, pp_mid_test, att_6mo, test_sens_e, "Stage 6: 14 weeks - 6 months post-partum")
update_stage(delivery + pp_mid_visit, delivery + pp_9months_visit - 1, inc_pp_late, pp_9months_test, att_9mo, test_sens_e, "Stage 7: 6-9 months post-partum")
update_stage(delivery + pp_9months_visit, delivery + end_of_model, inc_pp_late, pp_9months_test, att_9mo, test_sens_e, "Stage 8: 9-12 months post-partum")

