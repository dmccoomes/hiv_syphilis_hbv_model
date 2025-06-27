
#### Title: HIV syphilis retesting model
#### Purpose: Create HIV retesting model based on Excel models
#### Author: Akash Malhotra, David Coomes
#### Date: January 27, 2025

rm(list = ls())

# Load libraries
pacman::p_load(
  "EpiModel", 
  "tidyverse", 
  "here",
  "data.table",
  "dplyr"
)

parameters <- read.csv("inp_par.csv")
# DMC adding to include pathway to updated parameters
parameters <- read.csv(here("parameters", "inp_par.csv"))

PrEP <- 0


# only for the selected combinations
selected_params <- parameters %>%
  filter(model == "HIV_Syphilis_Dual_Preg", #will make a loop later
         country == "South Africa", #will make a loop later
         target_pop == "Pregnant population") #will make a loop later


for (i in 1:nrow(selected_params)) {
  param_name <- as.character(selected_params$var_name[i])
  
  assign(param_name, selected_params$value[i])
  assign(paste0(param_name, "_lb"), selected_params$value_lb[i])  # Assign lower bound
  assign(paste0(param_name, "_ub"), selected_params$value_ub[i])  # Assign upper bound
}

# Starting states
neg <- (1-hiv_prev) - (inc_preg_early*duration_of_recent_infection)
recent <- duration_of_recent_infection*inc_preg_early
established_uk_nart <- hiv_prev*(1-hiv_status_known)
established_k_nart <- hiv_prev*hiv_status_known*(1-p_ART)
established_art <- hiv_prev*hiv_status_known*p_ART

# set the test sensitivity by selecting e3 or e4
test_sens_e <- sens_e3 # replace or put an if test_sens_e <- ifelse(test_gen==3, sens_e3, sens_e4)
test_sens_c <- sens_c3
test_spec <- spec3

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
maternal_df[1, ] <- c(1, popsize * neg, popsize * recent, 0, 
                      popsize * established_uk_nart, popsize * established_k_nart, 
                      popsize * established_art, 0, 0)

maternal_df$hiv_neg <- as.numeric(maternal_df$hiv_neg)
maternal_df$ahiv_pos_nart <- as.numeric(maternal_df$ahiv_pos_nart)
maternal_df$ahiv_pos_art <- as.numeric(maternal_df$ahiv_pos_art)
maternal_df$chiv_pos_uk <- as.numeric(maternal_df$chiv_pos_uk)
maternal_df$chiv_pos_nart <- as.numeric(maternal_df$chiv_pos_nart)
maternal_df$chiv_pos_art <- as.numeric(maternal_df$chiv_pos_art)
maternal_df$deaths <- as.numeric(maternal_df$deaths)
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

first_ANC <- weekofanc1
second_ANC <- weekofanc2
delivery <- deliveryweek

early_pp_visit <- 6
pp_14weeks_visit <- 14
pp_mid_visit <- 26
pp_9months_visit <- 39
end_of_model <- 52

# update_stage <- function(start_week, end_week, inc_rate, test_prob, att_prob, test_sens, stage_number)
#first_ANC_test and late_ANC_test and delivery_test_prob probability set to 1 (the 4th parameter in the function)
# Also pp_early_test set to 1
# Update different pregnancy stages, pp_14weeks_test, pp_mid_test, pp_9months_test, all set to 1
update_stage(2, first_ANC - 1, inc_preg_early, 1, att_firstANC, test_sens_c, 0)
update_stage(first_ANC, second_ANC - 1, inc_preg_late, 1, att_firstANC, test_sens_c, 1)
update_stage(second_ANC, delivery - 1, inc_preg_late, 1, att_lategest, test_sens_c, 2)
update_stage(delivery, delivery, inc_preg_late, 1, att_delivery, test_sens_e, 3)
update_stage(delivery + 1, delivery + early_pp_visit - 1, inc_pp_early, 0, att_6wk, test_sens_e, 4)
update_stage(delivery + early_pp_visit, delivery + pp_14weeks_visit - 1, inc_pp_mid, 0, att_14wk, test_sens_e, 5)
update_stage(delivery + pp_14weeks_visit, delivery + pp_mid_visit - 1, inc_pp_mid, 0, att_6mo, test_sens_e, 6)
update_stage(delivery + pp_mid_visit, delivery + pp_9months_visit - 1, inc_pp_late, 0, att_9mo, test_sens_e, 7)
update_stage(delivery + pp_9months_visit, delivery + end_of_model, inc_pp_late, 0, att_9mo, test_sens_e, 8)

