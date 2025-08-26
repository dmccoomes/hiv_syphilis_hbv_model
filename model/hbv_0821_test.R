rm(list=ls())
# Load required libraries
library(tidyverse)
library(dplyr)
library(tidyr)

# Read parameter file
params <- read_csv("hbv_model_parameter_table.csv", show_col_types = FALSE) %>%
  filter(Population == "Pregnant", !is.na(`Value PE`)) %>%
  group_by(Symbol) %>%
  slice_tail(n = 1) %>%         # if duplicates per Symbol, keep last row
  ungroup()
# Define testing coverage by arm

prop_tested_list <- list(
  SOC = .31,    # Standard of care: 0% or 31% tested
  INT = 0.749   # Intervention: 80% tested 0.749
)


# Convert to named vector for easy access
param_vals <- params %>%
  filter(!is.na(`Value PE`)) %>% ## will pull in upper and lower bounds later to account for uncertainty
  select(Symbol, `Value PE`) %>%
  deframe()

# Define HBV health states for adults and infants
adult_states <- c("HBV_susceptible", "Acute_symptomatic_HBV", "HBV_immune",
                  "Chronic_HBV", "Compensated_cirrhosis", "Decompensated_cirrhosis",
                  "Liver_cancer", "Death")

infant_states <- c("HBV_negative_infant", "HBV_positive_infant", "Infant_immune",
                   "Infant_chronic_HBV", "Infant_compensated_cirrhosis",
                   "Infant_decompensated_cirrhosis", "Infant_liver_cancer", "Infant_death")

# Replace with realistic number of pregnant women else fine to leave it as below
N0 <- 1631470

HBV_susceptible         <- N0 * param_vals[["hbv_p_sus"]]
Acute_symptomatic_HBV   <- N0 * param_vals[["hbv_p_preg_acute"]]
Chronic_HBV             <- N0 * param_vals[["hbv_p_preg_chronic"]]
HBV_immune              <- N0 - (HBV_susceptible + Acute_symptomatic_HBV + Chronic_HBV)

state_vector <- c(
  HBV_susceptible       = HBV_susceptible,
  Acute_symptomatic_HBV = Acute_symptomatic_HBV,
  Chronic_HBV           = Chronic_HBV,
  HBV_immune            = HBV_immune
)



# Define all health states (adult + infant + coinfected)
adult_coinfection_states <- c("Adult_HIV_HBV_coinf", "Adult_Syphilis_HBV_coinf", "Adult_HIV_Syphilis_HBV_coinf")
infant_coinfection_states <- c("Infant_HIV_HBV_coinf", "Infant_Syphilis_HBV_coinf", "Infant_HIV_Syphilis_HBV_coinf")

coinfection_states <- c(adult_coinfection_states, infant_coinfection_states)

all_states <- unique(c(adult_states, infant_states, coinfection_states))


# Load mortality file
mort_df <- read_csv("mort_nonhiv.csv") %>%
  rename(
    age_group = `Age Group`,
    mort_both = `Both sexes`,
    mort_female = Female
  ) %>%
  mutate(
    age_start = as.numeric(str_extract(age_group, "^\\d+")),
    age_start = if_else(str_detect(age_group, "<1"), 0, age_start)
  )

# Helper function to pull mortality by age and sex
get_mortality <- function(age, sex = "female") {
  start_age <- if (age < 1) 0 else {
    possible <- sort(unique(na.omit(mort_df$age_start)))
    max(possible[possible <= age])
  }
  row <- mort_df %>% filter(age_start == start_age)
  if (sex == "female") return(row$mort_female)
  else return(row$mort_both)
}

# Initialize transition matrix list
transition_matrices <- list()

# Loop over time horizon (since BG mortlaity rates change each year)
for (time in 0:20) {
  # Initialize transition matrix for this year
  tm <- matrix(0, nrow = length(all_states), ncol = length(all_states),
               dimnames = list(from = all_states, to = all_states))
  
  # ====================
  
  # ADULT transitions 
  tm["HBV_susceptible", "Acute_symptomatic_HBV"] <- param_vals["hbv_r_inf"]
  tm["HBV_susceptible", "HBV_susceptible"] <- 1 - param_vals["hbv_r_inf"]
  
  tm["Acute_symptomatic_HBV", "HBV_immune"]  <- param_vals["hbv_r_acute_to_immune"]
  tm["Acute_symptomatic_HBV", "Chronic_HBV"] <- param_vals["hbv_r_acute_to_chronic"]
  tm["Acute_symptomatic_HBV", "Acute_symptomatic_HBV"] <- 1 - (
    param_vals["hbv_r_acute_to_immune"] + param_vals["hbv_r_acute_to_chronic"]
  )
  
  tm["Chronic_HBV", "Compensated_cirrhosis"] <- param_vals["hbv_r_chronic_to_comp_cirr"]
  tm["Chronic_HBV", "Death"] <- param_vals["hbv_r_chronic_death"]
  tm["Chronic_HBV", "Chronic_HBV"] <- 1 - (
    param_vals["hbv_r_chronic_to_comp_cirr"] + param_vals["hbv_r_chronic_death"]
  )
  
  tm["Compensated_cirrhosis", "Decompensated_cirrhosis"] <- param_vals["hbv_r_comp_to_decomp"]
  tm["Compensated_cirrhosis", "Liver_cancer"] <- param_vals["hbv_r_comp_to_cancer"]
  tm["Compensated_cirrhosis", "Compensated_cirrhosis"] <- 1 - (
    param_vals["hbv_r_comp_to_decomp"] + param_vals["hbv_r_comp_to_cancer"]
  )
  
  tm["Decompensated_cirrhosis", "Liver_cancer"] <- param_vals["hbv_r_decomp_to_cancer"]
  tm["Decompensated_cirrhosis", "Death"] <- param_vals["hbv_r_decomp_death"]
  tm["Decompensated_cirrhosis", "Decompensated_cirrhosis"] <- 1 - (
    param_vals["hbv_r_decomp_to_cancer"] + param_vals["hbv_r_decomp_death"]
  )
  
  tm["Liver_cancer", "Death"] <- param_vals["hbv_r_cancer_death"]
  tm["Liver_cancer", "Liver_cancer"] <- 1 - param_vals["hbv_r_cancer_death"]
  
  # INFANT transitions 
  tm["HBV_positive_infant", "Infant_chronic_HBV"] <- param_vals["hbv_r_infant_sym_to_chronic"]
  tm["HBV_positive_infant", "Infant_immune"] <- param_vals["hbv_r_infant_sym_to_immune"]
  tm["HBV_positive_infant", "HBV_positive_infant"] <- 1 - (
    param_vals["hbv_r_infant_sym_to_chronic"] + param_vals["hbv_r_infant_sym_to_immune"]
  )
  
  tm["Infant_chronic_HBV", "Infant_compensated_cirrhosis"] <- param_vals["hbv_r_chronic_to_comp_cirr_infant"]
  tm["Infant_chronic_HBV", "Infant_chronic_HBV"] <- 1 - param_vals["hbv_r_chronic_to_comp_cirr_infant"]
  
  tm["Infant_compensated_cirrhosis", "Infant_decompensated_cirrhosis"] <- param_vals["hbv_r_comp_to_decomp_infant"]
  tm["Infant_compensated_cirrhosis", "Infant_liver_cancer"] <- param_vals["hbv_r_comp_to_cancer_infant"]
  tm["Infant_compensated_cirrhosis", "Infant_compensated_cirrhosis"] <- 1 - (
    param_vals["hbv_r_comp_to_decomp_infant"] + param_vals["hbv_r_comp_to_cancer_infant"]
  )
  
  tm["Infant_decompensated_cirrhosis", "Infant_liver_cancer"] <- param_vals["hbv_r_decomp_to_cancer_infant"]
  tm["Infant_decompensated_cirrhosis", "Infant_death"] <- param_vals["hbv_r_decomp_death_infant"]
  tm["Infant_decompensated_cirrhosis", "Infant_decompensated_cirrhosis"] <- 1 - (
    param_vals["hbv_r_decomp_to_cancer_infant"] + param_vals["hbv_r_decomp_death_infant"]
  )
  
  tm["Infant_liver_cancer", "Infant_death"] <- param_vals["hbv_r_cancer_death_infant"]
  tm["Infant_liver_cancer", "Infant_liver_cancer"] <- 1 - param_vals["hbv_r_cancer_death_infant"]
  
  # --- Acute must clear within 1 year (no self-loops) ---
  tm["Acute_symptomatic_HBV", "Acute_symptomatic_HBV"] <- 0
  tm["HBV_positive_infant",  "HBV_positive_infant"]    <- 0
  
  # Fallbacks in case the parameter file had 0 prob. to both destinations
  if ((param_vals["hbv_r_acute_to_immune"] + param_vals["hbv_r_acute_to_chronic"]) == 0) {
    tm["Acute_symptomatic_HBV", "HBV_immune"] <- 1
  }
  if ((param_vals["hbv_r_infant_sym_to_chronic"] + param_vals["hbv_r_infant_sym_to_immune"]) == 0) {
    tm["HBV_positive_infant", "Infant_immune"] <- 1
  }
  
  # -------- Row-normalized background mortality (REPLACEMENT) --------
  # Make coinfection states absorbing (we handle their evolution externally)
  tm[coinfection_states, ] <- 0
  tm[cbind(coinfection_states, coinfection_states)] <- 1
  
  # Helper: rescale a row so that non-death transitions sum to survival_mass, and death gets the rest
  renorm_row <- function(mat, st, death_col, survival_mass) {
    non_death_cols <- setdiff(colnames(mat), c("Death","Infant_death"))
    row_vals <- mat[st, non_death_cols]
    s <- sum(row_vals)
    if (s > 0) {
      mat[st, non_death_cols] <- (row_vals / s) * survival_mass
    } else {
      # if no defined transitions, keep state with full survival mass
      mat[st, st] <- survival_mass
    }
    mat[st, death_col] <- 1 - survival_mass
    mat
  }
  
  # Background mortality for this year
  age_female <- 22 + time
  age_infant <- 0  + time
  bg_mort_female <- get_mortality(age_female, "female")
  bg_mort_infant <- get_mortality(age_infant - 1, "both")
  
  # Adults: normalize each adult row to (1 - bg_mort_female) and give bg_mort_female to Death
  for (st in setdiff(adult_states, "Death")) {
    tm <- renorm_row(tm, st, "Death", 1 - bg_mort_female)
  }
  
  # Infants: normalize each infant row to (1 - bg_mort_infant) and give bg_mort_infant to Infant_death
  for (st in setdiff(infant_states, "Infant_death")) {
    tm <- renorm_row(tm, st, "Infant_death", 1 - bg_mort_infant)
  }
  
  # Ensure death states are absorbing
  tm["Death", ] <- 0;         tm["Death", "Death"] <- 1
  tm["Infant_death", ] <- 0;  tm["Infant_death", "Infant_death"] <- 1
  
  transition_matrices[[time + 1]] <- tm
}



arm_results_list <- list()


## Main loop is below
for (arm in names(prop_tested_list)) {
  
  time_horizon <- 20
  
  # Re-seed state vector for each arm
  state_vector <- setNames(rep(0, length(all_states)), all_states)
  state_vector["HBV_susceptible"] <- N0 * param_vals["hbv_p_sus"]
  state_vector["Acute_symptomatic_HBV"] <- N0 * param_vals["hbv_p_preg_acute"]
  state_vector["Chronic_HBV"] <- N0 * param_vals["hbv_p_preg_chronic"]
  state_vector["HBV_immune"] <- N0 - (
    state_vector["HBV_susceptible"] +
      state_vector["Acute_symptomatic_HBV"] +
      state_vector["Chronic_HBV"]
  )
  
  # Set infant and coinfection states to zero
  missing_states <- setdiff(all_states, names(state_vector))
  for (s in missing_states) {
    state_vector[s] <- 0
  }
  state_vector <- state_vector[all_states]  # consistent order
  
  # Re-initialize results matrix
  results <- matrix(0, nrow = time_horizon + 1, ncol = length(all_states),
                    dimnames = list(time = 0:time_horizon, state = all_states))
  
  results[1, ] <- state_vector
  
  prop_tested <- prop_tested_list[[arm]]


  
  # Initialize full state vector
  full_state_vector <- setNames(rep(0, length(all_states)), all_states)
  full_state_vector[names(state_vector)] <- state_vector
  
  # Store state over time
  results <- matrix(0, nrow = time_horizon + 1, ncol = length(all_states),
                    dimnames = list(time = 0:time_horizon, state = all_states))
  
  results[1, ] <- full_state_vector
  
  # Track deaths 
  deaths_by_agegroup <- data.frame(
    time = 0:time_horizon,
    infant_deaths = numeric(time_horizon+1),
    adult_deaths = numeric(time_horizon+1),
    hiv_syph_infant_deaths = numeric(time_horizon+1),
    hiv_syph_adult_deaths = numeric(time_horizon+1)
    
  )
  

  
  # Simulate MTCT once at time 0 (single cohort model)
  #full_state_vector <- simulate_mtct(full_state_vector, param_vals, N0)
  
  state_vector <- state_vector[all_states]  # Force correct order
  
  # Initialize results matrix
  results[1, ] <- state_vector  # t = 0
  
  # Read this file
  # Read and pivot hiv_syph data (for coinfection adjustments)
  hiv_syph_raw <- read_csv("hiv_syph_raw_triplex.csv")
  
  
  # Compute death increments for both infants and adults (if cumulative deaths are provided)
  hiv_syph_inc <- hiv_syph_raw %>%
    arrange(arm, age) %>%
    group_by(arm) %>%
    mutate(
      HIV_infant_death_inc = HIV_infant_death - lag(HIV_infant_death, default = 0),
      Syphilis_infant_death_inc = Syphilis_infant_death - lag(Syphilis_infant_death, default = 0),
      HIV_adult_death_inc = HIV_adult_death - lag(HIV_adult_death, default = 0),
      Syphilis_adult_death_inc = Syphilis_adult_death - lag(Syphilis_adult_death, default = 0)
    ) %>%
    ungroup()
  
  # Pivot the long format for use in disaggregating coinfections
  hiv_syph_long <- hiv_syph_raw %>%
    select(age, arm,
           hiv_infected_infants, syphilis_infected_infants, coinfected_infants,
           HIV_infant_death, Syphilis_infant_death,
           hiv_infected_adults, syphilis_infected_adults, coinfected_adults,
           HIV_adult_death, Syphilis_adult_death,
           lbw_infants  # optional
    ) %>%
    pivot_longer(
      cols = -c(age, arm),
      names_to = "state",
      values_to = "count"
    ) %>%
    rename(time = age)
  
  
  # Define infant states that represent living infants (non-death states)
  infant_living_states <- c(
    "HBV_negative_infant", "HBV_positive_infant", "Infant_immune",
    "Infant_chronic_HBV", "Infant_compensated_cirrhosis",
    "Infant_decompensated_cirrhosis", "Infant_liver_cancer"
  )
  

  for (t in 0:time_horizon) {
    
    #  Transition ALL existing states (adults + infants)
    state_vector[is.na(state_vector)] <- 0
    next_state_vector <- as.numeric(state_vector[all_states]) %*% transition_matrices[[t+1]]
    next_state_vector <- as.vector(next_state_vector)
    names(next_state_vector) <- all_states
    
    if (t %in% c(0,1,2)) {
      cat("[", arm, "] t=", t,
          " pre-MTCT HBV_pos_inf=", next_state_vector["HBV_positive_infant"],
          " Inf_chronic=", next_state_vector["Infant_chronic_HBV"], "\n", sep = "")
    }
    
    
    # state before MTCT is applied at t = 1
    if (t == 0) {
      results[t + 1, ] <- state_vector  # true pre-birth values
    }
    
    # Accuracy-aware MTCT (no DNA confirm; Spec only affects cost, not infections)
    simulate_mtct <- function(state_vector, param_vals, n_births, arm, prop_tested) {
      # params
      births_chronic <- n_births * param_vals["hbv_p_preg_chronic"]
      births_acute   <- n_births * param_vals["hbv_p_preg_acute"]
      
      # choose sensitivity by arm
      triplex_sens <- ifelse(is.na(param_vals["triplex_sens"]), 1, param_vals["triplex_sens"])
      hbsag_sens   <- ifelse(is.na(param_vals["hbsag_sens"]),   1, param_vals["hbsag_sens"])
      sens <- if (arm == "INT") triplex_sens else hbsag_sens
      
      # who is "managed" among infected moms (known + tested*positive)
      p_known   <- param_vals["hbv_p_known_hbsag"]
      m         <- p_known + (1 - p_known) * prop_tested * sens
      m         <- pmin(pmax(m, 0), 1)
      
      # maternal TDF (applies only to managed infected mothers)
      p_tdf     <- param_vals["hbv_p_tenofovir_preg"]
      eff_tdf   <- param_vals["hbv_eff_tenofovir"]
      
      # infant prophylaxis (HBIG + birth dose)
      p_hbig    <- param_vals["hbv_p_infant_hbig"]
      p_vax     <- param_vals["hbv_p_infant_vaccine_birth"]
      eff_hbig  <- param_vals["hbv_eff_hbig"]
      eff_vax   <- param_vals["hbv_eff_vaccine_birth"]
      eff_inf   <- 1 - ((1 - eff_hbig * p_hbig) * (1 - eff_vax * p_vax))  # 0..1
      
      # base MTCT risks
      rC <- param_vals["hbv_mtct_chronic"]
      rA <- param_vals["hbv_mtct_acute"]
      
      # FINAL infections after maternal TDF (managed only) and infant prophylaxis (managed only)
      inf_chronic_managed   <- births_chronic * m       * rC * (1 - p_tdf * eff_tdf) * (1 - eff_inf)
      inf_chronic_unmanaged <- births_chronic * (1 - m) * rC
      inf_acute_managed     <- births_acute   * m       * rA * (1 - eff_inf)          # no TDF for acute
      inf_acute_unmanaged   <- births_acute   * (1 - m) * rA
      
      infected_final <- inf_chronic_managed + inf_chronic_unmanaged + inf_acute_managed + inf_acute_unmanaged
      
      # infants protected by prophylaxis (for the Immune bin)
      protected_by_proph <- births_chronic * m * rC * (1 - p_tdf * eff_tdf) * eff_inf +
        births_acute   * m * rA * eff_inf
      
      # allocate
      immune_add <- protected_by_proph
      hbv_pos_add <- infected_final
      hbv_neg_add <- max(0, n_births - immune_add - hbv_pos_add)  # remainder
      
      state_vector["Infant_immune"]         <- state_vector["Infant_immune"]       + immune_add
      state_vector["HBV_positive_infant"]   <- state_vector["HBV_positive_infant"] + hbv_pos_add
      state_vector["HBV_negative_infant"]   <- state_vector["HBV_negative_infant"] + hbv_neg_add
      
      state_vector
    }
    
    
    
    cat("[", arm, "] t=1 after MTCT: HBV_pos_inf=",
        next_state_vector["HBV_positive_infant"],
        " Infant_immune=", next_state_vector["Infant_immune"],
        " HBV_neg_inf=", next_state_vector["HBV_negative_infant"], "\n", sep = "")
    
    
    if (t == 1) {
      births <- N0
      next_state_vector <- simulate_mtct(next_state_vector, param_vals, births, arm, prop_tested)
    }
    
    
# --- Allocate HIV/Syph coinfections within HBV pools (every year from t >= 1) ---
    if (t >= 1) {
      # HBV disease pools to split
      adult_hbv_states_for_split  <- c("Chronic_HBV","Compensated_cirrhosis","Decompensated_cirrhosis","Liver_cancer")
      infant_hbv_states_for_split <- c("HBV_positive_infant","Infant_chronic_HBV",
                                       "Infant_compensated_cirrhosis","Infant_decompensated_cirrhosis","Infant_liver_cancer")
      
      # full living denominators for prevalence
      adult_all_living  <- c("HBV_susceptible","Acute_symptomatic_HBV","HBV_immune", adult_hbv_states_for_split)
      infant_all_living <- c("HBV_negative_infant","Infant_immune", infant_hbv_states_for_split)
      
      adult_hbv_pool   <- sum(next_state_vector[adult_hbv_states_for_split])
      infant_hbv_pool  <- sum(next_state_vector[infant_hbv_states_for_split])
      
      adult_total_pool  <- sum(next_state_vector[adult_all_living])
      infant_total_pool <- sum(next_state_vector[infant_all_living])
      
      get_count <- function(df, tt, arm_name, nm) {
        v <- df %>%
          dplyr::filter(time == tt, arm == !!arm_name, state == nm) %>%
          dplyr::summarise(v = sum(count, na.rm = TRUE)) %>%
          dplyr::pull(v)
        if (length(v) == 0 || is.na(v)) 0 else v
      }
      
      # external prevalences among all infants/adults
      n_adult_HIV   <- get_count(hiv_syph_long, t, arm, "hiv_infected_adults")
      n_adult_Syph  <- get_count(hiv_syph_long, t, arm, "syphilis_infected_adults")
      n_adult_pair  <- get_count(hiv_syph_long, t, arm, "coinfected_adults")       # HIV+Syph (no HBV)
      
      n_infant_HIV  <- get_count(hiv_syph_long, t, arm, "hiv_infected_infants")
      n_infant_Syph <- get_count(hiv_syph_long, t, arm, "syphilis_infected_infants")
      n_infant_pair <- get_count(hiv_syph_long, t, arm, "coinfected_infants")      # HIV+Syph (no HBV)
      
      clamp01 <- function(x) pmin(pmax(x, 0), 1)
      
      p_adult_HIV   <- ifelse(adult_total_pool  > 0, n_adult_HIV  / adult_total_pool, 0)
      p_adult_Syph  <- ifelse(adult_total_pool  > 0, n_adult_Syph / adult_total_pool, 0)
      p_adult_pair  <- ifelse(adult_total_pool  > 0, n_adult_pair / adult_total_pool, 0)
      
      p_infant_HIV  <- ifelse(infant_total_pool > 0, n_infant_HIV / infant_total_pool, 0)
      p_infant_Syph <- ifelse(infant_total_pool > 0, n_infant_Syph / infant_total_pool, 0)
      p_infant_pair <- ifelse(infant_total_pool > 0, n_infant_pair / infant_total_pool, 0)
      
      p_adult_HIV   <- clamp01(p_adult_HIV)
      p_adult_Syph  <- clamp01(p_adult_Syph)
      p_adult_pair  <- clamp01(pmin(p_adult_pair,  p_adult_HIV,  p_adult_Syph))
      
      p_infant_HIV  <- clamp01(p_infant_HIV)
      p_infant_Syph <- clamp01(p_infant_Syph)
      p_infant_pair <- clamp01(pmin(p_infant_pair, p_infant_HIV, p_infant_Syph))
      
      # shares (HBV-only remainder is 1 - (HIV + Syph - pair))
      share_none_adult  <- clamp01(1 - (p_adult_HIV + p_adult_Syph - p_adult_pair))
      share_none_infant <- clamp01(1 - (p_infant_HIV + p_infant_Syph - p_infant_pair))
      
      n_adult_triple    <- p_adult_pair                  * adult_hbv_pool
      n_adult_hiv_only  <- (p_adult_HIV  - p_adult_pair) * adult_hbv_pool
      n_adult_syph_only <- (p_adult_Syph - p_adult_pair) * adult_hbv_pool
      n_adult_hbv_only  <- share_none_adult              * adult_hbv_pool
      
      n_infant_triple    <- p_infant_pair                   * infant_hbv_pool
      n_infant_hiv_only  <- (p_infant_HIV  - p_infant_pair) * infant_hbv_pool
      n_infant_syph_only <- (p_infant_Syph - p_infant_pair) * infant_hbv_pool
      n_infant_hbv_only  <- share_none_infant              * infant_hbv_pool
      
      # numeric safety
      nz <- function(x) ifelse(is.finite(x), x, 0)
      n_adult_triple    <- nz(n_adult_triple)
      n_adult_hiv_only  <- nz(n_adult_hiv_only)
      n_adult_syph_only <- nz(n_adult_syph_only)
      n_adult_hbv_only  <- nz(n_adult_hbv_only)
      
      n_infant_triple    <- nz(n_infant_triple)
      n_infant_hiv_only  <- nz(n_infant_hiv_only)
      n_infant_syph_only <- nz(n_infant_syph_only)
      n_infant_hbv_only  <- nz(n_infant_hbv_only)
      
      # write to HBV coinfection compartments
     # next_state_vector["Adult_HIV_HBV_coinf"]          <- n_adult_hiv_only
    #  next_state_vector["Adult_Syphilis_HBV_coinf"]     <- n_adult_syph_only
     # next_state_vector["Adult_HIV_Syphilis_HBV_coinf"] <- n_adult_triple
      
    #  next_state_vector["Infant_HIV_HBV_coinf"]         <- n_infant_hiv_only
     # next_state_vector["Infant_Syphilis_HBV_coinf"]    <- n_infant_syph_only
    #  next_state_vector["Infant_HIV_Syphilis_HBV_coinf"]<- n_infant_triple
      
      # scale the base HBV disease states down to HBV-only remainder
      scale_adult  <- ifelse(adult_hbv_pool  > 0, n_adult_hbv_only  / adult_hbv_pool,  0)
      scale_infant <- ifelse(infant_hbv_pool > 0, n_infant_hbv_only / infant_hbv_pool, 0)
      
    #  next_state_vector[adult_hbv_states_for_split]  <- next_state_vector[adult_hbv_states_for_split]  * scale_adult
     # next_state_vector[infant_hbv_states_for_split] <- next_state_vector[infant_hbv_states_for_split] * scale_infant
    }
    

    if (t > 0) {    
    #  Save updated state
    results[t + 1, ] <- next_state_vector
    }
    
    # Adjust for HIV/Syphilis deaths from external CSV (infants + adults)
    if (t < time_horizon) {
      
      # --- INFANT DEATH ADJUSTMENT ---
      hiv_syph_deaths_infant_t <- hiv_syph_inc %>%
        filter(time == t, arm == !!arm) %>%
        summarise(
          total_deaths = sum(HIV_infant_death_inc + Syphilis_infant_death_inc, na.rm = TRUE)
        ) %>%
        pull(total_deaths)
      
      infant_living_states <- c("HBV_negative_infant", "HBV_positive_infant", "Infant_immune",
                                "Infant_chronic_HBV", "Infant_compensated_cirrhosis",
                                "Infant_decompensated_cirrhosis", "Infant_liver_cancer",
                                "Infant_HIV_HBV_coinf", "Infant_Syphilis_HBV_coinf", "Infant_HIV_Syphilis_HBV_coinf")
      
      infant_counts <- next_state_vector[infant_living_states]
      total_infant_alive <- sum(infant_counts)
      
      if (isTRUE(total_infant_alive > 0) && isTRUE(hiv_syph_deaths_infant_t > 0)) {
        prop_in_each_state <- infant_counts / total_infant_alive
        subtract_vector <- prop_in_each_state * hiv_syph_deaths_infant_t
        
        next_state_vector[infant_living_states] <- pmax(
          0, next_state_vector[infant_living_states] - subtract_vector
        )
      }
      
      # --- ADULT DEATH ADJUSTMENT ---
      hiv_syph_deaths_adult_t <- hiv_syph_inc %>%
        filter(time == t, arm == !!arm) %>%
        summarise(
          total_deaths = sum(HIV_adult_death_inc + Syphilis_adult_death_inc, na.rm = TRUE)
        ) %>%
        pull(total_deaths)
      
      adult_living_states <- c("Adult_HIV_HBV_coinf", "Adult_Syphilis_HBV_coinf", "Adult_HIV_Syphilis_HBV_coinf")
      
      adult_counts <- next_state_vector[adult_living_states]
      total_adult_alive <- sum(adult_counts)
      
      if (isTRUE(total_adult_alive > 0) && isTRUE(hiv_syph_deaths_adult_t > 0)) {
        prop_in_each_state <- adult_counts / total_adult_alive
        subtract_vector <- prop_in_each_state * hiv_syph_deaths_adult_t
        
        next_state_vector[adult_living_states] <- pmax(
          0, next_state_vector[adult_living_states] - subtract_vector
        )
      }
    }
    
    
    # Save HIV/Syph deaths for tracking
    deaths_by_agegroup$hiv_syph_infant_deaths[t + 1] <- ifelse(exists("hiv_syph_deaths_infant_t"), hiv_syph_deaths_infant_t, 0)
    deaths_by_agegroup$hiv_syph_adult_deaths[t + 1]  <- ifelse(exists("hiv_syph_deaths_adult_t"), hiv_syph_deaths_adult_t, 0)
    
    
    #  Update state for next loop
    state_vector <- next_state_vector
    
    discount_rate <- 0.03
    discount_factors <- 1 / ((1 + discount_rate) ^ (0:time_horizon))
    
    births <- ifelse(t == 1, N0, 0)
    
    
  }
  
  # Convert matrix to tidy format
 results_df <- as.data.frame(results) %>%
    rownames_to_column("time") %>%
    mutate(time = as.integer(time),
           arm = arm) %>%
     pivot_longer(-c(time, arm), names_to = "state", values_to = "count")
  
# MTCT function using dynamic effectiveness-based risk


# Simulation settings
years <- 0:20
n_years <- length(years)

# Define discount rate and discount factor vector
discount_rate <- 0.03
discount_factors <- 1 / (1 + discount_rate) ^ years


# Reorder safely
state_vector <- state_vector[all_states]

#Ensure Death is exactly zero
state_vector["Death"] <- 0
state_vector["Infant_death"] <- 0

# Simulation settings
time_horizon <- 20


# ✅ Pivot both infant and adult coinfection data to match results_df structure
hiv_syph_long <- hiv_syph_raw %>%
  select(
    age, arm,
    hiv_infected_infants, syphilis_infected_infants, coinfected_infants,
    HIV_infant_death, Syphilis_infant_death,
    hiv_infected_adults, syphilis_infected_adults, coinfected_adults,
    HIV_adult_death, Syphilis_adult_death,
    lbw_infants  # optional
  ) %>%
  pivot_longer(
    cols = -c(age, arm),
    names_to = "state",
    values_to = "count"
  ) %>%
  rename(time = age)


# ppend to results_df
results_df <- bind_rows(results_df, hiv_syph_long %>% filter(arm == !!arm))

# Save result
arm_results_list[[arm]] <- results_df

  }



  final_results <- bind_rows(arm_results_list)
  
  # Build an all-cause infant death series: background + HIV + Syph
  all_infant_deaths <- final_results %>%
    dplyr::filter(state %in% c("Infant_death","HIV_infant_death","Syphilis_infant_death")) %>%
    dplyr::group_by(arm, time) %>%
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(state = "All_infant_deaths")
  
  final_results <- dplyr::bind_rows(final_results, all_infant_deaths)
  
  # All-cause adult deaths: background + HIV + Syph
  all_adult_deaths <- final_results %>%
    dplyr::filter(state %in% c("Death","HIV_adult_death","Syphilis_adult_death")) %>%
    dplyr::group_by(arm, time) %>%
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(state = "All_adult_deaths")
  
  final_results <- dplyr::bind_rows(final_results, all_adult_deaths)
  

# adult plot
  plot_states_adult <- c(setdiff(adult_states, "Death"), "All_adult_deaths")
  
  ggplot(final_results %>% dplyr::filter(state %in% plot_states_adult),
         aes(x = time, y = count, color = arm, linetype = arm)) +
    geom_line(size = 1) +
    facet_wrap(~state, scales = "free_y") +
    scale_linetype_manual(values = c("solid","dashed")) +
    labs(title = "Adult HBV States Over Time (with All-Cause Deaths): SoC vs Intervention",
         x = "Time (years)", y = "Number of People") +
    theme_minimal()
  
  
# infant plot
  plot_states <- c(setdiff(infant_states, "Infant_death"), "All_infant_deaths")
  
  ggplot(final_results %>% dplyr::filter(state %in% plot_states),
         aes(x = time, y = count, color = arm, linetype = arm)) +
    geom_line(size = 1) +
    facet_wrap(~state, scales = "free_y") +
    scale_linetype_manual(values = c("solid","dashed")) +
    labs(title = "Infant HBV States Over Time (with Death): SoC vs Intervention",
         x = "Time (years)", y = "Number of People") +
    theme_minimal()
  



### DALYs

# Define duration adjustment for acute HBV states
acute_duration <- 1 / 3
acute_states <- c("Acute_symptomatic_HBV", "HBV_positive_infant")

# Create disability weight mapping (adult + infant states)
dw_lookup <- tibble::tribble(
  ~state,                         ~symbol,
  "Chronic_HBV",                  "hbv_dw_chronic",
  "Compensated_cirrhosis",        "hbv_dw_comp_cirrhosis",
  "Decompensated_cirrhosis",      "hbv_dw_decomp_cirrhosis",
  "Liver_cancer",                 "hbv_dw_liver_cancer",
  "Acute_symptomatic_HBV",        "hbv_dw_acute_sym_adult",
  "Infant_chronic_HBV",           "hbv_dw_chronic",
  "Infant_compensated_cirrhosis", "hbv_dw_comp_cirrhosis",
  "Infant_decompensated_cirrhosis", "hbv_dw_decomp_cirrhosis",
  "Infant_liver_cancer",          "hbv_dw_liver_cancer",
  "HBV_positive_infant",          "hbv_dw_acute_sym_adult",  # acute proxy
  "coinfected_infants",           "hiv_syph_dw_inf",
  "syphilis_infected_infants",    "syph_dw_inf",
  "lbw_infants",                 "syph_lbw_dw_inf",
  "hiv_infected_infants",         "hiv_dw_inf",
  "coinfected_adults",           "hiv_syph_dw_inf",
  "syphilis_infected_adults",    "syph_dw_inf",
 # "hiv_infected_adults",         "hiv_dw_inf",
  "HBV_immune",                  "hbv_dw_imm",
  "Infant_immune",               "hbv_dw_imm",
  "HBV_susceptible",             "hbv_dw_sus"
  
) %>%
  mutate(dw = param_vals[symbol])

#  Extract HBV-related states and their base DWs
hbv_state_dws <- tibble::tibble(
  base_state = c(
    "HBV_susceptible", "HBV_immune", "Chronic_HBV", "Compensated_cirrhosis", "Decompensated_cirrhosis", "Liver_cancer",
    "Acute_symptomatic_HBV",
    "Infant_immune", "Infant_chronic_HBV", "Infant_compensated_cirrhosis",
    "Infant_decompensated_cirrhosis", "Infant_liver_cancer",
    "HBV_positive_infant"  # if you're still using this as acute proxy
  ),
  dw_hbv = param_vals[c(
    "hbv_dw_sus", "hbv_dw_imm", "hbv_dw_chronic", "hbv_dw_comp_cirrhosis", "hbv_dw_decomp_cirrhosis", "hbv_dw_liver_cancer",
    "hbv_dw_acute_sym_adult",
    "hbv_dw_imm", "hbv_dw_chronic", "hbv_dw_comp_cirrhosis",
    "hbv_dw_decomp_cirrhosis", "hbv_dw_liver_cancer",
    "hbv_dw_acute_sym_adult"
  )]
)

# Define single-infection DWs
dw_hiv <- param_vals["hiv_dw_inf"]
dw_syph <- param_vals["syph_dw_inf"]

# Create coinfection DWs per HBV state
hbv_dw_coinfections <- hbv_state_dws %>%
  mutate(
    state_HIV   = paste0(base_state, "_HIV"),
    state_Syph  = paste0(base_state, "_Syph"),
    state_Triple= paste0(base_state, "_Triple"),
    dw_HIV      = 1 - (1 - dw_hbv) * (1 - dw_hiv),
    dw_Syph     = 1 - (1 - dw_hbv) * (1 - dw_syph),
    dw_Triple   = 1 - (1 - dw_hbv) * (1 - dw_hiv) * (1 - dw_syph)
  ) %>%
  select(
    HIV_state = state_HIV, HIV_dw = dw_HIV,
    Syph_state = state_Syph, Syph_dw = dw_Syph,
    Triple_state = state_Triple, Triple_dw = dw_Triple
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("type", ".value"),
    names_pattern = "(.*)_(state|dw)"
  ) %>%
  rename(state = state)



dw_lookup <- dw_lookup %>%
  bind_rows(hbv_dw_coinfections %>% rename(symbol = type))



# Compute person-years
final_df <- final_results %>%
  group_by(arm, state, time) %>%
  summarise(person_years = sum(count), .groups = "drop")

# Define HBV states for adults and infants
adult_hbv_states <- c("HBV_susceptible", "Acute_symptomatic_HBV", "HBV_immune",
                      "Chronic_HBV", "Compensated_cirrhosis", "Decompensated_cirrhosis", "Liver_cancer")

infant_hbv_states <- c("Infant_chronic_HBV", "Infant_compensated_cirrhosis",
                       "Infant_decompensated_cirrhosis", "Infant_liver_cancer", "HBV_positive_infant")


# Use these ONLY for the post-sim coinfection split (not for anything else)
hbv_infected_adult_split  <- c(
  "Chronic_HBV","Compensated_cirrhosis","Decompensated_cirrhosis","Liver_cancer"
)

hbv_infected_infant_split <- c(
  "HBV_positive_infant","Infant_chronic_HBV","Infant_compensated_cirrhosis",
  "Infant_decompensated_cirrhosis","Infant_liver_cancer"
)


# ----- Build denominators from the actual model states -----
infant_living_states_full <- c(
  "HBV_negative_infant","HBV_positive_infant","Infant_immune",
  "Infant_chronic_HBV","Infant_compensated_cirrhosis",
  "Infant_decompensated_cirrhosis","Infant_liver_cancer",
  "Infant_HIV_HBV_coinf","Infant_Syphilis_HBV_coinf","Infant_HIV_Syphilis_HBV_coinf"
)

adult_living_states_full <- c(
  "HBV_susceptible","Acute_symptomatic_HBV","HBV_immune",
  "Chronic_HBV","Compensated_cirrhosis","Decompensated_cirrhosis","Liver_cancer",
  "Adult_HIV_HBV_coinf","Adult_Syphilis_HBV_coinf","Adult_HIV_Syphilis_HBV_coinf"
)

pool_infant <- final_results %>%
  dplyr::filter(state %in% infant_living_states_full) %>%
  dplyr::group_by(arm, time) %>%
  dplyr::summarise(denom_infant = sum(count, na.rm = TRUE), .groups = "drop")

pool_adult <- final_results %>%
  dplyr::filter(state %in% adult_living_states_full) %>%
  dplyr::group_by(arm, time) %>%
  dplyr::summarise(denom_adult = sum(count, na.rm = TRUE), .groups = "drop")

# Person-years for each state (as before)
final_df <- final_results %>%
  dplyr::group_by(arm, state, time) %>%
  dplyr::summarise(person_years = sum(count), .groups = "drop")

# ========== INFANT SPLIT (uses full infant denominator) ==========
infant_df <- final_df %>%
  dplyr::filter(state %in% hbv_infected_infant_split) %>%
  dplyr::left_join(
    hiv_syph_raw %>%
      dplyr::select(age, arm, hiv_infected_infants, syphilis_infected_infants, coinfected_infants) %>%
      dplyr::rename(time = age),
    by = c("arm", "time")
  ) %>%
  dplyr::left_join(pool_infant, by = c("arm", "time")) %>%
  dplyr::group_by(arm, time) %>%
  dplyr::mutate(
    p_HIV   = if_else(denom_infant > 0, hiv_infected_infants      / denom_infant, 0),
    p_Syph  = if_else(denom_infant > 0, syphilis_infected_infants / denom_infant, 0),
    p_pair0 = if_else(denom_infant > 0, coinfected_infants        / denom_infant, 0),
    # inclusion–exclusion coherence
    p_pair  = pmin(pmax(p_pair0, p_HIV + p_Syph - 1, 0), pmin(p_HIV, p_Syph)),
    prop_triple = p_pair,
    prop_HIV    = pmax(p_HIV  - p_pair, 0),
    prop_Syph   = pmax(p_Syph - p_pair, 0),
    sum_coinf   = prop_triple + prop_HIV + prop_Syph,
    scale_down  = if_else(sum_coinf > 1, 1 / sum_coinf, 1),
    prop_triple = prop_triple * scale_down,
    prop_HIV    = prop_HIV    * scale_down,
    prop_Syph   = prop_Syph   * scale_down,
    prop_none   = 1 - (prop_triple + prop_HIV + prop_Syph)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    yrs_HBVonly    = person_years * prop_none,
    yrs_HBV_HIV    = person_years * prop_HIV,
    yrs_HBV_Syph   = person_years * prop_Syph,
    yrs_HBV_Triple = person_years * prop_triple
  ) %>%
  dplyr::select(arm, time, state, yrs_HBVonly, yrs_HBV_HIV, yrs_HBV_Syph, yrs_HBV_Triple)

# ========== ADULT SPLIT (fixed) ==========
# ========== ADULT SPLIT (uses full adult denominator) ==========
adult_df <- final_df %>%
  dplyr::filter(state %in% hbv_infected_adult_split) %>%
  dplyr::left_join(
    hiv_syph_raw %>%
      dplyr::select(age, arm, hiv_infected_adults, syphilis_infected_adults, coinfected_adults) %>%
      dplyr::rename(time = age),
    by = c("arm", "time")
  ) %>%
  dplyr::left_join(pool_adult, by = c("arm", "time")) %>%
  dplyr::group_by(arm, time) %>%
  dplyr::mutate(
    p_HIV   = if_else(denom_adult > 0, hiv_infected_adults      / denom_adult, 0),
    p_Syph  = if_else(denom_adult > 0, syphilis_infected_adults / denom_adult, 0),
    p_pair0 = if_else(denom_adult > 0, coinfected_adults        / denom_adult, 0),
    p_pair  = pmin(pmax(p_pair0, p_HIV + p_Syph - 1, 0), pmin(p_HIV, p_Syph)),
    prop_triple = p_pair,
    prop_HIV    = pmax(p_HIV  - p_pair, 0),
    prop_Syph   = pmax(p_Syph - p_pair, 0),
    sum_coinf   = prop_triple + prop_HIV + prop_Syph,
    scale_down  = if_else(sum_coinf > 1, 1 / sum_coinf, 1),
    prop_triple = prop_triple * scale_down,
    prop_HIV    = prop_HIV    * scale_down,
    prop_Syph   = prop_Syph   * scale_down,
    prop_none   = 1 - (prop_triple + prop_HIV + prop_Syph)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    yrs_HBVonly    = person_years * prop_none,
    yrs_HBV_HIV    = person_years * prop_HIV,
    yrs_HBV_Syph   = person_years * prop_Syph,
    yrs_HBV_Triple = person_years * prop_triple
  ) %>%
  dplyr::select(arm, time, state, yrs_HBVonly, yrs_HBV_HIV, yrs_HBV_Syph, yrs_HBV_Triple)


# Combine infant + adult splits
hbv_split <- bind_rows(infant_df, adult_df)

# Overlap person-years that carry HIV and/or Syph INSIDE the HBV-infected pools
overlap_by_year <- hbv_split %>%
  mutate(group = case_when(
    state %in% hbv_infected_adult_split  ~ "adult",
    state %in% hbv_infected_infant_split ~ "infant",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  group_by(arm, time, group) %>%
  summarise(
    hiv_with_hbv    = sum(yrs_HBV_HIV + yrs_HBV_Triple, na.rm = TRUE),
    syph_with_hbv   = sum(yrs_HBV_Syph + yrs_HBV_Triple, na.rm = TRUE),
    triple_with_hbv = sum(yrs_HBV_Triple,               na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = c(hiv_with_hbv, syph_with_hbv, triple_with_hbv),
    values_fill = 0
  )


# Pivot long format
hbv_long <- hbv_split %>%
  pivot_longer(cols = starts_with("yrs_"), names_to = "coinf_type", values_to = "person_years") %>%
  mutate(
    state = case_when(
      coinf_type == "yrs_HBVonly" ~ state,
      coinf_type == "yrs_HBV_HIV" ~ paste0(state, "_HIV"),
      coinf_type == "yrs_HBV_Syph" ~ paste0(state, "_Syph"),
      coinf_type == "yrs_HBV_Triple" ~ paste0(state, "_Triple")
    )
  ) %>%
  select(arm, time, state, person_years)

# since hiv adult dw weight is weighted
hbv_long <- hbv_long %>%
  mutate(
    symbol = case_when(
      state == "hiv_infected_adults" & arm == "SOC" ~ "hiv_dw_adt_soc",
      state == "hiv_infected_adults" & arm == "INT" ~ "hiv_dw_adt_int",
      TRUE ~ state  # use default state for other rows
    )
  )

final_df <- hbv_long %>%
  left_join(dw_lookup, by = "state") %>%
  mutate(
    adjusted_years = if_else(state %in% acute_states, person_years * acute_duration, person_years),
    yld = adjusted_years * dw,
    discount_factor = discount_factors[time + 1],
    yld_discounted = yld * discount_factor
  )


# Remove original HBV states
final_df <- final_df %>%
  filter(!(state %in% c(adult_hbv_states, infant_hbv_states)))

# Append coinfection-adjusted HBV states
final_df_x <- bind_rows(final_df, hbv_long)

# Join DWs and calculate YLDs
final_df <- hbv_long %>% ## CHECK
  left_join(dw_lookup, by = "state") %>%
  mutate(
    adjusted_years = if_else(state %in% acute_states, person_years * acute_duration, person_years),
    yld = adjusted_years * dw
  ) %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yld_discounted = yld * discount_factor
  )

## Syphilis related adjustments to disability weight
## Time = 3 onwards, lbw dw is 0, a
## Time = 5 onwards, in kids, all syphilis related weights are 0 (so wil ladjust coinfection weights too)
## Time = 1 onwards, in adults, all syphilis related weights are 0 (so wil ladjust coinfection weights too)
# Apply time- and state-dependent adjustments to DW and recalculate YLDs
final_df <- final_df %>%
  mutate(
    dw = case_when(
      # Syphilis-infected INFANTS: from time 5 onward, use DW from the base HBV state (strip "_Syph")
      grepl("_Syph$", state) & grepl("^Infant_", state) & time >= 5 ~ dw[match(gsub("_Syph$", "", state), state)],
      
      # Syphilis-infected ADULTS: from time 1 onward, use DW from the base HBV state
      grepl("_Syph$", state) & !grepl("^Infant_", state) & time >= 1 ~ dw[match(gsub("_Syph$", "", state), state)],
      
      # Triple-infected INFANTS: from time 5 onward, use DW from corresponding HIV state (replace "_Triple" with "_HIV")
      grepl("_Triple$", state) & grepl("^Infant_", state) & time >= 5 ~ dw[match(gsub("_Triple$", "_HIV", state), state)],
      
      # Triple-infected ADULTS: from time 1 onward, use DW from corresponding HIV state
      grepl("_Triple$", state) & !grepl("^Infant_", state) & time >= 1 ~ dw[match(gsub("_Triple$", "_HIV", state), state)],
      
      TRUE ~ dw
    )
  ) %>%
  mutate(
    yld = adjusted_years * dw,
    yld_discounted = yld * discount_factor
  )


## Account for lbw disability only for ages 1 and 2
lbw_df <- hiv_syph_raw %>%
  select(time = age, person_years = lbw_infants, arm) %>%
  mutate(
    state = "lbw_infants",
    dw = if_else(time >= 3, 0, param_vals["syph_lbw_dw_inf"]),
    adjusted_years = person_years,
    yld = adjusted_years * dw,
    discount_factor = discount_factors[time + 1],
    yld_discounted = yld * discount_factor
  )



final_df <- bind_rows(final_df, lbw_df)

# Adjust external series so totals = (external-only remainder) + (overlap sitting in HBV pools)
ext_adj <- hiv_syph_raw %>%
  dplyr::rename(time = age) %>%
  dplyr::left_join(overlap_by_year, by = c("arm","time")) %>%
  dplyr::mutate(
    # adults
    hiv_infected_adults       = pmax(0, hiv_infected_adults       - coalesce(hiv_with_hbv_adult,    0)),
    syphilis_infected_adults  = pmax(0, syphilis_infected_adults  - coalesce(syph_with_hbv_adult,   0)),
    coinfected_adults         = pmax(0, coinfected_adults         - coalesce(triple_with_hbv_adult, 0)),
    # infants
    hiv_infected_infants      = pmax(0, hiv_infected_infants      - coalesce(hiv_with_hbv_infant,    0)),
    syphilis_infected_infants = pmax(0, syphilis_infected_infants - coalesce(syph_with_hbv_infant,   0)),
    coinfected_infants        = pmax(0, coinfected_infants        - coalesce(triple_with_hbv_infant, 0))
  )

# Rebuild the external-only pools for YLD/costs
extra_states_df <- ext_adj %>%
  dplyr::select(
    time, arm,
    hiv_infected_adults,
    syphilis_infected_adults,
    coinfected_adults,
    hiv_infected_infants,
    syphilis_infected_infants,
    coinfected_infants
  ) %>%
  tidyr::pivot_longer(
    cols = -c(time, arm),
    names_to  = "state",
    values_to = "person_years"
  ) %>%
  dplyr::mutate(
    dw = case_when(
      # HIV-only adults vary by arm
      state == "hiv_infected_adults" & arm == "SOC" ~ param_vals["hiv_dw_adt_soc"],
      state == "hiv_infected_adults" & arm == "INT" ~ param_vals["hiv_dw_adt_int"],
      # HIV-only infants
      state == "hiv_infected_infants" ~ param_vals["hiv_dw_inf"],
      # Syph-only rules (as you had)
      state == "syphilis_infected_adults"  ~ if_else(time >= 1, 0, param_vals["syph_dw_inf"]),
      state == "syphilis_infected_infants" ~ if_else(time >= 5, 0, param_vals["syph_dw_inf"]),
      # HIV+Syph w/o HBV (the remainder after subtracting HBV triple)
      state == "coinfected_adults"  & arm == "SOC" & time >= 1 ~ param_vals["hiv_dw_adt_soc"],
      state == "coinfected_adults"  & arm == "INT" & time >= 1 ~ param_vals["hiv_dw_adt_int"],
      state == "coinfected_adults"  & time < 1                 ~ param_vals["hiv_syph_dw_inf"],
      state == "coinfected_infants" & time >= 5                ~ param_vals["hiv_dw_inf"],
      state == "coinfected_infants" & time < 5                 ~ param_vals["hiv_syph_dw_inf"],
      TRUE ~ 0
    ),
    adjusted_years   = person_years,
    yld              = adjusted_years * dw,
    discount_factor  = discount_factors[time + 1],
    yld_discounted   = yld * discount_factor
  )

# Then append as before:
final_df <- bind_rows(final_df, extra_states_df)

final_df <- final_df %>%
  select(-symbol.x, -symbol.y)

# === Add cost columns directly to final_df ===

# choose cost symbol given state & time with the remapping rules
choose_cost_symbol <- function(st, tt) {
  is_infant <- grepl("(?i)infant", st)           # infant anywhere in the name
  thr <- ifelse(is_infant, 5, 2)                 # infants >=5, adults >=2
  sym <- paste0("cost_", st)                     # default
  
  if (tt >= thr) {
    if (grepl("_Syph$", st)) {
      sym <- paste0("cost_", sub("_Syph$", "", st))                # drop Syph -> base
    } else if (grepl("_Triple$", st)) {
      sym <- paste0("cost_", sub("_Triple$", "_HIV", st))          # Triple -> HIV
    } else if (st %in% c("coinfected_infants", "coinfected_adults")) {
      sym <- if (is_infant) "cost_hiv_infected_infants" else "cost_hiv_infected_adults"
    }
  }
  sym
}

final_df <- final_df %>%
  dplyr::mutate(
    # symbol to look up in params
    cost_symbol = mapply(choose_cost_symbol, state, time),
    
    # unit cost from params (fall back to 0 if missing)
    unit_cost = as.numeric(param_vals[cost_symbol]),
    unit_cost = dplyr::coalesce(unit_cost, 0),
    
    # pure syph-only pools get 0 cost after the threshold (adults>=2, infants>=5)
    is_infant = grepl("(?i)infant", state),
    thr = dplyr::if_else(is_infant, 5L, 2L),
    unit_cost = dplyr::if_else(time >= thr & grepl("^syphilis_", state), 0, unit_cost),
    
    # totals
    total_cost = unit_cost * person_years,
    total_cost_discounted = total_cost * discount_factor
  ) %>%
  dplyr::select(-is_infant, -thr)



# ---------- PROGRAM COSTS (year 0; sensitivity & specificity matter here) ----------
# pull a year-0 discount factor from final_df (it’s already there)
df0 <- unique(final_df$discount_factor[final_df$time == 0])[1]
if (is.na(df0)) df0 <- 1

# helper: get sens/spec by arm (default 1 if missing)
get_sens <- function(arm) {
  s_trip <- ifelse(is.na(param_vals["triplex_sens"]), 1, param_vals["triplex_sens"])
  s_hbs  <- ifelse(is.na(param_vals["hbsag_sens"]),   1, param_vals["hbsag_sens"])
  if (arm == "INT") s_trip else s_hbs
}
get_spec <- function(arm) {
  c_trip <- ifelse(is.na(param_vals["triplex_spec"]), 1, param_vals["triplex_spec"])
  c_hbs  <- ifelse(is.na(param_vals["hbsag_spec"]),   1, param_vals["hbsag_spec"])
  if (arm == "INT") c_trip else c_hbs
}

prev_inf <- param_vals["hbv_p_preg_chronic"] + param_vals["hbv_p_preg_acute"]
p_known  <- param_vals["hbv_p_known_hbsag"]

prog_cost_list <- lapply(names(prop_tested_list), function(arm_name) {
  prop_tested <- prop_tested_list[[arm_name]]
  sens_arm <- get_sens(arm_name)
  spec_arm <- get_spec(arm_name)
  
  # testing volumes
  n_known   <- N0 * p_known
  n_tested  <- N0 * prop_tested
  n_inf     <- N0 * prev_inf
  
  n_tested_inf        <- n_tested * prev_inf
  n_tested_uninf      <- n_tested - n_tested_inf
  
  tp <- n_tested_inf * sens_arm
  fn <- n_tested_inf - tp
  fp <- n_tested_uninf * (1 - spec_arm)
  tn <- n_tested_uninf - fp
  
  # Managed (treated as positive) = known + TP + FP  (FP gets costs but no effects)
  managed_pos <- n_known + tp + fp
  
  # component costs
  test_cost <- if (arm_name == "INT") param_vals["triplex_cost"] else param_vals["hbv_c_test_hbsag"]
  cost_tests <- n_tested * test_cost
  
  # downstream products for anyone managed positive
  p_tdf   <- param_vals["hbv_p_tenofovir_preg"]
  p_hbig  <- param_vals["hbv_p_infant_hbig"]
  p_vaxb  <- param_vals["hbv_p_infant_vaccine_birth"]
  
  cost_tdf  <- managed_pos * p_tdf  * param_vals["hbv_c_tenofovir"]
  cost_hbig <- managed_pos * p_hbig * param_vals["hbv_c_hbig"]
  cost_vaxb <- managed_pos * p_vaxb * param_vals["hbv_c_vaccine_birth"]
  
  program_cost_y0 <- (cost_tests + cost_tdf + cost_hbig + cost_vaxb) * df0
  
  tibble(
    arm = arm_name,
    time = 0:time_horizon,
    program_costs_discounted = if_else(0:time_horizon == 0, program_cost_y0, 0)
  )
})

program_costs_by_year <- bind_rows(prog_cost_list)





# other costs

# ---------- EXTERNAL DISCOUNTED COSTS FROM hiv_syph_raw (keep separate) ----------
extra_costs_by_year <- hiv_syph_raw %>%
  rename(time = age) %>%
  group_by(arm, time) %>%
  summarise(
    extra_costs_discounted = sum(coalesce(disc_costs_infants, 0) + coalesce(disc_costs_adults, 0),
                                 na.rm = TRUE),
    .groups = "drop"
  )


# Quick getter for cost symbols (returns 0 if missing)
get_cost <- function(sym, default = 0) {
  v <- suppressWarnings(param_vals[[sym]])
  if (is.na(v)) default else v
}

# Pull costs from params (your symbols)
cost_triplex   <- get_cost("triplex_cost")                 # 2.9
cost_hbsag     <- get_cost("hbv_c_test_hbsag")             # 0.87
cost_tdf       <- get_cost("hbv_c_tenofovir")              # 32.4
cost_hbig      <- get_cost("hbv_c_hbig")                   # 45
cost_vax_birth <- get_cost("hbv_c_vaccine_birth")          # 1.35
cost_vax_full  <- get_cost("hbv_c_vaccine_full")           # 1.35 (unused here unless you add coverage)

# Parameters already in your model
prop_known      <- param_vals[["hbv_p_known_hbsag"]]
prop_tenofovir  <- param_vals[["hbv_p_tenofovir_preg"]]
p_infant_hbig   <- param_vals[["hbv_p_infant_hbig"]]
p_infant_vax    <- param_vals[["hbv_p_infant_vaccine_birth"]]
p_preg_chronic  <- param_vals[["hbv_p_preg_chronic"]]
p_preg_acute    <- param_vals[["hbv_p_preg_acute"]]


# Pull discount factors directly from final_df
df_by_time <- final_df %>%
  dplyr::distinct(time, discount_factor) %>%
  dplyr::arrange(time)

# Year-0 discount factor (should be 1.0, but make it robust)
df0 <- df_by_time %>%
  dplyr::filter(time == 0) %>%
  dplyr::pull(discount_factor)

if (length(df0) == 0 || is.na(df0)) df0 <- 1
# Your model's state-level discounted costs per year
state_costs_by_year <- final_df %>%
  group_by(arm, time) %>%
  summarise(state_costs_discounted = sum(total_cost_discounted, na.rm = TRUE),
            .groups = "drop")

# Combine state + external + program costs
costs_by_year <- state_costs_by_year %>%
  left_join(extra_costs_by_year,  by = c("arm","time")) %>%
  mutate(extra_costs_discounted = replace_na(extra_costs_discounted, 0)) %>%
  left_join(program_costs_by_year, by = c("arm","time")) %>%
  mutate(program_costs_discounted = replace_na(program_costs_discounted, 0),
         total_costs_discounted   = state_costs_discounted +
           extra_costs_discounted +
           program_costs_discounted) %>%
  arrange(arm, time)







# Summarize total YLDs by arm
total_ylds <- final_df %>%
  group_by(arm) %>%
  summarise(total_yld = sum(yld_discounted, na.rm = TRUE))


#YLL
#  Read life table and assign Kenya
life_exp_df <- read_csv("life_exp.csv") %>%
  rename(
    age_group = `Age Group`,
    life_exp_both = `Both sexes`,
    life_exp_female = Female
  ) %>%
  mutate(
    country = "Kenya",  # Force Kenya
    age_start = as.numeric(str_extract(age_group, "^\\d+")),
    age_start = if_else(str_detect(age_group, "<1"), 0, age_start)
  )

# Define helper to assign age group
get_age_group_start <- function(age) {
  if (age < 1) return(0)
  possible_starts <- sort(unique(na.omit(life_exp_df$age_start)))
  max(possible_starts[possible_starts <= age])
}

death_trace_adult <- final_results %>%
  filter(state == "Death") %>%
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_death = count - lag(count, default = 0)) %>%
  ungroup() %>%
  filter(new_death > 0) %>%
  mutate(
    country = "Kenya",
    age_at_death = 22 + time,
    age_bin = vapply(age_at_death, get_age_group_start, numeric(1))
    
  ) %>%
  left_join(life_exp_df, by = c("age_bin" = "age_start", "country")) %>%
  mutate(
    le_used = life_exp_female,
    yll = new_death * le_used
  )

death_trace_infant <- final_results %>%
  filter(state == "Infant_death", time > 0) %>%  # filter out t = 0
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_death = count - lag(count, default = 0)) %>%
  ungroup() %>%
  filter(new_death > 0) %>%
  mutate(
    country = "Kenya",
    age_at_death = time - 1,  # infants born in t = 1, die at end of same cycle
    age_bin = vapply(age_at_death, get_age_group_start, numeric(1))
  ) %>%
  left_join(life_exp_df, by = c("age_bin" = "age_start", "country")) %>%
  mutate(
    le_used = life_exp_both,
    yll = new_death * le_used
  )

# HIV infant deaths
death_trace_hiv_infant <- final_results %>%
  filter(state == "HIV_infant_death", time > 0) %>%
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_death = count - lag(count, default = 0)) %>%
  ungroup() %>%
  filter(new_death > 0) %>%
  mutate(
    country = "Kenya",
    age_at_death = time - 1,  # consistent with infant death timing
    age_bin = vapply(age_at_death, get_age_group_start, numeric(1))
  ) %>%
  left_join(life_exp_df, by = c("age_bin" = "age_start", "country")) %>%
  mutate(
    le_used = life_exp_both,
    yll = new_death * le_used
  )

# Syphilis infant deaths
death_trace_syph_infant <- final_results %>%
  filter(state == "Syphilis_infant_death", time > 0) %>%
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_death = count - lag(count, default = 0)) %>%
  ungroup() %>%
  filter(new_death > 0) %>%
  mutate(
    country = "Kenya",
    age_at_death = time - 1,
    age_bin = vapply(age_at_death, get_age_group_start, numeric(1))
    
  ) %>%
  left_join(life_exp_df, by = c("age_bin" = "age_start", "country")) %>%
  mutate(
    le_used = life_exp_both,
    yll = new_death * le_used
  )

# HIV adult deaths
death_trace_hiv_adult <- final_results %>%
  filter(state == "HIV_adult_death", time > 0) %>%
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_death = count - lag(count, default = 0)) %>%
  ungroup() %>%
  filter(new_death > 0) %>%
  mutate(
    country = "Kenya",
    age_at_death = 22 + time,
    age_bin = vapply(age_at_death, get_age_group_start, numeric(1))
  ) %>%
  left_join(life_exp_df, by = c("age_bin" = "age_start", "country")) %>%
  mutate(
    le_used = life_exp_female,
    yll = new_death * le_used
  )

# Syphilis adult deaths
death_trace_syph_adult <- final_results %>%
  filter(state == "Syphilis_adult_death", time > 0) %>%
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_death = count - lag(count, default = 0)) %>%
  ungroup() %>%
  filter(new_death > 0) %>%
  mutate(
    country = "Kenya",
    age_at_death = 22 + time,
    age_bin = vapply(age_at_death, get_age_group_start, numeric(1))
    
  ) %>%
  left_join(life_exp_df, by = c("age_bin" = "age_start", "country")) %>%
  mutate(
    le_used = life_exp_female,
    yll = new_death * le_used
  )



#discount
death_trace_adult <- death_trace_adult %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yll_discounted = yll * discount_factor
  )

death_trace_infant <- death_trace_infant %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yll_discounted = yll * discount_factor
  )

death_trace_hiv_infant <- death_trace_hiv_infant %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yll_discounted = yll * discount_factor
  )

death_trace_syph_infant <- death_trace_syph_infant %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yll_discounted = yll * discount_factor
  )

death_trace_hiv_adult <- death_trace_hiv_adult %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yll_discounted = yll * discount_factor
  )

death_trace_syph_adult <- death_trace_syph_adult %>%
  mutate(
    discount_factor = discount_factors[time + 1],
    yll_discounted = yll * discount_factor
  )

total_yll <- bind_rows(
  death_trace_adult,
  death_trace_infant,
  death_trace_hiv_infant,
  death_trace_syph_infant,
  death_trace_hiv_adult,
  death_trace_syph_adult
) %>%
  group_by(arm) %>%
  summarise(total_yll = sum(yll_discounted, na.rm = TRUE), .groups = "drop")



#  Merge with YLDs to compute DALYs
dalys_df <- final_df %>%
  group_by(arm) %>%
  summarise(total_yld = sum(yld_discounted, na.rm = TRUE), .groups = "drop") %>%
  left_join(total_yll, by = "arm") %>%
  mutate(total_daly = total_yld + total_yll)


yld_by_year <- final_df %>%
  group_by(arm, time, state) %>%
  summarise(yld_discounted = sum(yld_discounted, na.rm = TRUE), .groups = "drop") %>%
  arrange(arm, time, state)

yld_by_year_total <- yld_by_year %>%
  group_by(arm, time) %>%
  summarise(total_yld = sum(yld_discounted, na.rm = TRUE), .groups = "drop")

yll_by_year <- bind_rows(death_trace_adult, death_trace_infant, death_trace_hiv_infant, death_trace_syph_infant) %>%
  group_by(arm, time) %>%
  summarise(total_yll = sum(yll_discounted, na.rm = TRUE), .groups = "drop") %>%
  arrange(arm, time)

dalys_by_year <- full_join(yld_by_year_total, yll_by_year, by = c("arm", "time")) %>%
  mutate(
    total_yld = replace_na(total_yld, 0),
    total_yll = replace_na(total_yll, 0),
    total_dalys = total_yld + total_yll
  ) %>%
  arrange(arm, time)

discount_lookup <- tibble(
  time = 0:time_horizon,
  discount_factor = 1 / ((1 + discount_rate) ^ (0:time_horizon))
)


# Step 2: Get cumulative DALYs per arm
total_dalys <- dalys_by_year %>%
  group_by(arm) %>%
  summarise(total_daly = sum(total_dalys, na.rm = TRUE), .groups = "drop")



# Yearly ICER base table
icer_df_time <- dalys_by_year %>%           # has total_yld, total_yll, total_dalys by arm/time
  left_join(costs_by_year, by = c("arm","time")) %>%
  arrange(arm, time)


# Cumulative 20y totals by arm
icer_df <- icer_df_time %>%
  group_by(arm) %>%
  summarise(
    total_costs_discounted = sum(total_costs_discounted, na.rm = TRUE),
    total_dalys            = sum(total_dalys,            na.rm = TRUE),
    .groups = "drop"
  )

# INT vs SOC ICER
icer_summary <- icer_df %>%
  select(arm, total_costs_discounted, total_dalys) %>%
  tidyr::pivot_wider(names_from = arm,
                     values_from = c(total_costs_discounted, total_dalys)) %>%
  mutate(
    delta_cost    = total_costs_discounted_INT - total_costs_discounted_SOC,
    dalys_averted = total_dalys_SOC - total_dalys_INT,
    icer          = ifelse(dalys_averted > 0, delta_cost / dalys_averted, NA_real_)
  )


























## Optional calibration check

# --- define which single-infection states we want to keep ---
hbv_single <- c(
  "Acute_symptomatic_HBV","Chronic_HBV","Compensated_cirrhosis",
  "Decompensated_cirrhosis","Liver_cancer",
  "HBV_positive_infant","Infant_chronic_HBV","Infant_compensated_cirrhosis",
  "Infant_decompensated_cirrhosis","Infant_liver_cancer"
)

hiv_only   <- c("hiv_infected_adults","hiv_infected_infants")
syph_only  <- c("syphilis_infected_adults","syphilis_infected_infants")

drop_states <- c(
  "Death","Infant_death","All_infant_deaths","All_adult_deaths",
  "HIV_adult_death","Syphilis_adult_death","HIV_infant_death","Syphilis_infant_death",
  "HBV_susceptible","HBV_immune","HBV_negative_infant","Infant_immune"
)

# --- build calibration_check (long format) ---
calibration_check <- final_df %>%
  # drop any co-infected variants and other non-target states
  dplyr::filter(
    !grepl("(_HIV|_Syph|_Triple)$", state),
    !grepl("^coinfected_", state),
    !state %in% c("Adult_HIV_HBV_coinf","Adult_Syphilis_HBV_coinf","Adult_HIV_Syphilis_HBV_coinf",
                  "Infant_HIV_HBV_coinf","Infant_Syphilis_HBV_coinf","Infant_HIV_Syphilis_HBV_coinf",
                  drop_states)
  ) %>%
  # keep only the single-infection sets we care about
  dplyr::filter(state %in% c(hbv_single, hiv_only, syph_only)) %>%
  dplyr::group_by(arm, time, state) %>%
  dplyr::summarise(person_years = sum(person_years, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(arm, time, state)

# (optional) wide version for quick eyeballing vs IHME tables
calibration_check_wide <- calibration_check %>%
  tidyr::pivot_wider(names_from = state, values_from = person_years, values_fill = 0) %>%
  dplyr::arrange(arm, time)























# Calibration check 

# Define the states of interest
adult_focus_states <- c("Acute_symptomatic_HBV", "Chronic_HBV",
                        "Compensated_cirrhosis", "Decompensated_cirrhosis",
                        "Liver_cancer")

infant_focus_states <- c("HBV_positive_infant", "Infant_chronic_HBV",
                         "Infant_compensated_cirrhosis", "Infant_decompensated_cirrhosis",
                         "Infant_liver_cancer")

focus_states <- c(adult_focus_states, infant_focus_states)

# Filter for SOC arm, years 0, 5, 20, and focus states
soc_export <- final_results %>%
  filter(arm == "SOC",
         time %in% c(0, 1, 5, 20),
         state %in% focus_states)

# Optional: make it wider by year for cleaner view
soc_export_wide <- soc_export %>%
  pivot_wider(names_from = time, values_from = count)

# Save to CSV
#write_csv(soc_export, "soc_states_0_5_20_long.csv")
write_csv(soc_export_wide, "soc_states_0_5_20_wide.csv")
