
#### Title: HIV retesting model parameters
#### Purpose: Input parameters for HIV retesting model
#### Author: David Coomes
#### Date: January 9, 2025


# Kenya HIV retesting parameters

# Sensitivity parameters
#prev_reduction <- 0
#PrEP <- 0
#incidence_reduction <- 0

test_gen <- 3
disc_rate <- 0.03

# Miscellaneous
transition_ag <- 2.3
transition_ab <- 0.7
transition_r <- 6

# Calculated
duration_infection <- sum(transition_ag, transition_ab, transition_r)


# Test sensitivity and specificity
sens_Ag_neg <- 0
transition_ag <- 2.3
transition_ab <- 0.7
transition_r <- 6
sens_a3 <- 0
sens_a4 <- 0.5615
sens_c3 <- 1
sens_c4 <- 0.988
sens_e3 <- ((sens_Ag_neg*transition_ag) + (sens_a3*transition_ab) + (sens_c3*transition_r)) / duration_infection
sens_e4 <- ((sens_Ag_neg*transition_ag) + (sens_a4*transition_ab) + (sens_c4*transition_r)) / duration_infection
spec3 <- 0.9893
spec4 <- 0.991


# set the test sensitivity by selecting e3 or e4
test_sens_e <- ifelse(test_gen==4, sens_e4, sens_e3)
test_sens_c <- ifelse(test_gen==4, sens_c4, sens_c3)
test_spec <- ifelse(test_gen==4, spec4, spec3)

# Infant transmission parameters 
inutero_recent <- 1-exp(-0.1495/9)
inutero_established <- 1-exp(-0.0746/40)
early_pp_established <- 1-exp(-0.16/7)
early_pp_recent <- early_pp_established*1.3
mid_pp_recent <- 1-exp(-0.358/(16.5*4))
late_pp_recent <- 1-exp(-0.358/(16.5*4))
mid_pp_established <- 1-exp(-0.04/20)
late_pp_established <- 0.0005
viral_load_suppression <- 0.95




# maternal_input_func <- function(population) {

if (country == "Kenya") {

inf_ART_coverage <- 0.61          # Can change this 
inf_ART_coverage_base <- 0.61     # Don't change this
adult_ART_coverage <- 0.69 



red_prep <- 0.7049
avg_ART_adherence <- 0.84





# Testing probabilities
att_firstANC <- 0.96
att_lategest <- 1-0.033-0.04
att_delivery <- 0.618
att_6wk <- 0.96
att_14wk <- 0.88
att_6mo <- 0.865
att_9mo <- 0.85

test_accept <- 340/405
stockout <- 0.05
p_results <- 0.978




# ART parameters
p_ART <- 0.91
ART_ret <- 0.738
ART_drop_factor <- exp(log(ART_ret)/91)
p_VL <- 0.881

# Weekly mortality rates
mort_adult_females <- (0.001+0.002+0.003+0.003+0.004+0.005)/6/52
mort_maternal_perinatal <- (342/100000)/6




## Input parameters
# Gestational age
first_ANC <- 22
second_ANC <- 33
delivery <- 39
# infant weeks at MCH visit
early_pp_visit <- 6
pp_14weeks_visit <- 14
pp_mid_visit <- 26
pp_9months_visit <- 39
end_of_model <- 53


# Testing parameters (testing offered at indicated visit)
# We can set these up in the model output to run multiple models at a time
# first_ANC_test <- 1
# late_ANC_test <- 1
# pp_early_test <- 1
# pp_14weeks_test <- 1
# pp_mid_test <- 1
# pp_9months_test <- 1


# Testing probabilities
second_ANC_test_prob <- late_ANC_test * att_lategest * test_accept * (1-stockout)
delivery_test_prob <- (1-second_ANC_test_prob) * late_ANC_test * att_delivery * test_accept * (1-stockout)



# Costs
cpi_adj <- 1.23196

cost_rapid_screening_3gen <- 2.64*cpi_adj
cost_true_pos_screening <- 3.68*cpi_adj
cost_false_pos_screening <- 26.39*cpi_adj
cost_maternal_ART <- (234/52)*(112.317/104.047)*cpi_adj
cost_infant_prophylaxis <- (2.19)*(112.317/105.873)*cpi_adj
cost_maternal_PrEP <- 26.52/(30/7)*cpi_adj
cost_infant_ART_first_2weeks <- (349.92/52)*cpi_adj
cost_infant_ART_after_2weeks <- (349.92/52)*cpi_adj

test_cost <- ifelse(test_gen == 3, cost_rapid_screening_3gen, cost_rapid_screening_4gen)


# Infant incidence reduction parameters
p_arv <- 0.57*(13/188)+0.963*(175/188)
arv_red <- 0.675
nbf <- 1
prob_nbf_neg_early_pp <- 0.001
prob_nbf_pos_early_pp <- 0.025
prob_nbf_neg_mid_pp <- 0.0058
prob_nbf_pos_mid_pp <- 0.212
prob_nbf_neg_late_pp <- 0.009
prob_nbf_pos_late_pp <- 0.334


# Infant mortality parameters
mort_n <- 1-exp(-7*(7*79871.59+21*7534.39)/(100000*28*365))
mort_i <- 1-exp(-7*1609.69/100000/365)

# DALYs
dw_HIVtx <- 0.078
dw_HIVnotx <- 0.274
dw_AIDSnotx <- 0.582
propLE_HIV_notx <- 0.9





#############
## Population specific parameters ##
#############

if (population == "General") {
  # Population parameters
  population_size <- 1631470
  
  ## prevalence and incidence rates
  hiv_prev <- 0.061*(1-prev_reduction)
  status_known <- 0.57
  
  red_prep <- 0.7049
  
  # Weekly incidence parameters
  inc_preg_early <- (((3.44/100)/52)/2)*(1-incidence_reduction)
  inc_preg_late <- (((3.44/100)/52)/2)*(1-incidence_reduction)
  inc_pp_early <- 0
  inc_pp_mid <- ((1.4/100)/52)*(1-incidence_reduction)
  inc_pp_late <- ((1.4/100)/52)*(1-incidence_reduction)
  
} else if (population == "FSW") {
  
  # Population parameters
  population_size <- 19348
  
  ## prevalence and incidence rates
  hiv_prev <- 0.293*(1-prev_reduction)
  status_known <- 0.68
  
  # Weekly incidence parameters
  inc_preg_early <- ((4.7/100)/52)*(1-incidence_reduction)
  inc_preg_late <- ((4.7/100)/52)*(1-incidence_reduction)
  inc_pp_early <- 0
  inc_pp_mid <- ((4.7/100)/52)*(1-incidence_reduction)
  inc_pp_late <- ((4.7/100)/52)*(1-incidence_reduction)
  
  # ART parameters
  p_ART <- 0.73
  ART_ret <- 0.6
  
} else if (population == "PWID") {
  
  # Population parameters
  population_size <- 201
  
  ## prevalence and incidence rates
  hiv_prev <- 0.183*(1-prev_reduction)
  status_known <- 0.94
  
  # Weekly incidence parameters
  inc_preg_early <- 0.00109615384615385*(1-incidence_reduction)
  inc_preg_late <- 0.00109615384615385*(1-incidence_reduction)
  inc_pp_early <- 0
  inc_pp_mid <- 0.00109615384615385*(1-incidence_reduction)
  inc_pp_late <- 0.00109615384615385*(1-incidence_reduction)
  
  # ART parameters
  ART_ret <- 0.6
  
} else if (population == "serodiscordant") {
  
  # Population parameters
  population_size <- 39180
  
  ## prevalence and incidence rates
  hiv_prev <- 0*(1-prev_reduction)
  status_known <- 0.96
  
  # Weekly incidence parameters
  inc_preg_early <- 0.000721153846153846*(1-incidence_reduction)
  inc_preg_late <- 0.00135*(1-incidence_reduction)
  inc_pp_early <- 0
  inc_pp_mid <- 0.0009*(1-incidence_reduction)
  inc_pp_late <- 0.000240384615384615*(1-incidence_reduction)
  
  # ART parameters
  ART_ret <- 0.6
  ART_drop_factor <- exp(log(ART_ret)/91)
  

} else {
  
  print("No population selected")

}


##### 
# South Africa parameters
#####

} else if (country=="South Africa") {
  
  inf_ART_coverage <- 0.63          # Can change this 
  inf_ART_coverage_base <- 0.63     # Don't change this
  adult_ART_coverage <- 0.62 

  red_prep <- 0.7049
  avg_ART_adherence <- 0.84
  
  # Miscellaneous
  transition_ag <- 2.3
  transition_ab <- 0.7
  transition_r <- 6
  
  # Calculated
  duration_infection <- sum(transition_ag, transition_ab, transition_r)
  
  
  # Testing probabilities
  att_firstANC <- 0.94
  att_lategest <- 0.78
  att_delivery <- 0.96
  att_6wk <- 0.902
  att_14wk <- 0.727
  att_6mo <- 0.86
  att_9mo <- 0.619
  
  test_accept <- 0.98
  stockout <- 0.05
  p_results <- 81.8/(81.8+1.7)

  
  # ART parameters
  p_ART <- 0.87
  ART_ret <- 0.6
  ART_drop_factor <- exp(log(ART_ret)/91)
  p_VL <- 0.72
  
  # Weekly mortality rates
  mort_adult_females <- (0.001+0.002+0.003+0.005+0.008+0.007+0.008)/7/52
  mort_maternal_perinatal <- (117/100000)/6
  
  
  
  
  ## Input parameters
  # Gestational age
  first_ANC <- 18
  second_ANC <- 36
  delivery <- 39
  # infant weeks at MCH visit
  early_pp_visit <- 6
  pp_14weeks_visit <- 14
  pp_mid_visit <- 26
  pp_9months_visit <- 39
  end_of_model <- 53
  
  
  # Testing probabilities
  second_ANC_test_prob <- late_ANC_test * att_lategest * test_accept * (1-stockout)
  delivery_test_prob <- (1-second_ANC_test_prob) * late_ANC_test * att_delivery * test_accept * (1-stockout)
  

  # Costs
  cpi_adj <- 1.23196
  
  cost_rapid_screening_3gen <- 7.72*cpi_adj
  cost_true_pos_screening <- 11.39*cpi_adj
  cost_false_pos_screening <- 34.17*cpi_adj
  cost_maternal_ART <- (249.15/52)*cpi_adj
  cost_infant_prophylaxis <- (3.6*(112.317/105.873))*cpi_adj
  cost_maternal_PrEP <- (26.52/(30/7))*cpi_adj
  cost_infant_ART_first_2weeks <- (284.1/52)*cpi_adj
  cost_infant_ART_after_2weeks <- (284.1/52)*cpi_adj
  
  test_cost <- ifelse(test_gen == 3, cost_rapid_screening_3gen, cost_rapid_screening_4gen)
  
  
  # Infant incidence reduction parameters
  p_arv <- 0.987
  arv_red <- 0.675
  nbf <- 1
  prob_nbf_neg_early_pp <- 0.07
  prob_nbf_pos_early_pp <- 0.34
  prob_nbf_neg_mid_pp <- 0.19
  prob_nbf_pos_mid_pp <- 0.45
  prob_nbf_neg_late_pp <- 0.42
  prob_nbf_pos_late_pp <- 0.63
  
  
  # Infant mortality parameters
  mort_n <- 1-exp(-7*(7*44823.05+21*5306.55)/(100000*28*365))
  mort_i <- 1-exp(-7*1323.86/100000/365)
  
  # DALYs
  dw_HIVtx <- 0.078
  dw_HIVnotx <- 0.274
  dw_AIDSnotx <- 0.582
  propLE_HIV_notx <- 0.9
  
  
  if (population == "FSW") {
    
    # Population parameters
    population_size <- 9906
    first_ANC <- 19
    
    ## prevalence and incidence rates
    hiv_prev <- 0.577*(1-prev_reduction)
    status_known <- 0.81
    
    # Weekly incidence parameters
    inc_preg_early <- ((6.4/100)/52)*(1-incidence_reduction)
    inc_preg_late <- ((6.4/100)/52)*(1-incidence_reduction)
    inc_pp_early <- 0
    inc_pp_mid <- ((6.4/100)/52)*(1-incidence_reduction)
    inc_pp_late <- ((6.4/100)/52)*(1-incidence_reduction)
    
    att_firstANC <- 0.9
    
    # ART parameters
    p_ART <- 0.94
    ART_ret <- 0.6
    
  } else if (population == "PWID") {
    
    # Population parameters
    population_size <- 14309
    
    ## prevalence and incidence rates
    hiv_prev <- 0.577*(1-prev_reduction)
    status_known <- 0.80
    
    # Weekly incidence parameters
    inc_preg_early <- 0.00109615384615385*(1-incidence_reduction)
    inc_preg_late <- 0.00109615384615385*(1-incidence_reduction)
    inc_pp_early <- 0
    inc_pp_mid <- 0.00109615384615385*(1-incidence_reduction)
    inc_pp_late <- 0.00109615384615385*(1-incidence_reduction)
    
    # ART parameters
    p_ART <- 0.87
    ART_ret <- 0.6
    
  } else if (population == "serodiscordant") {
    
    # Population parameters
    population_size <- 55035
    
    ## prevalence and incidence rates
    hiv_prev <- 0*(1-prev_reduction)
    status_known <- 0.02
    
    # Weekly incidence parameters
    inc_preg_early <- 0.000721153846153846*(1-incidence_reduction)
    inc_preg_late <- 0.00135*(1-incidence_reduction)
    inc_pp_early <- 0
    inc_pp_mid <- 0.0009*(1-incidence_reduction)
    inc_pp_late <- 0.000240384615384615*(1-incidence_reduction)
    
    # ART parameters
    ART_ret <- 0.6
    ART_drop_factor <- exp(log(ART_ret)/91)
    
    
  } else {
    
    print("No population selected")
    
  }
  
  
  
} else if (country == "Ukraine") {
  
  inf_ART_coverage <- 0.95          # Can change this 
  inf_ART_coverage_base <- 0.95     # Don't change this
  adult_ART_coverage <- 0.52 
  
  red_prep <- 0.7049
  avg_ART_adherence <- 0.84
  
  # Miscellaneous
  transition_ag <- 2.3
  transition_ab <- 0.7
  transition_r <- 6
  
  # Calculated
  duration_infection <- sum(transition_ag, transition_ab, transition_r)
  
  
  # Testing probabilities
  att_firstANC <- 0.998
  att_lategest <- 0.9
  att_delivery <- 0.99
  att_6wk <- 0.99
  att_14wk <- 0.93
  att_6mo <- 0.96
  att_9mo <- 0.95
  
  test_accept <- 0.972
  stockout <- 0
  p_results <- 1
  
  
  # ART parameters
  p_ART <- 0.95
  ART_ret <- 0.6
  ART_drop_factor <- exp(log(ART_ret)/91)
  p_VL <- 0.881
  
  # Weekly mortality rates
  mort_adult_females <- (0+0+0.001*2+0.002*2+0.003)/7/52
  mort_maternal_perinatal <- (19/100000)/6
  
  
  
  
  ## Input parameters
  # Gestational age
  first_ANC <- 10
  second_ANC <- 28
  delivery <- 39
  # infant weeks at MCH visit
  early_pp_visit <- 6
  pp_14weeks_visit <- 14
  pp_mid_visit <- 26
  pp_9months_visit <- 39
  end_of_model <- 53
  
  
  # Testing probabilities
  second_ANC_test_prob <- late_ANC_test * att_lategest * test_accept * (1-stockout)
  delivery_test_prob <- (1-second_ANC_test_prob) * late_ANC_test * att_delivery * test_accept * (1-stockout)
  
  
  # Costs
  cpi_adj <- 1.23196
  
  cost_rapid_screening_3gen <- 3.99*cpi_adj
  cost_true_pos_screening <- 4.18*cpi_adj
  cost_false_pos_screening <- 19.8*cpi_adj
  cost_maternal_ART <- 32.8412037902989*cpi_adj
  cost_infant_prophylaxis <- 4*cpi_adj
  cost_maternal_PrEP <- 19.3811943772982*cpi_adj
  cost_infant_ART_first_2weeks <- cost_maternal_ART
  cost_infant_ART_after_2weeks <- cost_maternal_ART
  
  test_cost <- ifelse(test_gen == 3, cost_rapid_screening_3gen, cost_rapid_screening_4gen)
  
  
  # Infant incidence reduction parameters
  p_arv <- 0.984
  arv_red <- 0.675
  nbf <- 1
  prob_nbf_neg_early_pp <- 0.046
  prob_nbf_pos_early_pp <- 0.95
  prob_nbf_neg_mid_pp <- 0.023
  prob_nbf_pos_mid_pp <- 0.99
  prob_nbf_neg_late_pp <- 0.115
  prob_nbf_pos_late_pp <- 0.99
  
  
  # Infant mortality parameters
  mort_n <- 1-exp(-7*(7*17588.48+21*1979.87)/(100000*28*365))
  mort_i <- 1-exp(-7*300.72/100000/365)
  
  # DALYs
  dw_HIVtx <- 0.078
  dw_HIVnotx <- 0.274
  dw_AIDSnotx <- 0.582
  propLE_HIV_notx <- 0.9
  
  
  if (population == "FSW") {
    
    # Population parameters
    population_size <- 2413

    ## prevalence and incidence rates
    hiv_prev <- 0.052*(1-prev_reduction)
    status_known <- 0.58
    
    # Weekly incidence parameters
    inc_preg_early <- 0.00005*(1-incidence_reduction)
    inc_preg_late <- 0.00005*(1-incidence_reduction)
    inc_pp_early <- 0
    inc_pp_mid <- 0.00005*(1-incidence_reduction)
    inc_pp_late <- 0.00005*(1-incidence_reduction)
    
    # ART parameters
    ART_ret <- 0.6
    
  } else if (population == "PWID") {
    
    # Population parameters
    population_size <- 2882
    
    ## prevalence and incidence rates
    hiv_prev <- 0.314*(1-prev_reduction)
    status_known <- 0.43
    
    # Weekly incidence parameters
    inc_preg_early <- 0.0000224615384615385*(1-incidence_reduction)
    inc_preg_late <- 0.0000449230769230769*(1-incidence_reduction)
    inc_pp_early <- 0
    inc_pp_mid <- 0.0000336923076923077*(1-incidence_reduction)
    inc_pp_late <- 0.0000336923076923077*(1-incidence_reduction)
    
    # ART parameters
    p_ART <- 0.878
    ART_ret <- 0.6
    
    # Weekly mortality rates
    mort_adult_females <- 0.000251923076923077
    mort_maternal_perinatal <- 0.0000624
    
    
    
  } else {
    
    print("No population selected")
    
  }
  
  
  
} else {
  
  print("No country selected")
  
}






## Putting all variables into a data frame to add to parameters
# Get all objects in the global environment
vars <- ls()

# Filter only those that are length 1 and atomic (e.g., numeric, character, logical)
single_length_vars <- Filter(function(x) {
  obj <- get(x)
  is.atomic(obj) && length(obj) == 1
}, vars)

# Create a named list of those values
single_vals <- sapply(single_length_vars, get, simplify = FALSE)

# Convert to a data frame (single row)
df <- as.data.frame(single_vals)

# Optional: View the data frame
print(df)

df2 <- df %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "name")

# download variable worksheet
kenya_input <- readxl::read_xlsx(path = here("parameters", "inp_par_dmc.xlsx"), sheet = "pop_characteristics") %>%
  filter(country == "Kenya")

kenya_input <- kenya_input %>%
  left_join(df2, by = c("var_name" = "name")) %>%
  mutate(value = as.numeric(V1)) %>%
  select(-V1)


write.csv(kenya_input, here("parameters", "kenya_params.csv"))

# Starting states
# neg <- (1-hiv_prev) - (inc_preg_early*duration_infection)
# recent <- duration_infection*inc_preg_early
# established_uk_nart <- hiv_prev*(1-status_known)
# established_k_nart <- hiv_prev*status_known*(1-p_ART)
# established_art <- hiv_prev*status_known*p_ART


# parameter_list <- c(paste0("population_size=", population_size), paste0("neg=", neg))
# 
# return(parameter_list)
# 
# 
# }
# 
# parameter_list <- maternal_input_func(population="Kenya")
# maternal_func(parameter_list)

