
#### Title: HIV retesting model
#### Purpose: Create HIV retesting model based on Excel models
#### Author: David Coomes
#### Date: January 9, 2025


# Load libraries
pacman::p_load(
  "EpiModel", 
  "tidyverse", 
  "here",
  "data.table"
)


# Needs:
  # add to function warning/stop saying that if test_type = "dual" then syph_test must equal TRUE



###############
## HIV model ##
###############

hiv_maternal_func <- function(first_ANC_test, late_ANC_test,                       # first two inputs are for ANC testing (1=test administered, 0=no test)
                              test_type_1 = "rdt", test_type_2 = "rdt",            # these parameters indicate which test is used (default is "rdt")
                              syph_test = FALSE,                                   # include as TRUE if we want to include syphilis infant mortality
                              syph_infant_deaths, 
                              # these inputs are for sensitivity
                              prev_reduction=0, PrEP=0, incidence_reduction=0,   
                              # pp testing set to no-testing
                              pp_early_test=0, pp_14weeks_test=0, pp_mid_test=0, pp_9months_test=0,      # Post-partum testing set to none by default
                              pp_test_type = "rdt") {
  
  ### This function uses population parameter inputs, testing frequency, and sensitivity inputs 
  ### to create a data frame of population sizes based on HIV status, testing, and treatment
  ### Inputs: first and second ANC tests, reduction in prevalence and incidence as well as PrEP updake
    ### You can adjust post-partum testing but we won't likely change these for this work
  ### Output: dataframe of maternal population by HIV status, testing, and treatment

  syph_infant_deaths <- syph_infant_deaths
  
  # Starting states
  neg <- (1-hiv_prev) - (inc_preg_early*duration_infection)
  recent <- duration_infection*inc_preg_early
  established_uk_nart <- hiv_prev*(1-hiv_status_known)
  established_k_nart <- hiv_prev*hiv_status_known*(1-p_ART)
  established_art <- hiv_prev*hiv_status_known*p_ART
  
   
  # Setting up empty matrix
  maternal_matrix <- matrix(nrow = 92, ncol = 9, 
                       dimnames = list(NULL, 
                                       c("weeks", "hiv_neg", "ahiv_pos_nart", 
                                         "ahiv_pos_art", "chiv_pos_uk", 
                                         "chiv_pos_nart", "chiv_pos_art", "deaths",
                                         "stage")))
  
  maternal_df <- as.data.frame(maternal_matrix)
  
  maternal_df[1, ] <- c(1, population_size*neg, population_size*recent, 0, population_size*established_uk_nart,
                        population_size*established_k_nart, 0, 0, 0)

  maternal_df$stage <- as.character(maternal_df$stage) 
  maternal_df$stage[1] <- "Stage 0: Onset of pregnancy to first ANC"
  
  
  # Testing probabilities
  second_ANC_test_prob <- late_ANC_test * att_lategest * test_accept * (1-stockout)
  delivery_test_prob <- (1-second_ANC_test_prob) * late_ANC_test * att_delivery * test_accept * (1-stockout)
  
  
  ###########
  # STAGE 0 #
  ###########
  
  
  # Getting estimates for stage 0: onset of pregnancy to first ANC
  for (x in 2:(weekofanc1-1)) {
    
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_early * (1-(PrEP*red_prep))))
    
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_early * (1-(PrEP*red_prep)) +
                               maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) +
                               maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
    
    maternal_df$ahiv_pos_art[x] <- 0
    
    maternal_df$chiv_pos_uk[x] <- maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) + 
                             maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection)
    
    maternal_df$chiv_pos_nart[x] <- maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                               maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
    
    maternal_df$chiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                              maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                     maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 0: Onset of pregnancy to first ANC" 
    
  }
  
  
  ###########
  # STAGE 1 #
  ###########
  
  
  # Getting estimates for stage 1: First ANC to second ANC
  for (x in weekofanc1:(weekofanc2-1)) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(test_type_1 == "rdt" ~ rdt_sens_hiv_e,
                             test_type_1 == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(test_type_1 == "rdt" ~ rdt_sens_hiv_c,
                             test_type_1 == "dual" ~ dual.sens.hiv.c)
    
    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    # Calculating populations
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc1-1), 
                                maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_early * (1-(PrEP*red_prep)))),
                                maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_late * (1-(PrEP*red_prep))))
                                )
    
    maternal_df$ahiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc1-1),
                                      # HIV negative who become infected
                                      maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_early * (1-(PrEP*red_prep)) +
                                        # aHIV+ who remain positive but do not receive ART
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * 
                                          (1-first_ANC_test*att_firstANC*(1-stockout)*test_accept* 
                                             test_sens_combo_e *   # test_sens_e: DMC adding in HIV RDT  
                                           p_results*p_ART) +
                                        # aHIV+ who drop ART
                                        maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor),
                                      maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_late * (1-(PrEP*red_prep)) +
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) +
                                        maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
                                      )
      
    maternal_df$ahiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc1-1),
                                     maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * 
                                      first_ANC_test * att_firstANC * test_accept * (1-stockout) * test_sens_combo_e * p_results * p_ART,
                                     maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor
                                     )
      
    maternal_df$chiv_pos_uk[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc1-1),
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                           (1-(first_ANC_test * att_firstANC * (1-stockout) * test_accept * test_sens_combo_c * p_results)) + 
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                           (1-(first_ANC_test * att_firstANC * (1-stockout) * test_accept * test_sens_combo_c * p_results)),
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females)
                                         )
  
    maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc1-1),
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * first_ANC_test * 
                                             att_firstANC * test_accept * (1-stockout) * p_results * test_sens_combo_c * (1-p_ART) +
                                             maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * first_ANC_test * att_firstANC * test_accept *
                                             (1-stockout) * p_results * test_sens_combo_c * (1-p_ART) +
                                             maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                             maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor),
                                           maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                             maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
                                           )
    
    maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc1-1),
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * first_ANC_test *
                                            att_firstANC * test_sens_combo_c * test_accept * (1-stockout) * p_results * p_ART +
                                            maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                            maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * first_ANC_test * att_firstANC * test_sens_combo_c *
                                            test_accept * (1-stockout) * p_results * p_ART +
                                            maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor,
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                            maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
                                          )
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                 maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
  
    maternal_df$stage[x] <- "Stage 1: First ANC to second ANC" 
  
  }
    
  
  
  
  ###########
  # STAGE 2 #
  ###########
  
  
  # Getting estimates for Stage 2: Second ANC to deliveryweek
  for (x in weekofanc2:(deliveryweek-1)) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(test_type_2 == "rdt" ~ rdt_sens_hiv_e,
                             test_type_2 == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(test_type_2 == "rdt" ~ rdt_sens_hiv_c,
                             test_type_2 == "dual" ~ dual.sens.hiv.c)

    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    # Calculating population
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_late * (1-(PrEP*red_prep))))
  
    maternal_df$ahiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc2 - 1),
                                           maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_late * (1-(PrEP*red_prep)) +
                                             maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * 
                                             (1-late_ANC_test*att_lategest*(1-stockout)*test_accept*test_sens_combo_e*p_results*p_ART) +
                                             maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor),
                                           maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_late * (1-(PrEP*red_prep)) +
                                             maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) +
                                             maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
                                           )
    
      maternal_df$ahiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc2 - 1),
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * late_ANC_test *
                                            att_lategest * test_accept * (1-stockout) * test_sens_combo_e * p_results * p_ART +
                                            maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor,
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor
                                          )
  
    maternal_df$chiv_pos_uk[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc2 - 1),
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                           (1-(late_ANC_test * att_lategest * (1-stockout) * test_accept * test_sens_combo_c * p_results)) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                           (1-(late_ANC_test * att_lategest * (1-stockout) * test_accept * test_sens_combo_c * p_results)),
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females)
                                         )
    
    maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc2 - 1),
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * late_ANC_test *
                                             att_lategest * test_accept * (1-stockout) * p_results * test_sens_combo_c * (1-p_ART) +
                                             maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * late_ANC_test * att_lategest *
                                             test_accept * (1-stockout) * p_results * test_sens_combo_c * (1-p_ART) +
                                             maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) + 
                                             maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor),
                                           maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                             maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
                                           )
    
    maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (weekofanc2 - 1),
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * late_ANC_test *
                                            att_lategest * test_sens_combo_c * test_accept * (1-stockout) * p_results * p_ART +
                                            maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                            maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * late_ANC_test * att_lategest * test_sens_combo_c *
                                            test_accept * (1-stockout) * p_results * p_ART +
                                            maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor,
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                            maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
                                          )
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 2: Second ANC to deliveryweek" 
    
  }
  
  
  
  ###########
  # STAGE 3 #
  ###########
  
  
  
  # Getting estimates for Stage 3: deliveryweek
  for (x in deliveryweek:deliveryweek) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(test_type_2 == "rdt" ~ rdt_sens_hiv_e,
                             test_type_2 == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(test_type_2 == "rdt" ~ rdt_sens_hiv_c,
                             test_type_2 == "dual" ~ dual.sens.hiv.c)
    
    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    
    # Calculating populations
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * (1-(inc_preg_late * (1-(PrEP*red_prep)))) -
                                             # Adding in syphilis mortality only if syphilis testing is included
                                             ifelse(syph_test == 1, 
                                                    syph_infant_deaths * (maternal_df$hiv_neg[x-1] / sum(maternal_df[x-1, 2:7], na.rm=TRUE)),
                                                    0)
    
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * inc_preg_late * (1-(PrEP*red_prep)) +
                                             maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1-1/duration_infection) * 
                                             (1-late_ANC_test*att_delivery*(1-stockout)*test_accept*(1-second_ANC_test_prob)*test_sens_combo_e*p_results*p_ART) +
                                             maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-(1/duration_infection)) * (1-ART_drop_factor) -
                                             # Adding in syphilis mortality only if syphilis testing is included
                                             ifelse(syph_test == 1, 
                                                    syph_infant_deaths * (maternal_df$ahiv_pos_nart[x-1] / sum(maternal_df[x-1, 2:7], na.rm=TRUE)),
                                                    0)
  
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1-1/duration_infection) * late_ANC_test *
                                            att_delivery * test_accept * (1-stockout) * (1-second_ANC_test_prob) * test_sens_combo_e * p_results * p_ART +
                                            maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-(1/duration_infection)) * ART_drop_factor -
                                             # Adding in syphilis mortality only if syphilis testing is included
                                             ifelse(syph_test == 1, 
                                                    syph_infant_deaths*maternal_df$ahiv_pos_art[x-1] / sum(maternal_df[x-1, 2:7], na.rm=TRUE),
                                                    0)
  
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) * 
                                           (1-(late_ANC_test * att_delivery * (1-stockout) * test_accept * (1-second_ANC_test_prob) * 
                                           test_sens_combo_c * p_results)) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_maternal_perinatal) * 
                                           (1-(late_ANC_test * att_delivery * (1-stockout) * test_accept * (1-second_ANC_test_prob) *
                                           test_sens_combo_c * p_results)) -
                                           # Adding in syphilis mortality only if syphilis testing is included
                                           ifelse(syph_test == 1, 
                                                  syph_infant_deaths*maternal_df$chiv_pos_uk[x-1] / sum(maternal_df[x-1, 2:7], na.rm=TRUE),
                                                  0)

  
    maternal_df$chiv_pos_nart[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) * late_ANC_test *
                                             att_delivery * test_accept * (1-stockout) * p_results * (1-second_ANC_test_prob) * 
                                             test_sens_combo_c * (1-p_ART) +
                                             maternal_df$chiv_pos_uk[x-1] * (1-mort_maternal_perinatal) * late_ANC_test * att_delivery *
                                             test_accept * (1-stockout) * p_results * (1-second_ANC_test_prob) * test_sens_combo_c * (1-p_ART) +
                                             maternal_df$chiv_pos_nart[x-1] * (1-mort_maternal_perinatal) + 
                                             maternal_df$chiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-ART_drop_factor) -
                                             # Adding in syphilis mortality only if syphilis testing is included
                                             ifelse(syph_test == 1, 
                                                    syph_infant_deaths*maternal_df$chiv_pos_nart[x-1] / sum(maternal_df[x-1, 2:7], na.rm=TRUE),
                                                    0)
      
  
    maternal_df$chiv_pos_art[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) * late_ANC_test *
                                            att_delivery * test_sens_combo_c * test_accept * (1-stockout) * p_results * p_ART * (1-second_ANC_test_prob) +
                                            maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) +
                                            maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * late_ANC_test * att_delivery * test_sens_combo_c *
                                            test_accept * (1-stockout) * p_results * p_ART * (1-second_ANC_test_prob) +
                                            maternal_df$chiv_pos_art[x-1] * (1-mort_maternal_perinatal) * ART_drop_factor -
                                            # Adding in syphilis mortality only if syphilis testing is included
                                            ifelse(syph_test == 1, 
                                                   syph_infant_deaths*maternal_df$chiv_pos_art[x-1] / sum(maternal_df[x-1, 2:7], na.rm=TRUE),
                                                   0)

    
    maternal_df$deaths[x] <- mort_maternal_perinatal * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 3: deliveryweek" 
    
  }
  
  
  
  ###########
  # STAGE 4 #
  ###########
  
  
  # Getting estimates for Stage 4: 1-6 weeks post-partum
  for (x in (deliveryweek+1):(deliveryweek+early_pp_visit-1)) {
    
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * 
      (1-(ifelse(maternal_df$weeks[x-1] == deliveryweek, inc_preg_late, inc_pp_early) * (1-(PrEP*red_prep)))) 
    
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * 
                                      ifelse(maternal_df$weeks[x-1] == deliveryweek, inc_preg_late, inc_pp_early) * (1-(PrEP*red_prep)) +
                                    maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1-1/duration_infection)  +
                                    maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-(1/duration_infection)) * (1-ART_drop_factor)
    
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-(1/duration_infection)) * ART_drop_factor
    
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) +
                                  maternal_df$chiv_pos_uk[x-1] * (1-mort_maternal_perinatal)
    
    maternal_df$chiv_pos_nart[x] <- maternal_df$chiv_pos_nart[x-1] * (1-mort_maternal_perinatal) +
                                    maternal_df$chiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-ART_drop_factor)
  
    maternal_df$chiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) +
                                   maternal_df$chiv_pos_art[x-1] * (1-mort_maternal_perinatal) * ART_drop_factor
    
    maternal_df$deaths[x] <- mort_maternal_perinatal * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                           maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 4: 1-6 weeks post-partum" 
    
  }
  
  
  
  ###########
  # STAGE 5 #
  ###########
  
  
  # Getting estimates for Stage 5: 6-14 weeks post-partum
  for (x in (deliveryweek+early_pp_visit):(deliveryweek+pp_14weeks_visit-1)) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_e,
                             pp_test_type == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_c,
                             pp_test_type == "dual" ~ dual.sens.hiv.c)
    
    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    # Calculating populations
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * 
      (1-(ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, inc_pp_early, inc_pp_mid) * (1-(PrEP*red_prep)))) 
    
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * 
                                      ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, inc_pp_early, inc_pp_mid) * (1-(PrEP*red_prep)) +
                                    maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                      ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, (1 - pp_early_test * att_6wk * (1-stockout) * test_accept * (1- (second_ANC_test_prob + delivery_test_prob)) *
                                      test_sens_combo_e * p_results * p_ART), 1) +
                                    maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
    
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_early_test * 
                                          att_6wk * test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) * test_sens_combo_e *
                                          p_results * p_ART, 0)
                                   
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                           (1-pp_early_test * att_6wk * (1-stockout) * test_accept * (1-(second_ANC_test_prob + delivery_test_prob)) *
                                           test_sens_combo_c * p_results), 1) +
                                  maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                           (1-pp_early_test * att_6wk * (1-stockout) * test_accept * (1-(second_ANC_test_prob + delivery_test_prob)) * 
                                              test_sens_combo_c * p_results), 1)
   
    maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_early_test *
                                             att_6wk * test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) *
                                             p_results * test_sens_combo_c * (1-p_ART), 0) +
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_early_test * att_6wk * test_accept * (1-stockout) *
                                             (1-(second_ANC_test_prob + delivery_test_prob)) * p_results * test_sens_combo_c * (1-p_ART), 0) +      
                                    maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                    maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
   
    maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_early_test *
                                            att_6wk * test_sens_combo_c * test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) *
                                            p_results * p_ART, 0) +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+early_pp_visit)-1, 
                                          maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_early_test * att_6wk * test_sens_combo_c *
                                            test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) * p_results * p_ART, 0) +
                                   maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                   maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                           maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 5: 6-14 weeks post partum" 
    
  }
  
  
  
  
  ###########
  # STAGE 6 #
  ###########
  
  
  # Getting estimates for Stage 6:14 weeks - 6 months post-partum
  for (x in (deliveryweek+pp_14weeks_visit):(deliveryweek+pp_mid_visit-1)) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_e,
                             pp_test_type == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_c,
                             pp_test_type == "dual" ~ dual.sens.hiv.c)
    
    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    # Calculating populations
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_mid * (1-(PrEP*red_prep)))) 
    
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_mid * (1-(PrEP*red_prep)) +
                                    maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                      ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                             (1 - (pp_14weeks_test * att_14wk * (1-stockout) * test_accept * 
                                             test_sens_combo_e * p_results * p_ART)), 1) +
                                    maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
  
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_14weeks_test * 
                                          att_14wk * test_accept * (1-stockout) * test_sens_combo_e * p_results * p_ART, 0)
  
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                  ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                         (1 - pp_14weeks_test * att_14wk * (1-stockout) * test_accept * test_sens_combo_c * p_results), 1) +
                                  maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                  ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                         (1 - pp_14weeks_test * att_14wk * (1-stockout) * test_accept * test_sens_combo_c * p_results), 1)
    
    
    
    maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_14weeks_test *
                                             att_14wk * test_accept * (1-stockout) * p_results * test_sens_combo_c * (1-p_ART), 0) +
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_14weeks_test * att_14wk * test_accept * 
                                             (1-stockout) * p_results * test_sens_combo_c * (1-p_ART), 0) +      
                                    maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                    maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
    
    maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_14weeks_test *
                                            att_14wk * test_sens_combo_c * test_accept * (1-stockout) * p_results * p_ART, 0) +
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_14weeks_visit)-1, 
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_14weeks_test * att_14wk * test_sens_combo_c *
                                             test_accept * (1-stockout) * p_results * p_ART, 0) +
                                    maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                    maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 6:14 weeks - 6 months post-partum" 
    
  }
  
  
  
  
  ###########
  # STAGE 7 #
  ###########
  
  
  # Getting estimates for Stage 7: 6-9 months post-partum
  for (x in (deliveryweek+pp_mid_visit):(deliveryweek+pp_9months_visit-1)) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_e,
                             pp_test_type == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_c,
                             pp_test_type == "dual" ~ dual.sens.hiv.c)

    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    # Calculating populations
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit-1),
                                     maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_mid * (1-(PrEP*red_prep)))),
                                     maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_late * (1-(PrEP*red_prep))))) 
    
    
    maternal_df$ahiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit-1),
                                           maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_mid * (1-(PrEP*red_prep)),
                                           maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_late * (1-(PrEP*red_prep))) +
                                    maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                           (1 - (pp_mid_test * att_6mo * (1-stockout) * test_accept * 
                                                   test_sens_combo_e * p_results * p_ART)), 1) +
                                    maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
    
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_mid_test * 
                                            att_6mo * test_accept * (1-stockout) * test_sens_combo_e * p_results * p_ART, 0)
    
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                  ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                         (1 - pp_mid_test * att_6mo * (1-stockout) * test_accept * test_sens_combo_c * p_results), 1) +
                                  maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                  ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                         (1 - pp_mid_test * att_6mo * (1-stockout) * test_accept * test_sens_combo_c * p_results), 1)
    
    
    
    maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_mid_test *
                                             att_6mo * test_accept * (1-stockout) * p_results * test_sens_combo_c * (1-p_ART), 0) +
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_mid_test * att_6mo * test_accept * 
                                             (1-stockout) * p_results * test_sens_combo_c * (1-p_ART), 0) +      
                                    maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                    maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
    
    maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_mid_test *
                                            att_6mo * test_sens_combo_c * test_accept * (1-stockout) * p_results * p_ART, 0) +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_mid_visit)-1, 
                                          maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_mid_test * att_6mo * test_sens_combo_c *
                                             test_accept * (1-stockout) * p_results * p_ART, 0) +
                                   maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                   maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 7: 6-9 months post-partum" 
    
  }
  
  
  
  ###########
  # STAGE 8 #
  ###########
  
  
  # Getting estimates for Stage 8: 9-12 months post-partum
  for (x in (deliveryweek+pp_9months_visit):(deliveryweek+end_of_model)) {
    
    # Setting test sensitivity based on test type and stage
    test_sens_e <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_e,
                             pp_test_type == "dual" ~ dual.sens.hiv.e) 
    test_sens_c <- case_when(pp_test_type == "rdt" ~ rdt_sens_hiv_c,
                             pp_test_type == "dual" ~ dual.sens.hiv.c)
    
    # Accounting for those that only get an HIV RDT instead of dual test because they are treponemal positive
    test_sens_combo_e <- (test_sens_e + test_sens_e*syp_trep_pos - rdt_sens_hiv_e*syp_trep_pos)
    test_sens_combo_c <- (test_sens_c + test_sens_c*syp_trep_pos - rdt_sens_hiv_c*syp_trep_pos)
    
    # Calculating populations
    maternal_df$weeks[x] <- x
    
    maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_late * (1-(PrEP*red_prep))))
    
    
    maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_late * (1-(PrEP*red_prep)) +
                                    maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                           (1 - (pp_9months_test * att_9mo * (1-stockout) * test_accept * 
                                                   test_sens_combo_e * p_results * p_ART)), 1) +
                                    maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
    
    maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_9months_test * 
                                            att_9mo * test_accept * (1-stockout) * test_sens_combo_e * p_results * p_ART, 0)
    
    maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                  ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                         (1 - pp_9months_test * att_9mo * (1-stockout) * test_accept * test_sens_combo_c * p_results), 1) +
                                  maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                  ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                         (1 - pp_9months_test * att_9mo * (1-stockout) * test_accept * test_sens_combo_c * p_results), 1)
    
    maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_9months_test *
                                             att_9mo * test_accept * (1-stockout) * p_results * test_sens_combo_c * (1-p_ART), 0) +
                                    ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_9months_test * att_9mo * test_accept * 
                                             (1-stockout) * p_results * test_sens_combo_c * (1-p_ART), 0) +      
                                    maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                    maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
    
    maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                          maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_9months_test *
                                            att_9mo * test_sens_combo_c * test_accept * (1-stockout) * p_results * p_ART, 0) +
                                   ifelse(maternal_df$weeks[x-1] == (deliveryweek+pp_9months_visit)-1, 
                                          maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_9months_test * att_9mo * test_sens_combo_c *
                                            test_accept * (1-stockout) * p_results * p_ART, 0) +
                                   maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                   maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
    
    maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                      maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
    
    maternal_df$stage[x] <- "Stage 8: 9-12 months post-partum" 
    
  }

return(maternal_df)

}

#################








###########################
## Person-week dataframe ##
###########################


##### Getting output from model

# Getting total person weeks from model
# person_weeks <- maternal_df %>%
#   group_by(stage) %>%
#   summarize(across(hiv_neg:chiv_pos_art, ~sum(.)))


pw_func <- function(maternal_df) {
  out <- maternal_df %>%
    group_by(stage) %>%
    summarize(across(hiv_neg:chiv_pos_art, ~sum(.)))
  
  return(out)
}



# test_type_1 <- "dual"
# test_type_2 <- "dual"
# test_type_pp <- "dual"
# first_ANC_test <- 1
# late_ANC_test <- 0
# pp_early_test <- 0
# pp_mid_test <- 0
# pp_late_test <- 0
# pp_14weeks_test <- 0
# pp_9months_test <- 0
# pp_6months_test <-0
# 
# maternal_hiv_df <- hiv_maternal_func(1, 0, "dual", "dual", syph_infant_deaths = syph_infant_deaths)
# person_weeks <- pw_func(maternal_hiv_df)
# mat_costs <- mat_cost_func(person_weeks, maternal_hiv_df, 1, 1, "rdt", "dual")


##############################
## Maternal costs dataframe ##
##############################
 
mat_cost_func <- function(person_weeks, maternal_hiv_df,                        # requires person-week and maternal dataframes
                          first_ANC_test, late_ANC_test,                        # two inputs are for ANC testing (1=test administered, 0=no test)
                          test_type_1="rdt", test_type_2="rdt",                 # insert test type for first and second tests
                          # these inputs are for sensitivity
                          prev_reduction=0, PrEP=0, incidence_reduction=0,   
                          # pp testing set to no-testing
                          pp_early_test=0, pp_14weeks_test=0, pp_mid_test=0, pp_9months_test=0,
                          test_type_pp="rdt") {
  
  maternal_df <- maternal_hiv_df
  
  # Testing probabilities
  second_ANC_test_prob <- late_ANC_test * att_lategest * test_accept * (1-stockout)
  delivery_test_prob <- (1-second_ANC_test_prob) * late_ANC_test * att_delivery * test_accept * (1-stockout)
  
  ## DMC adding 8.21.25 - splitting costs by test, treatment, PrEP
  
  ####### Testing costs (HIV RDT or HIV-syphilis dual test) ############
    
  # Split by HIV-syphilis dual testing and HIV RDT testing costs
  # Starting with dual testing costs
  mat_costs_td0 <- c(rep(0, 6))
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_1 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_1 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("First test type not defined")
  }
  
  if (test_type_1 == "dual") {
    mat_costs_td1 <- c(maternal_df[[weekofanc1, "hiv_neg"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[weekofanc1, "ahiv_pos_nart"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                       (test_cost + test_sens_e * p_results * cost_true_pos_screening) * (1-syp_trep_pos),
                     0,
                     maternal_df[[weekofanc1, "chiv_pos_uk"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                       (test_cost + test_sens_c * p_results * cost_true_pos_screening) * (1-syp_trep_pos),
                     0,
                     0)
  } else if (test_type_1 == "rdt") {
    mat_costs_td1 <- c(rep(0, 6))
  } else {
    print("First test type not defined")
  }
  
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_2 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_2 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Second test type not defined")
  }
  
  if (test_type_2 == "dual") {
    mat_costs_td2 <- c(maternal_df[[weekofanc2, "hiv_neg"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[weekofanc2, "ahiv_pos_nart"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                       (test_cost + test_sens_e * p_results * cost_true_pos_screening) * (1-syp_trep_pos),
                     0,
                     maternal_df[[weekofanc2, "chiv_pos_uk"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                       (test_cost + test_sens_c * p_results * cost_true_pos_screening) * (1-syp_trep_pos),
                     0,
                     0)
    
    mat_costs_td3 <- c(maternal_df[[deliveryweek, "hiv_neg"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                       (1-second_ANC_test_prob) *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[deliveryweek, "ahiv_pos_nart"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                       (1-second_ANC_test_prob) * 
                       (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     maternal_df[[deliveryweek, "chiv_pos_uk"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                       (1-second_ANC_test_prob) *
                       (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     0)
    
  } else if (test_type_2 == "rdt") {
    mat_costs_td2 <- c(rep(0, 6))
    mat_costs_td3 <- c(rep(0, 6))
  } else {
    print("Second test type not defined")
  }
  
  mat_costs_td4 <- c(rep(0, 6))
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_pp == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_pp == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Postpartum test type not defined")
  }
  
  if (test_type_pp == "dual") {
    
    mat_costs_td5 <- c(maternal_df[[deliveryweek+early_pp_visit, "hiv_neg"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                       (1-(second_ANC_test_prob+delivery_test_prob)) *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[deliveryweek+early_pp_visit, "ahiv_pos_nart"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                       (1-(second_ANC_test_prob+delivery_test_prob)) * 
                       (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     maternal_df[[deliveryweek+early_pp_visit, "chiv_pos_uk"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                       (1-(second_ANC_test_prob+delivery_test_prob)) *
                       (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     0)
    
    
    mat_costs_td6 <- c(maternal_df[[deliveryweek+pp_14weeks_visit, "hiv_neg"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[deliveryweek+pp_14weeks_visit, "ahiv_pos_nart"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                       (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     maternal_df[[deliveryweek+pp_14weeks_visit, "chiv_pos_uk"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                       (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     0)
    
    # There is one cost here different than the Excel models
    # the excel models have the total maternal costs for acute HIV- women as $0 but I think there's a linking error
    mat_costs_td7 <- c(maternal_df[[deliveryweek+pp_mid_visit, "hiv_neg"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[deliveryweek+pp_mid_visit, "ahiv_pos_nart"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                       (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     maternal_df[[deliveryweek+pp_mid_visit, "chiv_pos_uk"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                       (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     0)
    
    
    mat_costs_td8 <- c(maternal_df[[deliveryweek+pp_9months_visit, "hiv_neg"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                       (test_cost + (1-test_spec) * cost_false_pos_screening) * (1-syp_trep_pos),
                     maternal_df[[deliveryweek+pp_9months_visit, "ahiv_pos_nart"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                       (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     maternal_df[[deliveryweek+pp_9months_visit, "chiv_pos_uk"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                       (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (1-syp_trep_pos),
                     0,
                     0)
  
  } else if (test_type_pp == "rdt") {
    
    mat_costs_td5 <- c(rep(0, 6))
    mat_costs_td6 <- c(rep(0, 6))
    mat_costs_td7 <- c(rep(0, 6))
    mat_costs_td8 <- c(rep(0, 6))

  } else {
    print("PP test type not defined")
  }
  
  
  
  # HIV RDT testing costs
    # HIV RDT testing costs are the full testing costs when using HIV RDT, 
    # and are testing costs for those with trep pos when using dual HIV-syphilis
  mat_costs_tr0 <- c(rep(0, 6))
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_1 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_1 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("First test type not defined")
  }
  
  if (test_type_1 == "rdt") {
    mat_costs_tr1 <- c(maternal_df[[weekofanc1, "hiv_neg"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening),
                       maternal_df[[weekofanc1, "ahiv_pos_nart"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                         (test_cost + test_sens_e * p_results * cost_true_pos_screening),
                       0,
                       maternal_df[[weekofanc1, "chiv_pos_uk"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                         (test_cost + test_sens_c * p_results * cost_true_pos_screening),
                       0,
                       0)
  } else if (test_type_1 == "dual") {
    
    mat_costs_tr1 <- c(maternal_df[[weekofanc1, "hiv_neg"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[weekofanc1, "ahiv_pos_nart"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                         (test_cost + test_sens_e * p_results * cost_true_pos_screening) * (syp_trep_pos),
                       0,
                       maternal_df[[weekofanc1, "chiv_pos_uk"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                         (test_cost + test_sens_c * p_results * cost_true_pos_screening) * (syp_trep_pos),
                       0,
                       0)
    
  } else {
    print("First test type not defined")
  }
  
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_2 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_2 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Second test type not defined")
  }
  
  if (test_type_2 == "rdt") {
    
    mat_costs_tr2 <- c(maternal_df[[weekofanc2, "hiv_neg"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening),
                       maternal_df[[weekofanc2, "ahiv_pos_nart"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                         (test_cost + test_sens_e * p_results * cost_true_pos_screening),
                       0,
                       maternal_df[[weekofanc2, "chiv_pos_uk"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                         (test_cost + test_sens_c * p_results * cost_true_pos_screening),
                       0,
                       0)
    
    mat_costs_tr3 <- c(maternal_df[[deliveryweek, "hiv_neg"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                         (1-second_ANC_test_prob) *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[deliveryweek, "ahiv_pos_nart"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                         (1-second_ANC_test_prob) * 
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       maternal_df[[deliveryweek, "chiv_pos_uk"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                         (1-second_ANC_test_prob) *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       0)
    
  } else if (test_type_2 == "dual") {

    mat_costs_tr2 <- c(maternal_df[[weekofanc2, "hiv_neg"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[weekofanc2, "ahiv_pos_nart"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                         (test_cost + test_sens_e * p_results * cost_true_pos_screening) * (syp_trep_pos),
                       0,
                       maternal_df[[weekofanc2, "chiv_pos_uk"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                         (test_cost + test_sens_c * p_results * cost_true_pos_screening) * (syp_trep_pos),
                       0,
                       0)
    
    mat_costs_tr3 <- c(maternal_df[[deliveryweek, "hiv_neg"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                         (1-second_ANC_test_prob) *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[deliveryweek, "ahiv_pos_nart"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                         (1-second_ANC_test_prob) * 
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       maternal_df[[deliveryweek, "chiv_pos_uk"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                         (1-second_ANC_test_prob) *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       0)
  } else {
    print("Second test type not defined")
  }
  
  mat_costs_tr4 <- c(rep(0, 6))
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_pp == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_pp == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Postpartum test type not defined")
  }
  
  if (test_type_pp == "rdt") {
    
    mat_costs_tr5 <- c(maternal_df[[deliveryweek+early_pp_visit, "hiv_neg"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                         (1-(second_ANC_test_prob+delivery_test_prob)) *
                         (test_cost + (1-test_spec) * cost_false_pos_screening),
                       maternal_df[[deliveryweek+early_pp_visit, "ahiv_pos_nart"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                         (1-(second_ANC_test_prob+delivery_test_prob)) * 
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)),
                       0,
                       maternal_df[[deliveryweek+early_pp_visit, "chiv_pos_uk"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                         (1-(second_ANC_test_prob+delivery_test_prob)) *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)),
                       0,
                       0)
    
    
    mat_costs_tr6 <- c(maternal_df[[deliveryweek+pp_14weeks_visit, "hiv_neg"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening),
                       maternal_df[[deliveryweek+pp_14weeks_visit, "ahiv_pos_nart"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)),
                       0,
                       maternal_df[[deliveryweek+pp_14weeks_visit, "chiv_pos_uk"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)),
                       0,
                       0)
    
    mat_costs_tr7 <- c(maternal_df[[deliveryweek+pp_mid_visit, "hiv_neg"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening),
                       maternal_df[[deliveryweek+pp_mid_visit, "ahiv_pos_nart"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)),
                       0,
                       maternal_df[[deliveryweek+pp_mid_visit, "chiv_pos_uk"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)),
                       0,
                       0)
    
    
    mat_costs_tr8 <- c(maternal_df[[deliveryweek+pp_9months_visit, "hiv_neg"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening),
                       maternal_df[[deliveryweek+pp_9months_visit, "ahiv_pos_nart"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)),
                       0,
                       maternal_df[[deliveryweek+pp_9months_visit, "chiv_pos_uk"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)),
                       0,
                       0)
    
  } else if (test_type_pp == "dual") {
    
    mat_costs_tr5 <- c(maternal_df[[deliveryweek+early_pp_visit, "hiv_neg"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                         (1-(second_ANC_test_prob+delivery_test_prob)) *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[deliveryweek+early_pp_visit, "ahiv_pos_nart"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                         (1-(second_ANC_test_prob+delivery_test_prob)) * 
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       maternal_df[[deliveryweek+early_pp_visit, "chiv_pos_uk"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                         (1-(second_ANC_test_prob+delivery_test_prob)) *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       0)
    
    
    mat_costs_tr6 <- c(maternal_df[[deliveryweek+pp_14weeks_visit, "hiv_neg"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[deliveryweek+pp_14weeks_visit, "ahiv_pos_nart"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       maternal_df[[deliveryweek+pp_14weeks_visit, "chiv_pos_uk"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       0)
    
    mat_costs_tr7 <- c(maternal_df[[deliveryweek+pp_mid_visit, "hiv_neg"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[deliveryweek+pp_mid_visit, "ahiv_pos_nart"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       maternal_df[[deliveryweek+pp_mid_visit, "chiv_pos_uk"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       0)
    
    
    mat_costs_tr8 <- c(maternal_df[[deliveryweek+pp_9months_visit, "hiv_neg"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                         (test_cost + (1-test_spec) * cost_false_pos_screening) * (syp_trep_pos),
                       maternal_df[[deliveryweek+pp_9months_visit, "ahiv_pos_nart"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                         (test_cost + test_sens_e * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       maternal_df[[deliveryweek+pp_9months_visit, "chiv_pos_uk"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                         (test_cost + test_sens_c * p_results * (cost_true_pos_screening)) * (syp_trep_pos),
                       0,
                       0)
    
  } else {
    print("PP test type not defined")
  }
  
  
  
  # Treatment costs
  # Getting maternal costs for model by stage
  mat_costs_a0 <- c(rep(0, 6))
  
  mat_costs_a1 <- c(0,
                   0,
                   person_weeks[[2,"ahiv_pos_art"]] * cost_maternal_ART,
                   0,
                   0,
                   person_weeks[[2,"chiv_pos_art"]] * cost_maternal_ART)
  
  mat_costs_a2 <- c(0,
                   0,
                   person_weeks[[3,"ahiv_pos_art"]] * cost_maternal_ART,
                   0,
                   0,
                   person_weeks[[3,"chiv_pos_art"]] * cost_maternal_ART)
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_2 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_2 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Second test type not defined")
  }

  mat_costs_a3 <- c(0,
                   maternal_df[[deliveryweek, "ahiv_pos_nart"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                     (1-second_ANC_test_prob) * 
                     (test_sens_e * p_results * (cost_infant_prophylaxis*p_arv)),
                   person_weeks[[4,"ahiv_pos_art"]] * (cost_infant_prophylaxis * p_arv + cost_maternal_ART),
                   maternal_df[[deliveryweek, "chiv_pos_uk"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                     (1-second_ANC_test_prob) *
                     (test_sens_c * p_results * (cost_infant_prophylaxis*p_arv)),
                   person_weeks[[4, "chiv_pos_nart"]] * p_arv * cost_infant_prophylaxis,
                   person_weeks[[4,"chiv_pos_art"]] * (cost_maternal_ART + p_arv * cost_infant_prophylaxis))
  
  mat_costs_a4 <- c(0,
                   0,
                   person_weeks[[5,"ahiv_pos_art"]] * cost_maternal_ART,
                   0,
                   0,
                   person_weeks[[5,"chiv_pos_art"]] * cost_maternal_ART)
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_pp == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_pp == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Postpartum test type not defined")
  }
  
  mat_costs_a5 <- c(0,
                   0,
                   person_weeks[[6,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+early_pp_visit, "chiv_pos_uk"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                     (1-(second_ANC_test_prob+delivery_test_prob)) *
                     (test_sens_c * p_results * (cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[6,"chiv_pos_art"]] * cost_maternal_ART)
  
  
  mat_costs_a6 <- c(0,
                   maternal_df[[deliveryweek+pp_14weeks_visit, "ahiv_pos_nart"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                     (test_sens_e * p_results * (cost_infant_prophylaxis*p_arv)),
                   person_weeks[[7,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+pp_14weeks_visit, "chiv_pos_uk"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                     (test_sens_c * p_results * (cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[7,"chiv_pos_art"]] * cost_maternal_ART)
  
  mat_costs_a7 <- c(0,
                   maternal_df[[deliveryweek+pp_mid_visit, "ahiv_pos_nart"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                     (test_sens_e * p_results * (cost_infant_prophylaxis*p_arv)),
                   person_weeks[[8,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+pp_mid_visit, "chiv_pos_uk"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                     (test_sens_c * p_results * (cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[8,"chiv_pos_art"]] * cost_maternal_ART)
  
  
  mat_costs_a8 <- c(0,
                   maternal_df[[deliveryweek+pp_9months_visit, "ahiv_pos_nart"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                     (test_sens_e * p_results * (cost_infant_prophylaxis*p_arv)),
                   person_weeks[[9,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+pp_9months_visit, "chiv_pos_uk"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                     (test_sens_c * p_results * (cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[9,"chiv_pos_art"]] * cost_maternal_ART)
  
  
  
  
  # PrEP costs
  # Getting maternal costs for model by stage
  mat_costs_p0 <- c(cost_maternal_PrEP * PrEP * person_weeks[[1,2]], 
                   cost_maternal_PrEP * PrEP * person_weeks[[1,3]], 
                   rep(0, 4))
  
  mat_costs_p1 <- c(cost_maternal_PrEP * PrEP * person_weeks[[2,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[2,3]],
                   0,
                   0,
                   0,
                   0)
  
  mat_costs_p2 <- c(cost_maternal_PrEP * PrEP * person_weeks[[3,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[3,3]],
                   0,
                   0,
                   0,
                   0)
  
  mat_costs_p3 <- c(cost_maternal_PrEP * PrEP * person_weeks[[4,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[4,3]],
                   0,
                   0,
                   0,
                   0)
  
  mat_costs_p4 <- c(person_weeks[[5, "hiv_neg"]] * cost_maternal_PrEP * PrEP,
                   person_weeks[[5, "ahiv_pos_nart"]] * cost_maternal_PrEP * PrEP,
                   0,
                   0,
                   0,
                   0)
  
  mat_costs_p5 <- c(cost_maternal_PrEP * PrEP * person_weeks[[6,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[6,3]],
                   0,
                   0,
                   0,
                   0)
  
  mat_costs_p6 <- c(cost_maternal_PrEP * PrEP * person_weeks[[7,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[7,3]],
                   0,
                   0,
                   0,
                   0)
  
  mat_costs_p7 <- c(cost_maternal_PrEP * PrEP * person_weeks[[8,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[8,3]],
                   0,
                   0,
                   0,
                   0)

  mat_costs_p8 <- c(cost_maternal_PrEP * PrEP * person_weeks[[9,2]],
                   cost_maternal_PrEP * PrEP * person_weeks[[9,3]],
                   0,
                   0,
                   0,
                   0)
  

  
  # Full costs - all in one (not using in this model)
  # Getting maternal costs for model by stage
  mat_costs_0 <- c(cost_maternal_PrEP * PrEP * person_weeks[[1,2]], 
                   cost_maternal_PrEP * PrEP * person_weeks[[1,3]], 
                   rep(0, 4))
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_1 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_1 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("First test type not defined")
  }
  

  mat_costs_1 <- c(maternal_df[[weekofanc1, "hiv_neg"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[2,2]],
                   maternal_df[[weekofanc1, "ahiv_pos_nart"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                     (test_cost + test_sens_e * p_results * cost_true_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[2,3]],
                   person_weeks[[2,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[weekofanc1, "chiv_pos_uk"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                     (test_cost + test_sens_c * p_results * cost_true_pos_screening),
                   0,
                   person_weeks[[2,"chiv_pos_art"]] * cost_maternal_ART)
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_2 == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_2 == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Second test type not defined")
  }
  
  mat_costs_2 <- c(maternal_df[[weekofanc2, "hiv_neg"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[3,2]],
                   maternal_df[[weekofanc2, "ahiv_pos_nart"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                     (test_cost + test_sens_e * p_results * cost_true_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[3,3]],
                   person_weeks[[3,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[weekofanc2, "chiv_pos_uk"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                     (test_cost + test_sens_c * p_results * cost_true_pos_screening),
                   0,
                   person_weeks[[3,"chiv_pos_art"]] * cost_maternal_ART)
  
  mat_costs_3 <- c(maternal_df[[deliveryweek, "hiv_neg"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                     (1-second_ANC_test_prob) *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[4,2]],
                   maternal_df[[deliveryweek, "ahiv_pos_nart"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                     (1-second_ANC_test_prob) * 
                     (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                     cost_maternal_PrEP * PrEP * person_weeks[[4,3]],
                   person_weeks[[4,"ahiv_pos_art"]] * (cost_infant_prophylaxis * p_arv + cost_maternal_ART),
                   maternal_df[[deliveryweek, "chiv_pos_uk"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                     (1-second_ANC_test_prob) *
                     (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                   person_weeks[[4, "chiv_pos_nart"]] * p_arv * cost_infant_prophylaxis,
                   person_weeks[[4,"chiv_pos_art"]] * (cost_maternal_ART + p_arv * cost_infant_prophylaxis))
  
  mat_costs_4 <- c(person_weeks[[5, "hiv_neg"]] * cost_maternal_PrEP * PrEP,
                   person_weeks[[5, "ahiv_pos_nart"]] * cost_maternal_PrEP * PrEP,
                   person_weeks[[5,"ahiv_pos_art"]] * cost_maternal_ART,
                   0,
                   0,
                   person_weeks[[5,"chiv_pos_art"]] * cost_maternal_ART)
  
  # Setting test sensitivity and cost based on test type and stage
  if (test_type_pp == "rdt") {
    test_sens_e <- rdt_sens_hiv_e
    test_sens_c <- rdt_sens_hiv_c
    test_cost <- cost_3gen
    test_spec <- spec3
  } else if (test_type_pp == "dual") {
    test_sens_e <- dual.sens.hiv.e
    test_sens_c <- dual.sens.hiv.c
    test_cost <- cost_dual
    test_spec <- hiv.spec.dual
  } else {
    print("Postpartum test type not defined")
  }
  
  mat_costs_5 <- c(maternal_df[[deliveryweek+early_pp_visit, "hiv_neg"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                     (1-(second_ANC_test_prob+delivery_test_prob)) *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[6,2]],
                   maternal_df[[deliveryweek+early_pp_visit, "ahiv_pos_nart"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                     (1-(second_ANC_test_prob+delivery_test_prob)) * 
                     (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                     cost_maternal_PrEP * PrEP * person_weeks[[6,3]],
                   person_weeks[[6,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+early_pp_visit, "chiv_pos_uk"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                     (1-(second_ANC_test_prob+delivery_test_prob)) *
                     (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[6,"chiv_pos_art"]] * cost_maternal_ART)
  
  
  mat_costs_6 <- c(maternal_df[[deliveryweek+pp_14weeks_visit, "hiv_neg"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[7,2]],
                   maternal_df[[deliveryweek+pp_14weeks_visit, "ahiv_pos_nart"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                     (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                     cost_maternal_PrEP * PrEP * person_weeks[[7,3]],
                   person_weeks[[7,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+pp_14weeks_visit, "chiv_pos_uk"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                     (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[7,"chiv_pos_art"]] * cost_maternal_ART)
  
  # There is one cost here different than the Excel models
    # the excel models have the total maternal costs for acute HIV- women as $0 but I think there's a linking error
  mat_costs_7 <- c(maternal_df[[deliveryweek+pp_mid_visit, "hiv_neg"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[8,2]],
                   maternal_df[[deliveryweek+pp_mid_visit, "ahiv_pos_nart"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                     (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                     cost_maternal_PrEP * PrEP * person_weeks[[8,3]],
                   person_weeks[[8,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+pp_mid_visit, "chiv_pos_uk"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                     (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[8,"chiv_pos_art"]] * cost_maternal_ART)
  
  
  mat_costs_8 <- c(maternal_df[[deliveryweek+pp_9months_visit, "hiv_neg"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                     (test_cost + (1-test_spec) * cost_false_pos_screening) +
                     cost_maternal_PrEP * PrEP * person_weeks[[9,2]],
                   maternal_df[[deliveryweek+pp_9months_visit, "ahiv_pos_nart"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                     (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                     cost_maternal_PrEP * PrEP * person_weeks[[9,3]],
                   person_weeks[[9,"ahiv_pos_art"]] * cost_maternal_ART,
                   maternal_df[[deliveryweek+pp_9months_visit, "chiv_pos_uk"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                     (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                   0,
                   person_weeks[[9,"chiv_pos_art"]] * cost_maternal_ART)
  
  # Putting all maternal costs into a df 
  maternal_costs <- data.frame(mat_costs_0, mat_costs_1, mat_costs_2, mat_costs_3, mat_costs_4, 
                               mat_costs_5, mat_costs_6, mat_costs_7, mat_costs_8) %>%
    t() %>%
    bind_cols(data.frame("stage" = person_weeks[,1])) %>%
    relocate("stage")
  
  
  # Putting all testing costs into a df
  maternal_costs_td <- data.frame(mat_costs_td0, mat_costs_td1, mat_costs_td2, mat_costs_td3, mat_costs_td4, 
               mat_costs_td5, mat_costs_td6, mat_costs_td7, mat_costs_td8) %>%
    t() %>%
    bind_cols(data.frame("stage" = person_weeks[,1])) %>%
    relocate("stage")
  
  maternal_costs_tr <- data.frame(mat_costs_tr0, mat_costs_tr1, mat_costs_tr2, mat_costs_tr3, mat_costs_tr4, 
                                  mat_costs_tr5, mat_costs_tr6, mat_costs_tr7, mat_costs_tr8) %>%
    t() %>%
    bind_cols(data.frame("stage" = person_weeks[,1])) %>%
    relocate("stage")
  
  # total testing costs
  maternal_costs_t <- (maternal_costs_td[ , 2:ncol(maternal_costs_td)] + maternal_costs_tr[ , 2:ncol(maternal_costs_tr)]) %>%
    bind_cols(maternal_costs_td[,1]) %>%
    relocate(last_col())
  
  
  # Putting treatment costs into df
  maternal_costs_a <- data.frame(mat_costs_a0, mat_costs_a1, mat_costs_a2, mat_costs_a3, mat_costs_a4, 
                                  mat_costs_a5, mat_costs_a6, mat_costs_a7, mat_costs_a8) %>%
    t() %>%
    bind_cols(data.frame("stage" = person_weeks[,1])) %>%
    relocate("stage")
  
  # Putting PrEP costs into df 
  maternal_costs_p <- data.frame(mat_costs_p0, mat_costs_p1, mat_costs_p2, mat_costs_p3, mat_costs_p4, 
                                  mat_costs_p5, mat_costs_p6, mat_costs_p7, mat_costs_p8) %>%
    t() %>%
    bind_cols(data.frame("stage" = person_weeks[,1])) %>%
    relocate("stage")

  colnames(maternal_costs) <- colnames(person_weeks)
  colnames(maternal_costs_td) <- colnames(person_weeks)
  colnames(maternal_costs_tr) <- colnames(person_weeks)
  colnames(maternal_costs_t) <- colnames(person_weeks)
  colnames(maternal_costs_a) <- colnames(person_weeks)
  colnames(maternal_costs_p) <- colnames(person_weeks)

  #return(maternal_costs) 
  return(list("maternal_costs" = maternal_costs, "dual_testing_costs" = maternal_costs_td,
              "rdt_testing_costs" = maternal_costs_tr, "total_testing_costs" = maternal_costs_t, 
              "treatment_costs" = maternal_costs_a, "prep_costs" = maternal_costs_p))

}

#total_maternal_costs <- sum(maternal_costs[ ,2:ncol(maternal_costs)])

## Getting total maternal costs
tot_mat_costs_func <- function(maternal_costs_list) {
  maternal_costs_df <- maternal_costs_list$maternal_costs
  out <- sum(maternal_costs_df[ ,2:ncol(maternal_costs_df)])
}














################################
## Number of infected infants ##
################################


inf_infect_func <- function(person_weeks_df) {
  
  person_weeks <- person_weeks_df
  
  # aligning input parameter names
  viral_load_suppression <- Vred

  # Getting infected infants from model
  inf_infants_0 <- c(0,
                     ifelse((1-exp(-inutero_recent*(weekofanc1-1)))*(person_weeks[[1,3]]/(weekofanc1-1))<0, 0, 
                            (1-exp(-inutero_recent*(weekofanc1-1)))*(person_weeks[[1,3]]/(weekofanc1-1))),
                     ifelse((1-exp(-inutero_recent*(weekofanc1-1)))*(person_weeks[[1,4]]/(weekofanc1-1))<0, 0, 
                            (1-exp(-inutero_recent*(weekofanc1-1)))*(person_weeks[[1,4]]/(weekofanc1-1))),
                     ifelse((1-exp(-inutero_established*(weekofanc1-1)))*(person_weeks[[1,5]]/(weekofanc1-1))<0, 0, 
                            (1-exp(-inutero_established*(weekofanc1-1)))*(person_weeks[[1,5]]/(weekofanc1-1))),
                     ifelse((1-exp(-inutero_established*(weekofanc1-1)))*(person_weeks[[1,6]]/(weekofanc1-1))<0, 0, 
                            (1-exp(-inutero_established*(weekofanc1-1)))*(person_weeks[[1,6]]/(weekofanc1-1))),
                     ifelse((1-exp(-(inutero_established*(1-viral_load_suppression*p_VL)*(weekofanc1-1))))*
                               (person_weeks[[1,7]]/(weekofanc1-1))<0, 0, 
                            (1-exp(-(inutero_established*(1-viral_load_suppression*p_VL)*(weekofanc1-1))))*
                              (person_weeks[[1,7]]/(weekofanc1-1)))
                     )
  
  
  inf_infants_1 <- c(0,
                     
                     ifelse((1-exp(-inutero_recent*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,3]]/(weekofanc2-weekofanc1))-inf_infants_0[2])<0, 0, 
                            (1-exp(-inutero_recent*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,3]]/(weekofanc2-weekofanc1))-inf_infants_0[2])),
                     
                     ifelse((1-exp(-inutero_recent*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,4]]/(weekofanc2-weekofanc1))-inf_infants_0[3])<0, 0, 
                            (1-exp(-inutero_recent*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,4]]/(weekofanc2-weekofanc1))-inf_infants_0[3])),
                     
                     ifelse((1-exp(-inutero_established*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,5]]/(weekofanc2-weekofanc1))-inf_infants_0[4])<0, 0, 
                            (1-exp(-inutero_established*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,5]]/(weekofanc2-weekofanc1))-inf_infants_0[4])),
                     
                     ifelse((1-exp(-inutero_established*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,6]]/(weekofanc2-weekofanc1))-inf_infants_0[5])<0, 0, 
                            (1-exp(-inutero_established*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,6]]/(weekofanc2-weekofanc1))-inf_infants_0[5])),
                     
                     ifelse((1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,7]]/(weekofanc2-weekofanc1))-inf_infants_0[6])<0, 0, 
                            (1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(weekofanc2-weekofanc1)))*
                              ((person_weeks[[2,7]]/(weekofanc2-weekofanc1))-inf_infants_0[6]))
                      )
  
  
  inf_infants_2 <- c(0,
                     
                     ifelse((1-exp(-inutero_recent*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,3]]/(deliveryweek-weekofanc2))-inf_infants_1[2])<0, 0, 
                            (1-exp(-inutero_recent*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,3]]/(deliveryweek-weekofanc2))-inf_infants_1[2])),
                     
                     ifelse((1-exp(-inutero_recent*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,4]]/(deliveryweek-weekofanc2))-inf_infants_1[3])<0, 0, 
                            (1-exp(-inutero_recent*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,4]]/(deliveryweek-weekofanc2))-inf_infants_1[3])),
                     
                     ifelse((1-exp(-inutero_established*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,5]]/(deliveryweek-weekofanc2))-inf_infants_1[4])<0, 0, 
                            (1-exp(-inutero_established*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,5]]/(deliveryweek-weekofanc2))-inf_infants_1[4])),
                     
                     ifelse((1-exp(-inutero_established*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,6]]/(deliveryweek-weekofanc2))-inf_infants_1[5])<0, 0, 
                            (1-exp(-inutero_established*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,6]]/(deliveryweek-weekofanc2))-inf_infants_1[5])),
                     
                     ifelse((1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,7]]/(deliveryweek-weekofanc2))-inf_infants_1[6])<0, 0, 
                            (1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(deliveryweek-weekofanc2)))*
                              ((person_weeks[[3,7]]/(deliveryweek-weekofanc2))-inf_infants_1[6]))
                      )
  
  
  inf_infants_3 <- c(0,
                     
                     ifelse((1-exp(-early_pp_recent))*((person_weeks[[4,3]])-inf_infants_2[2])*(1-mort_n)<0, 0, 
                            (1-exp(-early_pp_recent))*((person_weeks[[4,3]])-inf_infants_2[2])*(1-mort_n)),
                     
                     ifelse((1-exp(-early_pp_recent*(1-arv_red*p_arv)))*
                              (person_weeks[[4,4]]-inf_infants_2[3])*(1-mort_n)<0, 0, 
                            (1-exp(-early_pp_recent*(1-arv_red*p_arv)))*
                              (person_weeks[[4,4]]-inf_infants_2[3])*(1-mort_n)),
                     
                     ifelse((1-exp(-early_pp_established))*((person_weeks[[4,5]])-inf_infants_2[4])*(1-mort_n)<0, 0, 
                            (1-exp(-early_pp_established))*((person_weeks[[4,5]])-inf_infants_2[4])*(1-mort_n)),
                     
                     ifelse((1-exp(-early_pp_established*(1-arv_red*p_arv)))*
                              (person_weeks[[4,6]]-inf_infants_2[5])*(1-mort_n)<0, 0, 
                            (1-exp(-early_pp_established*(1-arv_red*p_arv)))*
                              (person_weeks[[4,6]]-inf_infants_2[5])*(1-mort_n)),
                     
                     ifelse((1-exp(-(early_pp_established*(1-(viral_load_suppression*p_VL)))*(1-arv_red*p_arv)))*
                              ((person_weeks[[4,7]])-inf_infants_2[6])*(1-mort_n)<0, 0, 
                            (1-exp(-(early_pp_established*(1-(viral_load_suppression*p_VL)))*(1-arv_red*p_arv)))*
                              ((person_weeks[[4,7]])-inf_infants_2[6])*(1-mort_n))
                      )
  
  
  
  inf_infants_4 <- c(0,
                     
                     ifelse((1-exp(-early_pp_recent*(1-nbf*prob_nbf_neg_early_pp)*(early_pp_visit-1)))*
                                (person_weeks[[5,3]]/(early_pp_visit-1)-inf_infants_3[2])*(1-mort_n) < 0, 0, 
                            (1-exp(-early_pp_recent*(1-nbf*prob_nbf_neg_early_pp)*(early_pp_visit-1)))*
                              (person_weeks[[5,3]]/(early_pp_visit-1)-inf_infants_3[2])*(1-mort_n)),
                     
                     ifelse((1-exp(-early_pp_recent*(1-nbf*prob_nbf_pos_early_pp)*(1-arv_red*p_arv)*(early_pp_visit-1)))*
                              (person_weeks[[5,4]]/(early_pp_visit-1)-inf_infants_3[3])*(1-mort_n) < 0, 0, 
                            (1-exp(-early_pp_recent*(1-nbf*prob_nbf_pos_early_pp)*(1-arv_red*p_arv)*(early_pp_visit-1)))*
                              (person_weeks[[5,4]]/(early_pp_visit-1)-inf_infants_3[3])*(1-mort_n)),
                     
                     ifelse((1-exp(-early_pp_established*(1-nbf*prob_nbf_neg_early_pp)*(early_pp_visit-1)))*
                               (person_weeks[[5,5]]/(early_pp_visit-1)-inf_infants_3[4])*(1-mort_n) < 0, 0, 
                             (1-exp(-early_pp_established*(1-nbf*prob_nbf_neg_early_pp)*(early_pp_visit-1)))*
                               (person_weeks[[5,5]]/(early_pp_visit-1)-inf_infants_3[4])*(1-mort_n)),
                     
                     ifelse( (1-exp(-early_pp_established*(1-arv_red*p_arv)*(1-nbf*prob_nbf_pos_early_pp)*(early_pp_visit-1))) *
                              (person_weeks[[5,6]]/(early_pp_visit-1)-inf_infants_3[5]) * (1-mort_n)  < 0, 0, 
                             (1-exp(-early_pp_established*(1-arv_red*p_arv)*(1-nbf*prob_nbf_pos_early_pp)*(early_pp_visit-1))) *
                               (person_weeks[[5,6]]/(early_pp_visit-1)-inf_infants_3[5]) * (1-mort_n)),
  
                     ifelse( (1-exp(-(early_pp_established*(1-arv_red*p_arv)*(1-nbf*prob_nbf_pos_early_pp)*
                                        (1-viral_load_suppression*p_VL)*(early_pp_visit-1))))*
                               (person_weeks[[5,7]]/(early_pp_visit-1)-inf_infants_3[6]) * (1-mort_n)< 0, 0, 
                             (1-exp(-(early_pp_established*(1-arv_red*p_arv)*(1-nbf*prob_nbf_pos_early_pp)*
                                        (1-viral_load_suppression*p_VL)*(early_pp_visit-1))))*
                               (person_weeks[[5,7]]/(early_pp_visit-1)-inf_infants_3[6]) * (1-mort_n))
  )
  
  
  num_weeks <- pp_14weeks_visit-early_pp_visit
  risk_reduction_pos <- (1-nbf*prob_nbf_pos_mid_pp)
  risk_reduction_neg <- (1-nbf*prob_nbf_neg_mid_pp)
  transm_risk <- c(mid_pp_recent*risk_reduction_neg,
                  mid_pp_recent*risk_reduction_pos,
                  mid_pp_established*risk_reduction_neg,
                  rep(mid_pp_established*risk_reduction_pos, 2))
  person_weeks_stage <- person_weeks[6,3:7]
  prev_inf_infants <- inf_infants_4[2:6]
  mort <- mort_i
  
  inf_infants_5 <- c(0,
                     
                     max((1-exp(-transm_risk[1]*num_weeks))*
                              (person_weeks_stage[[1]]/(num_weeks)-prev_inf_infants[1])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[2]*(1-arv_red*p_arv)*num_weeks))*
                              (person_weeks_stage[[2]]/num_weeks-prev_inf_infants[2])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[3]*num_weeks))*
                              (person_weeks_stage[[3]]/num_weeks-prev_inf_infants[3])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[4]*(1-arv_red*p_arv)*num_weeks)) *
                               (person_weeks_stage[[4]]/num_weeks-prev_inf_infants[4]) * (1-mort), 0),
                     
                     max((1-exp(-(transm_risk[5]*(1-arv_red*p_arv)*(1-viral_load_suppression*p_VL)*num_weeks)))*
                               (person_weeks_stage[[5]]/num_weeks-prev_inf_infants[5]) * (1-mort), 0)
  )
  
  
  
  num_weeks <- pp_mid_visit - pp_14weeks_visit
  risk_reduction_pos <- (1-nbf*prob_nbf_pos_mid_pp)
  risk_reduction_neg <- (1-nbf*prob_nbf_neg_mid_pp)
  transm_risk <- c(mid_pp_recent*risk_reduction_neg,
                   mid_pp_recent*risk_reduction_pos,
                   mid_pp_established*risk_reduction_neg,
                   rep(mid_pp_established*risk_reduction_pos, 2))
  person_weeks_stage <- person_weeks[7,3:7]
  prev_inf_infants <- inf_infants_5[2:6]
  
  inf_infants_6 <- c(0,
                     
                     max((1-exp(-transm_risk[1]*num_weeks))*
                           (person_weeks_stage[[1]]/(num_weeks)-prev_inf_infants[1])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[2]*(1-arv_red*p_arv)*num_weeks))*
                           (person_weeks_stage[[2]]/num_weeks-prev_inf_infants[2])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[3]*num_weeks))*
                           (person_weeks_stage[[3]]/num_weeks-prev_inf_infants[3])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[4]*(1-arv_red*p_arv)*num_weeks)) *
                           (person_weeks_stage[[4]]/num_weeks-prev_inf_infants[4]) * (1-mort), 0),
                     
                     max((1-exp(-(transm_risk[5]*(1-arv_red*p_arv)*(1-viral_load_suppression*p_VL)*num_weeks)))*
                           (person_weeks_stage[[5]]/num_weeks-prev_inf_infants[5]) * (1-mort), 0)
  )
  
  
  num_weeks <- pp_9months_visit - pp_mid_visit 
  risk_reduction_pos <- (1-nbf*prob_nbf_pos_late_pp)
  risk_reduction_neg <- (1-nbf*prob_nbf_neg_late_pp)
  transm_risk <- c(late_pp_recent*risk_reduction_neg,
                   late_pp_recent*risk_reduction_pos,
                   late_pp_established*risk_reduction_neg,
                   rep(late_pp_established*risk_reduction_pos, 2))
  person_weeks_stage <- person_weeks[8,3:7]
  prev_inf_infants <- inf_infants_6[2:6]
  mort <- mort_i
  
  inf_infants_7 <- c(0,
                     
                     max((1-exp(-transm_risk[1]*num_weeks))*
                           (person_weeks_stage[[1]]/(num_weeks)-prev_inf_infants[1])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[2]*(1-arv_red*p_arv)*num_weeks))*
                           (person_weeks_stage[[2]]/num_weeks-prev_inf_infants[2])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[3]*num_weeks))*
                           (person_weeks_stage[[3]]/num_weeks-prev_inf_infants[3])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[4]*(1-arv_red*p_arv)*num_weeks)) *
                           (person_weeks_stage[[4]]/num_weeks-prev_inf_infants[4]) * (1-mort), 0),
                     
                     max((1-exp(-(transm_risk[5]*(1-arv_red*p_arv)*(1-viral_load_suppression*p_VL)*num_weeks)))*
                           (person_weeks_stage[[5]]/num_weeks-prev_inf_infants[5]) * (1-mort), 0)
  )
  
  
  num_weeks <- end_of_model - pp_9months_visit + 1
  risk_reduction_pos <- (1-nbf*prob_nbf_pos_late_pp)
  risk_reduction_neg <- (1-nbf*prob_nbf_neg_late_pp)
  transm_risk <- c(late_pp_recent*risk_reduction_neg,
                   late_pp_recent*risk_reduction_pos,
                   late_pp_established*risk_reduction_neg,
                   rep(late_pp_established*risk_reduction_pos, 2))
  person_weeks_stage <- person_weeks[9,3:7]
  prev_inf_infants <- inf_infants_7[2:6]
  mort <- mort_i
  
  inf_infants_8 <- c(0,
                     
                     max((1-exp(-transm_risk[1]*num_weeks))*
                           (person_weeks_stage[[1]]/(num_weeks)-prev_inf_infants[1])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[2]*(1-arv_red*p_arv)*num_weeks))*
                           (person_weeks_stage[[2]]/num_weeks-prev_inf_infants[2])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[3]*num_weeks))*
                           (person_weeks_stage[[3]]/num_weeks-prev_inf_infants[3])*(1-mort), 0),
                     
                     max((1-exp(-transm_risk[4]*(1-arv_red*p_arv)*num_weeks)) *
                           (person_weeks_stage[[4]]/num_weeks-prev_inf_infants[4]) * (1-mort), 0),
                     
                     max((1-exp(-(transm_risk[5]*(1-arv_red*p_arv)*(1-viral_load_suppression*p_VL)*num_weeks)))*
                           (person_weeks_stage[[5]]/num_weeks-prev_inf_infants[5]) * (1-mort), 0)
  )
  
  
  infected_infants <- data.frame(inf_infants_0, inf_infants_1, inf_infants_2, inf_infants_3, inf_infants_4,
                                 inf_infants_5, inf_infants_6, inf_infants_7, inf_infants_8)
  
  return(infected_infants)

}












#####################
## Infant outcomes ## 
#####################



infant_outcomes_func <- function(infected_infants_df, maternal_df) {
  
  life_table <- lifetable_output()
  
  # Total number of infected infants
  total_infected_infants <- sum(infected_infants_df)
  
  ## Infant treatment costs
  # Getting infant treatment costs for Stage 0 - 3: conception to deliveryweek
  infant_treatment_costs_0_3 <- infected_infants_df[1:4] * inf_ART_coverage * (cost_infant_ART_after_2weeks*52) / (1+disc_rate)
  
  # Getting infant treatment costs for Stage 4-8: post-partum period
  infant_treatment_costs_4_8 <- matrix(nrow = 6, ncol = 5,
                                       dimnames = list(NULL,
                                                       c("inf_infants_4", 
                                                         "inf_infants_5",
                                                         "inf_infants_6",
                                                         "inf_infants_7",
                                                         "inf_infants_8")))
  
  infant_treatment_costs_4_8 <- as.data.frame(infant_treatment_costs_4_8)
  
  visit_time <- c(early_pp_visit, pp_14weeks_visit, pp_mid_visit, pp_9months_visit, end_of_model)
  
  for (x in 5:9) {
    infant_treatment_costs_4_8[,(x-4)] <- infected_infants_df[x] * inf_ART_coverage * 
      (cost_infant_ART_after_2weeks*52*(1-(visit_time[x-4]/end_of_model))) / (1+disc_rate) 
  }
  
  # Combining all infant treatment costs into a single data frame
  infant_treatment_costs <- cbind(infant_treatment_costs_0_3, infant_treatment_costs_4_8)
  total_infant_treatment_costs <- sum(infant_treatment_costs)
  
  surv_nhiv <- 1 - life_table$annual_cond_prob_nonhiv_death[1]
  
  ## Getting total number of infant deaths
  infant_deaths <- total_infected_infants * life_table[1,4] + 
    (population_size - sum(maternal_df$deaths[1:(deliveryweek-1)]) - total_infected_infants)*(1-surv_nhiv)
  
  return(list(total_infected_infants = total_infected_infants,
              total_infant_treatment_costs = total_infant_treatment_costs,
              infant_deaths = infant_deaths))

}
































