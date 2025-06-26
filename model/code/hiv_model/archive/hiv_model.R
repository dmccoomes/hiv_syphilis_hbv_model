
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


# Starting states
neg <- (1-hiv_prev) - (inc_preg_early*duration_infection)
recent <- duration_infection*inc_preg_early
established_uk_nart <- hiv_prev*(1-status_known)
established_k_nart <- hiv_prev*status_known*(1-p_ART)
established_art <- hiv_prev*status_known*p_ART

 
# Setting up empty matrix
maternal_matrix <- matrix(nrow = 92, ncol = 9, 
                     dimnames = list(NULL, 
                                     c("weeks", "hiv_neg", "ahiv_pos_nart", 
                                       "ahiv_pos_art", "chiv_pos_uk", 
                                       "chiv_pos_nart", "chiv_pos_art", "deaths",
                                       "stage")))

maternal_df <- as.data.frame(maternal_matrix)

# Setting up first week parameters
maternal_df[1, ] <- c(1, population_size*neg, population_size*recent, 0, population_size*established_uk_nart, 
                 population_size*established_k_nart, population_size*established_art, 0, 0)

maternal_df$stage <- as.character(maternal_df$stage) 
maternal_df$stage[1] <- "Stage 0: Onset of pregnancy to first ANC"


###########
# STAGE 0 #
###########


# Getting estimates for stage 0: onset of pregnancy to first ANC
for (x in 2:(first_ANC-1)) {
  
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
for (x in first_ANC:(second_ANC-1)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- ifelse(maternal_df$weeks[x-1] == (first_ANC-1), 
                              maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_early * (1-(PrEP*red_prep)))),
                              maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_late * (1-(PrEP*red_prep))))
                              )
  
  maternal_df$ahiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (first_ANC-1),
                                    maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_early * (1-(PrEP*red_prep)) +
                                      maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * 
                                        (1-first_ANC_test*att_firstANC*(1-stockout)*test_accept*test_sens_e*p_results*p_ART) +
                                      maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor),
                                    maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_late * (1-(PrEP*red_prep)) +
                                      maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) +
                                      maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
                                    )
    
  maternal_df$ahiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (first_ANC-1),
                                   maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * 
                                    first_ANC_test * att_firstANC * test_accept * (1-stockout) * test_sens_e * p_results * p_ART,
                                   maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor
                                   )
    
  maternal_df$chiv_pos_uk[x] <- ifelse(maternal_df$weeks[x-1] == (first_ANC-1),
                                       maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                         (1-(first_ANC_test * att_firstANC * (1-stockout) * test_accept * test_sens_c * p_results)) + 
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                         (1-(first_ANC_test * att_firstANC * (1-stockout) * test_accept * test_sens_c * p_results)),
                                       maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females)
                                       )

  maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (first_ANC-1),
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * first_ANC_test * 
                                           att_firstANC * test_accept * (1-stockout) * p_results * test_sens_c * (1-p_ART) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * first_ANC_test * att_firstANC * test_accept *
                                           (1-stockout) * p_results * test_sens_c * (1-p_ART) +
                                           maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                           maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor),
                                         maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                           maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
                                         )
  
  maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (first_ANC-1),
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * first_ANC_test *
                                          att_firstANC * test_sens_c * test_accept * (1-stockout) * p_results * p_ART +
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                          maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * first_ANC_test * att_firstANC * test_sens_c *
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


# Getting estimates for Stage 2: Second ANC to delivery
for (x in second_ANC:(delivery-1)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_preg_late * (1-(PrEP*red_prep))))

  maternal_df$ahiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (second_ANC - 1),
                                         maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_late * (1-(PrEP*red_prep)) +
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * 
                                           (1-late_ANC_test*att_lategest*(1-stockout)*test_accept*test_sens_e*p_results*p_ART) +
                                           maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor),
                                         maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_preg_late * (1-(PrEP*red_prep)) +
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) +
                                           maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
                                         )
  
    maternal_df$ahiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (second_ANC - 1),
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) * late_ANC_test *
                                          att_lategest * test_accept * (1-stockout) * test_sens_e * p_results * p_ART +
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor,
                                        maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor
                                        )

  maternal_df$chiv_pos_uk[x] <- ifelse(maternal_df$weeks[x-1] == (second_ANC - 1),
                                       maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                         (1-(late_ANC_test * att_lategest * (1-stockout) * test_accept * test_sens_c * p_results)) +
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                         (1-(late_ANC_test * att_lategest * (1-stockout) * test_accept * test_sens_c * p_results)),
                                       maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females)
                                       )
  
  maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (second_ANC - 1),
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * late_ANC_test *
                                           att_lategest * test_accept * (1-stockout) * p_results * test_sens_c * (1-p_ART) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * late_ANC_test * att_lategest *
                                           test_accept * (1-stockout) * p_results * test_sens_c * (1-p_ART) +
                                           maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) + 
                                           maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor),
                                         maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                           maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
                                         )
  
  maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (second_ANC - 1),
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * late_ANC_test *
                                          att_lategest * test_sens_c * test_accept * (1-stockout) * p_results * p_ART +
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                          maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * late_ANC_test * att_lategest * test_sens_c *
                                          test_accept * (1-stockout) * p_results * p_ART +
                                          maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor,
                                        maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                          maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
                                        )
  
  maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                    maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
  
  maternal_df$stage[x] <- "Stage 2: Second ANC to delivery" 
  
}



###########
# STAGE 3 #
###########


# Getting estimates for Stage 3: Delivery
for (x in delivery:delivery) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * (1-(inc_preg_late * (1-(PrEP*red_prep)))) 
  
  maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * inc_preg_late * (1-(PrEP*red_prep)) +
                                           maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1-1/duration_infection) * 
                                           (1-late_ANC_test*att_delivery*(1-stockout)*test_accept*(1-second_ANC_test_prob)*test_sens_e*p_results*p_ART) +
                                           maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-(1/duration_infection)) * (1-ART_drop_factor)

  maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1-1/duration_infection) * late_ANC_test *
                                          att_delivery * test_accept * (1-stockout) * (1-second_ANC_test_prob) * test_sens_e * p_results * p_ART +
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-(1/duration_infection)) * ART_drop_factor

  maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) * 
                                         (1-(late_ANC_test * att_delivery * (1-stockout) * test_accept * (1-second_ANC_test_prob) * 
                                         test_sens_c * p_results)) +
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_maternal_perinatal) * 
                                         (1-(late_ANC_test * att_delivery * (1-stockout) * test_accept * (1-second_ANC_test_prob) *
                                         test_sens_c * p_results))

  maternal_df$chiv_pos_nart[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) * late_ANC_test *
                                           att_delivery * test_accept * (1-stockout) * p_results * (1-second_ANC_test_prob) * 
                                           test_sens_c * (1-p_ART) +
                                           maternal_df$chiv_pos_uk[x-1] * (1-mort_maternal_perinatal) * late_ANC_test * att_delivery *
                                           test_accept * (1-stockout) * p_results * (1-second_ANC_test_prob) * test_sens_c * (1-p_ART) +
                                           maternal_df$chiv_pos_nart[x-1] * (1-mort_maternal_perinatal) + 
                                           maternal_df$chiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1-ART_drop_factor)

  maternal_df$chiv_pos_art[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) * late_ANC_test *
                                          att_delivery * test_sens_c * test_accept * (1-stockout) * p_results * p_ART * (1-second_ANC_test_prob) +
                                          maternal_df$ahiv_pos_art[x-1] * (1-mort_maternal_perinatal) * (1/duration_infection) +
                                          maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * late_ANC_test * att_delivery * test_sens_c *
                                          test_accept * (1-stockout) * p_results * p_ART * (1-second_ANC_test_prob) +
                                          maternal_df$chiv_pos_art[x-1] * (1-mort_maternal_perinatal) * ART_drop_factor
  
  maternal_df$deaths[x] <- mort_maternal_perinatal * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                    maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
  
  maternal_df$stage[x] <- "Stage 3: Delivery" 
  
}



###########
# STAGE 4 #
###########


# Getting estimates for Stage 4: 1-6 weeks post-partum
for (x in (delivery+1):(delivery+early_pp_visit-1)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * 
    (1-(ifelse(maternal_df$weeks[x-1] == delivery, inc_preg_late, inc_pp_early) * (1-(PrEP*red_prep)))) 
  
  maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_maternal_perinatal) * 
                                    ifelse(maternal_df$weeks[x-1] == delivery, inc_preg_late, inc_pp_early) * (1-(PrEP*red_prep)) +
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
for (x in (delivery+early_pp_visit):(delivery+pp_14weeks_visit-1)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * 
    (1-(ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, inc_pp_early, inc_pp_mid) * (1-(PrEP*red_prep)))) 
  
  maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * 
                                    ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, inc_pp_early, inc_pp_mid) * (1-(PrEP*red_prep)) +
                                  maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                    ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, (1 - pp_early_test * att_6wk * (1-stockout) * test_accept * (1- (second_ANC_test_prob + delivery_test_prob)) *
                                    test_sens_e * p_results * p_ART), 1) +
                                  maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
  
  maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_early_test * 
                                        att_6wk * test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) * test_sens_e *
                                        p_results * p_ART, 0)
                                 
  maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                  ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                         (1-pp_early_test * att_6wk * (1-stockout) * test_accept * (1-(second_ANC_test_prob + delivery_test_prob)) *
                                         test_sens_c * p_results), 1) +
                                maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                  ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                         (1-pp_early_test * att_6wk * (1-stockout) * test_accept * (1-(second_ANC_test_prob + delivery_test_prob)) * 
                                            test_sens_c * p_results), 1)
 
  maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_early_test *
                                           att_6wk * test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) *
                                           p_results * test_sens_c * (1-p_ART), 0) +
                                  ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_early_test * att_6wk * test_accept * (1-stockout) *
                                           (1-(second_ANC_test_prob + delivery_test_prob)) * p_results * test_sens_c * (1-p_ART), 0) +      
                                  maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                  maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
 
  maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_early_test *
                                          att_6wk * test_sens_c * test_accept * (1-stockout) * (1-(second_ANC_test_prob + delivery_test_prob)) *
                                          p_results * p_ART, 0) +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+early_pp_visit)-1, 
                                        maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_early_test * att_6wk * test_sens_c *
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
for (x in (delivery+pp_14weeks_visit):(delivery+pp_mid_visit-1)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_mid * (1-(PrEP*red_prep)))) 
  
  maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_mid * (1-(PrEP*red_prep)) +
                                  maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                    ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                           (1 - (pp_14weeks_test * att_14wk * (1-stockout) * test_accept * 
                                           test_sens_e * p_results * p_ART)), 1) +
                                  maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)

  maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_14weeks_test * 
                                        att_14wk * test_accept * (1-stockout) * test_sens_e * p_results * p_ART, 0)

  maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                       (1 - pp_14weeks_test * att_14wk * (1-stockout) * test_accept * test_sens_c * p_results), 1) +
                                maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                       (1 - pp_14weeks_test * att_14wk * (1-stockout) * test_accept * test_sens_c * p_results), 1)
  
  
  
  maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_14weeks_test *
                                           att_14wk * test_accept * (1-stockout) * p_results * test_sens_c * (1-p_ART), 0) +
                                  ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_14weeks_test * att_14wk * test_accept * 
                                           (1-stockout) * p_results * test_sens_c * (1-p_ART), 0) +      
                                  maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                  maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
  
  maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_14weeks_test *
                                          att_14wk * test_sens_c * test_accept * (1-stockout) * p_results * p_ART, 0) +
                                  ifelse(maternal_df$weeks[x-1] == (delivery+pp_14weeks_visit)-1, 
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_14weeks_test * att_14wk * test_sens_c *
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
for (x in (delivery+pp_mid_visit):(delivery+pp_9months_visit-1)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit-1),
                                   maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_mid * (1-(PrEP*red_prep)))),
                                   maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_late * (1-(PrEP*red_prep))))) 
  
  
  maternal_df$ahiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit-1),
                                         maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_mid * (1-(PrEP*red_prep)),
                                         maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_late * (1-(PrEP*red_prep))) +
                                  maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                  ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                         (1 - (pp_mid_test * att_6mo * (1-stockout) * test_accept * 
                                                 test_sens_e * p_results * p_ART)), 1) +
                                  maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
  
  maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_mid_test * 
                                          att_6mo * test_accept * (1-stockout) * test_sens_e * p_results * p_ART, 0)
  
  maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                       (1 - pp_mid_test * att_6mo * (1-stockout) * test_accept * test_sens_c * p_results), 1) +
                                maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                       (1 - pp_mid_test * att_6mo * (1-stockout) * test_accept * test_sens_c * p_results), 1)
  
  
  
  maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_mid_test *
                                           att_6mo * test_accept * (1-stockout) * p_results * test_sens_c * (1-p_ART), 0) +
                                  ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_mid_test * att_6mo * test_accept * 
                                           (1-stockout) * p_results * test_sens_c * (1-p_ART), 0) +      
                                  maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                  maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
  
  maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_mid_test *
                                          att_6mo * test_sens_c * test_accept * (1-stockout) * p_results * p_ART, 0) +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+pp_mid_visit)-1, 
                                        maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_mid_test * att_6mo * test_sens_c *
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
for (x in (delivery+pp_9months_visit):(delivery+end_of_model)) {
  
  maternal_df$weeks[x] <- x
  
  maternal_df$hiv_neg[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * (1-(inc_pp_late * (1-(PrEP*red_prep))))
  
  
  maternal_df$ahiv_pos_nart[x] <- maternal_df$hiv_neg[x-1] * (1-mort_adult_females) * inc_pp_late * (1-(PrEP*red_prep)) +
                                  maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-1/duration_infection) *
                                  ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                         (1 - (pp_9months_test * att_9mo * (1-stockout) * test_accept * 
                                                 test_sens_e * p_results * p_ART)), 1) +
                                  maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * (1-ART_drop_factor)
  
  maternal_df$ahiv_pos_art[x] <- maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * ART_drop_factor +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1-(1/duration_infection)) * pp_9months_test * 
                                          att_9mo * test_accept * (1-stockout) * test_sens_e * p_results * p_ART, 0)
  
  maternal_df$chiv_pos_uk[x] <- maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * 
                                ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                       (1 - pp_9months_test * att_9mo * (1-stockout) * test_accept * test_sens_c * p_results), 1) +
                                maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * 
                                ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                       (1 - pp_9months_test * att_9mo * (1-stockout) * test_accept * test_sens_c * p_results), 1)
  
  maternal_df$chiv_pos_nart[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                         maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_9months_test *
                                           att_9mo * test_accept * (1-stockout) * p_results * test_sens_c * (1-p_ART), 0) +
                                  ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                         maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_9months_test * att_9mo * test_accept * 
                                           (1-stockout) * p_results * test_sens_c * (1-p_ART), 0) +      
                                  maternal_df$chiv_pos_nart[x-1] * (1-mort_adult_females) +
                                  maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * (1-ART_drop_factor)
  
  maternal_df$chiv_pos_art[x] <- ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                        maternal_df$ahiv_pos_nart[x-1] * (1-mort_adult_females) * (1/duration_infection) * pp_9months_test *
                                          att_9mo * test_sens_c * test_accept * (1-stockout) * p_results * p_ART, 0) +
                                 ifelse(maternal_df$weeks[x-1] == (delivery+pp_9months_visit)-1, 
                                        maternal_df$chiv_pos_uk[x-1] * (1-mort_adult_females) * pp_9months_test * att_9mo * test_sens_c *
                                          test_accept * (1-stockout) * p_results * p_ART, 0) +
                                 maternal_df$ahiv_pos_art[x-1] * (1-mort_adult_females) * (1/duration_infection) +
                                 maternal_df$chiv_pos_art[x-1] * (1-mort_adult_females) * ART_drop_factor
  
  maternal_df$deaths[x] <- mort_adult_females * sum(maternal_df$hiv_neg[x-1], maternal_df$ahiv_pos_nart[x-1], maternal_df$ahiv_pos_art[x-1],
                                                    maternal_df$chiv_pos_uk[x-1], maternal_df$chiv_pos_nart[x-1], maternal_df$chiv_pos_art[x-1])
  
  maternal_df$stage[x] <- "Stage 8: 9-12 months post-partum" 
  
}

# return(maternal_df)
# 
# }

#################











##### Getting output from model

# Getting total person weeks from model
person_weeks <- maternal_df %>%
  group_by(stage) %>%
  summarize(across(hiv_neg:chiv_pos_art, ~sum(.)))


# Getting maternal costs for model by stage
mat_costs_0 <- c(cost_maternal_PrEP * PrEP * person_weeks[[1,2]], 
                 cost_maternal_PrEP * PrEP * person_weeks[[1,3]], 
                 rep(0, 4))

mat_costs_1 <- c(maternal_df[[first_ANC, "hiv_neg"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[2,2]],
                 maternal_df[[first_ANC, "ahiv_pos_nart"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                   (test_cost + test_sens_e * p_results * cost_true_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[2,3]],
                 person_weeks[[2,"ahiv_pos_art"]] * cost_maternal_ART,
                 maternal_df[[first_ANC, "chiv_pos_uk"]] * att_firstANC * test_accept * (1-stockout) * first_ANC_test *
                   (test_cost + test_sens_c * p_results * cost_true_pos_screening),
                 0,
                 person_weeks[[2,"chiv_pos_art"]] * cost_maternal_ART)

mat_costs_2 <- c(maternal_df[[second_ANC, "hiv_neg"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[3,2]],
                 maternal_df[[second_ANC, "ahiv_pos_nart"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                   (test_cost + test_sens_e * p_results * cost_true_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[3,3]],
                 person_weeks[[3,"ahiv_pos_art"]] * cost_maternal_ART,
                 maternal_df[[second_ANC, "chiv_pos_uk"]] * att_lategest * test_accept * (1-stockout) * late_ANC_test *
                   (test_cost + test_sens_c * p_results * cost_true_pos_screening),
                 0,
                 person_weeks[[3,"chiv_pos_art"]] * cost_maternal_ART)

mat_costs_3 <- c(maternal_df[[delivery, "hiv_neg"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                   (1-second_ANC_test_prob) *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[4,2]],
                 maternal_df[[delivery, "ahiv_pos_nart"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
                   (1-second_ANC_test_prob) * 
                   (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                   cost_maternal_PrEP * PrEP * person_weeks[[4,3]],
                 person_weeks[[4,"ahiv_pos_art"]] * (cost_infant_prophylaxis * p_arv + cost_maternal_ART),
                 maternal_df[[delivery, "chiv_pos_uk"]] * att_delivery * test_accept * (1-stockout) * late_ANC_test *
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


mat_costs_5 <- c(maternal_df[[delivery+early_pp_visit, "hiv_neg"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                   (1-(second_ANC_test_prob+delivery_test_prob)) *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[6,2]],
                 maternal_df[[delivery+early_pp_visit, "ahiv_pos_nart"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                   (1-(second_ANC_test_prob+delivery_test_prob)) * 
                   (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                   cost_maternal_PrEP * PrEP * person_weeks[[6,3]],
                 person_weeks[[6,"ahiv_pos_art"]] * cost_maternal_ART,
                 maternal_df[[delivery+early_pp_visit, "chiv_pos_uk"]] * att_6wk * test_accept * (1-stockout) * pp_early_test *
                   (1-(second_ANC_test_prob+delivery_test_prob)) *
                   (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                 0,
                 person_weeks[[6,"chiv_pos_art"]] * cost_maternal_ART)


mat_costs_6 <- c(maternal_df[[delivery+pp_14weeks_visit, "hiv_neg"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[7,2]],
                 maternal_df[[delivery+pp_14weeks_visit, "ahiv_pos_nart"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                   (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                   cost_maternal_PrEP * PrEP * person_weeks[[7,3]],
                 person_weeks[[7,"ahiv_pos_art"]] * cost_maternal_ART,
                 maternal_df[[delivery+pp_14weeks_visit, "chiv_pos_uk"]] * att_14wk * test_accept * (1-stockout) * pp_14weeks_test *
                   (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                 0,
                 person_weeks[[7,"chiv_pos_art"]] * cost_maternal_ART)

# There is one cost here different than the Excel models
  # the excel models have the total maternal costs for acute HIV- women as $0 but I think there's a linking error
mat_costs_7 <- c(maternal_df[[delivery+pp_mid_visit, "hiv_neg"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[8,2]],
                 maternal_df[[delivery+pp_mid_visit, "ahiv_pos_nart"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                   (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                   cost_maternal_PrEP * PrEP * person_weeks[[8,3]],
                 person_weeks[[8,"ahiv_pos_art"]] * cost_maternal_ART,
                 maternal_df[[delivery+pp_mid_visit, "chiv_pos_uk"]] * att_6mo * test_accept * (1-stockout) * pp_mid_test *
                   (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                 0,
                 person_weeks[[8,"chiv_pos_art"]] * cost_maternal_ART)


mat_costs_8 <- c(maternal_df[[delivery+pp_9months_visit, "hiv_neg"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                   (test_cost + (1-test_spec) * cost_false_pos_screening) +
                   cost_maternal_PrEP * PrEP * person_weeks[[9,2]],
                 maternal_df[[delivery+pp_9months_visit, "ahiv_pos_nart"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                   (test_cost + test_sens_e * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)) +
                   cost_maternal_PrEP * PrEP * person_weeks[[9,3]],
                 person_weeks[[9,"ahiv_pos_art"]] * cost_maternal_ART,
                 maternal_df[[delivery+pp_9months_visit, "chiv_pos_uk"]] * att_9mo * test_accept * (1-stockout) * pp_9months_test *
                   (test_cost + test_sens_c * p_results * (cost_true_pos_screening + cost_infant_prophylaxis*p_arv)),
                 0,
                 person_weeks[[9,"chiv_pos_art"]] * cost_maternal_ART)

# Putting all maternal costs into a df 
maternal_costs <- data.frame(mat_costs_0, mat_costs_1, mat_costs_2, mat_costs_3, mat_costs_4, 
                             mat_costs_5, mat_costs_6, mat_costs_7, mat_costs_8) %>%
  t() %>%
  bind_cols(data.frame("stage" = person_weeks[,1])) %>%
  relocate("stage")

colnames(maternal_costs) <- colnames(person_weeks)

total_maternal_costs <- sum(maternal_costs[ ,2:ncol(maternal_costs)])



# Getting infected infants from model
inf_infants_0 <- c(0,
                   ifelse((1-exp(-inutero_recent*(first_ANC-1)))*(person_weeks[[1,3]]/(first_ANC-1))<0, 0, 
                          (1-exp(-inutero_recent*(first_ANC-1)))*(person_weeks[[1,3]]/(first_ANC-1))),
                   ifelse((1-exp(-inutero_recent*(first_ANC-1)))*(person_weeks[[1,4]]/(first_ANC-1))<0, 0, 
                          (1-exp(-inutero_recent*(first_ANC-1)))*(person_weeks[[1,4]]/(first_ANC-1))),
                   ifelse((1-exp(-inutero_established*(first_ANC-1)))*(person_weeks[[1,5]]/(first_ANC-1))<0, 0, 
                          (1-exp(-inutero_established*(first_ANC-1)))*(person_weeks[[1,5]]/(first_ANC-1))),
                   ifelse((1-exp(-inutero_established*(first_ANC-1)))*(person_weeks[[1,6]]/(first_ANC-1))<0, 0, 
                          (1-exp(-inutero_established*(first_ANC-1)))*(person_weeks[[1,6]]/(first_ANC-1))),
                   ifelse((1-exp(-(inutero_established*(1-viral_load_suppression*p_VL)*(first_ANC-1))))*
                             (person_weeks[[1,7]]/(first_ANC-1))<0, 0, 
                          (1-exp(-(inutero_established*(1-viral_load_suppression*p_VL)*(first_ANC-1))))*
                            (person_weeks[[1,7]]/(first_ANC-1)))
                   )


inf_infants_1 <- c(0,
                   
                   ifelse((1-exp(-inutero_recent*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,3]]/(second_ANC-first_ANC))-inf_infants_0[2])<0, 0, 
                          (1-exp(-inutero_recent*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,3]]/(second_ANC-first_ANC))-inf_infants_0[2])),
                   
                   ifelse((1-exp(-inutero_recent*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,4]]/(second_ANC-first_ANC))-inf_infants_0[3])<0, 0, 
                          (1-exp(-inutero_recent*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,4]]/(second_ANC-first_ANC))-inf_infants_0[3])),
                   
                   ifelse((1-exp(-inutero_established*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,5]]/(second_ANC-first_ANC))-inf_infants_0[4])<0, 0, 
                          (1-exp(-inutero_established*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,5]]/(second_ANC-first_ANC))-inf_infants_0[4])),
                   
                   ifelse((1-exp(-inutero_established*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,6]]/(second_ANC-first_ANC))-inf_infants_0[5])<0, 0, 
                          (1-exp(-inutero_established*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,6]]/(second_ANC-first_ANC))-inf_infants_0[5])),
                   
                   ifelse((1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,7]]/(second_ANC-first_ANC))-inf_infants_0[6])<0, 0, 
                          (1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(second_ANC-first_ANC)))*
                            ((person_weeks[[2,7]]/(second_ANC-first_ANC))-inf_infants_0[6]))
                    )


inf_infants_2 <- c(0,
                   
                   ifelse((1-exp(-inutero_recent*(delivery-second_ANC)))*
                            ((person_weeks[[3,3]]/(delivery-second_ANC))-inf_infants_1[2])<0, 0, 
                          (1-exp(-inutero_recent*(delivery-second_ANC)))*
                            ((person_weeks[[3,3]]/(delivery-second_ANC))-inf_infants_1[2])),
                   
                   ifelse((1-exp(-inutero_recent*(delivery-second_ANC)))*
                            ((person_weeks[[3,4]]/(delivery-second_ANC))-inf_infants_1[3])<0, 0, 
                          (1-exp(-inutero_recent*(delivery-second_ANC)))*
                            ((person_weeks[[3,4]]/(delivery-second_ANC))-inf_infants_1[3])),
                   
                   ifelse((1-exp(-inutero_established*(delivery-second_ANC)))*
                            ((person_weeks[[3,5]]/(delivery-second_ANC))-inf_infants_1[4])<0, 0, 
                          (1-exp(-inutero_established*(delivery-second_ANC)))*
                            ((person_weeks[[3,5]]/(delivery-second_ANC))-inf_infants_1[4])),
                   
                   ifelse((1-exp(-inutero_established*(delivery-second_ANC)))*
                            ((person_weeks[[3,6]]/(delivery-second_ANC))-inf_infants_1[5])<0, 0, 
                          (1-exp(-inutero_established*(delivery-second_ANC)))*
                            ((person_weeks[[3,6]]/(delivery-second_ANC))-inf_infants_1[5])),
                   
                   ifelse((1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(delivery-second_ANC)))*
                            ((person_weeks[[3,7]]/(delivery-second_ANC))-inf_infants_1[6])<0, 0, 
                          (1-exp(-(inutero_established*(1-(viral_load_suppression*p_VL)))*(delivery-second_ANC)))*
                            ((person_weeks[[3,7]]/(delivery-second_ANC))-inf_infants_1[6]))
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

total_infected_infants <- sum(infected_infants)


# Getting infant treatment costs for Stage 0 - 3: conception to delivery
infant_treatment_costs_0_3 <- infected_infants[1:4] * inf_ART_coverage * (cost_infant_ART_after_2weeks*52) / (1+disc_rate)

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
  infant_treatment_costs_4_8[,(x-4)] <- infected_infants[x] * inf_ART_coverage * 
    (cost_infant_ART_after_2weeks*52*(1-(visit_time[x-4]/end_of_model))) / (1+disc_rate) 
}

# Combining all infant treatment costs into a single data frame
infant_treatment_costs <- cbind(infant_treatment_costs_0_3, infant_treatment_costs_4_8)

total_infant_treatment_costs <- sum(infant_treatment_costs)

# Getting total number of infant deaths
infant_deaths <- total_infected_infants * life_table[1,4] + 
  (population_size - sum(maternal_df$deaths[1:(delivery-1)]) - total_infected_infants)*(1-surv_nhiv)


































