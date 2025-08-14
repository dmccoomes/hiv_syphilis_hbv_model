

#### Title: Syphilis retesting model
#### Purpose: Create syphilis retesting model based on Excel models
#### Author: David Coomes
#### Date: April 18, 2025


#### Purpose ####
### This model uses population level parameters to estimate syphilis prevalence among pregnant women
  ## Input: population parameters from spreadsheet
  ## Output: Matrix of syphilis and testing status for population by week - 
    ## From 1-39 weeks during pregnancy


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table"
)



### Setting up syphilis data frame using function

syph_maternal_func <- function(first_ANC_test, late_ANC_test, test_type_1, test_type_2) {
  
  ### This function estimates populations of pregnant people 
  ### based on syphilis testing type and frequency during ANC 
    ### Input: syphilis population parameters, test frequency (1=test, 0=no test), test type
    ### Output: dataframe of population based on syphilis infection, knowledge, and treatment
    ### Choices for test type are 'dual' (dual HIV/syphilis test) or 'rdt' (syphilis RDT)
    ### Choices for testing regimen are 1. early ANC or 2. early and late ANC
  
  
  ## Making sure that our inputs are correct
  # First ANC must be 0 (no test) or 1 (tested)
  if (is.null(first_ANC_test) || !length(first_ANC_test) || !all(first_ANC_test %in% c(0, 1))) {
    stop("'first_ANC_test' argument must be one of: '1' (test),'0' (no test)")
  }
  # Second ANC must be 0 (no test) or 1 (tested)
  if (is.null(late_ANC_test) || !length(late_ANC_test) || !all(late_ANC_test %in% c(0, 1))) {
    stop("'late_ANC_test' argument must be one of: '1' (test),'0' (no test)")
  }
  # First test type must be selected
  if (is.null(test_type_1) || !length(test_type_1) || !all(test_type_1 %in% c("rdt", "lab", "dual", "multi"))) {
    stop("'test_type_1' argument must be one of: 'rdt', 'lab', 'dual', 'multi'")
  }
  # Second test type must be selected
  if (is.null(test_type_2) || !length(test_type_2) || !all(test_type_2 %in% c("rdt", "lab", "dual", "multi"))) {
    stop("'test_type_2' argument must be one of: 'rdt', 'lab', 'dual', 'multi'")
  }
  
  # Setting up empty matrix - include all states
  syph_matrix <- matrix(nrow = 39, ncol = 20, 
                        dimnames = list(NULL, 
                                        c("week", "syph_test", "untest_seroneg", "testneg_seroneg", "tunnelpos_seroneg",
                                          "seroneg_testpos_notx", "seroneg_testpos_tx", 
                                          "untest_new_seropos_noprior", "untest_new_seropos_prior",
                                          "seropos_testneg", 
                                          "tunnelpos_seropos", "seropos_notx", "seropos_tx",
                                          "seropos_earlytx", "seropos_latetx", "deaths", "anc_attend", 
                                          "maternal_mort", "test_sens", "pop_size")))
  
  syph_df <- as.data.frame(syph_matrix)
  
  # Renaming variables to match function
  base_sens_syp <- base.sens
  base_spec_syp <- base.spec
  dual_sens_syp <- dual.sens.syp
  dual_spec_syp <- dual.spec.syp
  p_syp_tx_dual <- p.syp.tx.dual
  p_syp_tx_base <- p.syp.tx.base
  thirdtrimester_start <- early.preg.wks
  
  # Entering missing data for RDT sensitivity and specificity - can add later if we get it
  rdt_sens_syp <- NA
  rdt_spec_syp <- NA
  multi_sens_syp <- NA
  multi_spec_syp <- NA
  
  
  ## Set up input that do not depend on the time within the model (on the week of testing)
  # Setting up test sensitivity and specificity based on test choice - NEED TO UPDATE WITH RDT AND MULTIPLEX
  first_test_sens <- case_when(test_type_1 == "lab" ~ base_sens_syp,
                               test_type_1 == "dual" ~ dual_sens_syp,
                               test_type_1 == "rdt" ~ rdt_sens_syp,
                               test_type_1 == "multi" ~ multi_sens_syp)
  first_test_spec <- case_when(test_type_1 == "lab" ~ base_spec_syp,
                               test_type_1 == "dual" ~ dual_spec_syp,
                               test_type_1 == "rdt" ~ rdt_spec_syp,
                               test_type_1 == "multi" ~ multi_spec_syp)
  second_test_sens <- case_when(test_type_2 == "lab" ~ base_sens_syp,
                               test_type_2 == "dual" ~ dual_sens_syp,
                               test_type_2 == "rdt" ~ rdt_sens_syp,
                               test_type_2 == "multi" ~ multi_sens_syp)
  second_test_spec <- case_when(test_type_2 == "lab" ~ base_spec_syp,
                               test_type_2 == "dual" ~ dual_spec_syp,
                               test_type_2 == "rdt" ~ rdt_spec_syp,
                               test_type_2 == "multi" ~ multi_spec_syp)
  
  ### Building model
  # Setting up first week parameters
  # based on population size, baseline syphilis prevalence, and female adult mortality
  syph_df[1, ] <- c(1, 0, population_size*(1-syp_prev), rep(0,4), population_size*syp_prev, rep(0, 9), mort_adult_females, 0, sum(syph_df[1,3:13]))
  
  ## Run this to get full matrix of population counts by syphilis prevalence and testing status
  for (x in 2:deliveryweek) {
    ## Setting up variables that change depending on ANC attendance and testing
    # mortality - background female mortality except for delivery week
    maternal_mort <- ifelse(x == deliveryweek, mort_maternal_perinatal, mort_adult_females)
    # testing to be early and/or late ANC depending on testing algorithm - this will 
    syph_test <- ifelse((x == weekofanc1 & first_ANC_test == 1) | (x == weekofanc2 & late_ANC_test == 1), 1, 0)
    # ANC attendance for first and second ANC
    anc_attend <- case_when(x == weekofanc1~ att_firstANC,
                            x == weekofanc2 ~ att_lategest,
                            TRUE ~ 0)
    
    test_type <- case_when(test_type_1 %in% c("dual", "multi") & x == weekofanc1 ~ "dual",
                           test_type_1 %in% c("rdt", "lab") & x == weekofanc1 ~ "rdt",
                           test_type_2 %in% c("dual", "multi") & x == weekofanc2 ~ "dual",
                           test_type_2 %in% c("rdt", "lab") & x == weekofanc2 ~ "rdt",
                           TRUE ~ NA)
    
    # Test coverage and type - MAY NEED TO ADD MULTIPLEX COVERAGE HERE (assuming same as dual for now)
    test_coverage <- case_when(test_type_1 %in% c("dual", "multi") & x == weekofanc1 ~ testcoverage_dual,
                               test_type_1 %in% c("rdt", "lab") & x == weekofanc1 ~ testcoverage_syp,
                               test_type_2 %in% c("dual", "multi") & x == weekofanc2 ~ testcoverage_dual,
                               test_type_2 %in% c("rdt", "lab") & x == weekofanc2 ~ testcoverage_syp,
                               TRUE ~ 0)
    
    # Syphilis treatment probabilities - MAY NEED TO ADD MULTIPLEX COVERAGE HERE (assuming same as dual for now)
    prob_syph_treat <- case_when(test_type_1 %in% c("dual", "multi") & x == weekofanc1 ~ p_syp_tx_dual,
                                 test_type_1 %in% c("rdt", "lab") & x == weekofanc1 ~ p_syp_tx_base,
                                 test_type_2 %in% c("dual", "multi") & x == weekofanc2 ~ p_syp_tx_dual,
                                 test_type_2 %in% c("rdt", "lab") & x == weekofanc2 ~ p_syp_tx_base,
                                 TRUE ~ p_syp_tx_base)
    
    # setting up test sensitivity to be dependent on timing of test and specific test
    test_sens <- case_when(x == weekofanc1~ first_test_sens,
                           x == weekofanc2 ~ second_test_sens,
                           TRUE ~ 0)
    test_spec <- case_when(x == weekofanc1~ first_test_spec, 
                           x == weekofanc2 ~ second_test_spec,
                           TRUE ~ 0)
      
    # Setting up testing 
    week <- x
    syph_df$week[x] <- x 
    syph_df$syph_test[x] <- syph_test
    syph_df$anc_attend[x] <- anc_attend
    syph_df$maternal_mort[x] <- maternal_mort    # making sure this runs correctly - can delete later
    syph_df$test_sens[x] <- test_sens            # checking - may delete later
    
    
    ## Getting populations for each syphilis and test category
    
    # Untested and seronegative
    syph_df$untest_seroneg[x] <- if (syph_df$syph_test[x] == 0) {
      (1-maternal_mort) * (1-syp_inc_wk) * syph_df$untest_seroneg[x-1]
      } else { (1-maternal_mort) * (1-syp_inc_wk) * syph_df$untest_seroneg[x-1] * (1-syph_df$anc_attend[x] * test_coverage)
        }
    
    # Tested negative and seronegative
    syph_df$testneg_seroneg[x] <- if (syph_test == 0) {
      (1-maternal_mort) * (1-syp_inc_wk) * syph_df$testneg_seroneg[x-1]
      } else if (syph_test == 1 & week == weekofanc1) {
        (1-maternal_mort) * (1-syp_inc_wk) * syph_df$testneg_seroneg[x-1] + 
          (syph_df$untest_seroneg[x-1]) * (1-maternal_mort) * (1-syp_inc_wk) * anc_attend * test_coverage * test_spec
        } else if (syph_test == 1 & week == weekofanc2) {
          (1-maternal_mort) * (1-syp_inc_wk) * (1-anc_attend * test_coverage) * syph_df$testneg_seroneg[x-1] + 
            ((syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) * (1-maternal_mort) * (1-syp_inc_wk) * 
               anc_attend * test_coverage * test_spec)
          } else { NA }
    
    # Tunnel test - keeping at 0 for now, not sure when this would not be 0
    syph_df$tunnelpos_seroneg[x] <- 0
    
    # seroneg, test positive, no tx
    syph_df$seroneg_testpos_notx[x] <- syph_df$seroneg_testpos_notx[x-1] * (1-maternal_mort) * (1-syp_inc_wk) +
      if (syph_test == 1 & test_type == "dual") {
        (syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) * (1-maternal_mort) * (1-syp_inc_wk) * 
          (syph_df$anc_attend[x] *testcoverage_dual * (1-dual_spec_syp) * (1-p_syp_tx_dual))
        } else {
          syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-syp_inc_wk) * (1-p_syp_tx_base)
          }
    
    syph_df$seroneg_testpos_tx[x] <- syph_df$seroneg_testpos_tx[x-1] * (1-maternal_mort) * (1-syp_inc_wk) +
      if (syph_test == 1 & test_type == "dual") {
        (syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) * (1-maternal_mort) * (1-syp_inc_wk) * 
          (syph_df$anc_attend[x] * testcoverage_dual * (1-dual_spec_syp) * (p_syp_tx_dual))
      } else {
        syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-syp_inc_wk) * (p_syp_tx_base)
      }
    
    # untested and new seropos based on incidence rate
    syph_df$untest_new_seropos_noprior[x] <- (syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) * 
                                              (1-maternal_mort) * syp_inc_wk +
                                              ifelse(syph_df$syph_test[x] == 1,
                                                     syph_df$untest_new_seropos_noprior[x-1] * (1-maternal_mort) *
                                                       (1-anc_attend*test_coverage),
                                                     syph_df$untest_new_seropos_noprior[x-1] * (1-maternal_mort))
    
    # untested and new seropos based on incidence rate and test sensitivity (?)
    syph_df$untest_new_seropos_prior[x] <- (syph_df$seroneg_testpos_notx[x-1] + syph_df$seroneg_testpos_tx[x-1]) * 
                                            (1-maternal_mort) * (syp_inc_wk) +
                                            syph_df$untest_new_seropos_prior[x-1] * (1-maternal_mort)
    
    # seropositive but tested negative (those that falsely tested negative) - based on test specificity
    syph_df$seropos_testneg[x] <- ifelse(syph_df$syph_test[x] == 1,
                                         syph_df$untest_new_seropos_noprior[x-1] * (1-maternal_mort) * anc_attend*(1-test_sens)*test_coverage +
                                           syph_df$seropos_testneg[x-1] * (1-maternal_mort) * (1-anc_attend*test_sens*test_coverage), 
                                         syph_df$seropos_testneg[x-1] * (1-maternal_mort))
                                         
    # tunnel test positive and seropositive
    syph_df$tunnelpos_seropos[x] <- ifelse(syph_df$syph_test[x] == 1 & test_type == "rdt",
                                           (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) *
                                             anc_attend * test_coverage * test_sens,
                                           0)
    
    # seropositive with no treatment
    syph_df$seropos_notx[x] <- syph_df$seropos_notx[x-1] * (1 - maternal_mort) +   
      if (syph_test == 1 & test_type == "dual") {
        (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
          anc_attend * test_coverage * test_sens * (1-prob_syph_treat)
        } else { 
          syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
            syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * syp_inc_wk * (1-prob_syph_treat)
          } 
    
       
    # seropositive with treatment
    syph_df$seropos_tx[x] <- syph_df$seropos_tx[x-1] * (1 - maternal_mort) +
      if (syph_test == 1 & test_type == "dual") {
        (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
          anc_attend * test_coverage * test_sens * (prob_syph_treat)
      } else { 
        syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
          syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * syp_inc_wk * (prob_syph_treat)
        } 
    
    # seropositive with early treatment (before third trimester start)
    syph_df$seropos_earlytx[x] <- syph_df$seropos_earlytx[x-1] * (1-maternal_mort) + 
      if (week < thirdtrimester_start & syph_test == 1 & test_type == "dual") {
        (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
          anc_attend * test_coverage * test_sens * (prob_syph_treat)
      } else if (week < thirdtrimester_start) {
        syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
          syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * syp_inc_wk * (prob_syph_treat)
      } else { 0 }
    
    
    # seropositive with late treatment (after third trimester start)
    syph_df$seropos_latetx[x] <- syph_df$seropos_latetx[x-1] * (1-maternal_mort) + 
      if (week >= thirdtrimester_start & syph_test == 1 & test_type == "dual") {
        (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
          anc_attend * test_coverage * test_sens * (prob_syph_treat)
      } else if (week >= thirdtrimester_start) {
        syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
          syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * syp_inc_wk * (prob_syph_treat)
      } else { 0 }
    
    
    #Calculating number of deaths
    syph_df$deaths[x] <- sum(syph_df[x-1, 3:13]) * maternal_mort + syph_df$deaths[x-1] 
    
    # Checking to make sure population size is stable
    syph_df$pop_size[x] <- sum(syph_df[x, c(3:13,16)])
    
  }

  return(syph_df)
  
}




################################
# Adverse birth outcomes model #
################################


# syph_adv_out <- matrix(nrow = 6, ncol = 13, 
#                       dimnames = list(NULL, 
#                                       c("untest_seroneg", "testneg_seroneg", "tunnelpos_seroneg",
#                                         "seroneg_testpos_notx", "seroneg_testpos_tx", 
#                                         "untest_new_seropos_noprior", "untest_new_seropos_prior",
#                                         "seropos_testneg", 
#                                         "tunnelpos_seropos", "seropos_notx", "seropos_tx",
#                                         "seropos_earlytx", "seropos_latetx")))
# 
# syph_adv_df <- as.data.frame(syph_adv_out)
# 
# rownames(syph_adv_df) <- c("stillbirths", "neonatal_deaths", "premature_lbw", "congenital_syphilis", "asymptomatic", "healthy_lb")



# get early treatment number - turn this into a function with inputs in the main function
# may need an ifelse to make this 0 if the first test is not run - I think those stillbirths will be captured in the other population groups
early_tx_single_test_func <- function(first_ANC_test = 1, test_type, maternal_syph_df) {
  out <- if (test_type %in% c("dual", "rdt")) {
    first_ANC_test * maternal_syph_df[weekofanc1, which(names(maternal_syph_df) == "seropos_earlytx")] 
  } else if (test_type == "lab") {
    first_ANC_test * maternal_syph_df[(weekofanc1+1), which(names(maternal_syph_df) == "seropos_earlytx")] 
  }
  return(out)
}

early_tx_repeat_test_func <- function(late_ANC_test = 1, test_type, maternal_syph_df) {
  out <- if(test_type %in% c("dual", "rdt")) {
    late_ANC_test * (maternal_syph_df[weekofanc2, which(names(maternal_syph_df) == "seropos_earlytx")] - 
                       maternal_syph_df[(weekofanc2-1), which(names(maternal_syph_df) == "seropos_earlytx")])
  } else if (test_type == "lab") {
    late_ANC_test * (maternal_syph_df[weekofanc2+1, which(names(maternal_syph_df) == "seropos_earlytx")] - 
                       maternal_syph_df[(weekofanc2), which(names(maternal_syph_df) == "seropos_earlytx")])
  }
  return(out)
} 





syph_adverse_outcomes_func <- function(maternal_syph_df, test_type_1, test_type_2) {
  
  #### Function description
  #### Inputs:
  #### Outputs:
  
  
  
  # Setting delay in results based on lab or rapid test (1 week delay for lab tests)
  if (test_type_1 %in% c("dual", "rdt")) {
    x1 <- 0
  } else if (test_type_1 == "lab") {
    x1 <- 1
  } else {
    print("No test type defined")
  }

  if (test_type_2 %in% c("dual", "rdt")) {
    x2 <- 0
  } else if (test_type_2 == "lab") {
    x2 <- 1
  } else {
    print("No test type defined")
  }
  

  # Getting stillbirths for seronegative syphilis births
  stillbirths1 <- maternal_syph_df[deliveryweek, 
                                  which(names(maternal_syph_df) == "untest_seroneg"):which(names(maternal_syph_df) == "seroneg_testpos_tx")] * 
                                  nosyp.p.stillbirth
  
  # Getting stillbirths for seropositive syphilis births
  stillbirths2 <- maternal_syph_df[deliveryweek, 
                                   c(which(names(maternal_syph_df) == "untest_new_seropos_noprior"):which(names(maternal_syph_df) == "seropos_notx"),
                                     which(names(maternal_syph_df) == "seropos_latetx"))] *
                                   (syp_pct_live*syp.notx.p.stillbirth + (1-syp_pct_live)*nosyp.p.stillbirth)

  
  
  # Getting stillbirths for seropositive syphilis but with early treatment
  # early_tx_single_test <- early_tx_single_test_func(1, "dual", maternal_syph_df)
  # early_tx_rpt_test <- early_tx_repeat_test_func(1, "dual", maternal_syph_df)
  
  early_tx_single_test <- early_tx_single_test_func(1, test_type = test_type_1, maternal_syph_df)
  early_tx_rpt_test <- early_tx_repeat_test_func(1, test_type = test_type_2, maternal_syph_df)
  
  
  
  
  # Getting early treatment number
    # We want x = 0 if using rapid tests or x = 1 if using lab tests
  stillbirths3 <- maternal_syph_df[deliveryweek, which(names(maternal_syph_df) == "seropos_earlytx")] * syp_pct_live * 
    (early_tx_single_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc1+x1)^0.5)*syp.tx.p.stillbirth) +
    (early_tx_rpt_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc2+x2)^0.5)*syp.tx.p.stillbirth) 
  
  # Putting all stillbirths together
  stillbirths <- c(unlist(stillbirths1), unlist(stillbirths2[1:5]), 0, unlist(stillbirths3), unlist(stillbirths2[6]))
  
  
  
  # Getting neonatal deaths
  neonatal1 <- maternal_syph_df[deliveryweek, 
                                c(which(names(maternal_syph_df) == "untest_new_seropos_noprior"):which(names(maternal_syph_df) == "seropos_notx"),
                                  which(names(maternal_syph_df) == "seropos_latetx"))] * syp_pct_live * syp.notx.p.death
  
  neonatal2 <- maternal_syph_df[deliveryweek, which(names(maternal_syph_df) == "seropos_earlytx")] * syp_pct_live * 
    (early_tx_single_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc1+x1)^0.5)*syp.tx.p.death) +
    (early_tx_rpt_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc2+x2)^0.5)*syp.tx.p.death) 
  
  neonatal <- c(rep(0, 5), unlist(neonatal1[1:5]), 0, unlist(neonatal2), unlist(neonatal1[6]))
  
  
  # Getting premature/LBW 
  prem1 <- maternal_syph_df[deliveryweek, 
                            c(which(names(maternal_syph_df) == "untest_new_seropos_noprior"):which(names(maternal_syph_df) == "seropos_notx"),
                              which(names(maternal_syph_df) == "seropos_latetx"))] * syp_pct_live * syp.notx.p.lbw
  
  prem2 <- maternal_syph_df[deliveryweek, which(names(maternal_syph_df) == "seropos_earlytx")] * syp_pct_live * 
    (early_tx_single_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc1+x1)^0.5)*syp.tx.p.lbw) +
    (early_tx_rpt_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc2+x2)^0.5)*syp.tx.p.lbw) 
  
  premature <- c(rep(0, 5), unlist(prem1[1:5]), 0, unlist(prem2), unlist(prem1[6]))
  
  
  # Getting clinical congenital syphilis
  cong1 <- maternal_syph_df[deliveryweek, 
                            c(which(names(maternal_syph_df) == "untest_new_seropos_noprior"):which(names(maternal_syph_df) == "seropos_notx"),
                              which(names(maternal_syph_df) == "seropos_latetx"))] * syp_pct_live * syp.notx.p.congsyp.symp
  
  cong2 <- maternal_syph_df[deliveryweek, which(names(maternal_syph_df) == "seropos_earlytx")] * syp_pct_live * 
    (early_tx_single_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc1+x1)^0.5)*syp.tx.p.congsyp.symp) +
    (early_tx_rpt_test / (early_tx_single_test + early_tx_rpt_test)) *
    (0.0011*exp((weekofanc2+x2)^0.5)*syp.tx.p.congsyp.symp) 
  
  congenital <- c(rep(0, 5), unlist(cong1[1:5]), 0, unlist(cong2), unlist(cong1[6]))
  
  
  # Asymptomatic - mother untreated
  asympt <- rep(0, 13)
  for (i in c(6:10, 13)) {
    asympt[i] <- maternal_syph_df[deliveryweek, (i+2)] - (stillbirths[i] + neonatal[i] + premature[i] + congenital[i])
  }
  
  
  # Healthy infants
  healthy <- rep(0, 13) 
  for (i in c(1:5, 12)) {
    healthy[i] <- maternal_syph_df[deliveryweek, (i+2)] - (stillbirths[i] + neonatal[i] + premature[i] + congenital[i])
  }
  
  
  # Put adverse outcomes together
  adverse_df <- as.data.frame(stillbirths) %>%
    bind_cols(as.data.frame(neonatal)) %>%
    bind_cols(as.data.frame(premature)) %>%
    bind_cols(as.data.frame(congenital)) %>%
    bind_cols(as.data.frame(asympt)) %>%
    bind_cols(as.data.frame(healthy)) %>%
    t() %>% as.data.frame() 
  
  colnames(adverse_df) <- colnames(maternal_syph_df)[3:15]
  
  return(adverse_df)
  #out <- list("adverse" = adverse_df, "asympt" = asympt, "healthy" = healthy)

}




# Getting infant deaths due to syphilis (stillbirths and neonatal deaths)
syph_infant_deaths_func <- function(adverse_df) {
  
  ### Function description
  
  out <- sum(adverse_df[1:2, ], na.rm=TRUE)
  return(out)
}







# Getting infant outcomes at one year
syph_infant_outcomes_1yr_func <- function(adverse_df) {
  
  ### This function sums up the number of infants in each outcome category by the end of the first year
  
  out <- NULL
  out$stillbirths <- sum(adverse_df[1, ], na.rm=TRUE)
  out$neonatal <- sum(adverse_df[2, ], na.rm=TRUE)
  out$premature <- sum(adverse_df[3, ], na.rm=TRUE) * surv_nhiv
  out$cong_syphilis <- sum(adverse_df[4, ], na.rm=TRUE) * surv_nhiv
  out$asympt <- sum(adverse_df[5, ], na.rm=TRUE) * surv_nhiv
  
  return(out)
}



##########################
## Syphilis costs model ##
##########################


# Wrap into a function that takes test type as an input

# Needs for function
  # test type
  # first and second ANC test - we don't need this because our costs depend on finished testing matrix


syph_cost_func <- function(maternal_syph_df, adverse_df,
                           test_type_1, test_type_2) {
  
  # Inputs are the test type: "dual", "rdt", or "lab", and maternal syphilis dataframe
    # Currently "rdt" is the same as "dual"
  

  # Setting up empty matrix - include all states
  syph_cost_matrix <- matrix(nrow = deliveryweek, ncol = 8, 
                        dimnames = list(NULL, 
                                        c("rpr_tests", "tpha_tests", "dual_tests", "treatments_fp", "treatments_tp", "test_cost",
                                          "treatment_cost_fp", "treatment_cost_tp")))
  
  syph_cost_df <- as.data.frame(syph_cost_matrix)
  
  # Starting states
  syph_cost_df[1, ] <- rep(0, 8)
  
  # Setting up testing costs matrix
  for (x in 2:deliveryweek) {
    
    # Getting attendance at ANC for first and second ANC
    att_anc <- case_when(x == weekofanc1 ~ att_firstANC,
                         x == weekofanc2 ~ att_lategest,
                         TRUE ~ as.numeric(0))
    
    # Getting test type
    test_type <- case_when(x <= weekofanc2 ~ test_type_1,
                           x >= weekofanc2 ~ test_type_2)
    
    # Update this if RDT test coverage is different than dual test coverage
    testcoverage_rdt <- testcoverage_dual
    testcoverage <- case_when(test_type == "dual" ~ testcoverage_dual,
                              test_type == "rdt" ~ testcoverage_rdt,
                              test_type == "lab" ~ testcoverage_syp)
    
    # Update if RDT test sensitivity is different than dual test sensitivity
    rdt.sens.syp <- dual.sens.syp
    sensitivity <- case_when(test_type == "dual" ~ dual.sens.syp,
                             test_type == "rdt" ~ rdt.sens.syp,
                             test_type == "lab" ~ base.sens)
    
    # Update if probability of treatment is different for rdt compared to dual
    p.syp.tx.rdt <- p.syp.tx.dual
    prob_treatment <- case_when(test_type == "dual" ~ p.syp.tx.dual,
                                test_type == "rdt" ~ p.syp.tx.rdt,
                                test_type == "lab" ~ p.syp.tx.base)
    
    # Update this when we get costs of syphilis RDT
    cost_rdt <- cost_dual
    cost_dual_rdt <- case_when(test_type == "dual" ~ cost_dual,
                               test_type == "rdt" ~ cost_rdt,
                               TRUE ~ as.numeric(0))
    
    
    # Number of RPR tests
    syph_cost_df$rpr_tests[x] <- if (test_type %in% c("dual", "rdt")) {
      0 
      } else if (test_type %in% c("lab")) {
        maternal_syph_df$syph_test[x] * (maternal_syph_df$untest_seroneg[x-1] + maternal_syph_df$testneg_seroneg[x-1]) *
      (1-mort_adult_females) * (1-syp_inc_wk) * att_anc * testcoverage + 
      (maternal_syph_df$untest_new_seropos_noprior[x-1] + maternal_syph_df$seropos_testneg[x-1]) * (1-mort_adult_females) *
      att_anc * testcoverage
      } else {
        print("No test type selected")
      }
    
    # Number of TPHA tests
    syph_cost_df$tpha_tests[x] <- if (test_type %in% c("dual", "rdt")) {
      0 
    } else if (test_type %in% c("lab")) {
      maternal_syph_df$syph_test[x] * (maternal_syph_df$testneg_seroneg[x-1] + maternal_syph_df$tunnelpos_seroneg[x-1]) *
      (1-mort_adult_females) * (1-syp_inc_wk) * (att_anc * testcoverage * (1-rpr.spec)) +
      (maternal_syph_df$untest_new_seropos_noprior[x-1] + maternal_syph_df$seropos_testneg[x-1]) * (1-mort_adult_females) *
      att_anc * testcoverage * rpr.sens
    } else {
      print("No test type selected")
    }
    
    # Number of dual or RDT tests
    syph_cost_df$dual_tests[x] <- if (test_type %in% c("dual", "rdt")) {
      maternal_syph_df$syph_test[x] * (maternal_syph_df$untest_seroneg[x-1] + maternal_syph_df$testneg_seroneg[x-1]) * 
      (1-mort_adult_females) * (1-syp_inc_wk) * att_anc * testcoverage +
      maternal_syph_df$syph_test[x] * (maternal_syph_df$untest_new_seropos_noprior[x-1] + maternal_syph_df$seropos_testneg[x-1]) * (1-mort_adult_females) *
      att_anc * testcoverage_dual
    } else if (test_type %in% c("lab")) {
      0
    } else {
      print("No test type selected")
    }
    
    # Number of FP treatments
    syph_cost_df$treatments_fp[x] <- if (test_type %in% c("dual", "rdt")) {
      maternal_syph_df$syph_test[x] * (maternal_syph_df$untest_seroneg[x-1] + maternal_syph_df$testneg_seroneg[x-1]) * 
        (1-mort_adult_females) * (1-syp_inc_wk) * att_anc * testcoverage_dual * (1-dual.spec.syp) * p.syp.tx.dual
    } else if (test_type %in% c("lab")) {
      maternal_syph_df$syph_test[x] * maternal_syph_df$tunnelpos_seroneg[x-1] * (1-mort_adult_females) * (1-syp_inc_wk) * p.syp.tx.base
    } else {
      print("No test type selected")
    }
    
    # Number of TP treatments
    syph_cost_df$treatments_tp[x] <- if (test_type %in% c("dual", "rdt")) {
      maternal_syph_df$syph_test[x] * (maternal_syph_df$untest_new_seropos_noprior[x-1] + maternal_syph_df$seropos_testneg[x-1]) *
        (1-mort_adult_females) * (att_anc) * testcoverage * sensitivity * prob_treatment 
    } else if (test_type %in% c("lab")) {
      maternal_syph_df$syph_test[x-1] * maternal_syph_df$tunnelpos_seropos[x-1] * (1-mort_adult_females) * prob_treatment +
        maternal_syph_df$tunnelpos_seroneg[x-1] * (1-mort_adult_females) * syp_inc_wk * prob_treatment
    } else {
      print("No test type selected")
    }
      
    # Test costs - we include the cost of dual testing here, but need to remove for CE model because it's also captured in HIV testing
    syph_cost_df$test_cost[x] <- syph_cost_df$rpr_tests[x]*cost_rpr + syph_cost_df$tpha_tests[x]*cost_tpha
    
    # Treatment costs
    syph_cost_df$treatment_cost_fp[x] <- syph_cost_df$treatments_fp[x] * cost_penicillin
    
    syph_cost_df$treatment_cost_tp[x] <- syph_cost_df$treatments_tp[x] * cost_penicillin
    
    
  }
  
  syph_total_costs <- colSums(syph_cost_df)
  syph_total_costs$cong_syph_costs <- sum(adverse_df[4, ]) * cost_cong_tx
  # DMC - Why does the excel model not include some costs associated with congenital syphilis? 

return(syph_total_costs)

}


















# Needs - 
  # update test coverage to be flexible. Can add dual or individual tests































##############
## OLD CODE ##
##############



# ### Setting up syphilis data frame
# ## Run this to get full matrix of population counts by syphilis prevalence and testing status
# for (x in 2:delivery) {
#   ## Setting up variables that change depending on ANC attendance and testing
#   # mortality - background female mortality except for delivery week
#   maternal_mort <- ifelse(x == delivery, mort_maternal_perinatal, mort_adult_females)
#   # testing to be early and/or late ANC depending on testing algorithm
#   syph_test <- ifelse((x == weekofanc1& first_ANC_test == 1) | (x == weekofanc2 & late_ANC_test == 1), 1, 0)
#   # ANC attendance for first and second ANC
#   anc_attend <- case_when(x == weekofanc1~ att_firstANC,
#                           x == weekofanc2 ~ att_lategest,
#                           TRUE ~ 0)
#   
#   # Test coverage and type - turn into function to use syphilis and dual test coverage
#   test_coverage <- testcoverage_dual
#   test_type <- "dual"
#   
#   # Syphilis treatment probabilities - depends on test type
#   prob_syph_treat <- case_when(test_type == "dual" ~ p_syp_tx_dual,
#                                test_type == "rdt" ~ p_syp_tx_base,
#                                TRUE ~ NA)
#   
#   # Setting up test sensitivity - may need to adjust based on testing scenario
#   first_test_sens <- dual_sens_syp
#   second_test_sens <- dual_sens_syp
#   first_test_spec <- dual_spec_syp
#   second_test_spec <- dual_spec_syp
#   # setting up test sensitivity to be dependent on timing of test and specific test
#   test_sens <- case_when(x == weekofanc1~ first_test_sens,
#                          x == weekofanc2 ~ second_test_sens,
#                          TRUE ~ 0)
#   test_spec <- case_when(x == weekofanc1~ first_test_spec, 
#                          x == weekofanc2 ~ second_test_spec,
#                          TRUE ~ 0)
#   
#   # Setting up testing 
#   syph_df$week[x] <- x 
#   syph_df$syph_test[x] <- syph_test
#   syph_df$anc_attend[x] <- anc_attend
#   syph_df$maternal_mort[x] <- maternal_mort    # making sure this runs correctly - can delete later
#   syph_df$test_sens[x] <- test_sens            # checking - may delete later
#   
#   
#   ## Getting populations for each syphilis and test category
#   # Untested and seronegative
#   syph_df$untest_seroneg[x] <- ifelse(syph_df$syph_test[x] == 0,
#                                       (1-maternal_mort) * (1-syp_inc_wk) * syph_df$untest_seroneg[x-1],
#                                       (1-maternal_mort) * (1-syp_inc_wk) * syph_df$untest_seroneg[x-1] * 
#                                         (1-syph_df$anc_attend[x] * test_coverage))
#   
#   # Tested negative and seronegative
#   syph_df$testneg_seroneg[x] <- if (syph_df$syph_test[x] == 0) {
#     (1-maternal_mort) * (1-syp_inc_wk) * syph_df$testneg_seroneg[x-1]
#   } else if (syph_df$syph_test[x] == 1 & syph_df$week[x] == first_ANC) {
#     (1-maternal_mort) * (1-syp_inc_wk) * syph_df$testneg_seroneg[x-1] +
#       (syph_df$untest_seroneg[x-1] * (1-maternal_mort) * (1-syp_inc_wk) * 
#          anc_attend * test_coverage * test_spec)
#   } else if (syph_df$syph_test[x] == 1 & syph_df$week[x] == second_ANC) {
#     (1-maternal_mort) * (1-syp_inc_wk) * (1-anc_attend * test_coverage) * syph_df$testneg_seroneg[x-1] +
#       ((syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) * (1-maternal_mort) * (1-syp_inc_wk) * 
#          anc_attend * test_coverage * test_spec)
#   } else {
#     NA
#   }
#   
#   # Tunnel test - keeping at 0 for now, not sure when this would not be 0
#   syph_df$tunnelpos_seroneg[x] <- 0
#   
#   # seroneg, test positive, no tx
#   syph_df$seroneg_testpos_notx[x] <- syph_df$seroneg_testpos_notx[x-1] * (1-maternal_mort) * (1-syp_inc_wk) +
#     ifelse(syph_df$syph_test[x] == 1,
#            (syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) *
#              (1-maternal_mort) * (1-syp_inc_wk) * (syph_df$anc_attend[x] *
#                                                      testcoverage_dual *
#                                                      (1-dual_spec_syp) *
#                                                      (1-p_syp_tx_dual)),
#            syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-syp_inc_wk) *
#              (1-p_syp_tx_base))
#   
#   # seroneg, test positive, tx
#   syph_df$seroneg_testpos_tx[x] <- syph_df$seroneg_testpos_tx[x-1] * (1-maternal_mort) * (1-syp_inc_wk) +
#     ifelse(syph_df$syph_test[x] == 1,
#            (syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) *
#              (1-maternal_mort) * (1-syp_inc_wk) * (syph_df$anc_attend[x] *
#                                                      testcoverage_dual *
#                                                      (1-dual_spec_syp) *
#                                                      (p_syp_tx_dual)),
#            syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-syp_inc_wk) *
#              (p_syp_tx_base))
#   
#   # untested and new seropos based on incidence rate
#   syph_df$untest_new_seropos_noprior[x] <- (syph_df$untest_seroneg[x-1] + syph_df$testneg_seroneg[x-1]) * 
#     (1-maternal_mort) * syp_inc_wk +
#     ifelse(syph_df$syph_test[x] == 1,
#            syph_df$untest_new_seropos_noprior[x-1] * (1-maternal_mort) *
#              (1-anc_attend*test_coverage),
#            syph_df$untest_new_seropos_noprior[x-1] * (1-maternal_mort))
#   
#   # untested and new seropos based on incidence rate and test sensitivity (?)
#   syph_df$untest_new_seropos_prior[x] <- (syph_df$seroneg_testpos_notx[x-1] + syph_df$seroneg_testpos_tx[x-1]) * 
#     (1-maternal_mort) * (syp_inc_wk) +
#     syph_df$untest_new_seropos_prior[x-1] * (1-maternal_mort)
#   
#   # seropositive but tested negative (those that falsely tested negative) - based on test specificity
#   syph_df$seropos_testneg[x] <- ifelse(syph_df$syph_test[x] == 1,
#                                        syph_df$untest_new_seropos_noprior[x-1] * (1-maternal_mort) * anc_attend*(1-test_sens)*test_coverage +
#                                          syph_df$seropos_testneg[x-1] * (1-maternal_mort) * (1-anc_attend*test_sens*test_coverage), 
#                                        syph_df$seropos_testneg[x-1] * (1-maternal_mort))
#   
#   # tunnel test positive and seropositive
#   syph_df$tunnelpos_seropos[x] <- ifelse(syph_df$syph_test[x] == 1 & test_type == "rdt",
#                                          (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) *
#                                            anc_attend * test_coverage * test_sens,
#                                          0)
#   
#   # seropositive with no treatment
#   syph_df$seropos_notx[x]<- syph_df$seropos_notx[x-1] * (1 - maternal_mort) +
#     ifelse(syph_test == 1 & test_type == "dual",
#            (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
#              anc_attend * test_coverage * test_sens * (1-prob_syph_treat), 0) +
#     ifelse(test_type == "rdt",
#            syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
#              syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * syp_inc_wk * (1-prob_syph_treat), 0)
#   
#   # seropositive with treatment
#   syph_df$seropos_tx[x] <- syph_df$seropos_tx[x-1] * (1 - maternal_mort) +
#     ifelse(syph_test == 1 & test_type == "dual",
#            (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
#              anc_attend * test_coverage * test_sens * (prob_syph_treat), 0) +
#     ifelse(test_type == "rdt",
#            syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
#              syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * syp_inc_wk * (prob_syph_treat), 0)
#   
#   # seropositive with early treatment (before third trimester start)
#   syph_df$seropos_earlytx[x] <- syph_df$seropos_earlytx[x-1] * (1-maternal_mort) +
#     ifelse(syph_df$week[x] < thirdtrimester_start,
#            ifelse(syph_test == 1 & test_type == "dual",
#                   (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
#                     anc_attend * test_coverage * test_sens * (prob_syph_treat), 0) +
#              ifelse(test_type == "rdt",
#                     syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
#                       syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * syp_inc_wk * (prob_syph_treat), 0), 0)
#   
#   # seropositive with late treatment (after third trimester start)
#   syph_df$seropos_latetx[x] <- syph_df$seropos_latetx[x-1] * (1-maternal_mort) +
#     ifelse(syph_df$week[x] >= thirdtrimester_start,
#            ifelse(syph_test == 1 & test_type == "dual",
#                   (syph_df$untest_new_seropos_noprior[x-1] + syph_df$seropos_testneg[x-1]) * (1-maternal_mort) * 
#                     anc_attend * test_coverage * test_sens * (prob_syph_treat), 0) +
#              ifelse(test_type == "rdt",
#                     syph_df$tunnelpos_seroneg[x-1] * (1-maternal_mort) * (1-prob_syph_treat) +
#                       syph_df$tunnelpos_seropos[x-1] * (1-maternal_mort) * syp_inc_wk * (prob_syph_treat), 0), 0)
#   
#   #Calculating number of deaths
#   syph_df$deaths[x] <- sum(syph_df[x-1, 3:13]) * maternal_mort + syph_df$deaths[x-1] 
#   
#   # Checking to make sure population size is stable
#   syph_df$pop_size[x] <- sum(syph_df[x, c(3:13,16)])
#   
# }







