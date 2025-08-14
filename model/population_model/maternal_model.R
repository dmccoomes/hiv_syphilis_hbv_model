
#### Title: HIV-syphilis testing maternal model 
#### Purpose: Create HIV retesting maternal model
#### Author: David Coomes
#### Date: July 8, 2025


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table"
)

# rm(list = ls())
# 
# source(here("model", "supplemental", "life_table.R"))
# source(here("model", "syphilis_model", "syphilis_model.R"))
# source(here("model", "hiv_model", "hiv_model_final.R"))
# source(here("model", "model_output", "hiv_syphilis_hbv_model_output.R"))


## Inputs: maternal df
## Outputs: Dataframe of maternal HIV and syphilis outcomes for 20 years

## We want to run the maternal dataframe to get the number of women by HIV status
  ## Then we'll run this forward, making assumptions on testing and HIV incidence 
  ## Needs:
    # Testing frequency among non-pregnant women
    # HIV incidence among non-pregnant women
    # the numbers don't quite add up for the maternal population, but it's close (we lose < 0.0001%)


# Load mortality file
# mort_df <- read_csv(here("parameters", "kenya_mort_ihme.csv"))

# Helper to pull mortality by age and sex - not using this here at the moment
get_mortality <- function(age, sex = "female") {
  start_age <- if (age < 1) 0 else {
    possible <- sort(unique(na.omit(mort_df$age_start)))
    max(possible[possible <= age])
  }
  
  row <- mort_df %>% filter(age_start == start_age)
  if (sex == "female") return(row$mort_female)
  else return(row$mort_both)
  
}


# # Read in parameters 
# pop_param <- pop_param_func(here("parameters", "inp_par_dmc.xlsx"), "HIV_Syphilis_Dual_Preg",
#                             "Kenya", "Pregnant population")
# 
# list2env(pop_param, envir = parent.frame())
# 
# # Adding test 
# first_ANC_test <- 1
# late_ANC_test <- 0
# first_syphilis_test <- "lab"
# late_syphilis_test <- "dual"
# first_hiv_test <- "rdt"
# late_hiv_test <- "dual"
# syph_test <- 1
# pct_ART_female <- 0.9
# 
# test_type_1 <- "rdt"
# test_type_2 <- "rdt"

## Run syphilis models
## Run syphilis model - choices are "lab" or "rdt"
# maternal_syph_df <- syph_maternal_func(first_ANC_test = first_ANC_test, 
#                                        late_ANC_test = late_ANC_test,
#                                        test_type_1 = first_syphilis_test, 
#                                        test_type_2 = late_syphilis_test)
# 
# adverse_df <- syph_adverse_outcomes_func(maternal_syph_df, 
#                                          test_type_1 = first_syphilis_test,
#                                          test_type_2 = late_syphilis_test)
# 
# syph_infant_deaths <- syph_infant_deaths_func(adverse_df)
# 
# syph_infant_outcomes_1yr <- syph_infant_outcomes_1yr_func(adverse_df)
# 
# syph_cost_df <- syph_cost_func(test_type_1 = first_syphilis_test,
#                                test_type_2 = late_syphilis_test,
#                                maternal_syph_df = maternal_syph_df, 
#                                adverse_df = adverse_df)
# 
# 
# # We don't use this - why not? 
# maternal_hiv_syph_df <- hiv_maternal_func(first_ANC_test = 1, late_ANC_test = 1,
#                                           test_type_1 = "dual", test_type_2 = "dual",
#                                           syph_test = TRUE, 
#                                           syph_infant_deaths = syph_infant_deaths)
# 
# 
# 
# 
# 
# 
# hiv_inc <- 0.00053
# mort_adult_females_hiv <- mort_adult_females
# starting_age <- 22
# adult_screen_hiv <- 0.50




## Getting maternal HIV population
# Start with maternal HIV-syphilis dataframe
  # this dataframe tracks categories of hiv and syphilis infections but only reports weekly deaths
maternal_output_func <- function(first_ANC_test=1, late_ANC_test, 
                                 test_type_1, test_type_2,
                                 syph_test = FALSE,
                                 syph_infant_deaths = syph_infant_deaths, 
                                 starting_age = 22) {
  
  # Run syphilis models
  # Run syphilis model - choices are "lab" or "rdt"
  maternal_syph_df <- syph_maternal_func(first_ANC_test = first_ANC_test,
                                         late_ANC_test = late_ANC_test,
                                         test_type_1 = test_type_1,
                                         test_type_2 = test_type_2)

  adverse_df <- syph_adverse_outcomes_func(maternal_syph_df,
                                           test_type_1 = test_type_1,
                                           test_type_2 = test_type_2)

  syph_infant_deaths <- syph_infant_deaths_func(adverse_df)

  syph_infant_outcomes_1yr <- syph_infant_outcomes_1yr_func(adverse_df)

  syph_cost_df <- syph_cost_func(test_type_1 = test_type_1,
                                 test_type_2 = test_type_2,
                                 maternal_syph_df = maternal_syph_df,
                                 adverse_df = adverse_df)
  
  

  
  
  # Get maternal hiv syphilis cases with no syphilis deaths
  maternal_hiv_syph_nodeath <- hiv_maternal_func(first_ANC_test = first_ANC_test, late_ANC_test = late_ANC_test,
                                                 test_type_1 = test_type_1, test_type_2 = test_type_2,
                                                 syph_test = FALSE, 
                                                 syph_infant_deaths = syph_infant_deaths)
  

  mat_hiv_df <- maternal_hiv_syph_nodeath[39, 2:7] %>%                       # before deliveryweek
    mutate(dead = sum(maternal_hiv_syph_nodeath[1:39, 8], na.rm=TRUE)) %>%  # Add in sum of deaths before this time
    # Combine chronic and acute HIV categories into known status and unknown status
    mutate(hiv_pos_nart = ahiv_pos_nart + chiv_pos_nart,
           hiv_pos_uk = chiv_pos_uk,
           hiv_pos_art = ahiv_pos_art + chiv_pos_art) %>%
    bind_rows(
    maternal_hiv_syph_nodeath[92, 2:7] %>%
    mutate(dead = sum(maternal_hiv_syph_nodeath[ , 8], na.rm=TRUE)) %>%
    mutate(hiv_pos_nart = ahiv_pos_nart + chiv_pos_nart,
           hiv_pos_uk = chiv_pos_uk,
           hiv_pos_art = ahiv_pos_art + chiv_pos_art)
    ) %>%
    select(hiv_neg, hiv_pos_nart, hiv_pos_uk, hiv_pos_art, dead) 
  
  rownames(mat_hiv_df) <- NULL
  mat_hiv_df$year <- c(0, 1)
  
  ## Getting maternal syphilis population
  mat_syph_df <- maternal_syph_df[39, ] %>%
    mutate(syph_neg = untest_seroneg + testneg_seroneg + tunnelpos_seroneg + 
             seroneg_testpos_notx + seroneg_testpos_tx + seropos_tx,
           syph_pos = untest_new_seropos_noprior + untest_new_seropos_prior +
             seropos_testneg + tunnelpos_seropos + seropos_notx) %>%
    select(syph_neg, syph_pos) %>%
    rbind(c(syph_neg = NA, syph_pos = NA))
  
  rownames(mat_syph_df) <- NULL
  
  
  mat_df <- mat_hiv_df %>%
    bind_cols(mat_syph_df) %>%
    select(-syph_neg)
  
  ### Maternal data frame with year 0 and year 1
  mat_df <- mat_df %>%
    relocate(year, .before = "hiv_neg") %>%
    relocate(dead, .after = "syph_pos")
  
  
  
  
  # Adding starting age
  mat_df$age <- c(starting_age, starting_age+1)
  mat_df$hiv_pos_k <- mat_df$hiv_pos_art + mat_df$hiv_pos_nart
  mat_df <- mat_df %>%
    mutate(total = hiv_neg + hiv_pos_uk + hiv_pos_k + dead)
  
  mat_df <- mat_df %>%
    select(year, age, hiv_neg, hiv_pos_nart, hiv_pos_art, hiv_pos_k, hiv_pos_uk, syph_pos, dead, total)
  
  # Add hiv and non-hiv deaths
  mat_df$dead_nhiv <- mat_df$dead
  mat_df$dead_hiv <- 0
  
  # Get adult mortality rate for HIV pos and HIV neg women
  mort_adult_df <- mort_df %>%
    filter(`Age Group` %in% c("20-24 years", "25-29 years", "30-34 years", "35-39 years", 
                              "40-44 years", "45-49 years")) 
  mort_adult_df <- mort_adult_df[rep(seq_len(nrow(mort_adult_df)), each = 5), ]
  mort_adult_df$age <- rep(seq(20, 49, 1), 2)
  mort_adult_df <- mort_adult_df %>%
    select(Indicator, Female, age) %>%
    pivot_wider(id_cols = "age", names_from = "Indicator", values_from = "Female")
    
  
  ### Running maternal model
  for (x in 3:21) {
    # Adding blank row
    blank_row <- mat_df[1, ]
    blank_row[] <- NA
    mat_df <- rbind(mat_df, blank_row)
    
    # Year
    mat_df$year[x] <- (x-1)
    # Get age 
    mat_df$age[x] <- mat_df$age[x-1] + 1
    
    # Get mortality rate
    mort_rate_hiv <- mort_adult_df$`mortality_HIV+_total`[mort_adult_df$age == mat_df$age[x]]
    mort_rate_nhiv <- mort_adult_df$`mortality_HIV-_total`[mort_adult_df$age == mat_df$age[x]]
    
    
    # Getting non-hiv population
    mat_df$hiv_neg[x] <- mat_df$hiv_neg[x-1] * (1-hiv_inc) * (1-mort_rate_nhiv)
    #mat_df$hiv_neg[x] <- mat_df$hiv_neg[x-1] * (1-hiv_inc) * (1-mort_adult_females)
    
    # Getting HIV pos with unk status
    mat_df$hiv_pos_uk[x] <- mat_df$hiv_pos_uk[x-1]*(1-adult_screen_hiv)*(1-mort_rate_hiv) +    # old infection but not screened
      mat_df$hiv_neg[x]*hiv_inc*(1-adult_screen_hiv)*(1-mort_rate_hiv)     # newly infected but not screened
    
    # Getting HIV pos with known status
    mat_df$hiv_pos_k[x] <- mat_df$hiv_pos_k[x-1]*(1-mort_rate_hiv) +   # status already known
      mat_df$hiv_pos_uk[x-1]*adult_screen_hiv*(1-mort_rate_hiv) +        # old infection
      mat_df$hiv_neg[x-1]*hiv_inc*adult_screen_hiv*(1-mort_rate_hiv)     # new infection
    
    # Deaths
    mat_df$dead[x] <- mat_df$dead[x-1] + 
      mat_df$hiv_neg[x-1]*mort_rate_nhiv + 
      mat_df$hiv_pos_k[x-1]*mort_rate_hiv +
      mat_df$hiv_pos_uk[x-1]*mort_rate_hiv
    
    # splitting into HIV and non-HIV deaths
    mat_df$dead_hiv[x] <- mat_df$dead_hiv[x-1] +
      mat_df$hiv_pos_k[x-1]*mort_rate_hiv +
      mat_df$hiv_pos_uk[x-1]*mort_rate_hiv
    
    mat_df$dead_nhiv[x] <- mat_df$dead_nhiv[x-1] +
      mat_df$hiv_neg[x-1]*mort_rate_nhiv
    
    
    # Total
    mat_df$total[x] <- mat_df$hiv_neg[x] + mat_df$hiv_pos_k[x] + mat_df$hiv_pos_uk[x] + mat_df$dead[x]
    
  }
  
  
  dw_HIV_uk <- dw_HIVnotx*propLE_HIV_notx + dw_AIDSnotx*(1-propLE_HIV_notx)
  dw_HIV_k <- pct_ART_female*dw_HIVtx + (1-pct_ART_female)*(dw_HIVnotx*propLE_HIV_notx + dw_AIDSnotx * 
                                                              (1-propLE_HIV_notx))
  
  
  ## Getting DALYs and costs for maternal model
  mat_df <- mat_df %>%
    select(-c(hiv_pos_art, hiv_pos_nart, total)) %>%
    mutate(maternal_hiv_treatment_costs = (hiv_pos_k * pct_ART_female * cost_maternal_ART*52) / 
             (1+disc_rate)^(mat_df$year[x]+1),
           syph_dalys = syph_pos * 0.04,
           dw_HIV_k = dw_HIV_k,
           dw_HIV_uk = dw_HIV_uk)
  
  return(mat_df)
  
}
  

# mat_df <- maternal_output_func(1,1,"rdt", "rdt", FALSE, syph_infant_deaths)

#write.csv(mat_df, here("output", "maternal_model_base.csv"), row.names = FALSE)



## Needs
  # HIV and non-HIV mortality rates (keeping them the same here)














