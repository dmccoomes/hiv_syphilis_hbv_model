

#### Title: Adult outcomes model
#### Purpose: Adult outcomes with time horizon
#### Author: David Coomes
#### Date: July 17, 2025


## Purpose: this script will 


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table"
)

rm(list = ls())


# We want to take the output of the hiv-syphilis model and run the adults through a time horizon (here 20 years)
  # Populations:
    # HIV infected - on treatment, not on treatment, unknown status
    # syphilis infected 
    # tertiary syphilis
    # HIV-syphilis coinfected - on treatment, not on treatment, unknown status
    # Uninfected
    # STARTING states from HIV-syphilis maternal model after 1 year, year 0 is birth
      # Keep maternal deaths - subtract out infant syphilis deaths
  # Transition probabilities
    # treatment <--> no treatment
    # testing rates
    # testing times
    # transition to tertiary syphilis
    # mortality
      # HIV
      # non-HIV
      # tertiary syphilis
      # 


## Steps:
  # Extract 0 year and 1 year population estimates from the HIV-syphilis maternal models
    # load parameters
    # Run HIV and syphilis maternal models and combine
  # Advance these over 20 years based on our transition probabilities
  # Calculate DALYs based on disability weights


# This is the adult model - 
  # where we take an adult population and estimate the impacts of 10 years of testing on 20 years of outcomes
  # we will use a dynamic Markov model where we change the force of infection based on the number of infectious individuals



# Questions:
  # Can we incorporate interactions between key populations and general populations? 


# Needs:
  # Incorporate ART treatment in model
    # Incidence rate is a function of ART use. Assuming x% reduction in the probability of transmission 
    # due to ART use, we can apportion the incidence rate into from those on ART and those not on ART
  # Incorportate testing strategies in model
  # Incorporate some kind of weighting for mortality based on age groups? 


# Getting some population parameters - delete this after uploaded to excel file
hiv_inc <- 0.00053
adult_pop <- 33971376
hiv_prev_adult <- 0.03
hiv_red_ART <- 0.93


# Load mortality file
mort_df <- read_csv(here("parameters", "kenya_mort_ihme.csv"))
# Getting average adult mortality with and without HIV for now
mort_est <- mort_df %>%
  filter(!"Age Group" %in% c("<1 year", "1-4 years", "5-9 years", "10-14 years")) %>%
  group_by(Indicator) %>%
  summarize(avg_mort_both = mean(`Both sexes`, na.rm = TRUE),
            avg_mort_male = mean(Male, na.rm = TRUE),
            avg_mort_female = mean(Female, na.rm = TRUE))

hiv_mort <- mort_est$avg_mort_both[mort_est$Indicator == "mortality_HIV+_total"]
nhiv_mort <- mort_est$avg_mort_both[mort_est$Indicator == "mortality_HIV-_total"]


hiv_inc_apport <- hiv_inc*


# Starting with very simple dynamic Markov model
adult_matrix <- matrix(nrow = 21, ncol = 6, 
                          dimnames = list(NULL, 
                                          c("year", "hiv_neg", "hiv_pos", 
                                            "dead_hiv", "dead_nhiv", "total")))

adult_df <- as.data.frame(adult_matrix)

# Set up starting states
adult_df[1, ] <- c(0, adult_pop*(1-hiv_prev_adult), adult_pop*hiv_prev_adult, rep(0, 2), adult_pop)

# Working on no intervention model
for (x in 1:20) {
  
  adult_df$year[x+1] <- x 
  
  # No HIV 
  adult_df$hiv_neg[x+1] <- adult_df$hiv_neg[x] - adult_df$hiv_neg[x]*nhiv_mort - adult_df$hiv_neg[x]*hiv_inc
  
  # HIV positive
  adult_df$hiv_pos[x+1] <- adult_df$hiv_pos[x] + adult_df$hiv_neg[x]*hiv_inc - adult_df$hiv_pos[x]*hiv_mort
  
  # HIV deaths
  adult_df$dead_hiv[x+1] <- adult_df$dead_hiv[x] + adult_df$hiv_pos[x]*hiv_mort
  
  # nHIV deaths
  adult_df$dead_nhiv[x+1] <- adult_df$dead_nhiv[x] + adult_df$hiv_neg[x]*nhiv_mort
  
  # Get total numbers 
  adult_df$total[x+1] <- sum(adult_df[(x+1), 2:5])
  

}



































