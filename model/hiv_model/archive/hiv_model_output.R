
#### Title: HIV retesting model output
#### Purpose: Run HIV retesting models and collect output
#### Author: David Coomes
#### Date: February 25, 2025


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table"
)

rm(list = ls())
source(here("code", "test_frequency.R"))

prev_reduction <- 0
incidence_reduction <- 0
PrEP <- 0

# Set population 
country <- "Ukraine"
population <- "FSW"


# Get data frame of all needed output
output_matrix <- matrix(nrow = 8, ncol = 9, 
                          dimnames = list(NULL, 
                                          c("infant_infections", "tot_inf_infections_averted", "pct_inf_infections_averted", 
                                            "infant_deaths", "total_cost", 
                                            "incremental_cost", "total_dalys", "incremental_dalys_averted",
                                            "icer")))

output_df <- as.data.frame(output_matrix)

rownames(output_df) <- c("No retesting", "late_ANC", "late_ANC_14wk", "late_ANC_6mo", "late_ANC_9mo", 
                         "late_ANC_14wk_6mo", "late_ANC_14wk_9mo", "late_ANC_every3mo")




## Get number of infant infections
test_freq_list <- c("No retesting", "Late ANC", "Late ANC + 14 weeks", "Late ANC + 6 months", "Late ANC + 9 months",
                    "Late ANC + 14 weeks + 6 months", "Late ANC + 14 weeks + 9 months", "Late ANC + every 3 months")

infant_infection_output <- vector()
infant_death_output <- vector()
total_cost_output <- vector()
total_dalys_output <- vector()

for (x in 1:length(test_freq_list)){

  test_list <- test_freq_function(test_freq_list[x])

  first_ANC_test <- test_list[1]
  late_ANC_test <- test_list[2]
  pp_early_test <- test_list[3]
  pp_14weeks_test <- test_list[4]
  pp_mid_test <- test_list[5]
  pp_9months_test <- test_list[6]

  # running model
  # source(here("code", "input_parameters.R"))
  # source(here("code", "lifetable.R"))
  # source(here("code", "hiv_model.R"))
  
  
  source(here("code", "input_parameters.R"), local=TRUE)
  source(here("code", "lifetable.R"), local=TRUE)
  source(here("code", "hiv_model.R"), local=TRUE)
  source(here("code", "infant_model.R"), local=TRUE)
  
  if (x == 1) {
    infant_infection_output <- total_infected_infants
    infant_death_output <- infant_deaths
    total_cost_output <- total_maternal_costs + total_infant_costs
    total_dalys_output <- total_infant_dalys
  } else {
    infant_infection_output <- c(infant_infection_output, total_infected_infants)
    infant_death_output <- c(infant_death_output, infant_deaths)
    total_cost_output <- c(total_cost_output, total_maternal_costs + total_infant_costs)
    total_dalys_output <- c(total_dalys_output, total_infant_dalys)
  }
}


output_df$infant_deaths <- infant_death_output
output_df$total_dalys <- total_dalys_output
output_df$infant_infections <- infant_infection_output
output_df$total_cost <- total_cost_output

# Arrange cost in order with 'No retesting at the top'
output_df <- output_df["No retesting", ] %>%
  rbind(output_df %>%
          filter(!row.names(output_df) == "No retesting") %>%
          arrange(total_cost)) %>%
  mutate(tot_inf_infections_averted = first(infant_infections) - infant_infections,
         pct_inf_infections_averted = tot_inf_infections_averted / first(infant_infections)*100,
         incremental_cost = total_cost - lag(total_cost),
         incremental_dalys_averted = lag(total_dalys) - total_dalys,
         icer = incremental_cost / incremental_dalys_averted)






