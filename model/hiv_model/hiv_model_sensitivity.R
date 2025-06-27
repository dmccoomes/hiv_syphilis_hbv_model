
####
####
####
####


pacman::p_load(
  "tidyverse", 
  "here",
  "data.table"
)

rm(list = ls())
source(here("code", "test_frequency.R"))

# Starting states
prev_reduction <- 0
incidence_reduction <- 0
PrEP <- 0


## Get number of infant infections
test_freq_list <- list("No retesting", "Late ANC")
#population_list <- c("Kenya FSW", "Kenya PWID", "Kenya serodiscordant")
#population_list <- c("Kenya serodiscordant")

sens_func <- function(test_freq_list, prep_pct, inc_red, prev_red) {
  
  test_list <- test_freq_function(test_freq_list)
  
  first_ANC_test <<- test_list[1]
  late_ANC_test <<- test_list[2]
  pp_early_test <<- test_list[3]
  pp_14weeks_test <<- test_list[4]
  pp_mid_test <<- test_list[5]
  pp_9months_test <<- test_list[6]
  
  print(late_ANC_test)
  
  PrEP <- prep_pct
  incidence_reduction <- inc_red
  
  prev_reduction <- prev_red
  #hiv_prev <- hiv_prev*prevalence_reduction
  test_out <- test_freq_list
  
  source(here("code", "input_parameters.R"), local=TRUE)
  source(here("code", "lifetable.R"), local=TRUE)
  source(here("code", "hiv_model.R"), local=TRUE)
  source(here("code", "infant_model.R"), local=TRUE)
  
  total_cost_output <- total_maternal_costs + total_infant_costs
  total_dalys_output <- total_infant_dalys
  test_freq <- test_out
  population_out <- population
  prep_out <- PrEP
  inc_red_out <- incidence_reduction
  prev_red_out <- prev_reduction
  
  df <- data.frame(total_cost = total_cost_output, 
                   total_dalys = total_dalys_output,
                   test_freq = test_freq, 
                   population = population_out,
                   prep = PrEP,
                   inc_red = inc_red_out, 
                   prev_red = prev_red_out) 
  
  return(df)
  
}

#sens_df1 <- sens_func(test_freq_list=test_freq_list[2], prep_pct=0, inc_red=0.5, prev_red=0.3)

# Set up list of prep, incidence, and prevalence reduction
prep_list <- as.list(seq(0, 0.5, 0.1))
inc_list <- as.list(seq(0, 0.5, 0.1))
prev_list <- as.list(seq(0, 0.5, 0.1))
# Get every combination of prep, incidence, and prevalence
args <- expand.grid(test_freq_list = test_freq_list, prep_list = prep_list, inc_list = inc_list, prev_list = prev_list)

# sens_df1 <- mapply(test_freq_list=test_freq_list, prep_pct=prep_list, inc_red=inc_list, prev_red=prev_list, sens_func,
#                    SIMPLIFY = TRUE)

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_df1 <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                   inc_red=args$inc_list, prev_red=args$prev_list,
                   SIMPLIFY = TRUE) %>%
  t()

sens_df2 <- sens_df1

# Getting ICER
sens_df2 <- sens_df2 %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


## Plotting
sens_df2 %>%
  mutate(inc_red = factor(inc_red),
         prev_red = factor(prev_red)) %>%
  ggplot() +
  geom_point(aes(x=icer, y=prep, color=inc_red, shape=prev_red)) +
  theme_bw()







# Run for each population
country <- "Kenya"
population <- "FSW"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_kenya_fsw <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                   inc_red=args$inc_list, prev_red=args$prev_list,
                   SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


population <- "PWID"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_kenya_pwid <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                         inc_red=args$inc_list, prev_red=args$prev_list,
                         SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


population <- "serodiscordant"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_kenya_sero <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                          inc_red=args$inc_list, prev_red=args$prev_list,
                          SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")

# Put all Kenya populations together
sens_kenya <- sens_kenya_fsw %>%
  bind_rows(sens_kenya_pwid) %>% 
  bind_rows(sens_kenya_sero)

# Plot all Kenya 
## Plotting
sens_kenya %>%
  mutate(inc_red = factor(inc_red),
         prev_red = factor(prev_red),
         population = as.character(population)) %>%
  ggplot() +
  geom_point(aes(x=icer, y=prep, color=inc_red, shape=prev_red)) +
  geom_line(aes(x=icer, y=prep, , group=interaction(inc_red, prev_red), color=inc_red), alpha=0.25) +
  geom_vline(xintercept=500, color="red", linetype="dotted") +
  theme_bw() +
  xlab("ICER") + ylab("PrEP use") +
  guides(color=guide_legend(title="Incidence \nreduction"),
         shape=guide_legend(title="Prevalence \nreduction")) +
  facet_wrap(~population, ncol=1) 
  

# sens_kenya_fsw %>%
#   mutate(inc_red = factor(inc_red),
#          prev_red = factor(prev_red),
#          population = as.character(population)) %>%
#   ggplot() +
#   geom_point(aes(x=icer, y=prep, color=inc_red, shape=prev_red)) +
#   theme_bw() 


##### 
## South Africa 
#####
country <- "South Africa"
population <- "FSW"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_sa_fsw <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                         inc_red=args$inc_list, prev_red=args$prev_list,
                         SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


population <- "PWID"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_sa_pwid <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                          inc_red=args$inc_list, prev_red=args$prev_list,
                          SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


population <- "serodiscordant"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_sa_sero <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                          inc_red=args$inc_list, prev_red=args$prev_list,
                          SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")

# Put all Kenya populations together
sens_sa <- sens_sa_fsw %>%
  bind_rows(sens_sa_pwid) %>% 
  bind_rows(sens_sa_sero)

# Plot all Kenya 
## Plotting
sens_sa %>%
  mutate(inc_red = factor(inc_red),
         prev_red = factor(prev_red),
         population = as.character(population)) %>%
  ggplot() +
  geom_point(aes(x=icer, y=prep, color=inc_red, shape=prev_red)) +
  geom_line(aes(x=icer, y=prep, , group=interaction(inc_red, prev_red), color=inc_red), alpha=0.25) +
  geom_vline(xintercept=750, color="red", linetype="dotted") +
  theme_bw() +
  xlim(0,1800) +
  xlab("ICER") + ylab("PrEP use") +
  guides(color=guide_legend(title="Incidence \nreduction"),
         shape=guide_legend(title="Prevalence \nreduction")) +
  facet_wrap(~population, ncol=1)



##### 
## Ukraine
#####
country <- "Ukraine"
population <- "FSW"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_ukraine_fsw <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                      inc_red=args$inc_list, prev_red=args$prev_list,
                      SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


population <- "PWID"

# Use function to go through every combination and output for no retesting and late ANC retesting
sens_ukraine_pwid <- mapply(sens_func, test_freq_list=test_freq_list, prep_pct=args$prep_list, 
                       inc_red=args$inc_list, prev_red=args$prev_list,
                       SIMPLIFY = TRUE) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(c(total_cost:total_dalys, prep:prev_red), ~as.numeric(unlist(.)))) %>%
  group_by(population, prep, inc_red, prev_red) %>%
  mutate(incr_cost = total_cost - lag(total_cost),
         incr_daly = lag(total_dalys) - total_dalys,
         icer = incr_cost / incr_daly) %>%
  filter(test_freq=="Late ANC")


# Put all Ukraine populations together
sens_ukraine <- sens_ukraine_fsw %>%
  bind_rows(sens_ukraine_pwid) 

# Plot all Kenya 
## Plotting
sens_ukraine %>%
  mutate(inc_red = factor(inc_red),
         prev_red = factor(prev_red),
         population = as.character(population)) %>%
  ggplot() +
  geom_point(aes(x=icer, y=prep, color=inc_red, shape=prev_red)) +
  geom_line(aes(x=icer, y=prep, , group=interaction(inc_red, prev_red), color=inc_red), alpha=0.25) +
  geom_vline(xintercept=1000, color="red", linetype="dotted") +
  theme_bw() +
  xlab("ICER") + ylab("PrEP use") +
  guides(color=guide_legend(title="Incidence \nreduction"),
         shape=guide_legend(title="Prevalence \nreduction")) +
  facet_wrap(~population, ncol=1)






















#########
# OLD CODE ##
############

#Get data frame of all needed output
output_matrix <- matrix(nrow = 8, ncol = 9,
                        dimnames = list(NULL,
                                        c("infant_infections", "tot_inf_infections_averted", "pct_inf_infections_averted",
                                          "infant_deaths", "total_cost",
                                          "incremental_cost", "total_dalys", "incremental_dalys_averted",
                                          "icer")))

output_df <- as.data.frame(output_matrix)

rownames(output_df) <- c("No retesting", "late_ANC", "late_ANC_14wk", "late_ANC_6mo", "late_ANC_9mo",
                         "late_ANC_14wk_6mo", "late_ANC_14wk_9mo", "late_ANC_every3mo")



total_cost_output <- vector()
total_dalys_output <- vector()
test_freq <- vector()
population_out <- vector()
prep_out <- vector()
inc_red_out <- vector()
prev_red_out <- vector()

for (i in 1:length(population_list)) {
  
  population <- population_list[i]
  
  for (x in 1:length(test_freq_list)){
    
    test_list <- test_freq_function(test_freq_list[x])
    
    print(population)
    print(test_freq_list[x])
    
    pop_out <- population
    test_out <- test_freq_list[x]
    
    first_ANC_test <- test_list[1]
    late_ANC_test <- test_list[2]
    pp_early_test <- test_list[3]
    pp_14weeks_test <- test_list[4]
    pp_mid_test <- test_list[5]
    pp_9months_test <- test_list[6]
    
    source(here("code", "input_parameters.R"), local=TRUE)
    source(here("code", "lifetable.R"), local=TRUE)
    
    for (prep_pct in seq(0, 0.2, 0.1)) {
      
      PrEP <- prep_pct
      print(PrEP)
      
      for (inc_red in seq(0, 0.2, 0.1)) {
        
        incidence_reduction <- inc_red
        print(incidence_reduction)
        
        inc_preg_early <- inc_preg_early*(1-incidence_reduction)
        inc_preg_late <- inc_preg_late*(1-incidence_reduction)
        inc_pp_early <- inc_pp_early*(1-incidence_reduction)
        inc_pp_mid <- inc_pp_mid*(1-incidence_reduction)
        inc_pp_late <- inc_pp_mid*(1-incidence_reduction)
        
        for (prev_red in seq(0, 0.2, 0.1)) {
          
          prevalence_reduction <- prev_red
          print(prevalence_reduction) 
          
          hiv_prev <- hiv_prev*prevalence_reduction
          
          source(here("code", "infant_model.R"), local=TRUE)
          
          if (i == 1) {
            total_cost_output <- total_maternal_costs + total_infant_costs
            total_dalys_output <- total_infant_dalys
            test_freq <- test_out
            population_out <- population
            prep_out <- PrEP
            inc_red_out <- incidence_reduction
            prev_red_out <- prevalence_reduction
          } else {
            total_cost_output <- c(total_cost_output, total_maternal_costs + total_infant_costs)
            total_dalys_output <- c(total_dalys_output, total_infant_dalys)
            test_freq <- c(test_freq, test_out)
            population_out <- c(population_out, population)
            prep_out <- c(prep_out, PrEP)
            inc_red_out <- c(inc_red_out, incidence_reduction)
            prev_red_out <- c(prev_red_out, prevalence_reduction)
          }
        }
      }
    }
  }
}



# Put together sensitivity output
sens_df <- data.frame(total_cost=total_cost_output, total_dalys=total_dalys_output,
                      test_freq=test_freq, population=population_out, prep=prep_out, 
                      inc_red=inc_red_out, prev_red=prev_red_out)

sens_df <- sens_df %>%
  group_by(population, prep, inc_red, prev_red) %>%
  arrange(population, prep, inc_red, prev_red) %>%
  mutate(incremental_cost = total_cost - lag(total_cost),
         incremental_dalys = lag(total_dalys) - total_dalys,
         icer = incremental_cost / incremental_dalys) %>%
  filter(test_freq == "Late ANC") %>%
  select(population, icer, prep, inc_red, prev_red)



for (i in 1:length(population_list)) {
  
  population <- population_list[i]
  
  for (x in 1:length(test_freq_list)){
    
    test_list <- test_freq_function(test_freq_list[x])
    
    print(population)
    print(test_freq_list[x])
    
    pop_out <- population
    test_out <- test_freq_list[x]
    
    first_ANC_test <- test_list[1]
    late_ANC_test <- test_list[2]
    pp_early_test <- test_list[3]
    pp_14weeks_test <- test_list[4]
    pp_mid_test <- test_list[5]
    pp_9months_test <- test_list[6]
    
    source(here("code", "input_parameters.R"))
    source(here("code", "lifetable.R"))
    
    for (prep_pct in seq(0, 0.5, 0.05)) {
      
      PrEP <- prep_pct
      print(PrEP)
      
      for (inc_red in seq(0, 0.5, 0.05)) {
        
        incidence_reduction <- inc_red
        print(incidence_reduction)
        
        inc_preg_early <- inc_preg_early*(1-incidence_reduction)
        inc_preg_late <- inc_preg_late*(1-incidence_reduction)
        inc_pp_early <- inc_pp_early*(1-incidence_reduction)
        inc_pp_mid <- inc_pp_mid*(1-incidence_reduction)
        inc_pp_late <- inc_pp_mid*(1-incidence_reduction)
        
        for (prev_red in seq(0, 0.5, 0.05)) {
          
          prevalence_reduction <- prev_red
          print(prevalence_reduction) 
          
          hiv_prev <- hiv_prev*prevalence_reduction
          
          source(here("code", "infant_model.R"))
          
          if (i == 1) {
            total_cost_output <- total_maternal_costs + total_infant_costs
            total_dalys_output <- total_infant_dalys
            test_freq <- test_out
            population_out <- population
            prep_out <- PrEP
            inc_red_out <- incidence_reduction
            prev_red_out <- prevalence_reduction
          } else {
            total_cost_output <- c(total_cost_output, total_maternal_costs + total_infant_costs)
            total_dalys_output <- c(total_dalys_output, total_infant_dalys)
            test_freq <- c(test_freq, test_out)
            population_out <- c(population_out, population)
            prep_out <- c(prep_out, PrEP)
            inc_red_out <- c(inc_red_out, incidence_reduction)
            prev_red_out <- c(prev_red_out, prevalence_reduction)
          }
        }
      }
    }
  }
}













