


test_freq_function <- function(test_frequency) {
  
  if(test_frequency=="No retesting") {
    
    first_ANC_test <- 1
    late_ANC_test <- 0
    pp_early_test <- 0
    pp_14weeks_test <- 0
    pp_mid_test <- 0
    pp_9months_test <- 0
    
  } else if(test_frequency=="Late ANC") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 0
    pp_mid_test <- 0
    pp_9months_test <- 0
    
  } else if(test_frequency=="Late ANC + 14 weeks") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 1
    pp_mid_test <- 0
    pp_9months_test <- 0
    
    
  } else if(test_frequency=="Late ANC + 6 months") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 0
    pp_mid_test <- 1
    pp_9months_test <- 0
    
  } else if(test_frequency=="Late ANC + 9 months") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 0
    pp_mid_test <- 0
    pp_9months_test <- 1
    
  } else if(test_frequency=="Late ANC + 14 weeks + 6 months") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 1
    pp_mid_test <- 1
    pp_9months_test <- 0
    
  } else if(test_frequency=="Late ANC + 14 weeks + 9 months") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 1
    pp_mid_test <- 0
    pp_9months_test <- 1
    
  } else if(test_frequency=="Late ANC + every 3 months") {
    
    first_ANC_test <- 1
    late_ANC_test <- 1
    pp_early_test <- 1
    pp_14weeks_test <- 1
    pp_mid_test <- 1
    pp_9months_test <- 1
    
  } else {
    
    print("Invalid test choice")
    
  }
  
  return(c(first_ANC_test, late_ANC_test, pp_early_test, pp_14weeks_test, pp_mid_test, pp_9months_test))
  
}
