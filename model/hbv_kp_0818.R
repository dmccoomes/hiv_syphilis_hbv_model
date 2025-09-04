rm(list=ls())
library(tidyverse); library(readr); library(stringr)


# 0) LOAD INPUTS (KP ONLY)
params_kp <- read_csv("hbv_model_parameter_table.csv", show_col_types = FALSE) %>%
  filter(Population == "KP", !is.na(`Value PE`)) %>%
  group_by(Symbol) %>%
  slice_tail(n = 1) %>%   # if duplicates, keep last
  ungroup()

if (nrow(params_kp) == 0) {
  stop("No rows with Population == 'KP' found in hbv_model_parameter_table.csv")
}

# Make it a LIST so [[ "name" ]] is safe
param_vals <- params_kp %>%
  select(Symbol, `Value PE`) %>%
  deframe() %>%
  as.list()

# Robust getter
p <- function(sym, default) {
  if (is.null(param_vals) || is.null(param_vals[[sym]]) || is.na(param_vals[[sym]])) {
    return(default)
  }
  as.numeric(param_vals[[sym]])
}


# Mortality table (same structure as pregnancy build)
mort_df <- read_csv("mort_nonhiv.csv", show_col_types = FALSE) %>%
  rename(age_group = `Age Group`, mort_both = `Both sexes`, mort_female = Female) %>%
  mutate(age_start = as.numeric(str_extract(age_group, "^\\d+")),
         age_start = if_else(str_detect(age_group, "<1"), 0, age_start))

# ---- Background mortality by age; default to BOTH sexes for KP ----
get_mortality <- function(age, sex = c("both","female")) {
  sex <- match.arg(sex)
  start_age <- if (age < 1) 0 else {
    possible <- sort(unique(na.omit(mort_df$age_start)))
    max(possible[possible <= age])
  }
  row <- mort_df %>% dplyr::filter(age_start == start_age)
  if (sex == "both") as.numeric(row$mort_both) else as.numeric(row$mort_female)
}




# 1) States

time_horizon <- 20
years <- 0:time_horizon

# Testing coverage (constants for Step 1; we can swap to vectors/ramps later)
prop_tested_list <- list(
  SOC = p("hbv_cov_soc",  0.31),   # fallback to 0.31 if missing in KP sheet
  INT = p("hbv_cov_int",  0.749)   # fallback to 0.749 if missing in KP sheet
)

coverage_by_year <- function(t, arm) {
  if (arm == "SOC") return(prop_tested_list$SOC)
  if (t >= 1 && t <= 10) prop_tested_list$INT else prop_tested_list$SOC
}


# --- Test profiles by arm and year ---
test_profile <- function(t, arm) {
  use_triplex <- (arm == "INT" && t >= 1 && t <= 10)
  
  list(
    type = if (use_triplex) "triplex" else "hbsag",
    sens = if (use_triplex) p("triplex_sens", 0.98) else p("hbsag_sens", 0.95),
    spec = if (use_triplex) p("triplex_spec", 0.99) else p("hbsag_spec", 0.98),
    cost = if (use_triplex) p("triplex_cost", 3.00)  else p("hbv_c_test_hbsag", 0.87)
  )
}

# Adult HBV states (single chronic bucket for reporting; internal on/off tx rows)
adult_states <- c(
  "HBV_susceptible",
  "Acute_symptomatic_HBV",
  "HBV_immune",
  "Chronic_HBV_nontx",
  "Chronic_HBV_tx",
  "Compensated_cirrhosis",
  "Decompensated_cirrhosis",
  "Liver_cancer",
  "Death"
)


# 2) INITIALIZATION

N0       <- p("hbv_kp_N0", 500000)
p_sus0   <- p("hbv_p_sus",      0.75)
p_acute0 <- p("hbv_p_acute",   0.001)
p_chron0 <- p("hbv_p_chronic", 0.05)
p_ctx0   <- p("hbv_p_C_tx0",    0.30)
p_CC0    <- p("hbv_p_CC0",      0.00)
p_DC0    <- p("hbv_p_DC0",      0.00)
p_HCC0   <- p("hbv_p_HCC0",     0.000)

init_vector <- function() {
  Ctot <- N0 * p_chron0
  Ctx  <- Ctot * p_ctx0
  Cntx <- Ctot - Ctx
  v <- setNames(numeric(length(adult_states)), adult_states)
  v["HBV_susceptible"]       <- N0 * p_sus0
  v["Acute_symptomatic_HBV"] <- N0 * p_acute0
  v["HBV_immune"]            <- N0 * (1 - (p_sus0 + p_acute0 + p_chron0))
  v["Chronic_HBV_nontx"]     <- Cntx
  v["Chronic_HBV_tx"]        <- Ctx
  v["Compensated_cirrhosis"] <- N0 * p_CC0
  v["Decompensated_cirrhosis"] <- N0 * p_DC0
  v["Liver_cancer"]          <- N0 * p_HCC0
  v["Death"] <- 0
  v
}

# 3) FOI & NATURAL HISTORY PARAMETERS

# Semi-dynamic FOI components
hbv_beta      <- p("hbv_beta",         1.6) #0.9
hbv_mix_kp    <- p("hbv_mix_kp",       1.2) # 1
r_beta_treat  <- p("hbv_r_beta_treat", 0.05)  #0.1 # infectiousness multiplier on treatment
foi_floor     <- p("hbv_foi_floor",    0.0)

# Acute outcomes (no waning; acute clears within the year)
r_acute_to_immune <- p("hbv_r_acute_to_immune", 0.90)
r_acute_to_chronic<- p("hbv_r_acute_to_chronic",0.09)

# Progression (annual)
r_C_to_CC   <- p("hbv_r_chronic_to_comp_cirr", 0.0029)
r_CC_to_DC  <- p("hbv_r_comp_to_decomp",       0.04)
r_CC_to_HCC <- p("hbv_r_comp_to_cancer",       0.01)

r_DC_to_D  <- p("hbv_r_decomp_death", 0.20)  # Decompensated cirrhosis → Death
r_HCC_to_D <- p("hbv_r_cancer_death",    0.35)  # Liver cancer → Death
r_DC_to_HCC <- p("hbv_r_decomp_to_cancer", 0.03)


# Step 1: no explicit disease death; background mortality only for all alive states.

# Treatment effect (slows C -> CC)
hr_treat_C_to_CC <- p("hbv_hr_prog_treat_C_to_CC", 0.50)

# Testing cascade + vaccination (no waning)
linkage <- p("hbv_linkage",    0.85)
tx_init <- p("hbv_tx_init",    0.85)
tx_pers <- p("hbv_tx_persist", 0.85)
vax_full<- p("hbv_vax_full",   0.60)

# KP starting age for background mortality
kp_age0 <- p("hbv_kp_age0", 25)


# 4) TRANSITION MATRIX FOR YEAR t (built from CURRENT state)

build_tm <- function(state_vec, t) {
  tm <- matrix(0, nrow = length(adult_states), ncol = length(adult_states),
               dimnames = list(from = adult_states, to = adult_states))
  
  alive_states <- setdiff(adult_states, "Death")
  alive <- sum(state_vec[alive_states])
  
  # --- FOI (hazard-consistent probability) ---
  C_tx  <- state_vec["Chronic_HBV_tx"]
  C_ntx <- state_vec["Chronic_HBV_nontx"]
  prev_eff   <- if (alive > 0) (C_ntx + r_beta_treat * C_tx) / alive else 0
  lambda_raw <- max(0, hbv_beta * hbv_mix_kp * prev_eff + foi_floor)
  p_inf      <- 1 - exp(-lambda_raw)  # annual infection probability
  
  # --- Natural history BEFORE background mortality (rows will be renormalized) ---
  # S -> A
  tm["HBV_susceptible", "Acute_symptomatic_HBV"] <- p_inf
  tm["HBV_susceptible", "HBV_susceptible"]       <- 1 - p_inf
  
  # A -> {Immune, Chronic_nontx} (no self-loop)
  sA <- r_acute_to_immune + r_acute_to_chronic
  rA_i <- r_acute_to_immune / sA
  rA_c <- r_acute_to_chronic / sA
  tm["Acute_symptomatic_HBV","HBV_immune"]        <- rA_i
  tm["Acute_symptomatic_HBV","Chronic_HBV_nontx"] <- rA_c
  
  # Chronic (on/off tx differ via HR)
  tm["Chronic_HBV_nontx","Compensated_cirrhosis"] <- r_C_to_CC
  tm["Chronic_HBV_nontx","Chronic_HBV_nontx"]     <- 1 - r_C_to_CC
  tm["Chronic_HBV_tx","Compensated_cirrhosis"]    <- r_C_to_CC * hr_treat_C_to_CC
  tm["Chronic_HBV_tx","Chronic_HBV_tx"]           <- 1 - r_C_to_CC * hr_treat_C_to_CC
  
  # CC -> {DC, HCC, stay}
  stay_CC <- max(0, 1 - (r_CC_to_DC + r_CC_to_HCC))
  tm["Compensated_cirrhosis","Decompensated_cirrhosis"] <- r_CC_to_DC
  tm["Compensated_cirrhosis","Liver_cancer"]            <- r_CC_to_HCC
  tm["Compensated_cirrhosis","Compensated_cirrhosis"]   <- stay_CC
  
  #  (progression + stay; death added later by competing-risks renorm)
  tm["Decompensated_cirrhosis","Liver_cancer"]            <- r_DC_to_HCC
  tm["Decompensated_cirrhosis","Decompensated_cirrhosis"] <- 1 - r_DC_to_HCC
  
  tm["Liver_cancer","Liver_cancer"] <- 1
  
  # Immune stays immune (no waning)
  tm["HBV_immune","HBV_immune"] <- 1
  
  # Death absorbing
  tm["Death","Death"] <- 1
  
  # --- Competing-risks row normalization (background + disease death) ---
  renorm_row_CR <- function(mat, st, p_bg, p_dis) {
    # Survive both background and disease death this cycle:
    surv_mass <- (1 - p_bg) * (1 - p_dis)      # ∈ [0,1]
    non_death <- setdiff(colnames(mat), "Death")
    s <- sum(mat[st, non_death])
    if (s > 0) {
      mat[st, non_death] <- mat[st, non_death] / s * surv_mass
    } else {
      mat[st, st] <- surv_mass
    }
    mat[st, "Death"] <- 1 - surv_mass          # total death prob (no double-count)
    mat
  }
  
  # --- Apply background (Both sexes) + disease-specific deaths by state ---
  age_any <- kp_age0 + t
  p_bg    <- get_mortality(age_any, "both")    # KP uses BOTH sexes
  
  for (st in setdiff(adult_states, "Death")) {
    p_dis <- switch(st,
                    "Decompensated_cirrhosis" = r_DC_to_D,
                    "Liver_cancer"            = r_HCC_to_D,
                    # "Compensated_cirrhosis" = r_CC_to_D,   # uncomment if you added CC→D
                    0
    )
    tm <- renorm_row_CR(tm, st, p_bg = p_bg, p_dis = p_dis)
  }
  
  # Ensure Death remains absorbing
  tm["Death", ] <- 0; tm["Death", "Death"] <- 1
  
  
  # Expected new infections (for incidence reporting)
  new_inf <- state_vec["HBV_susceptible"] * p_inf
  
  list(tm = tm, p_inf = p_inf, lambda_raw = lambda_raw, prev_eff = prev_eff, new_inf = new_inf)
}


# 5) TESTING & VACCINATION (post-TM; INT active in years 1–10)

apply_testing_vax <- function(vec, t, arm) {
  # totals
  alive_states <- setdiff(names(vec), "Death")
  alive <- sum(vec[alive_states])
  
  cov <- coverage_by_year(t, arm)
  n_tested <- alive * cov
  
  prof <- test_profile(t, arm)  # <-- sens/spec/cost switch here
  
  # Approx prevalence among those tested: HBsAg-positive ~ Chronic + (some Acute)
  C_tot <- vec["Chronic_HBV_nontx"] + vec["Chronic_HBV_tx"]
  A_t   <- vec["Acute_symptomatic_HBV"]
  prev_hbsag <- (C_tot + A_t) / max(alive, 1e-9)
  
  # Test outcomes
  n_inf_tested   <- n_tested * prev_hbsag
  n_uninf_tested <- n_tested - n_inf_tested
  
  tp <- n_inf_tested   * prof$sens
  fn <- n_inf_tested   - tp
  fp <- n_uninf_tested * (1 - prof$spec)
  tn <- n_uninf_tested - fp
  
  # Of infected positives, proportion that are chronic (not acute)
  share_chronic_among_inf <- ifelse((C_tot + A_t) > 0, C_tot / (C_tot + A_t), 1)
  detected_chronic <- tp * share_chronic_among_inf
  
  # Among detected chronic, share that are untreated at baseline of year t
  share_nontx <- ifelse(C_tot > 0, vec["Chronic_HBV_nontx"] / C_tot, 1)
  detected_nontx <- detected_chronic * share_nontx
  
  # Care cascade
  linkage <- p("hbv_linkage", 0.80)
  tx_init <- p("hbv_tx_init", 0.85)
  tx_pers <- p("hbv_tx_persist", 0.90)
  
  newly_tx <- detected_nontx * linkage * tx_init * tx_pers
  mov_tx <- min(newly_tx, vec["Chronic_HBV_nontx"])
  vec["Chronic_HBV_nontx"] <- vec["Chronic_HBV_nontx"] - mov_tx
  vec["Chronic_HBV_tx"]    <- vec["Chronic_HBV_tx"] + mov_tx
  
  # Vaccinate among HBsAg-negative tested (tn + fn). Allocate to susceptibles only.
  vax_full <- p("hbv_vax_full", 0.70)
  negs <- tn + fn
  newly_vax <- 0
  if (negs > 0) {
    non_hbsag_pool <- alive - (C_tot + A_t)
    sus_share <- ifelse(non_hbsag_pool > 0, vec["HBV_susceptible"] / non_hbsag_pool, 0)
    newly_vax <- negs * vax_full * sus_share
    mov_vax <- min(newly_vax, vec["HBV_susceptible"])
    vec["HBV_susceptible"] <- vec["HBV_susceptible"] - mov_vax
    vec["HBV_immune"]      <- vec["HBV_immune"] + mov_vax
  }
  
  list(
    vec = vec,
    n_tested = as.numeric(n_tested),
    tp = as.numeric(tp), fp = as.numeric(fp), tn = as.numeric(tn), fn = as.numeric(fn),
    newly_tx = as.numeric(mov_tx),
    newly_vax = as.numeric(newly_vax),
    test_cost = as.numeric(n_tested * prof$cost),
    test_type = prof$type,
    sens = prof$sens, spec = prof$spec
  )
}



# 6) RUNNER (one arm)

run_arm <- function(arm = "SOC") {
  # start-of-year 0 state
  sv <- init_vector()

  # results storage (t = 0..time_horizon)
  results <- matrix(0, nrow = time_horizon + 1, ncol = length(adult_states),
                    dimnames = list(time = 0:time_horizon, state = adult_states))
  results[1, ] <- sv

  # diagnostics (one row per year; row index = t + 1)
  diag <- tibble(
    time       = years,
    arm        = rep(arm, length(years)),
    p_inf      = NA_real_,
    lambda_raw = NA_real_,
    prev_eff   = NA_real_,
    incidence  = NA_real_,
    n_tested   = NA_real_,
    tp = NA_real_, fp = NA_real_, tn = NA_real_, fn = NA_real_,
    newly_tx   = NA_real_,
    newly_vax  = NA_real_,
    test_cost  = NA_real_,
    test_type  = NA_character_,
    sens = NA_real_, spec = NA_real_
  )

  # optional: keep transition matrices if you want to validate row sums later
  tm_list <- vector("list", time_horizon)

  # baseline (t = 0) diagnostics
  diag$p_inf[1]       <- 0
  diag$lambda_raw[1]  <- 0
  diag$prev_eff[1]    <- 0
  diag$incidence[1]   <- 0
  diag$n_tested[1]    <- 0
  diag$newly_tx[1]    <- 0
  diag$newly_vax[1]   <- 0
  diag$test_cost[1]   <- 0
  diag$test_type[1]   <- "hbsag"   # intervention starts at t=1, so t=0 uses SoC test profile

  # yearly loop (produces end-of-year t state at row t+1)
  for (t in 1:time_horizon) {
    # 1) Build TM using current state (semi-dynamic FOI)
    b <- build_tm(sv, t)         # must return list(tm, p_inf, lambda_raw, prev_eff, new_inf)
    tm_list[[t]] <- b$tm

    # 2) Evolve one Markov step
    next_sv <- as.numeric(sv) %*% b$tm
    next_sv <- setNames(as.vector(next_sv), names(sv))

    # 3) Apply testing/treatment/vaccination for this year & arm
    tv <- apply_testing_vax(next_sv, t, arm)   # must return list(vec, n_tested, tp, fp, tn, fn, newly_tx, newly_vax, test_cost, test_type, sens, spec)
    sv <- tv$vec

    # 4) Record end-of-year t state
    results[t + 1, ] <- sv

    # 5) Record diagnostics for year t (row = t+1)
    diag$p_inf[t + 1]       <- b$p_inf
    diag$lambda_raw[t + 1]  <- b$lambda_raw
    diag$prev_eff[t + 1]    <- b$prev_eff
    diag$incidence[t + 1]   <- b$new_inf

    diag$n_tested[t + 1]    <- tv$n_tested
    diag$tp[t + 1]          <- tv$tp
    diag$fp[t + 1]          <- tv$fp
    diag$tn[t + 1]          <- tv$tn
    diag$fn[t + 1]          <- tv$fn
    diag$newly_tx[t + 1]    <- tv$newly_tx
    diag$newly_vax[t + 1]   <- tv$newly_vax
    diag$test_cost[t + 1]   <- tv$test_cost
    diag$test_type[t + 1]   <- tv$test_type
    diag$sens[t + 1]        <- tv$sens
    diag$spec[t + 1]        <- tv$spec
  }

  # return
  list(results = results, diag = diag)  # (optionally add tm_list = tm_list)
}



# 4) CALIBRATION — tune hbv_beta (and optionally r_beta_treat)

# Metrics for SoC over a chosen window (e.g., years 1:3).
# You can calibrate using incidence per 100 PY with denominator either "alive" or "susceptible".
soc_metrics <- function(beta_override = NULL,
                        r_beta_treat_override = NULL,
                        years_window = 1:3,
                        denom = c("alive","susceptible")) {
  denom <- match.arg(denom)
  
  # save/restore globals
  old_beta <- hbv_beta
  old_rbt  <- r_beta_treat
  on.exit({ hbv_beta <<- old_beta; r_beta_treat <<- old_rbt }, add = TRUE)
  
  if (!is.null(beta_override))         hbv_beta      <<- beta_override
  if (!is.null(r_beta_treat_override)) r_beta_treat  <<- r_beta_treat_override
  
  res <- run_arm("SOC")
  
  # Start-of-year t composition (row "t-1" in results, because results[t+1] is end of year t)
  results_df <- as.data.frame(res$results)
  alive_states <- c("HBV_susceptible","Acute_symptomatic_HBV","HBV_immune",
                    "Chronic_HBV_nontx","Chronic_HBV_tx",
                    "Compensated_cirrhosis","Decompensated_cirrhosis","Liver_cancer")
  
  start_alive <- rowSums(results_df[as.character(0:time_horizon), alive_states], na.rm = TRUE)
  start_sus   <- results_df[as.character(0:time_horizon), "HBV_susceptible", drop = TRUE]
  start_chron <- rowSums(results_df[as.character(0:time_horizon), c("Chronic_HBV_nontx","Chronic_HBV_tx")], na.rm = TRUE)
  prev_start  <- ifelse(start_alive > 0, start_chron / start_alive, NA_real_)
  
  diag_df <- res$diag %>% dplyr::select(time, incidence)  # new infections during year t
  
  tset <- years_window[years_window >= 1 & years_window <= time_horizon]
  inc_counts <- diag_df$incidence[match(tset, diag_df$time)]
  den_series <- if (denom == "alive") start_alive else start_sus
  den_window <- den_series[tset]
  
  inc_per100 <- 100 * sum(inc_counts, na.rm = TRUE) / pmax(sum(den_window, na.rm = TRUE), 1e-12)
  prev_mean  <- mean(prev_start[tset], na.rm = TRUE)
  
  tibble::tibble(
    beta = hbv_beta,
    r_beta_treat = r_beta_treat,
    inc_per100py = inc_per100,
    prev_mean = prev_mean
  )
}

# Calibrate hbv_beta only (primary target: incidence per 100 PY; optional prevalence via weight w_prev)
calibrate_beta <- function(target_inc_per100py,
                           years_window = 1:3,
                           target_prev = NULL,
                           w_prev = 0,                 # 0 = ignore prevalence, else 0..1
                           denom = c("alive","susceptible"),
                           beta_bounds = c(1e-6, 10),  # widen if needed
                           tol = 1e-6,
                           maxit = 60,
                           apply = TRUE,
                           quiet = TRUE) {
  denom <- match.arg(denom)
  
  obj <- function(b) {
    m <- soc_metrics(beta_override = b, years_window = years_window, denom = denom)
    inc_err <- (m$inc_per100py - target_inc_per100py)
    if (!is.null(target_prev) && is.finite(target_prev) && w_prev > 0) {
      prev_err <- (m$prev_mean - target_prev)
      val <- (inc_err / max(target_inc_per100py, 1e-9))^2 +
        w_prev * (prev_err / max(target_prev, 1e-9))^2
    } else {
      val <- (inc_err / max(target_inc_per100py, 1e-9))^2
    }
    if (!quiet) message(sprintf("beta=%.6f -> inc=%.4f (tgt=%.4f), prev=%.4f",
                                b, m$inc_per100py, target_inc_per100py, m$prev_mean))
    val
  }
  
  opt <- optimize(obj, interval = beta_bounds, tol = tol, maximum = FALSE)
  beta_star <- opt$minimum
  fit <- soc_metrics(beta_override = beta_star, years_window = years_window, denom = denom)
  
  if (apply) hbv_beta <<- beta_star
  
  list(beta_star = beta_star,
       fit = fit,
       obj = opt$objective,
       years_window = years_window,
       denom = denom,
       applied = apply)
}

# Optional: grid search over r_beta_treat, calibrating beta for each candidate
calibrate_beta_and_rbt <- function(target_inc_per100py,
                                   years_window = 1:3,
                                   rbt_grid = c(0.05, 0.1, 0.2, 0.3),
                                   target_prev = NULL, w_prev = 0,
                                   denom = c("alive","susceptible"),
                                   beta_bounds = c(1e-6, 10),
                                   tol = 1e-6,
                                   quiet = TRUE,
                                   apply = TRUE) {
  denom <- match.arg(denom)
  
  grid <- purrr::map_dfr(rbt_grid, function(rbt) {
    old_rbt <- r_beta_treat
    r_beta_treat <<- rbt
    on.exit({ r_beta_treat <<- old_rbt }, add = TRUE)
    
    obj <- function(b) {
      m <- soc_metrics(beta_override = b, r_beta_treat_override = rbt,
                       years_window = years_window, denom = denom)
      inc_err <- (m$inc_per100py - target_inc_per100py)
      if (!is.null(target_prev) && is.finite(target_prev) && w_prev > 0) {
        prev_err <- (m$prev_mean - target_prev)
        (inc_err / max(target_inc_per100py, 1e-9))^2 +
          w_prev * (prev_err / max(target_prev, 1e-9))^2
      } else {
        (inc_err / max(target_inc_per100py, 1e-9))^2
      }
    }
    
    opt <- optimize(obj, interval = beta_bounds, tol = tol, maximum = FALSE)
    tibble::tibble(r_beta_treat = rbt, beta = opt$minimum, obj = opt$objective)
  })
  
  best <- grid %>% dplyr::arrange(obj) %>% dplyr::slice(1)
  
  if (apply) {
    r_beta_treat <<- best$r_beta_treat
    hbv_beta     <<- best$beta
  }
  
  list(grid = grid, best = best, years_window = years_window, denom = denom, applied = apply)
}

# 7) SCENARIOS & QUICK OUTPUTS

res_SOC <- run_arm("SOC")
res_INT <- run_arm("INT")


# 3.1 Discounting
disc_rate <- p("disc_rate", 0.03)
disc_vec  <- 1 / (1 + disc_rate) ^ (0:time_horizon)


# KP adult: HBV×(HIV,Syph,both) split WITH explicit per-state outputs


# 1) Long traces
to_long <- function(res, arm) {
  as.data.frame(res$results) %>%
    tibble::rownames_to_column("time") %>% 
    dplyr::mutate(time = as.integer(time), arm = arm) %>%
    tidyr::pivot_longer(-c(time, arm), names_to = "state", values_to = "count")
}
trace_all <- dplyr::bind_rows(to_long(res_SOC,"SOC"), to_long(res_INT,"INT"))

alive_states <- setdiff(adult_states, "Death")
alive_df <- trace_all %>%
  dplyr::filter(state %in% alive_states) %>%
  dplyr::group_by(arm, time) %>%
  dplyr::summarise(alive = sum(count, na.rm = TRUE), .groups = "drop")

# 2) External HIV/Syph series
hiv_syph_raw <- readr::read_csv("hiv_syph_raw_triplex.csv", show_col_types = FALSE)

if (!"time" %in% names(hiv_syph_raw) && "age" %in% names(hiv_syph_raw)) {
  hiv_syph_raw <- dplyr::rename(hiv_syph_raw, time = age)
}

prev_df <- hiv_syph_raw %>%
  dplyr::select(time, arm, hiv_infected_adults, syphilis_infected_adults, coinfected_adults) %>%
  dplyr::left_join(alive_df, by = c("arm","time")) %>%
  dplyr::mutate(
    p_HIV   = dplyr::if_else(alive > 0, hiv_infected_adults      / alive, 0),
    p_Syph  = dplyr::if_else(alive > 0, syphilis_infected_adults / alive, 0),
    p_pair0 = dplyr::if_else(alive > 0, coinfected_adults        / alive, 0),
    p_pair  = pmin(pmax(p_pair0, p_HIV + p_Syph - 1, 0), pmin(p_HIV, p_Syph)),
    prop_triple = p_pair,
    prop_HIV    = pmax(p_HIV  - p_pair, 0),
    prop_Syph   = pmax(p_Syph - p_pair, 0),
    sum_coinf   = prop_triple + prop_HIV + prop_Syph,
    scale_down  = dplyr::if_else(sum_coinf > 1, 1/sum_coinf, 1),
    prop_triple = prop_triple * scale_down,
    prop_HIV    = prop_HIV    * scale_down,
    prop_Syph   = prop_Syph   * scale_down,
    prop_none   = 1 - (prop_triple + prop_HIV + prop_Syph)
  ) %>%
  dplyr::select(arm, time, dplyr::starts_with("prop_"))

# 3) HBV states to split (INCLUDES CC, DC, HCC)
hbv_overlap_states <- c(
  "Acute_symptomatic_HBV",
  "Chronic_HBV_nontx","Chronic_HBV_tx",
  "Compensated_cirrhosis","Decompensated_cirrhosis","Liver_cancer"
)

hbv_base_py <- trace_all %>%
  dplyr::filter(state %in% hbv_overlap_states) %>%
  dplyr::group_by(arm, time, state) %>%
  dplyr::summarise(person_years = sum(count, na.rm = TRUE), .groups = "drop")
# ^ keep this untouched for epi & calibration

# 4) Split into explicit states (_HIV/_Syph/_Triple)
hbv_split_long <- hbv_base_py %>%
  dplyr::left_join(prev_df, by = c("arm","time")) %>%
  dplyr::mutate(
    yrs_HBVonly    = person_years * prop_none,
    yrs_HBV_HIV    = person_years * prop_HIV,
    yrs_HBV_Syph   = person_years * prop_Syph,
    yrs_HBV_Triple = person_years * prop_triple
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("yrs_"),
    names_to = "coinf",
    values_to = "py"
  ) %>%
  dplyr::mutate(
    state = dplyr::case_when(
      coinf == "yrs_HBVonly"    ~ state,
      coinf == "yrs_HBV_HIV"    ~ paste0(state, "_HIV"),
      coinf == "yrs_HBV_Syph"   ~ paste0(state, "_Syph"),
      coinf == "yrs_HBV_Triple" ~ paste0(state, "_Triple")
    )
  ) %>%
  dplyr::select(arm, time, state, person_years = py)

# 5) Disability weights (multiplicative combination)
dw_hbv_tbl <- tibble::tibble(
  base_state = hbv_overlap_states,
  dw_hbv = c(
    p("hbv_dw_acute_sym_adult", 0.20),
    p("hbv_dw_chronic",         0.10),
    p("hbv_dw_chronic",         0.10),
    p("hbv_dw_comp_cirrhosis",  0.20),
    p("hbv_dw_decomp_cirrhosis",0.50),
    p("hbv_dw_liver_cancer",    0.54)
  )
)

dw_hiv_adult  <- p("hiv_dw_adt",  0.15)
dw_syph_adult <- p("syph_dw_adt", 0.05)

dw_lookup <- dw_hbv_tbl %>%
  tidyr::expand_grid(coinf = c("HBVonly","HIV","Syph","Triple")) %>%
  dplyr::mutate(
    state = dplyr::case_when(
      coinf == "HBVonly" ~ base_state,
      coinf == "HIV"     ~ paste0(base_state, "_HIV"),
      coinf == "Syph"    ~ paste0(base_state, "_Syph"),
      coinf == "Triple"  ~ paste0(base_state, "_Triple")
    ),
    dw = dplyr::case_when(
      coinf == "HBVonly" ~ dw_hbv,
      coinf == "HIV"     ~ 1 - (1 - dw_hbv) * (1 - dw_hiv_adult),
      coinf == "Syph"    ~ 1 - (1 - dw_hbv) * (1 - dw_syph_adult),
      coinf == "Triple"  ~ 1 - (1 - dw_hbv) * (1 - dw_hiv_adult) * (1 - dw_syph_adult)
    )
  ) %>%
  dplyr::select(state, dw)

# Discount vector if missing
if (!exists("disc_vec")) {
  disc_rate <- p("disc_effects", 0.03)
  disc_vec  <- 1 / ((1 + disc_rate) ^ (0:time_horizon))
}

# --- Build layered YLD state from the split pools ---
hbv_layered_yld_state <- hbv_split_long %>%
  dplyr::left_join(dw_lookup, by = "state") %>%
  dplyr::mutate(
    person_years_adj = ifelse(grepl("^Acute_symptomatic_HBV", state),
                              person_years * p("acute_duration_years", 1/3),
                              person_years),
    yld      = person_years_adj * dw,
    yld_disc = yld * disc_vec[time + 1]
  )

# Adult Syph DW rule: from year >= 2, drop the Syph component
#  - _Syph   -> use base HBV DW (as if no syph)
#  - _Triple -> use the corresponding _HIV DW (drop syph piece)
hbv_layered_yld_state <- hbv_layered_yld_state %>%
  mutate(
    base_state = sub("_(HIV|Syph|Triple)$","", state),
    hiv_state  = paste0(base_state, "_HIV"),
    dw = case_when(
      grepl("_Syph$", state)   & time >= 2 ~ dw_lookup$dw[match(base_state, dw_lookup$state)],
      grepl("_Triple$", state) & time >= 2 ~ dw_lookup$dw[match(hiv_state,  dw_lookup$state)],
      TRUE ~ dw
    ),
    # recompute YLDs with adjusted DWs
    yld      = person_years_adj * dw,
    yld_disc = yld * disc_vec[time + 1]
  ) %>%
  select(-base_state, -hiv_state)


# 6) Costs (base HBV + adders; Chronic_HBV_tx gets tx cost)
cost_map <- tibble::tibble(
  base_state = hbv_overlap_states,
  cost_base = c(
    p("cost_acute",   0),
    p("cost_chronic", 0),
    p("cost_chronic", 0),   # base for Chronic_HBV_tx; tx added below
    p("cost_CC",      0),
    p("cost_DC",      0),
    p("cost_HCC",     0)
  ),
  is_tx = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
)

cost_tx_annual <- p("cost_tx_annual",   0)
cost_add_hiv   <- p("cost_add_hiv_adt", 0)
cost_add_syph  <- p("cost_add_syph_adt",0)

cost_lookup <- cost_map %>%
  tidyr::expand_grid(coinf = c("HBVonly","HIV","Syph","Triple")) %>%
  dplyr::mutate(
    state = dplyr::case_when(
      coinf == "HBVonly" ~ base_state,
      coinf == "HIV"     ~ paste0(base_state, "_HIV"),
      coinf == "Syph"    ~ paste0(base_state, "_Syph"),
      coinf == "Triple"  ~ paste0(base_state, "_Triple")
    ),
    cost_add = dplyr::case_when(
      coinf == "HBVonly" ~ 0,
      coinf == "HIV"     ~ cost_add_hiv,
      coinf == "Syph"    ~ cost_add_syph,
      coinf == "Triple"  ~ cost_add_hiv + cost_add_syph
    )
  ) %>%
  dplyr::select(state, cost_base, is_tx, cost_add)

# --- Build layered COST state from the split pools ---
hbv_layered_cost_state <- hbv_split_long %>%
  dplyr::left_join(cost_lookup, by = "state") %>%
  dplyr::mutate(
    c_base = person_years * cost_base,
    c_tx   = ifelse(grepl("^Chronic_HBV_tx", state),
                    person_years * p("cost_tx_annual", 0),
                    0),
    c_add  = person_years * cost_add,
    cost      = c_base + c_add + c_tx,
    cost_disc = cost * disc_vec[time + 1]
  )

# Adult Syph cost rule: from year >= 2, remove syph cost adders
#  - _Syph   -> 0 syph adder
#  - _Triple -> keep only HIV adder (subtract syph component)
hbv_layered_cost_state <- hbv_layered_cost_state %>%
  mutate(
    cost_add_eff = case_when(
      grepl("_Syph$", state)   & time >= 2 ~ 0,
      grepl("_Triple$", state) & time >= 2 ~ cost_add - cost_add_syph,
      TRUE ~ cost_add
    ),
    c_add    = person_years * cost_add_eff,
    cost     = c_base + c_add + c_tx,
    cost_disc= cost * disc_vec[time + 1]
  ) %>%
  select(-cost_add_eff)


# 7) If you need the compact per-year totals for CE:
yld_hbv_layered_disc_by_year <- hbv_layered_yld_state %>%
  dplyr::group_by(arm, time) %>%
  dplyr::summarise(yld_hbv_layered_disc = sum(yld_disc, na.rm = TRUE), .groups = "drop")

cost_hbv_layered_disc_by_year <- hbv_layered_cost_state %>%
  dplyr::group_by(arm, time) %>%
  dplyr::summarise(cost_hbv_layered_disc = sum(cost_disc, na.rm = TRUE), .groups = "drop")












# 3.2 Tidy traces (aggregate chronic for reporting)
tidy_trace <- function(res, arm) {
  as.data.frame(res$results) %>%
    rownames_to_column("time") %>%
    mutate(time = as.integer(time), arm = arm) %>%
    as_tibble() %>%
    mutate(
      Chronic_HBV = Chronic_HBV_nontx + Chronic_HBV_tx,
      alive = HBV_susceptible + Acute_symptomatic_HBV + HBV_immune +
        Chronic_HBV + Compensated_cirrhosis + Decompensated_cirrhosis + Liver_cancer
    )
}
tr_SOC <- tidy_trace(res_SOC, "SOC")
tr_INT <- tidy_trace(res_INT, "INT")
trace_all <- bind_rows(tr_SOC, tr_INT)



# 3.3 Incidence (we already track it in diag$new_inf). Join for convenience.
diag_all <- bind_rows(res_SOC$diag, res_INT$diag) %>%
  dplyr::select(
    arm, time,
    incidence, p_inf,
    n_tested, newly_tx, newly_vax,
    test_cost, test_type, sens, spec
  )

epi <- trace_all %>%
  select(arm, time, alive, HBV_susceptible, HBV_immune, Chronic_HBV,
         Compensated_cirrhosis, Decompensated_cirrhosis, Liver_cancer, Death) %>%
  left_join(diag_all, by = c("arm","time")) %>%
  mutate(chronic_prev = ifelse(alive > 0, Chronic_HBV / alive, NA_real_))

# 3.4 YLD scaffold (state DWs pulled from params; defaults 0)
dw_map <- tibble::tribble(
  ~state,                    ~sym,
  "Acute_symptomatic_HBV",   "hbv_dw_acute_sym_adult",
  "Chronic_HBV",             "dw_chronic", #0 
  "Compensated_cirrhosis",   "dw_cc", #0 
  "Decompensated_cirrhosis", "hbv_dw_decomp_cirrhosis",
  "Liver_cancer",            "hbv_dw_liver_cancer",
  "HBV_susceptible",         "dw_sus",    # usually 0
  "HBV_immune",              "dw_imm"     # usually 0
) %>%
  mutate(dw = purrr::map_dbl(sym, ~ p(.x, 0)))

# acute lasts part of a year (optional)
acute_duration <- p("acute_duration_years", 1/3)

yld_df <- trace_all %>%
  select(arm, time,
         HBV_susceptible, Acute_symptomatic_HBV, HBV_immune,
         Chronic_HBV, Compensated_cirrhosis, Decompensated_cirrhosis, Liver_cancer) %>%
  pivot_longer(-c(arm,time), names_to = "state", values_to = "py") %>%
  left_join(dw_map %>% select(state, dw), by = "state") %>%
  mutate(
    py_adj = ifelse(state == "Acute_symptomatic_HBV", py * acute_duration, py),
    yld = py_adj * coalesce(dw, 0),
    disc = disc_vec[time + 1],
    yld_disc = yld * disc
  )

yld_by_year <- yld_df %>%
  group_by(arm, time) %>%
  summarise(total_yld_disc = sum(yld_disc, na.rm = TRUE), .groups = "drop")

# 3.5 YLL scaffold (Both sexes life table)
life_exp_df <- read_csv("life_exp.csv", show_col_types = FALSE) %>%
  rename(age_group = `Age Group`, le_both = `Both sexes`, le_female = Female) %>%
  mutate(age_start = as.numeric(str_extract(age_group, "^\\d+")),
         age_start = if_else(str_detect(age_group, "<1"), 0, age_start))

get_agebin <- function(age) {
  if (age < 1) return(0)
  possible <- sort(unique(na.omit(life_exp_df$age_start)))
  max(possible[possible <= age])
}

death_increments <- trace_all %>%
  arrange(arm, time) %>%
  group_by(arm) %>%
  mutate(new_deaths = Death - lag(Death, default = 0)) %>%
  ungroup() %>%
  filter(new_deaths > 0) %>%
  mutate(
    age_at_death = kp_age0 + time,          # KP adult cohort age progression
    age_bin = vapply(age_at_death, get_agebin, numeric(1))
  ) %>%
  left_join(life_exp_df %>% select(age_start, le_both), by = c("age_bin" = "age_start")) %>%
  mutate(
    yll = new_deaths * le_both,
    disc = disc_vec[time + 1],
    yll_disc = yll * disc
  )

yll_by_year <- death_increments %>%
  group_by(arm, time) %>%
  summarise(total_yll_disc = sum(yll_disc, na.rm = TRUE), .groups = "drop")

# 3.6 DALYs by year and total
dalys_by_year <- full_join(yld_by_year, yll_by_year, by = c("arm","time")) %>%
  mutate(
    total_yld_disc = coalesce(total_yld_disc, 0),
    total_yll_disc = coalesce(total_yll_disc, 0),
    total_dalys_disc = total_yld_disc + total_yll_disc
  )

dalys_total <- dalys_by_year %>%
  group_by(arm) %>%
  summarise(total_dalys_disc = sum(total_dalys_disc), .groups = "drop")

# 3.7 Costs scaffold (program + states). All unit costs default 0 unless in KP sheet.
cost_test       <- p("cost_test", 2)
cost_vax_full   <- p("hbv_c_vaccine_full", 0)
cost_tx_annual  <- p("cost_tx_annual", 0)
cost_cc         <- p("cost_Compensated_cirrhosis", 0)
cost_dc         <- p("cost_Decompensated_cirrhosis", 0)
cost_hcc        <- p("cost_Liver_cancer", 0)

# Per-year volumes
vol_prog_by_year <- diag_all %>%
  mutate(disc = disc_vec[time + 1]) %>%
  mutate(
    cost_tests_disc = coalesce(test_cost, 0) * disc,                 # <- uses diag$test_cost
    cost_vax_disc   = coalesce(newly_vax, 0) * p("cost_vax_full", 0) * disc
  ) %>%
  group_by(arm, time) %>%
  summarise(
    cost_tests_disc = sum(cost_tests_disc, na.rm = TRUE),
    cost_vax_disc   = sum(cost_vax_disc,   na.rm = TRUE),
    .groups = "drop"
  )



state_costs_by_year <- trace_all %>%
  transmute(
    arm, time, disc = disc_vec[time + 1],
    cost_tx_disc  = Chronic_HBV_tx         * cost_tx_annual * disc,
    cost_cc_disc  = Compensated_cirrhosis  * cost_cc        * disc,
    cost_dc_disc  = Decompensated_cirrhosis* cost_dc        * disc,
    cost_hcc_disc = Liver_cancer           * cost_hcc       * disc
  ) %>%
  group_by(arm, time) %>%
  summarise(across(starts_with("cost_"), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

costs_by_year <- full_join(vol_prog_by_year, state_costs_by_year, by = c("arm","time")) %>%
  mutate(across(starts_with("cost_"), ~ coalesce(.x, 0)),
         total_costs_disc = cost_tests_disc + cost_vax_disc + cost_tx_disc + cost_cc_disc + cost_dc_disc + cost_hcc_disc)

costs_total <- costs_by_year %>%
  group_by(arm) %>%
  summarise(total_costs_disc = sum(total_costs_disc), .groups = "drop")

# 3.8 ICER snapshot (INT vs SOC)
icer <- full_join(costs_total, dalys_total, by = "arm") %>%
  pivot_wider(names_from = arm, values_from = c(total_costs_disc, total_dalys_disc)) %>%
  mutate(
    delta_cost    = total_costs_disc_INT - total_costs_disc_SOC,
    dalys_averted = total_dalys_disc_SOC - total_dalys_disc_INT,
    icer          = ifelse(dalys_averted > 0, delta_cost / dalys_averted, NA_real_)
  )

# 3.9 Quick prints (you can comment out later)
print(head(epi, 12))
print(head(dalys_by_year, 6))
print(icer)



## USE POST calibration 
# target_inc_per100py <- 2.0   # SoC incidence per 100 PY over years 1–3
# target_prev         <- 0.06  # SoC chronic prevalence (mean over years 1–3)
#win <- 1:3


#cal <- calibrate_beta(target_inc_per100py,
 #                     years_window = win,
  #                    target_prev = target_prev,
   #                   w_prev = 0.25,            # 0 to ignore prevalence
    #                  denom = "alive",          # or "susceptible"
     #                 beta_bounds = c(1e-6, 10),
      #                apply = TRUE, quiet = TRUE)
#cal$fit        # prints fitted inc/prev
#hbv_beta       # now updated in global scope

