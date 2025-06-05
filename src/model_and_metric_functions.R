# src/model_and_metric_functions.R

# --- Global Constants (relevant to model/metric calculation) ---
F_KON_CONVERSION <- 6.022e23 * 1e-15 # for conversion from M to molecule based k_on with the assumed voume of an E. coli cell 1e-15 L
PSEUDO_COUNT <- 1e-9
TIME_FOR_STEADY_STATE <- 50000 # Or pass as argument if it varies

# --- ODE Model Function ---
model_func_gsa <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dm <- k_elon * m_pre - beta_m * m - k_on * m * s + k_off * ms
    dms <- k_on * m * s + k_elon * m_pre_s - beta_ms * ms - k_off * ms
    ds <- alpha_s - beta_s * s - k_on * m * s - k_on * m_pre * s + k_off * ms + k_off * m_pre_s + k_off * m_trunc_s
    dm_pre <- k_init - k_elon * m_pre - k_on * s * m_pre + k_off * m_pre_s - beta_m_co * m_pre
    dm_pre_s <- k_on * s * m_pre - k_off * m_pre_s - k_elon * m_pre_s - beta_ms * m_pre_s - k_coreg * m_pre_s
    dm_trunc_s <- k_coreg * m_pre_s - k_off * m_trunc_s - beta_ms * m_trunc_s
    dp <- k_m * m + k_ms * ms - beta_p * p
    list(c(dm, dms, ds, dm_pre, dm_pre_s, dm_trunc_s, dp))
  })
}

# --- Output Metric Function ---
calculate_protein_regulation_strength <- function(current_params,
                                                  initial_states,
                                                  ode_model_func,
                                                  sim_times = c(0, TIME_FOR_STEADY_STATE), 
                                                  rtol = 1e-6, atol = 1e-10) {

  params_scenario <- current_params
 
  out_scenario_df <- tryCatch({
    out <- deSolve::ode(y = initial_states, times = sim_times, func = ode_model_func,
                        parms = params_scenario, method = "lsoda", rtol = rtol, atol = atol)
    as.data.frame(out)
  }, error = function(e) {
    # warning(paste("ODE solve error (scenario):", e$message));
    NULL
    })
  if (is.null(out_scenario_df) || nrow(out_scenario_df) == 0) return(NA)
  p_scenario <- out_scenario_df[nrow(out_scenario_df), "p"]

  params_no_sRNA <- params_scenario
  params_no_sRNA[["k_on"]] <- 0 # Simulate with no sRNA interaction
 
  out_no_sRNA_df <- tryCatch({
    out <- deSolve::ode(y = initial_states, times = sim_times, func = ode_model_func,
                        parms = params_no_sRNA, method = "lsoda", rtol = rtol, atol = atol)
    as.data.frame(out)
  }, error = function(e) {
    # warning(paste("ODE solve error (no sRNA):", e$message));
    NULL
    })
  if (is.null(out_no_sRNA_df) || nrow(out_no_sRNA_df) == 0) return(NA)
  p_no_sRNA_interaction <- out_no_sRNA_df[nrow(out_no_sRNA_df), "p"]

  if(is.na(p_scenario) || is.na(p_no_sRNA_interaction)) return(NA)

  p_scenario_adj <- p_scenario + PSEUDO_COUNT # Use global PSEUDO_COUNT
  p_no_sRNA_interaction_adj <- p_no_sRNA_interaction + PSEUDO_COUNT

  regulation_strength <- if (abs(p_no_sRNA_interaction_adj) < 1e-12) { # Avoid division by near-zero
                           ifelse(abs(p_scenario_adj) < 1e-12, 0, NA) # If both are zero, no regulation (log2(1)=0)
                         } else {
                           log2(p_scenario_adj / p_no_sRNA_interaction_adj)
                         }

  if(!is.finite(regulation_strength)) return(NA)
  return(regulation_strength)
}
