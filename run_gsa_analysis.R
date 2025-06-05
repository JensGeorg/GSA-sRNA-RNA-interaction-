# run_gsa_analysis.R

# Load Required Packages
library(sensitivity)
library(deSolve)
library(dplyr)
library(parallel)
library(ggplot2)
library(tidyr)

# Source helper functions
source("src/model_and_metric_functions.R")
source("src/gsa_sobol_functions.R")

source("model_and_metric_functions.R")
source("gsa_sobol_functions.R")


# --- 1. Define Baseline Parameters and Initial States ---
# F_KON_CONVERSION is defined in model_and_metric_functions.R
k_elon_default_rate <- 25/1000 # Example for a gene of length 1000nt and an elongation rate 0 25nt/s, adjust as needed


# baseline_params are a fallback if not all parameters are sampled for the analysis
baseline_params <- list(
  k_init    = 0.96,
  alpha_s   = 0.36,
  beta_s    = 1.5e-3,
  beta_m    = 0.002310491,
  beta_ms   = 0.002310491,
  beta_m_co = 0.002310491,
  beta_p    = 0.005,
  k_on      = (1.3e6) / F_KON_CONVERSION, # Molar to molecular conversion
  k_off     = 0.3,
  k_elon    = k_elon_default_rate,
  k_m       = 14.0,
  k_ms      = 14.0,
  k_coreg   = 1.0 
)

default_initial_states <- c(
  m = 0, ms = 0, s = baseline_params$alpha_s / baseline_params$beta_s,
  m_pre = 0, m_pre_s = 0, m_trunc_s = 0, p = 0
)

# --- 2. Define Parameters for GSA ---
# List of parameters to vary, their ranges, and distributions
gsa_parameters_definition <- list(
  k_init    = list(min=0.1,    max=2,    dist="uniform"),
  alpha_s   = list(min=0.1,    max=2,    dist="uniform"),
  k_on      = list(min=1.660578e-06, max=0.1660578 , dist="loguniform"),
  k_off     = list(min=1e-3,   max=100,  dist="loguniform"),
  beta_m    = list(min=0.0002, max=0.06, dist="loguniform"),
  beta_ms   = list(min=0.0002, max=0.06, dist="loguniform"),
  beta_s    = list(min=0.0002, max=0.06, dist="loguniform"),
  k_coreg   = list(min=0.01,   max=0.99, dist="uniform"),
  k_elon    = list(min=5/1000, max=60/1000, dist="uniform"),
  beta_m_co = list(min=0.0001, max=0.06, dist="loguniform"),
  beta_p    = list(min=0.0005, max=0.05, dist="loguniform"), 
  k_m       = list(min=0.1,    max=30,   dist="uniform"),
  k_ms      = list(min=0.1,    max=30,   dist="uniform")
)

gsa_param_names_to_vary <- names(gsa_parameters_definition)
k_gsa_varied <- length(gsa_param_names_to_vary)

min_ranges_for_gsa <- sapply(gsa_parameters_definition, function(p) p$min)
max_ranges_for_gsa <- sapply(gsa_parameters_definition, function(p) p$max)
param_distributions_for_gsa <- sapply(gsa_parameters_definition, function(p) p$dist)


# --- 3. Run Sobol GSA (in Parallel) ---
n_base_samples_sobol <- 1e5 # ADJUST AS NEEDED (e.g., 1e3 for quick, 1e4+ for robust)
cat("Preparing Sobol GSA with", n_base_samples_sobol, "base samples for", k_gsa_varied, "parameters.\n")

# Generate Sobol design matrices (scaled 0-1)
X1_sobol <- data.frame(matrix(runif(n_base_samples_sobol * k_gsa_varied), nrow = n_base_samples_sobol, byrow = TRUE))
X2_sobol <- data.frame(matrix(runif(n_base_samples_sobol * k_gsa_varied), nrow = n_base_samples_sobol, byrow = TRUE))

sobol_design_object <- sensitivity::sobolSalt(model = NULL, X1 = X1_sobol, X2 = X2_sobol, nboot = 100) 

# Setup for parallel execution
num_cores_to_use <- max(1, parallel::detectCores() - 1)
cat("Using", num_cores_to_use, "cores for parallel execution.\n")
gsa_cluster <- parallel::makeCluster(num_cores_to_use)

# Export necessary variables and functions to the cluster
parallel::clusterExport(gsa_cluster, c("single_gsa_run", "calculate_protein_regulation_strength", "model_func_gsa",
                                  "gsa_param_names_to_vary", "min_ranges_for_gsa", "max_ranges_for_gsa",
                                  "param_distributions_for_gsa", "baseline_params", "default_initial_states",
                                  "F_KON_CONVERSION", "PSEUDO_COUNT", "TIME_FOR_STEADY_STATE")) # Ensure all globals are exported
parallel::clusterEvalQ(gsa_cluster, { library(deSolve) }) # Load deSolve on each worker

cat("Starting parallel model evaluations for GSA...\n")
# sobol_design_object$X contains all parameter sets to run (N * (k+2) rows)
gsa_model_outputs <- NA
tryCatch({
  gsa_model_outputs <- parallel::parApply(gsa_cluster, sobol_design_object$X, 1, single_gsa_run,
                                    gsa_param_names_arg = gsa_param_names_to_vary,
                                    min_ranges_arg = min_ranges_for_gsa,
                                    max_ranges_arg = max_ranges_for_gsa,
                                    param_distributions_arg = param_distributions_for_gsa,
                                    baseline_params_arg = baseline_params,
                                    initial_states_arg = default_initial_states,
                                    ode_model_func_arg = model_func_gsa,
                                    calculate_metric_func_arg = calculate_protein_regulation_strength)
}, finally = {
  parallel::stopCluster(gsa_cluster) # Always stop the cluster
  cat("Parallel execution finished or stopped.\n")
})


if(all(is.na(gsa_model_outputs))) {
    stop("All GSA model outputs are NA. Cannot proceed. Check errors during single_gsa_run execution or parameter ranges.")
} else if (anyNA(gsa_model_outputs)) {
    num_na_outputs <- sum(is.na(gsa_model_outputs))
    warning(paste(num_na_outputs, "NA values produced in GSA output. Results might be unreliable. Consider checking parameter ranges or ODE stability. Imputing with median for Sobol calculation."))
    median_output_for_imputation <- median(gsa_model_outputs, na.rm=TRUE)
    if (is.na(median_output_for_imputation) && num_na_outputs < length(gsa_model_outputs)) { # if median is NA but not all are NA
        warning("Median for imputation is NA, attempting to use mean.")
        mean_output_for_imputation <- mean(gsa_model_outputs, na.rm=TRUE)
         if(is.na(mean_output_for_imputation)){
             stop("Cannot impute NAs as both median and mean are NA. Too many NAs or problematic output values.")
         }
        gsa_model_outputs[is.na(gsa_model_outputs)] <- mean_output_for_imputation
    } else if (is.na(median_output_for_imputation) && num_na_outputs == length(gsa_model_outputs)) {
         stop("All GSA outputs are NA even after trying to get median. Cannot proceed.")
    } else {
         gsa_model_outputs[is.na(gsa_model_outputs)] <- median_output_for_imputation
    }
    if(all(is.na(gsa_model_outputs))) stop("All GSA outputs are still NA after imputation. Cannot proceed.") 
}


# --- 4. Calculate and Print Sobol Indices ---
cat("Calculating Sobol indices...\n")
sensitivity::tell(sobol_design_object, gsa_model_outputs) # Assigns y to the sobol object
cat("\nSobol GSA Results Object Summary:\n")
print(sobol_design_object) # This will print S1, ST etc.

# --- 5. Plot GSA Results ---
cat("\nPlotting Sobol indices...\n")
sobol_plot <- plot_sobol_indices(sobol_design_object, gsa_param_names_to_vary)
if (!is.null(sobol_plot)) {
    print(sobol_plot)
    ggsave(file.path("results", "sobol_indices.png"), plot = sobol_plot, width = 10, height = 7)
    cat("Sobol indices plot saved to results/sobol_indices.png\n")
}

# --- 6. Unscale Parameters and Combine with Outputs for Further Analysis ---
unscaled_parameters_df_gsa <- unscale_sobol_matrix(
    sobol_design_object$X,
    gsa_param_names_to_vary,
    min_ranges_for_gsa,
    max_ranges_for_gsa,
    param_distributions_for_gsa
)

# Combine unscaled parameters with GSA output values
if (length(gsa_model_outputs) == nrow(unscaled_parameters_df_gsa)) {
  gsa_results_with_inputs_df <- cbind(unscaled_parameters_df_gsa, regulation_strength = gsa_model_outputs)
  print(head(gsa_results_with_inputs_df))

  # Save the combined data and the sobol object for further analysis
  gsa_output_data <- list(
      gsa_results_with_inputs_df = gsa_results_with_inputs_df, # Data for regression, plotting etc.
      sobol_indices_summary = sobol_design_object, # Contains S, T, etc.
      gsa_parameter_definitions = gsa_parameters_definition # For reference
  )
  save(gsa_output_data, file = file.path("results", "gsa_analysis_outputs.RData"))

} else {
  cat(paste("Error: 'gsa_model_outputs' length (", length(gsa_model_outputs),
            ") does not match the number of GSA samples in design matrix (",
            nrow(unscaled_parameters_df_gsa), ").\n"))
}
