# src/gsa_sobol_functions.R

# --- Wrapper function for a single GSA run (for parallel execution) ---
single_gsa_run <- function(scaled_sample_row, # A single row from Sobol design matrix (values 0-1)
                           gsa_param_names_arg,
                           min_ranges_arg, max_ranges_arg,
                           param_distributions_arg, # Named vector: param_name -> "uniform" or "loguniform"
                           baseline_params_arg,
                           initial_states_arg,
                           ode_model_func_arg,
                           calculate_metric_func_arg) {

  unscaled_params <- numeric(length(gsa_param_names_arg))
  names(unscaled_params) <- gsa_param_names_arg

  for (j in 1:length(gsa_param_names_arg)) {
    p_name <- gsa_param_names_arg[j]
    p_scaled_val <- scaled_sample_row[j]
    p_min <- min_ranges_arg[[p_name]] 
    p_max <- max_ranges_arg[[p_name]] 
    dist_type <- param_distributions_arg[[p_name]] 

    if (is.null(dist_type)) stop(paste("Distribution type not specified for parameter:", p_name))

    if (dist_type == "uniform") {
      unscaled_params[p_name] <- p_min + p_scaled_val * (p_max - p_min)
    } else if (dist_type == "loguniform") {
      if (p_min <= 0 || p_max <= 0) stop(paste("Log-uniform ranges for", p_name, "must be positive."))
      log_p_min <- log10(p_min)
      log_p_max <- log10(p_max)
      unscaled_params[p_name] <- 10^(log_p_min + p_scaled_val * (log_p_max - log_p_min))
    } else {
      stop(paste("Unknown distribution type:", dist_type, "for parameter:", p_name))
    }
  }

  params_for_run <- baseline_params_arg
  params_for_run[names(unscaled_params)] <- unscaled_params

  if (any(is.na(params_for_run[gsa_param_names_arg]))) {
     # warning(paste("NA found in parameters for run after unscaling for params:", paste(gsa_param_names_arg, collapse=", ")))
     return(NA)
  }

  return(calculate_metric_func_arg(
    current_params = params_for_run,
    initial_states = initial_states_arg,
    ode_model_func = ode_model_func_arg
  ))
}


# --- Function to Unscale Sobol Matrix ---
unscale_sobol_matrix <- function(sobol_X_matrix, gsa_param_names_arg, min_ranges_arg, max_ranges_arg, param_distributions_arg) {
    unscaled_matrix <- matrix(NA, nrow = nrow(sobol_X_matrix), ncol = ncol(sobol_X_matrix))
    colnames(unscaled_matrix) <- gsa_param_names_arg

    for (j in 1:ncol(sobol_X_matrix)) {
        param_name <- gsa_param_names_arg[j]
        scaled_column_values <- sobol_X_matrix[, j]
        p_min <- min_ranges_arg[param_name]
        p_max <- max_ranges_arg[param_name]
        dist_type <- param_distributions_arg[param_name]

        if (is.na(dist_type)) stop(paste("Distribution type not specified for parameter:", param_name))

        if (dist_type == "uniform") {
            unscaled_matrix[, j] <- p_min + scaled_column_values * (p_max - p_min)
        } else if (dist_type == "loguniform") {
            if (p_min <= 0 || p_max <= 0) stop(paste("Log-uniform ranges for '", param_name, "' must be positive."))
            log_p_min <- log10(p_min)
            log_p_max <- log10(p_max)
            unscaled_matrix[, j] <- 10^(log_p_min + scaled_column_values * (log_p_max - log_p_min))
        } else {
            stop(paste("Unknown distribution type:", dist_type, "for parameter:", param_name))
        }
    }
    return(as.data.frame(unscaled_matrix))
}


# --- Function to Plot Sobol Results ---
plot_sobol_indices <- function(sobol_results_object, varied_param_names) {
    if (!is.null(sobol_results_object$S) && !is.null(sobol_results_object$T)) {
        s1_original <- sobol_results_object$S$original
        s1_min_ci <- sobol_results_object$S$`min. c.i.` 
        s1_max_ci <- sobol_results_object$S$`max. c.i.`
        st_original <- sobol_results_object$T$original
        st_min_ci <- sobol_results_object$T$`min. c.i.`
        st_max_ci <- sobol_results_object$T$`max. c.i.`

        results_df_plot <- data.frame(
            Parameter = factor(varied_param_names, levels = varied_param_names), 
            S1 = s1_original,
            S1.conf.low = s1_min_ci,
            S1.conf.high = s1_max_ci,
            ST = st_original,
            ST.conf.low = st_min_ci,
            ST.conf.high = st_max_ci
        )
        results_long_df_plot <- tidyr::pivot_longer(results_df_plot,
                                                    cols = c(S1, ST),
                                                    names_to = "Index_Type",
                                                    values_to = "Value")
        results_long_df_plot$CI_low <- ifelse(results_long_df_plot$Index_Type == "S1",
                                              results_long_df_plot$S1.conf.low,
                                              results_long_df_plot$ST.conf.low)
        results_long_df_plot$CI_high <- ifelse(results_long_df_plot$Index_Type == "S1",
                                               results_long_df_plot$S1.conf.high,
                                               results_long_df_plot$ST.conf.high)

        p <- ggplot(results_long_df_plot, aes(x = Parameter, y = Value, fill = Index_Type)) +
            geom_bar(stat="identity", position=position_dodge(width=0.9)) +
            geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                          position = position_dodge(width=0.9), width=0.25, linewidth=0.3) +
            labs(title="Sobol Sensitivity Indices",
                 y="Sensitivity Index Value", x = "Parameter") +
            theme_minimal(base_size=12) +
            theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
                  plot.title = element_text(hjust=0.5, face="bold")) +
            scale_fill_brewer(palette="Set1", name="Index Type",
                              labels=c("S1 (First Order)", "ST (Total Order)"))
        return(p)
    } else {
        cat("Could not generate Sobol plot components (S or T missing from object).\n")
        return(NULL)
    }
}
