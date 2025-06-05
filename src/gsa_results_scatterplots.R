# gsa_results_scatterplots.R



# Ensure required packages are loaded
library(ggplot2)
library(dplyr)
library(viridis) # For scale_fill_viridis_c

# --- Define Generic Plotting Function ---
generate_gsa_visualization <- function(
    data_df,
    x_var_col,               # Character string: name of the column for X-axis
    y_var_col,               # Character string: name of the column for Y-axis
    plot_title = "",         # Character string: Title of the plot
    x_axis_label = x_var_col,# Character string or expression() for X-axis label
    y_axis_label = y_var_col,# Character string or expression() for Y-axis label
    use_log_x_scale = FALSE, # Boolean: TRUE to apply scale_x_log10()
    filter_x_positive = FALSE, # Boolean: TRUE to filter x_var_col > epsilon before plotting (esp. for log scale)
    vline_xintercept = NULL, # Numeric: x-value for a vertical line, or NULL
    hline_yintercepts = NULL,# Numeric vector: y-value(s) for horizontal line(s), or NULL
    num_bins_2d = 50,        # Integer: number of bins for geom_bin2d
    gam_color = "red",
    gam_linewidth = 1.5,
    point_alpha = 0.5,       # Alpha for points if geom_point is used instead of geom_bin2d
    base_font_size = 11,
    epsilon_val = 1e-9       # Epsilon for filtering positive x if filter_x_positive is TRUE
) {

  # --- Input Checks ---
  if (!is.data.frame(data_df)) {
    stop("Input 'data_df' must be a dataframe.")
  }
  if (!x_var_col %in% names(data_df)) {
    stop(paste("x_var_col '", x_var_col, "' not found in data_df."))
  }
  if (!y_var_col %in% names(data_df)) {
    stop(paste("y_var_col '", y_var_col, "' not found in data_df."))
  }
  if (!is.numeric(data_df[[x_var_col]])) {
    stop(paste("Column '", x_var_col, "' must be numeric for x-axis."))
  }
  if (!is.numeric(data_df[[y_var_col]])) {
    stop(paste("Column '", y_var_col, "' must be numeric for y-axis."))
  }

  # --- Data Preparation ---
  plot_data <- data_df %>%
    select(x_val = all_of(x_var_col), y_val = all_of(y_var_col)) %>%
    filter(is.finite(x_val) & is.finite(y_val))

  if (filter_x_positive) {
    plot_data <- plot_data %>% filter(x_val > epsilon_val)
  }
  
  if (use_log_x_scale && any(plot_data$x_val <= 0, na.rm = TRUE)) {
      original_rows <- nrow(plot_data)
      plot_data <- plot_data %>% filter(x_val > epsilon_val) # Ensure positive values for log scale
      cat(paste("Warning: For log x-scale, filtered out", original_rows - nrow(plot_data),
                "rows where '", x_var_col, "' was not positive. Remaining rows:", nrow(plot_data),"\n"))
  }


  if (nrow(plot_data) < 10) {
    cat(paste("Not enough data points (found", nrow(plot_data), ") for '", x_var_col, "' vs '", y_var_col, "' plot after filtering. Minimum 10 required.\n"))
    return(invisible(NULL)) # Return NULL if not enough data
  }

  # --- Create Plot ---
  p <- ggplot(plot_data, aes(x = x_val, y = y_val)) +
    geom_bin2d(bins = num_bins_2d, na.rm = TRUE) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), # cs for cubic regression splines
                color = gam_color, se = FALSE, na.rm = TRUE, linewidth = gam_linewidth) +
    scale_fill_viridis_c(name = "Count", option = "plasma", na.value = "transparent") +
    labs(title = plot_title, x = x_axis_label, y = y_axis_label) +
    theme_bw(base_size = base_font_size) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  if (use_log_x_scale) {
    p <- p + scale_x_log10()
  }

  if (!is.null(hline_yintercepts)) {
    for (y_int in hline_yintercepts) {
      p <- p + geom_hline(yintercept = y_int, linetype = "solid", color = "black", linewidth = 0.4)
    }
  }

  if (!is.null(vline_xintercept)) {
    p <- p + geom_vline(xintercept = vline_xintercept, linetype = "dashed", color = "darkgreen", linewidth = 1)
  }

  return(p)
}

