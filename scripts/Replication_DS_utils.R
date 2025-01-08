

###*************************************###
###*************************************###
##### ***REPLICATION DS: UTILITIES*** #####
###*************************************###
###*************************************###

###****************************###
##### ***DATA PROCESSING*** ######
###****************************###

#' @description Utility functions for the data processing pipeline, including:
#' 1. calculate_perpendicular_line()
#' 2. calculate_corridor_corners()
#' 3. filter_participants_in_corridor()
#' 4. check_corridor_membership()
#' 5. plot_corridor_check()
#' 6. plot_corridor_samples()
#' 7. calculate_distances()
#' 8. validate_distance_calculations()
#' 9. add_distance_validation


###*********************###
##### DEFINE CORRIDOR #####
###*********************###

###**********************************###
##### Calculate perpendicular line #####
###**********************************###

#' Calculate perpendicular line through confederate point
#'
#' @description Creates a perpendicular line through the confederate's position that
#' intersects with the regression line of participant movement
#'
#' @param lm_obj Linear model object from participant trajectory
#' @param conf_coords Tibble containing confederate x,z coordinates, derived from the respective median
#' @return List containing:
#'   - line_function: Function representing perpendicular line
#'   - intersection: Tibble with intersection point coordinates
#'   - conf_coords: Original confederate coordinates

calculate_perpendicular_line <- function(lm_obj, conf_coords) {
  # Get slope of regression line
  slope <- coef(lm_obj)[2]

  # Get intercept of regression line
  intercept <- coef(lm_obj)[1]

  # Calculate perpendicular slope (negative reciprocal)
  perp_slope <- -1 / slope

  # Calculate perpendicular line intercept through confederate point
  # y = mx + b -> b = y - mx
  perp_intercept <- conf_coords$conf_z - (perp_slope * conf_coords$conf_x)

  # Find intersection point
  x_intersect <- (perp_intercept - intercept) / (slope - perp_slope)
  z_intersect <- slope * x_intersect + intercept

  # Create and return function representing the perpendicular line
  list(
    line_function = function(x) perp_slope * x + perp_intercept,
    intersection = tibble(
      x_intersect = x_intersect,
      z_intersect = z_intersect
    ),
    conf_coords = conf_coords
  )
}

###*********************************###
##### Calculate corridor corners #####
###*********************************###

#' @description Defines a rectangular corridor centered on the perpendicular line
#' with specified width
#'
#' @param part_data Participant trajectory data
#' @param conf_data Confederate position data
#' @param lm_obj Linear model of participant trajectory
#' @param perp Perpendicular line object from calculate_perpendicular_line()
#' @return Tibble containing x,z coordinates of corridor corners

calculate_corridor_corners <- function(part_data, conf_data, lm_obj, perp) {
  # Get x range from confederate data
  x_min <- min(conf_data$centered_x) - 10
  x_max <- max(conf_data$centered_x) + 10

  # Width of corridor (1 unit on each side)
  corridor_width <- 4

  # Get perpendicular line function
  perp_line <- perp$line_function

  # Calculate corners of corridor
  corners <- tibble(
    x = c(x_min, x_min, x_max, x_max),
    z = c(
      perp_line(x_min) - corridor_width/2,  # bottom left
      perp_line(x_min) + corridor_width/2,  # top left
      perp_line(x_max) + corridor_width/2,  # top right
      perp_line(x_max) - corridor_width/2   # bottom right
    ),
    corner = c("bl", "tl", "tr", "br")
  )

  return(corners)
}

###*************************************###
##### Filter participants in corridor #####
###*************************************###

#' Filter participant positions to those within corridor
#'
#' @description Identifies which participant positions fall within the defined
#' analysis corridor
#'
#' @param part_data Participant position data
#' @param rect Tibble defining corridor corners from calculate_corridor_corners()
#' @return Filtered participant data containing only in-corridor positions

filter_participants_in_corridor <- function(part_data, rect) {
  # Make sure corridor points are in the correct order
  rect_ordered <- rect %>%
    arrange(match(corner, c("bl", "tl", "tr", "br")))

  # For each participant point, check if it falls within the polygon
  part_data %>%
    mutate(
      in_corridor = point.in.polygon(
        point.x = centered_x,
        point.y = centered_z,
        pol.x = rect_ordered$x,
        pol.y = rect_ordered$z
      ) > 0
    ) %>%
    filter(in_corridor) %>%
    select(-in_corridor)
}


###*******************************###
##### Check corridor membership #####
###*******************************###

##### Create column indicating corridor membership #####

#' Check which participant positions fall within corridor
#'
#' @description Takes participant position data and corridor definition, returns the
#' full dataset with a new column indicating whether each position falls within the corridor
#'
#' @param full_df Full dataset containing confederate and participant positions
#' @param part_df Participant position data
#' @param rect Tibble defining corridor corners (bl, tl, tr, br)
#' @return Original dataframe with new 'in_corridor' column ("yes"/"no"). Confederates are set to "no".

check_corridor_membership <- function(full_df, part_df, rect) {
  # Make sure corridor points are in the correct order
  rect_ordered <- rect %>%
    arrange(match(corner, c("bl", "tl", "tr", "br")))

  # Find which participant points are in the corridor
  part_in_corridor <- part_df %>%
    mutate(
      in_corridor = point.in.polygon(
        point.x = centered_x,
        point.y = centered_z,
        pol.x = rect_ordered$x,
        pol.y = rect_ordered$z
      ) > 0
    ) %>%
    filter(in_corridor) %>%
    select(frame, pid)  # Keep only identifying columns

  # Join with full dataset to mark points in corridor
  full_df %>%
    left_join(
      part_in_corridor %>% mutate(in_corridor = "yes"),
      by = c("frame", "pid")
    ) %>%
    mutate(in_corridor = coalesce(in_corridor, "no"))
}


##### Visual check of corridor #####

##### Inspect specific videos

#' Create visual check plot for corridor analysis
#'
#' @description Creates a plot showing participant trajectories, corridor boundaries,
#' regression line, and perpendicular line for a specific video file
#'
#' @param df_name String name of the video file to plot (without path)
#' @return ggplot object showing:
#'   - Corridor polygon (shaded area)
#'   - Participant points (colored by corridor membership)
#'   - Confederate points
#'   - Regression line (red)
#'   - Perpendicular line (blue)

plot_corridor_check <- function(df_name) {
  # Find the index matching the df_name
  index <- which(nested_dfs_check_corridor$df_name == df_name)

  if (length(index) == 0) {
    stop("df_name not found in nested_dfs_check_corridor")
  }

  # Get the data for this df_name
  plot_data <- nested_dfs_check_corridor$data[[index]]

  # Get confederate data for this index
  conf_data <- nested_dfs_corridor$confederates[[index]]

  # Calculate data bounds
  x_min <- min(plot_data$centered_x)
  x_max <- max(plot_data$centered_x)
  z_min <- min(plot_data$centered_z)
  z_max <- max(plot_data$centered_z)

  plot_data %>%
    mutate(in_corridor = ifelse(in_corridor == "yes", "In corridor", "Out of corridor")) %>%
    ggplot(aes(x = centered_x, y = centered_z)) +
    # Add corridor polygon
    geom_polygon(data = nested_dfs_corridor$corridor[[index]],
                 aes(x = x, y = z),
                 fill = "#2c444b",
                 alpha = 0.2) +
    # Add participant points colored by corridor status
    geom_point(aes(color = in_corridor), size = 2) +
    # Add confederate points
    geom_point(data = conf_data,
               aes(x = centered_x, y = centered_z, color = "Confederate"), size = 2) +
    # Add regression line
    geom_line(data = tibble(
      x = seq(x_min, x_max, length.out = 100),
      y = coef(nested_dfs_corridor$part_lm[[index]])[1] +
          coef(nested_dfs_corridor$part_lm[[index]])[2] * x
    ), aes(x = x, y = y), color = "#e84c4c", linewidth = 1) +
    # Add perpendicular line
    geom_line(data = tibble(
      x = seq(x_min, x_max, length.out = 100),
      y = map_dbl(x, nested_dfs_corridor$perpendicular[[index]]$line_function)
    ), aes(x = x, y = y), color = "#517eb9", linewidth = 1) +
    # Set axis limits to data range
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(z_min, z_max)
    ) +
      theme(
        axis.text = element_blank(),
        text = element_text(size = 25)
      ) +
      labs(
        x = "",
        y = "",
        color = "Data points",
        title = glue("Replication DS (check corridor): {df_name}")
      )

}


##### Sample four videos within each location

#' Create sample plots of corridor analysis
#'
#' @description Creates a multi-panel plot showing corridor analysis results for a
#' random sample of videos from each location. Uses plot_corridor_check() to create
#' the plots.
#'
#' @param nested_data Nested dataframe containing corridor analysis results
#' @param n_samples Number of videos to sample from each location (default = 4)
#' @return Saves a PDF file containing the multi-panel plot to the output directory
#' @details Creates a 4x4 grid of plots, each showing:
#'   - Corridor polygon
#'   - Participant points (colored by corridor membership)
#'   - Regression and perpendicular lines
#'   Plots are saved with consistent axis limits for comparison

plot_corridor_samples <- function(nested_data, n_samples = 4) {
  # Get sample indices for each location
  sampled_indices <- nested_data %>%
    mutate(row_number = row_number()) %>%
    group_by(location) %>%
    slice_sample(n = n_samples) %>%
    pull(row_number)

  # First pass to determine global axis limits
  all_bounds <- map_dfr(sampled_indices, function(index) {
    plot_data <- nested_dfs_check_corridor$data[[index]]
    tibble(
      x_min = min(plot_data$centered_x),
      x_max = max(plot_data$centered_x),
      z_min = min(plot_data$centered_z),
      z_max = max(plot_data$centered_z)
    )
  })

  global_x_min <- min(all_bounds$x_min)
  global_x_max <- max(all_bounds$x_max)
  global_z_min <- min(all_bounds$z_min)
  global_z_max <- max(all_bounds$z_max)

  # Create list of plots
  plots <- map(sampled_indices, function(index) {
    # Get the data for this index
    plot_data <- nested_dfs_check_corridor$data[[index]]
    df_name <- nested_dfs_check_corridor$df_name[index]

    plot_data %>%
      ggplot(aes(x = centered_x, y = centered_z)) +
      # Add corridor polygon
      geom_polygon(data = nested_dfs_corridor$corridor[[index]],
                   aes(x = x, y = z),
                   fill = "#2c444b",
                   alpha = 0.2) +
      # Add points colored by corridor status
      geom_point(aes(color = in_corridor)) +
      # Add regression line
      geom_line(data = tibble(
        x = seq(global_x_min, global_x_max, length.out = 100),
        y = coef(nested_dfs_corridor$part_lm[[index]])[1] +
            coef(nested_dfs_corridor$part_lm[[index]])[2] * x
      ), aes(x = x, y = y, linetype = "Regression Line"), color = "#e84c4c", linewidth = 1) +
      # Add perpendicular line
      geom_line(data = tibble(
        x = seq(global_x_min, global_x_max, length.out = 100),
        y = map_dbl(x, nested_dfs_corridor$perpendicular[[index]]$line_function)
      ), aes(x = x, y = y, linetype = "Perpendicular Line"), color = "#517eb9", linewidth = 1) +
      # Set consistent axis limits across all plots
      coord_cartesian(
        xlim = c(global_x_min, global_x_max),
        ylim = c(global_z_min, global_z_max)
      ) +
      theme(
        text = element_text(size = 25),
        axis.text = element_text(size = 12),
        legend.position = "none"  # Remove individual legends
      ) +
      labs(
        x = "",
        y = "",
        subtitle = df_name
      ) +
      scale_linetype_manual(
        name = "Lines",
        values = c("Regression Line" = "solid", "Perpendicular Line" = "solid")
      ) +
      scale_color_discrete(name = "In corridor")
  })

  # Combine plots using patchwork with a shared legend at the bottom
  combined_plot <- wrap_plots(plots, ncol = 4) +
    plot_annotation(
      title = "Check corridor",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(2, "lines")
      )
    ) &
    theme(legend.box = "horizontal")  # Arrange legends horizontally

  # Save plot
  ggsave("Replication_DS_corridor check sample.pdf",
         plot = combined_plot,
         path = path_output_figures,
         device = cairo_pdf,
         width = 18,
         height = 15)
}


###**********************************************###
##### CALCULATE DISTANCES FROM REFERENCE POINT #####
###**********************************************###

#' Calculate minimum distance between each participant-observation and the closer confederates
#' for the corresponding frame
#'
#' @description For each participant position within the corridor, calculates the
#' minimum distance to any confederate in the same frame
#'
#' @param data Dataframe containing participant positions with 'in_corridor' column
#' @param confederates Dataframe containing confederate positions
#' @return Dataframe containing frame, pid, and min_distance columns for in-corridor positions
#' @details
#'   - Only processes rows where in_corridor == "yes"
#'   - For each participant position, calculates distances to all confederates in that frame
#'   - Returns the minimum of these distances

calculate_distances <- function(data, confederates) {
  # Calculate distances for in_corridor rows
  data %>%
    filter(in_corridor == "yes") %>%
    # Group by frame and pid
    group_by(frame, pid) %>%
    mutate(
      min_distance = map2_dbl(centered_x, centered_z, function(x, z) {
        # Get confederate positions for this frame
        conf_positions <- confederates %>%
          filter(frame == frame[1])  # Using currentframe

        # Calculate distances to both confederates
        distances <- sqrt(
          (x - conf_positions$centered_x)^2 +
          (z - conf_positions$centered_z)^2
        )

        # Return the smaller distance
        min(distances)
      })
    ) %>%
    ungroup() %>%
    select(frame, pid, min_distance)  # Keep only necessary columns for join

}


###***************************###
##### CHECK REFERENCE POINT #####
###***************************###

##### Validate distance calculations #####

#' @description Performs detailed validation of distance calculations by computing
#' distances to each confederate separately and identifying which confederate is closer
#'
#' @param data Dataframe containing participant positions with 'in_corridor' column
#' @param confederates Dataframe containing confederate positions
#' @param n_samples Optional number of random samples to validate. If NULL, validates all rows
#' @return Dataframe containing detailed distance information:
#'   - frame: Frame number
#'   - pid: Participant ID
#'   - centered_x/z: Participant coordinates
#'   - closer_conf: ID of closer confederate (1 or 2)
#'   - dist_conf1/2: Distance to each confederate
#'   - min_dist: Minimum distance

validate_distance_calculations <- function(data, confederates, n_samples = NULL) {
  # Use all rows if n_samples is NULL, otherwise sample
  validation_data <- data %>%
    filter(in_corridor == "yes") %>%
    {if (!is.null(n_samples)) slice_sample(., n = n_samples) else .}

  # Calculate detailed distance information
  validation_results <- validation_data %>%
    group_by(frame, pid) %>%
    mutate(
      distance_details = map2(centered_x, centered_z, function(x, z) {
        # Get confederate positions for this frame
        conf_positions <- confederates %>%
          filter(frame == frame[1]) %>%
          mutate(conf_id = row_number())

        # Calculate distances to both confederates
        distances <- sqrt(
          (x - conf_positions$centered_x)^2 +
          (z - conf_positions$centered_z)^2
        )

        # Return detailed information
        list(
          distances = distances,
          closer_confederate = which.min(distances),
          min_distance = min(distances),
          confederate_positions = conf_positions
        )
      })
    ) %>%
    ungroup()

  # Unnest the results for easier viewing
  detailed_results <- validation_results %>%
    mutate(
      closer_conf = map_dbl(distance_details, ~.$closer_confederate),
      dist_conf1 = map_dbl(distance_details, ~.$distances[1]),
      dist_conf2 = map_dbl(distance_details, ~.$distances[2]),
      min_dist = map_dbl(distance_details, ~.$min_distance)
    ) %>%
    select(frame, pid, centered_x, centered_z,
           closer_conf, dist_conf1, dist_conf2, min_dist)

  return(detailed_results)
}


##### Add info on whether all distances are correct #####

#' @description Validates that the minimum distance matches the smaller of the two
#' individual confederate distances
#'
#' @param distances_check_df Dataframe containing distance validation results
#' @return Logical value indicating whether all distance calculations are correct
#' @details Uses a small tolerance (1e-10) to account for floating point arithmetic

add_distance_validation <- function(distances_check_df) {
  # Use exact equality check instead of all.equal
  all(abs(
    distances_check_df$min_dist -
    pmin(distances_check_df$dist_conf1, distances_check_df$dist_conf2)
  ) < 1e-10)  # Very small tolerance for floating point arithmetic
}



###****************************###
##### ***DATA ANALYSIS*** ########
###****************************###

#' @description Utility functions for the data analysis pipeline, including:
#' 1. create_dist_plot()
#' 2. identify_outliers()
#' 3. perform_ate_analysis_cluster_robust()
#' 4. perform_ate_analysis_wild_block_bootstrap()
#' 5. extract_regression_results()
#' 6. plot_estimate_distributions()
#' 7. create_specification_curve()


###*******************************###
##### DISTRIBUTION OF DISTANCES #####
###*******************************###

###******************************###
##### Create distribution plot #####
###******************************###

#' @description Creates a histogram showing the distribution of minimum distances
#' between participants and confederates
#'
#' @param data Dataframe containing minimum distance calculations
#' @return ggplot object showing histogram of distances with:

create_dist_plot <- function(data) {

  data %>%
    filter(!is.na(min_distance)) %>%   # Remove NA values
    ggplot(aes(x = min_distance)) +
      geom_histogram(color = "black") +
      theme(text = element_text(size = 20)) +
      labs(
        x = "Distance",
        y = "Count"
      )

}


###***********************###
##### IDENTIFY OUTLIERS #####
###***********************###

#' @description Identifies outliers in the minimum distances between participants and confederates
#' using the 1.5 * IQR rule for regular outliers and 3 * IQR rule for strong outliers
#'
#' @param data Dataframe containing minimum distance calculations
#' @return Filtered dataframe containing only outlier observations with columns:

identify_outliers <- function(data) {
  data %>%
    filter(!is.na(min_distance)) %>%
    mutate(
      Q1 = quantile(min_distance, 0.25, na.rm = TRUE),
      Q3 = quantile(min_distance, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower_bound = Q1 - 1.5 * IQR,
      upper_bound = Q3 + 1.5 * IQR,
      outlier_category = case_when(
        min_distance < (Q1 - 3 * IQR) ~ "Strong Outlier (Lower)",
        min_distance < lower_bound ~ "Outlier (Lower)",
        min_distance > upper_bound ~ "Outlier (Upper)",
        min_distance > (Q3 + 3 * IQR) ~ "Strong Outlier (Upper)",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(outlier_category)) %>%
    select(
      df_name,
      frame,
      pid,
      location,
      location_type,
      exp_condition,
      block_position,
      block_id,
      time_of_day,
      min_distance
    )
}


###************************###
##### PERFORM REGRESSION #####
###************************###

###************************************###
##### Cluster robust standard errors #####
###************************************###

#' @description Conducts regression analysis of treatment effects using cluster-robust
#' standard errors at the sub-block level. Runs the regression for both ATE (i.e., without covariates)
#' and CTE (i.e., with covariates).
#'
#' @param data Dataframe containing minimum distances and experimental conditions
#' @param include_covariates Logical indicating whether to include block_id as covariate (default = FALSE)
#' @return List containing:
#'   - model: Full lm object
#'   - coefficients: Treatment effects with cluster-robust statistics
#'   - summary: Model summary
#'   - cluster_vcov: Cluster-robust variance-covariance matrix
#'   - cluster_se: Cluster-robust standard errors
#' @details
#'   - Standardizes outcome variable (min_distance)
#'   - Uses "white" as reference category for exp_condition
#'   - Calculates cluster-robust SEs using HC1 correction
#'   - Returns 90% and 95% confidence intervals

perform_ate_analysis_cluster_robust <- function(data, include_covariates = FALSE) {
  # Prepare analysis data and standardize outcome
  analysis_data <- data %>%
    filter(!is.na(min_distance)) %>%
    mutate(
      outcome = scale(min_distance)[,1],
      exp_condition = relevel(factor(exp_condition), ref = "white"),
      block_id = factor(block_id),
      sub_block_id = factor(sub_block_id)  # Ensure sub_block_id is a factor
    )

  # Keep only complete cases for the model
  model_data <- analysis_data %>%
        filter(., !is.na(block_id), !is.na(sub_block_id), !is.na(outcome))

  # Create formula
  formula <- if (include_covariates) {
    as.formula("outcome ~ exp_condition + block_id")
  } else {
    as.formula("outcome ~ exp_condition")
  }

  # Fit model
  model <- lm(formula, data = model_data)

  # Get model summary
  model_summary <- summary(model)

  # Get cluster-robust standard errors
  cluster_vcov <- vcovCL(model, cluster = model_data$sub_block_id, type = "HC1")
  cluster_se <- sqrt(diag(cluster_vcov))

  # Extract treatment effects and CIs with cluster-robust SEs
  treatment_effects <- tidy(model) %>%
    filter(str_detect(term, "^exp_condition")) %>%
    mutate(
      std.error = cluster_se[str_detect(names(cluster_se), "^exp_condition")],
      statistic = estimate / std.error,
      p.value = 2 * pt(abs(statistic), df = df.residual(model), lower.tail = FALSE),
      lower_95 = estimate - qt(0.975, df.residual(model)) * std.error,
      upper_95 = estimate + qt(0.975, df.residual(model)) * std.error,
      lower_90 = estimate - qt(0.95, df.residual(model)) * std.error,
      upper_90 = estimate + qt(0.95, df.residual(model)) * std.error,
      # Clean up term names
      term = str_replace_all(term, "exp_condition", "") %>%
        str_replace_all("_", " ") %>%
        str_trim()
    )

  list(
    model = model,
    coefficients = treatment_effects,
    summary = model_summary,
    cluster_vcov = cluster_vcov,
    cluster_se = cluster_se
  )
}


###**************************###
##### Wild block bootstrap #####
###**************************###

#' @description Conducts regression analysis of treatment effects using wild block
#' bootstrap for standard error estimation
#'
#' @param data Dataframe containing minimum distances and experimental conditions
#' @param include_covariates Logical indicating whether to include block_id as covariate (default = FALSE)
#' @param R Number of bootstrap replications (default = 1000)
#' @return List containing:
#'   - model: Full lm object
#'   - coefficients: Treatment effects with bootstrapped statistics
#'   - boot_vcov: Bootstrap variance-covariance matrix
#'   - summary: Model summary
#' @details
#'   - Uses wild bootstrap for cluster-robust inference
#'   - Standardizes outcome variable (min_distance)
#'   - Uses "white" as reference category
#'   - Computes standard errors and confidence intervals using t-distribution

perform_ate_analysis_wild_block_bootstrap <- function(data, include_covariates = FALSE, R = 1000) {
  # Prepare analysis data and standardize outcome
  analysis_data <- data %>%
    filter(!is.na(min_distance)) %>%
    mutate(
      outcome = scale(min_distance)[,1],
      exp_condition = relevel(factor(exp_condition), ref = "white"),
      block_id = factor(block_id),
      sub_block_id = factor(sub_block_id)
    )

  # Initial model to get complete cases
  formula <- if (include_covariates) {
    as.formula("outcome ~ exp_condition + block_id")
  } else {
    as.formula("outcome ~ exp_condition")
  }

  initial_model <- lm(formula, data = analysis_data)

  # Keep only complete cases using the same approach as the template
  model_data <- analysis_data %>%
    slice(which(!is.na(residuals(initial_model))))

  # Refit model with complete cases
  model <- lm(formula, data = model_data)

  # Get wild bootstrap standard errors (matching template)
  boot_se <- cluster.boot(
    model,
    boot_type = "wild",
    cluster = model_data %>% pull(sub_block_id),
    parallel = TRUE,
    R = R
  )

  # Get indices for exp_condition coefficients
  exp_indices <- which(str_detect(names(coef(model)), "^exp_condition"))

  # Extract coefficients and compute confidence intervals
  treatment_effects <- tibble(
    term = names(coef(model))[exp_indices],
    coef = coef(model)[exp_indices],
    se = sqrt(diag(boot_se))[exp_indices]
  ) %>%
    mutate(
      t_stat = coef / se,
      p_value = 2 * pt(abs(t_stat), df = model$df.residual, lower.tail = FALSE),
      # Calculate standard CIs using t-distribution
      lower_95 = coef - qt(0.975, df = model$df.residual) * se,
      upper_95 = coef + qt(0.975, df = model$df.residual) * se,
      lower_90 = coef - qt(0.95, df = model$df.residual) * se,
      upper_90 = coef + qt(0.95, df = model$df.residual) * se,
      # Clean up term names
      term = str_replace_all(term, "exp_condition", "") %>%
        str_replace_all("_", " ") %>%
        str_trim()
    )

  list(
    model = model,
    coefficients = treatment_effects,
    boot_vcov = boot_se,
    summary = coef(summary(model))
  )
}


###********************************###
##### EXTRACT REGRESSION RESULTS #####
###********************************###

#' @description Extracts and formats regression results from different model specifications
#' (cluster-robust SE and wild bootstrap, with and without covariates)
#'
#' @param data Dataframe containing regression results from different model specifications
#' @return Long-format dataframe containing:
#'   - sample_no: Sample identifier
#'   - term: Treatment condition
#'   - model: Model specification type
#'   - estimate: Coefficient estimate
#'   - std.error: Standard error
#'   - p.value: P-value
#'   - lower/upper_90/95: Confidence interval bounds

extract_regression_results <- function(data) {
  data %>%
    mutate(
      ate_clustered_se_no_cov_stats = map(ate_clustered_se_no_cov, ~ .x$coefficients),
      ate_clustered_se_with_cov_stats = map(ate_clustered_se_with_cov, ~ .x$coefficients),
      ate_wild_no_cov_stats = map(ate_wild_no_cov, ~ .x$coefficients),
      ate_wild_with_cov_stats = map(ate_wild_with_cov, ~ .x$coefficients)
    ) %>%
    select(sample_no, ate_clustered_se_no_cov_stats:ate_wild_with_cov_stats) %>%
    mutate(
      ate_clustered_se_no_cov_stats = map(ate_clustered_se_no_cov_stats, ~ .x %>%
        select(term, estimate, std.error, p.value:upper_90)),  # Added std.error
      ate_clustered_se_with_cov_stats = map(ate_clustered_se_with_cov_stats, ~ .x %>%
        select(term, estimate, std.error, p.value:upper_90)),  # Added std.error
      ate_wild_no_cov_stats = map(ate_wild_no_cov_stats, ~ .x %>%
        rename(estimate = coef, p.value = p_value, std.error = se) %>%  # Renamed se to std.error
        select(term, estimate, std.error, p.value:upper_90)),  # Added std.error
      ate_wild_with_cov_stats = map(ate_wild_with_cov_stats, ~ .x %>%
        rename(estimate = coef, p.value = p_value, std.error = se) %>%  # Renamed se to std.error
        select(term, estimate, std.error, p.value:upper_90))   # Added std.error
    ) %>%
    unnest(cols = ends_with("_stats"), names_sep = "_") %>%
    select(-ate_clustered_se_with_cov_stats_term, -ate_wild_no_cov_stats_term, -ate_wild_with_cov_stats_term) %>%
    rename(term = ate_clustered_se_no_cov_stats_term) %>%
    pivot_longer(
      cols = !c(sample_no, term),
      names_to = "model",
      values_to = c("value")
    ) %>%
    separate(model, into = c("model", "stat"), sep = "_(?=[^_]+$)") %>%
    mutate(stat = if_else(str_detect(model, "_lower"), glue("lower_{stat}"), stat),
           stat = if_else(str_detect(model, "_upper"), glue("upper_{stat}"), stat),
           model = str_remove(model, "_lower|_upper")) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(model = case_when(
      model == "ate_clustered_se_no_cov_stats" ~ "Clustered SE (ATE)",
      model == "ate_clustered_se_with_cov_stats" ~ "Clustered SE (CTE)",
      model == "ate_wild_no_cov_stats" ~ "Wild Bootstrap (ATE)",
      model == "ate_wild_with_cov_stats" ~ "Wild Bootstrap (CTE)"
    ))
}


###********************************###
##### CREATE SPECIFICATION CURVE #####
###********************************###

#' @description Creates a specification curve plot showing effect size estimates
#' and confidence intervals across sub-samples anddifferent model specifications
#'
#' @param data Dataframe containing regression results
#' @return ggplot object showing:
#'   - Effect size estimates (points)
#'   - 90% confidence intervals (thick lines)
#'   - 95% confidence intervals (thin lines)
#'   - Horizontal reference line at zero
#'   - Faceted by treatment and model type
#' @details
#'   - Points and intervals are colored by statistical significance (p ≤ 0.05)
#'   - Estimates are ordered by p-value within each facet
#'   - Saves plot as PDF in output directory

create_specification_curve <- function(data) {
  # Automatically get and format the data name
  subtitle <- deparse(substitute(data)) %>%
  str_replace_all("_", " ") %>%
  str_to_title() %>%
  str_remove("Extracted Results ")

  data %>%
    group_by(term, model) %>%
    mutate(
      rank = rank(p.value),
      plot_order = dense_rank(p.value),
      sig_color = ifelse(p.value <= 0.05, "black", "grey70")
    ) %>%
  ungroup() %>%
  ggplot() +
  # Upper panel: Effect size estimates and CIs
  geom_hline(yintercept = 0, color = "grey80", linetype = 2, linewidth = 1.5) +
  # 90% CI (thicker)
  geom_linerange(aes(x = plot_order, y = estimate,
                     ymin = lower_90, ymax = upper_90,
                     color = sig_color),
                 linewidth = 1.5) +
  # 95% CI (thinner)
  geom_linerange(aes(x = plot_order, y = estimate,
                     ymin = lower_95, ymax = upper_95,
                     color = sig_color),
                 linewidth = 0.75) +
  geom_point(aes(x = plot_order, y = estimate, color = sig_color),
             size = 5,
             shape = 21,
             fill = "white") +
  scale_color_identity() +  # Use the actual colors specified in sig_color
  facet_grid(term ~ model, scales = "free") +
  theme(
    text = element_text(size = 20),
    panel.spacing = unit(2, "lines"),
    axis.text.x = element_blank(),
    plot.caption = element_text(hjust = 0, vjust = 1, lineheight = 1.2),
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10, unit = "pt")
  ) +
  labs(
      x = "Subsamples (ordered by p-value)",
      y = "Treatment Effect on Standardized Distance",
      title = "Coefficient Estimates and CIs for 50 Random Sub-Samples",
      subtitle = subtitle
    )

}
