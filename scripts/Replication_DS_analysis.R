

###************************************###
###************************************###
##### ***REPLICATION DS: ANALYSIS*** #####
###************************************###
###************************************###

###******************************###
##### ***SET UP WORK SPACE *** #####
###******************************###

###*******************###
##### LOAD PACKAGES #####
###*******************###

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,
  glue,
  ggthemr,
  tidyverse,
  patchwork,
  conflicted,
  broom,
  sandwich,
  lmtest,
  conflicted,
  flextable,
  multiwayvcov
)

###***************###
##### CONFLICTS #####
###***************###

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("lag", "dplyr")
conflict_prefer("first", "dplyr")
conflict_prefer("slice", "dplyr")


###*****************************###
##### SOURCE HELPER FUNCTIONS #####
###*****************************###

source(here("scripts", "Replication_DS_utils.R"))


###*********************###
##### PATHS & OUTPUTS #####
###*********************###

path_input <- here("data", "processed")
path_output_figures <- here("outputs", "figures")
path_output_tables <- here("outputs", "tables")


###***************###
##### LOAD DATA #####
###***************###

subsamples_df <- read_rds(file.path(path_input, "Replication_DS_subsamples.rds"))


###************###
##### THEMES #####
###************###

###******************###
##### ggplot theme #####
###******************###

fresh_palette <- c(
  "#72B5A4", # Green
  "#C1D7A3", # Light green
  "#F5A25D", # Orange
  "#FF4D4F", # Red
  "#FF79C6", # Pink
  "#6A5ACD", # Violet
  "#BF9DDA", # Lavender
  "#85B8D6", # Light blue
  "#58A3CB", # Blue
  "#000000", # Black
  "#FFDA79", # Yellow
  "#FFA07A"  # Coral
)

# Set ggthemr theme
ggthemr("fresh")

# Get the current theme settings from ggthemr
current_theme <- ggplot2::theme_get()

# Modify the theme to include custom text elements
modified_theme <- current_theme +
  theme(text = element_text(size = 15, family = "Montserrat"))

# Set the modified theme as the default
ggplot2::theme_set(modified_theme)

# Set the global discrete color and fill scales
scale_colour_discrete <- function(...) scale_colour_manual(values = fresh_palette)
scale_fill_discrete <- function(...) scale_fill_manual(values = fresh_palette)


###*********************###
##### flextable theme #####
###*********************###

set_flextable_defaults(
  font.family = "Montserrat",
  font.size = 10,
  border.color = "gray",
  big.mark = "",
  padding = 6
)


###********************###
##### ***ANALYSIS*** #####
###********************###

###*******************************###
##### DISTRIBUTION OF DISTANCES #####
###*******************************###

###***********************************###
##### Create all distribution plots #####
###***********************************###

#' @description Generates histograms of minimum distances for each subsample and
#' combines them into a multi-panel plot
#'
#' @details
#'   - Uses create_dist_plot() from Replication_DS_utils.R
#'   - Creates histogram for each subsample
#'   - Combines plots using patchwork::wrap_plots()
#'   - Saves as PDF

dist_plots <- map(subsamples_df$data, create_dist_plot)

wrap_plots(dist_plots, ncol = 5)

  # Save with appropriate filename
ggsave("Replication_DS_distance_distribution.pdf",
       path = path_output_figures,
       device = cairo_pdf,
       width = 16,
       height = 30)


###******************###
##### RUN ANALYSES #####
###******************###

###****************###
##### Fit models #####
###****************###

#' @description Estimates treatment effects using multiple model specifications
#' for different geographic subsets of the data
#'
#' @details
#'   - Creates separate results for:
#'     * All locations combined
#'     * Individual locations:
#'       - Neue Bahnhofstrasse
#'       - Frankfurter Allee
#'       - Gürtelstraße
#'     * Location types:
#'       - Residential areas
#'       - Commercial areas
#'   - For each subset, fits four model specifications:
#'     1. Cluster-robust SE without covariates (ATE)
#'     2. Cluster-robust SE with covariates (CTE)
#'     3. Wild bootstrap without covariates (ATE)
#'     4. Wild bootstrap with covariates (CTE)
#'   - Uses functions from Replication_DS_utils.R:
#'     * perform_ate_analysis_cluster_robust()
#'     * perform_ate_analysis_wild_block_bootstrap()
#'   - Results stored in nested dataframes for each location/subset


###*******************###
##### All locations #####
###*******************###

results_all_locations <- subsamples_df %>%
  mutate(ate_clustered_se_no_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = FALSE)),
         ate_clustered_se_with_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = TRUE)),
         ate_wild_no_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = FALSE)),
         ate_wild_with_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = TRUE)))


###*************************###
##### Neue Bahnhofstrasse #####
###*************************###

results_neue_bahnhofstrasse <- subsamples_df %>%
  mutate(data = map(data, ~ filter(., location == "Neue Bahnhofstraße"))) %>%
  mutate(ate_clustered_se_no_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = FALSE)),
         ate_clustered_se_with_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = TRUE)),
         ate_wild_no_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = FALSE)),
         ate_wild_with_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = TRUE)))


###***********************###
##### Frankfurter Allee #####
###***********************###

results_frankfurter_allee <- subsamples_df %>%
  mutate(data = map(data, ~ filter(., location == "Frankfurter Allee"))) %>%
  mutate(ate_clustered_se_no_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = FALSE)),
         ate_clustered_se_with_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = TRUE)),
         ate_wild_no_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = FALSE)),
         ate_wild_with_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = TRUE)))


###******************###
##### Gürtelstraße #####
###******************###

results_guertelstrasse <- subsamples_df %>%
  mutate(data = map(data, ~ filter(., location == "Gürtelstraße"))) %>%
  mutate(ate_clustered_se_no_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = FALSE)),
         ate_clustered_se_with_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = TRUE)),
         ate_wild_no_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = FALSE)),
         ate_wild_with_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = TRUE)))


###**********************###
##### Residential area #####
###**********************###

results_residential <- subsamples_df %>%
  mutate(data = map(data, ~ filter(., location_type == "Residential"))) %>%
  mutate(ate_clustered_se_no_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = FALSE)),
         ate_clustered_se_with_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = TRUE)),
         ate_wild_no_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = FALSE)),
         ate_wild_with_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = TRUE)))


###*********************###
##### Commercial area #####
###*********************###

results_commercial <- subsamples_df %>%
  mutate(data = map(data, ~ filter(., location_type == "Commercial"))) %>%
  mutate(
    ate_clustered_se_no_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = FALSE)),
    ate_clustered_se_with_cov = map(data, ~ perform_ate_analysis_cluster_robust(., include_covariates = TRUE)),
    ate_wild_no_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = FALSE)),
    ate_wild_with_cov = map(data, ~ perform_ate_analysis_wild_block_bootstrap(., include_covariates = TRUE))
  )

###********************************###
##### Extract degrees of freedom #####
###********************************###

#' @description Extracts and summarizes degrees of freedom from different model
#' specifications and location types. Table is not included in main report, but information is used to compile table 2
#'
#' @details
#'   - Creates summary for residential and commercial areas
#'   - For each area type:
#'     * Extracts df from ATE and CTE models
#'     * Calculates min, mean, and max df
#'   - Combines results into single dataframe

df_summary <- bind_rows(
  # Residential area
  tibble(
    area = "Residential",
    model = "ATE",
    df = map_dbl(results_residential$ate_wild_no_cov, ~.$model$df.residual)
  ),
  tibble(
    area = "Residential",
    model = "CTE",
    df = map_dbl(results_residential$ate_wild_with_cov, ~.$model$df.residual)
  ),
  # Commercial area
  tibble(
    area = "Commercial",
    model = "ATE",
    df = map_dbl(results_commercial$ate_wild_no_cov, ~.$model$df.residual)
  ),
  tibble(
    area = "Commercial",
    model = "CTE",
    df = map_dbl(results_commercial$ate_wild_with_cov, ~.$model$df.residual)
  )
) %>%
  group_by(area, model) %>%
  summarise(
    min_df = min(df),
    mean_df = mean(df),
    max_df = max(df),
    .groups = "drop"
  )


###*********************###
##### Extract results #####
###*********************###

#' @description Extracts and formats regression results from all model specifications
#' for each location subset
#'
#' @details
#'   - Uses extract_regression_results() from Replication_DS_utils.R
#'   - Creates formatted results for:
#'     * All locations combined
#'     * Individual locations:
#'       - Neue Bahnhofstrasse
#'       - Frankfurter Allee
#'       - Gürtelstraße
#'     * Location types:
#'       - Residential areas
#'       - Commercial areas
#'   - Each results object contains:
#'     * Coefficient estimates
#'     * Standard errors
#'     * P-values
#'     * 90% and 95% confidence intervals

extracted_results_all_locations <- extract_regression_results(results_all_locations)

extracted_results_neue_bahnhofstrasse <- extract_regression_results(results_neue_bahnhofstrasse)

extracted_results_frankfurter_allee <- extract_regression_results(results_frankfurter_allee)

extracted_results_guertelstrasse <- extract_regression_results(results_guertelstrasse)

extracted_results_residential <- extract_regression_results(results_residential)

extracted_results_commercial <- extract_regression_results(results_commercial)


###******************###
##### Plot results #####
###******************###

#' @description Generates specification curve plots for each location subset.
#' Produces figures 2, 3A, 3B, A.1A, A.1B, and A.1C
#'
#' @details
#'   - Uses create_specification_curve() from Replication_DS_utils.R
#'   - Creates separate plots for:
#'     * All locations
#'     * Each individual location
#'     * Each location type
#'   - Each plot shows:
#'     * Effect estimates ordered by p-value
#'     * 90% and 95% confidence intervals
#'     * Points/intervals colored by significance
#'   - Saves each plot as separate PDF

##### Figure 2: All locations #####

create_specification_curve(extracted_results_all_locations)

ggsave(
  filename = glue("Replication_DS_figure 2.pdf"),
  path = path_output_figures,
  device = cairo_pdf,
  width = 25,
  height = 12
)


##### Figure 3A: Residential area #####

create_specification_curve(extracted_results_residential)

ggsave(
  filename = glue("Replication_DS_figure 3A.pdf"),
  path = path_output_figures,
  device = cairo_pdf,
  width = 25,
  height = 12
)


##### Figure 3B: Commercial area #####

create_specification_curve(extracted_results_commercial)

ggsave(
  filename = glue("Replication_DS_figure 3B.pdf"),
  path = path_output_figures,
  device = cairo_pdf,
  width = 25,
  height = 12
)


##### Figure A.1A: Frankfurter Allee #####

create_specification_curve(extracted_results_frankfurter_allee)

ggsave(
  filename = glue("Replication_DS_figure A.1A.pdf"),
  path = path_output_figures,
  device = cairo_pdf,
  width = 25,
  height = 12
)


##### Figure A.1B: Gürtelstraße #####

create_specification_curve(extracted_results_guertelstrasse)

ggsave(
  filename = glue("Replication_DS_figure A.1B.pdf"),
  path = path_output_figures,
  device = cairo_pdf,
  width = 25,
  height = 12
)


##### Figure A.1C: Neue Bahnhofstrasse #####

create_specification_curve(extracted_results_neue_bahnhofstrasse)

ggsave(
  filename = glue("Replication_DS_figure A.1C.pdf"),
  path = path_output_figures,
  device = cairo_pdf,
  width = 25,
  height = 12
)


###*****************************************###
##### Create summary table for all models #####
###*****************************************###

#' @description Generates summary tables of model results across all locations
#' and model specifications
#'
#' @details
#'   - Creates list of data frames to analyze:
#'     * All Locations
#'     * Individual locations (Neue Bahnhofstraße, Frankfurter Allee, Gürtelstraße)
#'     * Location types (Residential, Commercial)
#'   - For each location/type:
#'     * Filters for "muslim garb" treatment
#'     * Groups by model specification
#'     * Calculates summary statistics:
#'       - Min/mean/max estimates
#'       - Min/mean/max standard errors
#'       - Share of significant results (p ≤ 0.05)
#'   - Creates two formatted tables:
#'     1. Area types (All, Residential, Commercial)
#'     2. Specific locations (All, individual streets)
#'   - Saves as Word documents

# Create list of data frames to analyze
summary_data_variants <- list(
  "All Locations" = extracted_results_all_locations,
  "Neue Bahnhofstraße" = extracted_results_neue_bahnhofstrasse,
  "Frankfurter Allee" = extracted_results_frankfurter_allee,
  "Gürtelstraße" = extracted_results_guertelstrasse,
  "Residential" = extracted_results_residential,
  "Commercial" = extracted_results_commercial
)

# Create nested dataframe with results
results_summary <- tibble(
  location = names(summary_data_variants),
  data = summary_data_variants
) %>%
  mutate(results = map(data, ~ .x %>%
    filter(term == "muslim garb") %>%
    group_by(model) %>%
    summarize(
      min_estimate = min(estimate),
      mean_estimate = mean(estimate),
      max_estimate = max(estimate),
      max_se = max(std.error),
      min_se = min(std.error),
      mean_se = mean(std.error),
      share_significant_05 = sum(p.value <= 0.05) / n()
    ))) %>%
  select(location, results)

##### Table 3: All locations and area types #####

results_summary %>%
  filter(location %in% c("All Locations", "Residential", "Commercial")) %>%
  unnest(results) %>%
  select(
    location, model,
    min_estimate, mean_estimate, max_estimate,
    min_se, mean_se, max_se,
    share_significant_05
  ) %>%
  flextable() %>%
  set_caption("Summary Statistics by Area Type") %>%
  set_header_labels(
    location = "Location",
    model = "Model",
    min_estimate = "Min",
    mean_estimate = "Mean",
    max_estimate = "Max",
    min_se = "Min",
    mean_se = "Mean",
    max_se = "Max",
    share_significant_05 = "Share Significant (p ≤ 0.05)"
  ) %>%
  add_header_row(
    values = c("", "Estimate", "Standard Error", ""),
    colwidths = c(2, 3, 3, 1)
  ) %>%
  colformat_double(
    j = c("min_estimate", "max_estimate", "mean_estimate",
          "min_se", "max_se", "mean_se"),
    digits = 3
  ) %>%
  colformat_double(
    j = "share_significant_05",
    digits = 2
  ) %>%
  bold(part = "header") %>%
  autofit() %>%
  save_as_docx(path = file.path(path_output_tables, "Replication_DS_table 3.docx"))

##### Table A1: All locations and specific locations #####

results_summary %>%
  filter(location %in% c("All Locations", "Frankfurter Allee", "Gürtelstraße", "Neue Bahnhofstraße")) %>%
  unnest(results) %>%
  select(
    location, model,
    min_estimate, mean_estimate, max_estimate,
    min_se, mean_se, max_se,
    share_significant_05
  ) %>%
  flextable() %>%
  set_caption("Summary Statistics by Location") %>%
  set_header_labels(
    location = "Location",
    model = "Model",
    min_estimate = "Min",
    mean_estimate = "Mean",
    max_estimate = "Max",
    min_se = "Min",
    mean_se = "Mean",
    max_se = "Max",
    share_significant_05 = "Share Significant (p ≤ 0.05)"
  ) %>%
  add_header_row(
    values = c("", "Estimate", "Standard Error", ""),
    colwidths = c(2, 3, 3, 1)
  ) %>%
  colformat_double(
    j = c("min_estimate", "max_estimate", "mean_estimate",
          "min_se", "max_se", "mean_se"),
    digits = 3
  ) %>%
  colformat_double(
    j = "share_significant_05",
    digits = 2
  ) %>%
  bold(part = "header") %>%
  autofit() %>%
  save_as_docx(path = file.path(path_output_tables, "Replication_DS_table A1.docx"))
