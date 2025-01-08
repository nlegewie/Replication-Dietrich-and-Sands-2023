

###*******************************************###
###*******************************************###
##### ***REPLICATION DS: DATA PROCESSING*** #####
###*******************************************###
###*******************************************###

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
  tidyverse,
  readxl,
  sp,
  conflicted,
  patchwork,
  ggthemr,
  flextable
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

path_input <- here("data", "raw")
path_output_data <- here("data", "processed")
path_output_figures <- here("outputs", "figures")
path_output_tables <- here("outputs", "tables")


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


###******************###
##### READ IN DATA #####
###******************###

###**********************###
##### Create file list #####
###**********************###

file_list <- list.files(path_input, pattern = "^202406.*\\.csv$", full.names = TRUE)


###***********###
##### Apply #####
###***********###

nested_dfs <- tibble(df_name = file_list,
                     data = map(file_list, read_csv))


###**********************###
##### Clean up df_name #####
###**********************###

nested_dfs_names_cleaned <- nested_dfs %>%
  mutate(df_name = basename(df_name)) %>%
  mutate(df_name = str_remove(df_name, "_centered_marked"))


###*************************###
##### ***ADD META DATA*** #####
###*************************###

###******************************###
##### ADD LOCATION INFORMATION #####
###******************************###

#' @description Adds location information to the dataframe
#' @param nested_dfs_names_cleaned Dataframe containing cleaned data
#' @return Dataframe with location information
#' @details
#'   - Videos from June 10 are from Neue Bahnhofstraße
#'   - Videos from June 11 and 20 are from Frankfurter Allee
#'   - Videos from June 18 are from Gürtelstraße
#'   - Neue Bahnhofstraße and Gürtelstraße are residential
#'   - Frankfurter Allee is commercial


nested_dfs_location <- nested_dfs_names_cleaned %>%
  mutate(location = case_when(
    str_detect(df_name, "20240610.*\\.csv$") ~ "Neue Bahnhofstraße",
    str_detect(df_name, "20240611.*\\.csv$") ~ "Frankfurter Allee",
    str_detect(df_name, "20240618.*\\.csv$") ~ "Gürtelstraße",
    str_detect(df_name, "20240620.*\\.csv$") ~ "Frankfurter Allee",
    TRUE ~ NA_character_)) %>%
  mutate(location_type = case_when(
    str_detect(df_name, "20240610.*\\.csv$") ~ "Residential",
    str_detect(df_name, "20240611.*\\.csv$") ~ "Commercial",
    str_detect(df_name, "20240618.*\\.csv$") ~ "Residential",
    str_detect(df_name, "20240620.*\\.csv$") ~ "Commercial",
    TRUE ~ NA_character_
  ))


###********************************###
##### ADD EXPERIMENTAL CONDITION #####
###********************************###

##### Read in experimenal condition data #####

#' @description Reads in experimental condition data from metadata file
#' @param path_input Path to input directory
#' @return Dataframe with experimental condition data
#' @details
#'   - Reads in experimental condition data from Excel file
#'   - Selects relevant columns: df_name, exp_condition, block_position, block_id, sub_block_id

condition_data <- read_excel(glue("{path_input}/Replication_DS_metadata.xlsx")) %>%
  select(df_name = CSV_Name, exp_condition = Experimental_Condition, block_position, block_id, sub_block_id)


##### Join with main data #####

nested_dfs_condition <- nested_dfs_location %>%
  left_join(condition_data)


###****************************###
##### ADD TIME OF DAY COLUMN #####
###****************************###

nested_dfs_time <- nested_dfs_condition %>%
  mutate(time_of_day = as.POSIXct(
    str_extract(df_name, "[0-9]{8}_[0-9]{6}") %>%
    str_replace("_", ""),
    format = "%Y%m%d%H%M%S"
  ))


###***************************###
##### ***CREATE CORRIDOR*** #####
###***************************###

#' @description First step of the main processing pipeline that:
#' 1. Splits data into participants and confederates
#' 2. Fits linear models to participant trajectories
#' 3. Creates analysis corridors using calculate_perpendicular_line() [lines 14-46]
#' 4. Calculates corridor corners using calculate_corridor_corners() [lines 48-86]
#' 5. Filters for in-corridor positions using filter_participants_in_corridor() [lines 88-112]


###***************************************************###
##### SPLIT DATA INTO PARTICIPANTS AND CONFEDERATES #####
###***************************************************###

#' @description Creates separate dataframes for participant and confederate trajectories
#' within a nested dataframe structure
#'
#' @details
#'   - Uses dplyr::filter() to separate rows based on 'role' column
#'   - Creates two new columns in nested structure:
#'     * participants: Contains only participant trajectory data
#'     * confederates: Contains only confederate trajectory data
#'   - Maintains original nested structure with df_name and other metadata

nested_dfs_split <- nested_dfs_time %>%
  mutate(
    # Get participant and confederate data
    participants = map(data, ~ filter(.x, role == "participant")),
    confederates = map(data, ~ filter(.x, role == "confederate"))
  )

###****************************************###
##### REMOVE TRIALS WITH NO PARTICIPANTS #####
###****************************************###

#' @description Removes any trials where no participant trajectories were recorded
#'
#' @details
#'   - Uses purrr::map_lgl() to check if participant data exists
#'   - Removes entire row from nested structure if participants column is empty

nested_dfs_filtered <- nested_dfs_split %>%
  filter(map_lgl(participants, ~ nrow(.x) > 0))


###**************************###
##### CREATE LINEAR MODELS #####
###**************************###

#' @description Creates linear regression models of participant movement paths
#' for use in corridor analysis
#'
#' @details
#'   - Fits linear model for each participant trajectory
#'   - Uses centered_z as dependent and centered_x as independent variable
#'   - Models represent the general direction of participant movement

nested_dfs_lm <- nested_dfs_filtered %>%
  mutate(
    part_lm = map(participants, ~ lm(centered_z ~ centered_x, data = .x))
  )


###************************************###
##### CREATE CONFEDERATE COORDINATES #####
###************************************###

#' @description Computes median x,z coordinates for confederates in each video
#'
#' @details
#'   - Uses median to reduce impact of possible outliers
#'   - Creates reference points for corridor construction

nested_dfs_conf_coords <- nested_dfs_lm %>%
  mutate(
    confederates_coords = map(confederates, ~ .x %>% summarize(
      conf_x = median(centered_x),
      conf_z = median(centered_z)
    ))
  )


###***********************************###
##### CALCULATE PERPENDICULAR LINES #####
###***********************************###

#' @description Creates perpendicular lines through confederate positions that
#' intersect with participant trajectory regression lines
#'
#' @details
#'   - Uses calculate_perpendicular_line() from Replication_DS_utils.R
#'   - For each video:
#'     * Takes participant trajectory linear model
#'     * Takes confederate median position
#'     * Returns perpendicular line function and intersection point

nested_dfs_perp_lines <- nested_dfs_conf_coords %>%
  mutate(
    perpendicular = map2(part_lm, confederates_coords, calculate_perpendicular_line)
  )


###********************************###
##### CALCULATE CORRIDOR CORNERS #####
###********************************###

#' @description Defines rectangular corridors for analyzing participant trajectories
#'
#' @details
#'   - Uses calculate_corridor_corners() from Replication_DS_utils.R
#'   - For each video:
#'     * Takes participant data, confederate data, linear model, and perpendicular line
#'     * Returns coordinates of corridor corners (bottom-left, top-left, top-right, bottom-right)
#'   - Creates standardized analysis region close to confederate positions

nested_dfs_corridor <- nested_dfs_perp_lines %>%
  mutate(
    corridor = pmap(list(participants, confederates, part_lm, perpendicular),
                   calculate_corridor_corners)
  )


###*******************************************************************###
##### FILTER PARTICIPANTS TO ONLY INCLUDE THOSE WITHIN THE CORRIDOR #####
###*******************************************************************###

#' @description Identifies and retains only participant positions that fall within
#' the defined analysis corridors. Creates a new nested column 'corridor_data'.
#'
#' @details
#'   - Uses filter_participants_in_corridor() from Replication_DS_utils.R
#'   - For each video:
#'     * Takes participant positions and corridor definition
#'     * Returns only positions that fall within corridor boundaries

nested_dfs_in_corridor <- nested_dfs_corridor %>%
  mutate(
    corridor_data = pmap(list(participants, corridor),
                        filter_participants_in_corridor)
  )


###*****************************###
##### CHECK CORRIDOR APPROACH #####
###*****************************###

#' @description Adds corridor membership indicator to full dataset
#'
#' @details
#'   - Uses check_corridor_membership() from Replication_DS_utils.R
#'   - For each video:
#'     * Takes full dataset, participant data, and corridor definition
#'     * Adds 'in_corridor' column ("yes"/"no") to all positions

nested_dfs_check_corridor <- nested_dfs_in_corridor %>%
  mutate(
    data = pmap(list(data, participants, corridor), check_corridor_membership)
  )


###******************************###
##### VISUAL CHECK OF CORRIDOR #####
###******************************###

#' @description Creates plots for visual inspection of corridor analysis. Produced figure is not included in main report.
#'
#' @details
#'   - Uses plot_corridor_check() and plot_corridor_samples() from Replication_DS_utils.R
#'     * Creates multi-panel plot of random samples
#'   - Shows:
#'     * Corridor boundaries
#'     * Participant positions (colored by corridor membership)
#'     * Confederate positions
#'     * Regression and perpendicular lines

#' Note: You can use plot_corridor_check() to inspect specific videos.
#' Insert any csv file name to visually inspect corridor, e.g.,
#' plot_corridor_check("20240620_104321_0_6100_1.csv")

plot_corridor_samples(nested_dfs_check_corridor)


###**********************###
##### CORRIDOR EXAMPLE #####
###**********************###

#' @description Generates and saves a visualization an example of the corridor analysis approach. Produces figure 1.
#'
#' @details
#'   - Uses plot_corridor_check() from Replication_DS_utils.R
#'   - Creates plot showing:
#'     * Corridor boundaries
#'     * Participant trajectories (colored by corridor membership)
#'     * Confederate positions
#'     * Regression and perpendicular lines
#'   - Saves plot in two formats:
#'     1. PDF
#'     2. PNG (for quick viewing and presentations)

plot_corridor_check("20240610_093009_0_839_1.csv") +
  labs(title = "")

ggsave("Replication_DS_figure 1.pdf",
  path = path_output_figures,
  device = cairo_pdf,
  width = 6,
  height = 8
)

ggsave("Replication_DS_figure 1.png",
       path = path_output_figures,
       width = 6,
       height = 8)


###*******************************###
##### ***CALCULATE DISTANCES*** #####
###*******************************###

#' @description Computes minimum distances for all in-corridor positions
#'
#' @details
#'   - Uses calculate_distances() from Replication_DS_utils.R
#'   - For each video:
#'     * Takes participant and confederate positions
#'     * Returns minimum distance to any confederate for each participant position within the corresponding frame
#'   - Unnests the distances column to create a flat dataframe

nested_dfs_corridor_distance <- nested_dfs_check_corridor %>%
  mutate(
    distances = map2(data, confederates, calculate_distances)
  ) %>%
  unnest(distances) %>%
  select(
    df_name,
    frame,
    pid,
    location,
    location_type,
    exp_condition,
    block_position,
    block_id,
    sub_block_id,
    time_of_day,
    min_distance
  )


###*********************************###
##### CHECK DISTANCE CALCULATIONS #####
###*********************************###

#' @description Performs detailed validation of distance calculations and checks results
#'
#' @details
#'   - Uses validate_distance_calculations() and add_distance_validation() from Replication_DS_utils.R
#'     * Computes distances to each confederate separately
#'     * Identifies which confederate is closer
#'     * Verifies that minimum distance matches closer confederate
#'   - Performs two validation checks:
#'     1. Identifies any videos where distances are incorrect
#'     2. Verifies that closer_conf matches actual minimum distance

nested_dfs_distance_check <- nested_dfs_check_corridor %>%
  mutate(
    distances_check = map2(data, confederates,
                          function(d, c) validate_distance_calculations(d, c, n_samples = NULL)),
    all_correct = unlist(map(distances_check, add_distance_validation))
  ) %>%
  select(df_name, distances_check, all_correct)

# Check if all distances are correct
nested_dfs_distance_check %>%
  filter(!all_correct)

# Check if the closer_conf matches which distance is smaller
nested_dfs_distance_check %>%
  unnest(distances_check) %>%
  mutate(
    actual_closer = if_else(dist_conf1 < dist_conf2, 1, 2)
  ) %>%
  filter(closer_conf != actual_closer)


###************************###
##### ***DRAW SAMPLES*** #####
###************************###

#' @description Creates 50 random sub-samples of 3000 observations each,
#' based on our estimations that the locations will result in a sample size around N = 3000.
#'
#' @details
#'   - Sets random seed for reproducibility
#'   - Creates tibble with 50 samples
#'   - Each sample contains 3000 randomly selected observations

set.seed(321)

subsamples_df <- tibble(sample_no = 1:50, data = map(1:50, ~ nested_dfs_corridor_distance %>%
  sample_n(3000)))


###********************************************###
##### TABLE: NUMBER OF OBS IN EACH CONDITION #####
###********************************************###

###*******************************###
##### By experimental condition #####
###*******************************###

#' @description Generates and saves a table showing the number of observations
#' in each experimental condition
#'
#' @details
#'   - Cleans condition names for output
#'   - Creates formatted flextable
#'   - Saves as Word document

nested_dfs_corridor_distance %>%
  group_by(exp_condition) %>%
  summarise(n = n()) %>%
  mutate(
    exp_condition = case_when(
      exp_condition == "muslim_garb" ~ "Muslim garb`",
      exp_condition == "muslim_no_garb" ~ "Muslim no garb",
      exp_condition == "white" ~ "White",
      TRUE ~ exp_condition
    )
  ) %>%
  flextable() %>%
  set_caption("Number of Observations by Condition") %>%
  set_header_labels(
    exp_condition = "Experimental Condition",
    n = "Number of Observations"
  ) %>%
  autofit() %>%
  bold(part = "header") %>%
  save_as_docx(path = file.path(path_output_tables, "condition_n.docx"))


###*****************************************###
##### By experimental condition and block #####
###*****************************************###

#' @description Generates and saves a table showing the number of observations
#' in each experimental condition and block

nested_dfs_corridor_distance %>%
  group_by(exp_condition, block_id) %>%
  summarise(n = n()) %>%
  mutate(
    exp_condition = case_when(
      exp_condition == "muslim_garb" ~ "Muslim garb`",
      exp_condition == "muslim_no_garb" ~ "Muslim no garb",
      exp_condition == "white" ~ "White",
      TRUE ~ exp_condition
    )
  ) %>%
  flextable() %>%
  set_caption("Number of Observations by Condition and Block") %>%
  set_header_labels(
    exp_condition = "Experimental Condition",
    block_id = "Block",
    n = "Number of Observations"
  ) %>%
  autofit() %>%
  bold(part = "header") %>%
  save_as_docx(path = file.path(path_output_tables, "condition_block_n.docx"))


###*****************###
##### By location #####
###*****************###

#' @description Generates and saves a table showing the number of observations
#' in each location

nested_dfs_corridor_distance %>%
  group_by(location) %>%
  summarise(n = n()) %>%
  flextable() %>%
  set_caption("Number of Observations by Location") %>%
  set_header_labels(
    location = "Location",
    n = "Number of Observations"
  ) %>%
  autofit() %>%
  bold(part = "header") %>%
  save_as_docx(path = file.path(path_output_tables, "location_n.docx"))


###*********************###
##### ***SAVE DATA*** #####
###*********************###

#' @description Saves to RDS containing processed data

write_rds(subsamples_df, file.path(path_output_data, "Replication_DS_subsamples.rds"))


