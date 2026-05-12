# ============================================================
# Train a Random Forest model
# ============================================================

# Load required libraries
library(sits)
library(ggplot2)
library(stringr)
library(randomForestExplainer, lib.loc = "/opt/r/R/x86_64-pc-linux-gnu-library/4.4")

# Define the parameters: These are user-defined variables
time_series_name  <- "TS-tiles_012014-012015-013014-013015_1y_2024-08-01_2025-07-31_all-samples-new-pol-avg-false_2026-02-24_20h01m.rds"

# Extract the tiles and date of the string separated by "_"
tiles      <- str_split(str_extract(time_series_name, "(?<=tiles_)[^_]+"), "-")[[1]]
start_date <- stringr::str_split_i(time_series_name, "_", 4)
end_date   <- stringr::str_split_i(time_series_name, "_", 5)

# Calculate the number of years in the training cube
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")
tiles_train <- paste(sort(tiles), collapse = "-")
no.cubes <- paste0(length(tiles_train), "t")

# Function to read class names and their colors::IMPORTANT
read_class_config <- function(config_file = "class_config.txt") {
  
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }
  
  lines <- readLines(config_file, encoding = "UTF-8", warn = FALSE)
  
  # Remove empty lines and comments
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0 & !startsWith(lines, "#")]
  
  # Identify sections and populate lists
  current_section  <- NULL
  class_trans_list <- list()
  colors_list      <- list()
  
  for (line in lines) {
    if (startsWith(line, "[") && endsWith(line, "]")) {
      current_section <- gsub("\\[|\\]", "", line)
      next
    }
    
    if (!is.null(current_section) && grepl("=", line)) {
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      key   <- trimws(parts[1])
      value <- trimws(paste(parts[-1], collapse = "=")) # preserves '=' in hex codes
      
      if (current_section == "CLASS_TRANSLATION") {
        class_trans_list[[key]] <- value
      } else if (current_section == "COLORS") {
        colors_list[[key]] <- value
      }
    }
  }
  
  class_translation <- unlist(class_trans_list)
  my_colors         <- unlist(colors_list)
  
  message(sprintf("Config loaded: %d class translations | %d colors",
                  length(class_translation), length(my_colors)))
  
  return(list(
    class_translation = class_translation,
    my_colors         = my_colors
  ))
}

# Date and time of the start of processing
date_process    <- format(Sys.Date(), "%Y-%m-%d_")
time_process    <- format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo")
process_version <- paste0(date_process, time_process)

# File and folder paths
time_series_path  <- file.path("data/rds/time_series/", time_series_name)
rds_path          <- "data/rds/"
plots_path        <- "data/plots/"
config_dir        <- ".."

# Identifier to distinguish this model run from previous versions
var <- stringr::str_split_i(time_series_name, "_", 6)

# Plots organized by var
plots_dir <- file.path(plots_path, var)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Cross-validation of training data
# ============================================================

# Step 1.1 -- Reading training samples
train_samples <- readRDS(time_series_path)

# Step 1.2 -- Load color palette from external config file
config     <- read_class_config(file.path(config_dir, "class_config.txt"))
my_colors  <- config$my_colors
my_colors  <- my_colors[names(my_colors) %in% unique(train_samples$label)]

# Step 1.3 -- Using k-fold validation
sits_kfold_validate_start <- Sys.time()
rfor_validate <- sits_kfold_validate(
  samples = train_samples,
  folds = 5, # how many times to split the data (default = 5)
  ml_method = sits_rfor(),
  multicores = 28,
  progress = TRUE) # adapt to your computer CPU core availability
sits_kfold_validate_end <- Sys.time()
sits_kfold_validate_time <- as.numeric(sits_kfold_validate_end - sits_kfold_validate_start, units = "secs")
sprintf("SITS kfold_validate process duration (HH:MM): %02d:%02d", 
        as.integer(sits_kfold_validate_time / 3600),
        as.integer((sits_kfold_validate_time %% 3600) / 60))

# Step 1.3.1 -- Plot the confusion matrix
plot(rfor_validate, type = "confusion_matrix")

# Step 1.3.2 -- Plot the metrics by class
plot(rfor_validate, type = "metrics")

# Step 1.4 -- Save confusion matrix plot
g_cm <- plot(rfor_validate, type = "confusion_matrix")
ggplot2::ggsave(
  filename = file.path(
    plots_dir,
    paste0(
      "Kfold-confusion-matrix_",
      tiles_train, "_",
      start_date, "_", end_date, "_",
      var, "_",
      format(Sys.Date(), "%Y-%m-%d"),
      ".png"
    )
  ),
  plot = g_cm,
  width = 1600,
  height = 1000,
  units = "px",
  dpi = 200
)

# Step 1.4.1 -- Save metrics plot
g_metrics <- plot(rfor_validate, type = "metrics")
ggplot2::ggsave(
  filename = file.path(
    plots_dir,
    paste0(
      "Kfold-metrics_",
      tiles_train, "_",
      start_date, "_", end_date, "_",
      var, "_",
      format(Sys.Date(), "%Y-%m-%d"),
      ".png"
    )
  ),
  plot = g_metrics,
  width = 1600,
  height = 1000,
  units = "px",
  dpi = 200
)

# ============================================================
# 2. Training and saving model
# ============================================================

# Step 2.1 -- Set a seed of random number generator (RNG) for reproducibility
set.seed(88)

# Step 2.2 -- Train the model
rf_model <- sits_train(
  samples   = train_samples,
  ml_method = sits_rfor(num_trees = 100)
)

# Step 2.3 -- Save the ML model to a R file
saveRDS(rf_model,
        paste0(rds_path, "model/random_forest/",
               paste("rf-model", no.cubes,
                     tiles_train, no.years,
                     start_date, end_date,
                     var, process_version, sep = "_"),
               ".rds"))

print("Model trained successfully!")

# ============================================================
# 3. Plotting Section
# ============================================================

# Step 3.1 -- Define the function to plot and save the most important variables of the model
save_rf_model_plot <- function(
    rf_model,
    plots_dir,
    tiles,
    no.years,
    start_date,
    end_date,
    var,
    width  = 1200,
    height = 800,
    res    = 150,
    scale  = 1
) {
  
  # Generate the native sits/randomForestExplainer plot
  g <- plot(rf_model)
  
  # Render in RStudio
  print(g)
  
  # Build file name
  tiles_str <- paste(tiles, collapse = "-")
  file_name <- paste0(
    "RF-minimal-tree-depth",
    "_", tiles_str,
    "_", no.years,
    "_", start_date,
    "_", end_date,
    "_", var,
    "_", format(Sys.Date(), "%Y-%m-%d"),
    ".png"
  )
  
  # Save
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  full_path <- file.path(plots_dir, file_name)
  
  ggplot2::ggsave(full_path, plot = g, width = width, height = height,
                  units = "px", dpi = res, scale = scale)
  
  message("Plot saved: ", full_path)
  invisible(full_path)
}

# Step 3.2 -- Run the function to plot and save the most important variables of the model
save_rf_model_plot(
  rf_model   = rf_model,
  plots_dir  = plots_dir,
  tiles      = tiles,
  no.years   = no.years,
  start_date = start_date,
  end_date   = end_date,
  var        = var,
  width      = 1600,   # width in pixels
  height     = 1000,   # height in pixels
  res        = 200,    # DPI
  scale      = 1       # increases all elements proportionally  
)

# Step 3.3 --  Define the function to plot and save Out of Box error by the number of trees
save_rf_oob_plot <- function(
    rf_model,
    plots_dir,
    tiles,
    no.years,
    start_date,
    end_date,
    var,
    width  = 1200,
    height = 800,
    res    = 150,
    scale  = 1
) {
  
  # Export the model object
  rf_model2 <- environment(rf_model)$model
  
  # Convert err.rate matrix to tidy data frame for ggplot
  err_df <- as.data.frame(rf_model2$err.rate)
  err_df$ntree <- seq_len(nrow(err_df))
  
  err_long <- tidyr::pivot_longer(
    err_df,
    cols      = -ntree,
    names_to  = "Class",
    values_to = "OOB_Error"
  )
  
  # Build ggplot
  g <- ggplot2::ggplot(err_long, ggplot2::aes(x = ntree, y = OOB_Error, color = Class)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(
      title    = "Out-of-Bag Error by Number of Trees",
      subtitle = paste0(paste(tiles, collapse = ", "), " | ", start_date, " to ", end_date),
      x        = "Number of Trees (ntree)",
      y        = "OOB Error",
      color    = "Class"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle      = ggplot2::element_text(color = "gray40", size = 9),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "right"
    )
  print(g)
  
  # Build file name
  tiles_str <- paste(tiles, collapse = "-")
  file_name <- paste0(
    "RF-oob-ntree-mde",
    "_", tiles_str,
    "_", no.years,
    "_", start_date,
    "_", end_date,
    "_", var,
    "_", format(Sys.Date(), "%Y-%m-%d"),
    ".png"
  )
  
  # Save
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  full_path <- file.path(plots_dir, file_name)
  
  ggplot2::ggsave(full_path, plot = g, width = width, height = height,
                  units = "px", dpi = res, scale = scale)
  
  message("Plot saved: ", full_path)
  invisible(full_path)
}

# Step 3.4 --  Define the function to plot and save Out of Box error by the number of trees
save_rf_oob_plot(
  rf_model   = rf_model,
  plots_dir  = plots_dir,
  tiles      = tiles,
  no.years   = no.years,
  start_date = start_date,
  end_date   = end_date,
  var        = var,
  width      = 1600,   # width in pixels
  height     = 1000,   # height in pixels
  res        = 200,    # DPI
  scale      = 1.5     # increases all elements proportionally  
)
