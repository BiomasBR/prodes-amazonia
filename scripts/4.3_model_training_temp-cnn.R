# ============================================================
# Train a Temporal Convolution Neural Network
# ============================================================

# Load required libraries
library(sits)
library(ggplot2)
library(torch)
library(luz)
library(stringr)

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

# Step 2.1 -- Set random seed to ensure reproducibility
set.seed(88)

# Step 2.2 -- Train TempCNN model using training samples
tempcnn_model <- sits_train(
  samples = train_samples,
  ml_method = sits_tempcnn(
    
    # Convolutional layers configuration
    cnn_layers = c(64, 64, 64),              # number of filters per layer
    cnn_kernels = c(3, 3, 3),                # kernel size for each convolution
    cnn_dropout_rates = c(0.2, 0.2, 0.2),    # dropout rate to reduce overfitting
    
    # Fully connected (dense) layer configuration
    dense_layer_nodes = 256,                 # number of neurons in dense layer
    dense_layer_dropout_rate = 0.5,          # dropout for dense layer
    
    # Training configuration
    epochs = 150,                            # number of training epochs
    batch_size = 1024,                       # batch size for training
    validation_split = 0.2,                  # proportion of data used for validation
    
    # Optimizer configuration (AdamW)
    optimizer = torch::optim_adamw,
    opt_hparams = list(
      lr = 5e-04,                            # learning rate
      eps = 1e-08,                           # numerical stability term
      weight_decay = 1e-06                   # L2 regularization
    ),
    
    # Learning rate scheduling
    lr_decay_epochs = 1,                     # frequency of decay (in epochs)
    lr_decay_rate = 0.95,                    # decay factor
    
    # Early stopping configuration
    patience = 20,                           # epochs to wait without improvement
    min_delta = 0.01,                        # minimum improvement threshold
    verbose = FALSE                          # disable training logs
  )
)

# Step 2.2.1 -- Plot model diagnostics (e.g., training history or feature relevance)
plot(tempcnn_model)

# Step 2.3 -- Save trained model to disk (RDS format)
saveRDS(
  tempcnn_model,
  paste0(
    rds_path, "model/temp_cnn/",
    paste(
      "tcnn-model",
      length(tiles),
      tiles_train, no.years,
      start_date, end_date,
      var, process_version,sep = "_"),
    ".rds"))

print("Model trained successfully!")

# ============================================================
# Step 3 -- Save model parameters to txt
# ============================================================

# Folder: same as the model RDS
model_dir <- file.path(rds_path, "model/temp_cnn/")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

# File name uses same pattern as model
params_filename <- paste0(
  "TCNN-parameters_",
  length(tiles), "-tiles-", tiles_train, "_",
  no.years, "-period-",
  start_date, "_", end_date,
  "_", var, "_", process_version, ".txt"
)

# Parameters used in training
params_lines <- c(
  "# ============================================================",
  "# TempCNN Model Parameters",
  "# ============================================================",
  "",
  "# --- General information ---",
  paste0("process_version         : ", process_version),
  paste0("var                     : ", var),
  paste0("tiles                   : ", tiles_train),
  paste0("period_years            : ", no.years),
  paste0("number_of_tiles         : ", length(tiles)),
  paste0("start_date              : ", start_date),
  paste0("end_date                : ", end_date),
  "",
  "# --- CNN architecture ---",
  paste0("cnn_layers              : ", paste(c(64, 64, 64), collapse = ", ")),
  paste0("cnn_kernels             : ", paste(c(3, 3, 3), collapse = ", ")),
  paste0("cnn_dropout_rates       : ", paste(c(0.2, 0.2, 0.2), collapse = ", ")),
  "",
  "# --- Dense layer ---",
  paste0("dense_layer_nodes       : ", 256),
  paste0("dense_dropout_rate      : ", 0.5),
  "",
  "# --- Training configuration ---",
  paste0("epochs                  : ", 150),
  paste0("batch_size              : ", 1024),
  paste0("validation_split        : ", 0.2),
  paste0("random_seed             : ", 88),
  "",
  "# --- Optimizer (AdamW) ---",
  paste0("optimizer               : AdamW"),
  paste0("learning_rate           : ", 5e-04),
  paste0("eps                     : ", 1e-08),
  paste0("weight_decay            : ", 1e-06),
  "",
  "# --- Learning rate decay ---",
  paste0("lr_decay_epochs         : ", 1),
  paste0("lr_decay_rate           : ", 0.95),
  "",
  "# --- Early stopping ---",
  paste0("patience                : ", 20),
  paste0("min_delta               : ", 0.01),
  "",
  "# --- Execution info ---",
  paste0("training_datetime       : ", Sys.time())
)
# Save txt file
writeLines(
  params_lines,
  file.path(model_dir, params_filename)
)
message("Parameters saved: ", params_filename)
