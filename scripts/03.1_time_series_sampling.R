# ============================================================
# Time series extraction from the training samples
# ============================================================

# Load required libraries
library(sits)
library(tibble)
library(dplyr)
library(ggplot2)

# Define the parameters: These are user-defined variables
tiles           <- c("012014","012015","013014","013015")
start_date      <- "2024-08-01"
end_date        <- "2025-07-31"
var             <- "all-samples-new-pol-avg-false" #ALWAYS SPACE THE WORDS WITH "-"
sampling_date   <- "2026-02-24"                    # Date of the sampling file (YYYY-MM-DD)

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

# File and folder paths
sample_path   <- "data/raw/samples"
ts_path       <- "data/rds/time_series/"
mixture_path  <- "data/raw/mixture_model"
plots_dir    <- "data/plots/time_series"
config_dir    <- ".."

# ============================================================
# 1. Define and Load Data Cubes
# ============================================================

# Step 1.1 -- Create a training cube from a collection
cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
  tiles       = tiles,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 1.2 -- Calculate the number of years in the training cube
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")

# Step 1.3 -- Concatenates all the names of the training tiles into a single string separated by '-'
tiles_train <- paste(sort(tiles), collapse = "-")

# Step 1.4 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source      = "BDC",
  collection  = "SENTINEL-2-16D",
  bands       = c("SOIL", "VEG", "WATER"),
  tiles       = tiles,
  data_dir    = mixture_path,
  start_date  = start_date,
  end_date    = end_date,
  progress    = TRUE)

# Step 1.5 -- Merge the Training Cube with Mixture Model Cube
cube_merge_lsmm_train <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Load and Explore Training Sample Data
# ============================================================

# Step 2.1 -- Read training samples (rewrite the name of your samples file)
samples_name    <- paste("training-samples", tiles_train, var, sampling_date, sep = "_")
samples_train   <- sf::st_read(file.path(sample_path, paste0(samples_name, ".gpkg")))

# Step 2.2 -- Load class translation from external config file
config     <- read_class_config(file.path(config_dir, "class_config.txt"))
class_translation <- config$class_translation

# Step 2.2.1 -- Apply translation: keep the original label if no translation is found
samples_train$label <- ifelse(
  samples_train$label %in% names(class_translation),
  class_translation[samples_train$label],
  samples_train$label
)

print(table(samples_train$label))

# ============================================================
# 3. Time Series Extraction
# ============================================================

# Step 3.1 -- Extract Time Series from samples_train and calculate the process duration
sits_get_data_start <- Sys.time()
samples <- sits_get_data(
  cube        = cube_merge_lsmm_train,
  samples     = samples_train,
  n_sam_pol   = 16,
  pol_avg     = FALSE,
  label       = "label",
  multicores  = 28,       # adapt to your computer CPU core availability
  progress    = TRUE)
sits_get_data_end <- Sys.time()
sits_get_data_time <- as.numeric(sits_get_data_end - sits_get_data_start, units = "secs")
sprintf("SITS get data process duration (HH:MM): %02d:%02d",
        as.integer(sits_get_data_time / 3600),
        as.integer((sits_get_data_time %% 3600) / 60))
print("Time series extracted successfully!")

# Step 3.2 -- Save the samples Time Series to a R file
saveRDS(samples, 
        paste0(ts_path,
               paste("TS-tiles", tiles_train,
                     no.years, start_date,
                     end_date, var, format(Sys.Date(), "%Y-%m-%d"), 
                     format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo"),
                     sep = "_"),
               ".rds"))

# ============================================================
# 4. Plotting Time Series Patterns
# ============================================================

# Step 4.1 -- Define function to visualize and save the plot of the time-spectral patterns of all features
save_sits_patterns_plot <- function(samples,
                                    start_date,
                                    end_date,
                                    plots_dir,
                                    tiles,
                                    var,
                                    labels          = NULL,
                                    bands           = NULL,
                                    vline_dates     = NULL,   # ex: c("08-01") → mouth/day in each year
                                    legend_text_size = NULL,  # font size band legend
                                    class_text_size  = NULL,  # font size class name
                                    line_width       = NULL,  # width of the band lines
                                    vline_width      = 0.6,   # width of the vertical line
                                    width            = 1200,
                                    height           = 800,
                                    res              = 150) {
  
  # Filter only if specified
  s <- samples
  if (!is.null(labels)) s <- sits_select(s, labels = labels)
  if (!is.null(bands))  s <- sits_select(s, bands  = bands)
  
  # It calculates the pattern only once
  p <- sits_patterns(s)
  
  # Captures the ggplot object without rendering it to the screen yet
  g <- plot(p)
  
  # --- Vertical dotted line indicating the year ---
  if (!is.null(vline_dates)) {
    years <- seq(
      as.integer(format(as.Date(start_date), "%Y")),
      as.integer(format(as.Date(end_date),   "%Y"))
    )
    vlines <- as.Date(paste0(years, "-", vline_dates))
    vlines <- vlines[vlines >= as.Date(start_date) & vlines <= as.Date(end_date)]
    
    g <- g + ggplot2::geom_vline(
      xintercept = as.numeric(vlines),
      linetype   = "dashed",
      color      = "gray40",
      linewidth  = vline_width
    )
  }
  
  # Optional theme customizations
  theme_args <- list()
  
  if (!is.null(legend_text_size))
    theme_args$legend.text <- ggplot2::element_text(size = legend_text_size)
  
  if (!is.null(class_text_size))
    theme_args$strip.text <- ggplot2::element_text(size = class_text_size)
  
  if (length(theme_args) > 0)
    g <- g + do.call(ggplot2::theme, theme_args)
  
  # --- Grossura das linhas das bandas ---
  if (!is.null(line_width))
    g <- g + ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(linewidth = line_width)
      )
    ) + ggplot2::geom_line(linewidth = line_width) +
    ggplot2::labs(color = "Bands")              # <-- força o título correto
  
  # Renders on screen
  print(g)
  
  # Set file name
  suffix <- ""
  
  if (!is.null(labels) && length(labels) == 1) {
    suffix <- paste0("_", labels[1])
    if (!is.null(bands))
      suffix <- paste0(suffix, "-", paste(tolower(bands), collapse = "-"))
  } else if (!is.null(labels) || !is.null(bands)) {
    if (!is.null(labels))
      suffix <- paste0("_", length(labels), "classes")
    if (!is.null(bands))
      suffix <- paste0(suffix, "-", paste(tolower(bands), collapse = "-"))
  }
  
  tiles_str <- paste(tiles, collapse = "-")
  file_name <- paste0("sits-patterns",
                      "_tiles-", tiles_str,
                      "_", start_date,
                      "_", end_date,
                      "_", var,
                      suffix,
                      ".png")
  
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  full_path <- file.path(plots_dir, file_name)
  ggplot2::ggsave(full_path, plot = g, width = width, height = height,
                  units = "px", dpi = res)
  
  message("Plot saved in: ", full_path)
  invisible(full_path)
}

# Step 4.2 -- Run function to visualize and save the plot of the time-spectral patterns of all features
save_sits_patterns_plot(
  samples          = samples,
  start_date       = unique(samples$start_date),
  end_date         = unique(samples$end_date),
  plots_dir        = plots_dir,
  tiles            = tiles,
  var              = var,
  bands            = c('B11','EVI','NDVI'), # NULL to plot and save all bands patterns
  #labels            = c('DESMAT_ARVORE_REMANESCE'), # NULL to plot and save all classes patterns
  vline_dates      = "08-01",   # vertical doted line on August 1st of each year
  legend_text_size = 10, 
  class_text_size  = 12,
  line_width       = 1.5, # line width of each spectral band
  vline_width      = 0.6  # vertical line width
)
