# ============================================================
#  Validation and accuracy measurement in Non-Forest Areas
# ============================================================

# Load required libraries
library(tibble)
library(sits)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(stringr)
library(segmetric)
library(openxlsx)

# Define the parameters: These are user-defined variables
tiles      = '000000'
model_name <- "rf-model_2t_014002-015002_2y_2023-07-28_2025-07-28_com-nuvens-cheias_2026-04-07_14h45m.rds"

# Extract the date of the string separated by "_"
start_date <- stringr::str_split_i(model_name, "_", 5)
end_date   <- stringr::str_split_i(model_name, "_", 6)

# File and folder paths
models <- c("rf"   = "random_forest",
            "xgb"  = "xgboost",
            "ltae" = "ltae",
            "tcnn" = "temp_cnn",
            "rnet" = "res_net",
            "lstm" = "ltsm")
model_type       <- stringr::str_split_i(model_name, "-", 1)
model_path       <- file.path("data/rds/model", models[model_type], model_name)
model            <- readRDS(model_path)
class_dir        <- "data/class"
samples_dir      <- "data/raw/samples/validation_samples"
plots_path       <- "data/plots"
mask_dir         <- "data/raw/auxiliary/masks"
version          <- paste(stringr::str_split_i(model_name, "-", 1),
                          stringr::str_split_i(model_name, "_", 4),
                          stringr::str_split_i(model_name, "_", 7),
                          sep = "-")

# Plots organized by version
plots_dir <- file.path(plots_path, version)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# List of validation sample files matching the version pattern in the samples directory
pattern <- paste0(".*", tiles, ".*", ".*", version, ".*\\.gpkg$")

samples_validation_list <- dir(
  samples_dir,
  pattern = pattern,
  full.names = TRUE
)

prodes_avaliation <- function(x){
  data.frame(
    t(
      c(
        Mínimo = format(min(x, na.rm = TRUE), digits = 3),
        Máximo = format(max(x, na.rm = TRUE), digits = 3),
        Média  = format(mean(x, na.rm = TRUE), digits = 3),
        Mediana = format(median(x, na.rm = TRUE), digits = 3),
        Desvio_Padrão = format(sd(x, na.rm = TRUE), digits = 3)
      )
    )
  )
}

# Plotting function
plot_accuracy <- function(acc, version, tile, plots_dir, prefix) {
  
  today     <- format(Sys.Date(), "%Y-%m-%d")
  save_plot <- function(filename, width = 3529, height = 1578) {
    ggsave(
      filename = paste0(filename, "-", prefix, "_", tile, "_", version, "_", today, ".png"),
      path     = plots_dir,
      scale    = 1,
      width    = width,
      height   = height,
      units    = "px",
      dpi      = 350
    )
  }
  
  label_other <- c("Other_Classes" = "Other Classes")
  
  # ── 1. Confusion Matrix ───────────────────────────────────────────────────────
  matriz_conf <- tibble(as.data.frame(acc$error_matrix))
  
  p1 <- ggplot(matriz_conf, aes(x = Var2, y = Var1, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 4) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "Cases") +
    labs(x = "Reference", y = "Predicted", title = "Confusion Matrix") +
    scale_y_discrete(limits = rev, labels = label_other) +
    scale_x_discrete(labels = label_other) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title  = element_text(hjust = 0.5)
    )
  
  print(p1)
  save_plot("confusion-matrix")
  
  # ── 2. Accuracies ─────────────────────────────────────────────────────────────
  acuracias <- data.frame(
    class         = names(acc$accuracy[[1]]),
    user_accuracy = round(as.numeric(acc$accuracy[[1]]), 3),
    prod_accuracy = round(as.numeric(acc$accuracy[[2]]), 3)
  ) %>%
    dplyr::mutate(
      f1_score = round(2 * (user_accuracy * prod_accuracy) / (user_accuracy + prod_accuracy), 3)
    ) %>%
    tidyr::pivot_longer(cols = -class, names_to = "tipo_acuracia", values_to = "acuracia")
  
  p2 <- ggplot(acuracias, aes(x = tipo_acuracia, y = class, fill = acuracia)) +
    geom_tile(color = "white") +
    geom_text(aes(label = acuracia), size = 4) +
    scale_fill_distiller(palette   = "RdYlGn",
                         direction = 1,
                         name      = "Value",
                         limits    = c(0, 1)) +
    labs(
      y       = "Class",
      x       = "Accuracy",
      title   = "Accuracies",
      caption = paste0("Global Accuracy: ", round(acc$accuracy[[3]], 3))
    ) +
    scale_y_discrete(limits = rev, labels = label_other) +
    scale_x_discrete(labels = c("prod_accuracy" = "Prod Acc",
                                "user_accuracy" = "User Acc",
                                "f1_score"      = "F1 Score")) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 1, size = 11, margin = margin(t = 10)),
      panel.grid   = element_blank()
    )
  
  print(p2)
  save_plot("metrics")
  
  # ── 3. Area Metrics ───────────────────────────────────────────────────────────
  class_areas <- data.frame(
    class             = names(acc$area_pixels),
    mapped_area_ha    = round(as.numeric(acc$area_pixels), 2),
    error_adj_area_ha = round(as.numeric(acc$error_ajusted_area), 2),
    conf_interval_ha  = round(as.numeric(acc$conf_interval), 2)
  ) %>%
    tidyr::pivot_longer(cols = -class, names_to = "tipo_area", values_to = "area")
  
  p3 <- ggplot(class_areas, aes(x = tipo_area, y = class, fill = area)) +
    geom_tile(color = NA, fill = "white") +
    geom_text(aes(label = area), size = 4) +
    geom_vline(xintercept = c(1.5, 2.5), color = "grey30") +
    labs(y = "Class", x = "Metrics", title = "Area Metrics") +
    scale_y_discrete(limits = rev, labels = label_other) +
    scale_x_discrete(
      limits = rev,
      labels = c(
        "mapped_area_ha"    = "Mapped Area (ha)",
        "error_adj_area_ha" = "Error-Adjusted Area (ha)",
        "conf_interval_ha"  = "Conf Interval (ha)"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid      = element_blank()
    )
  
  print(p3)
  save_plot("areas")
  
  message("Plots salvos em: ", plots_dir)
  
  invisible(list(confusion_matrix = p1, accuracies = p2, area_metrics = p3))
}

# ============================================================
# 1. Accuracy assessment of Full Map classified images
# ============================================================

# Step 1.1 -- Get labels associated to the trained model data set (Enumerate them in the order they appear according to "sits_labels(model)")
pattern <- paste0(".*", tiles, ".*", version, ".*\\.tif$")

cube_dirs <- list.dirs(class_dir, recursive = TRUE) |> 
  purrr::keep(~ length(list.files(.x, pattern = pattern)) > 0)

# Step 1.2 -- Retrieve local cube of Full Map classified
cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = c("1"  = "Burn_Scar", # List the classes according to each number sequence in your raster
             "2"  = "Lake",
             "3"  = "River",
             "4"  = "Clearing_For_Agriculture",
             "5"  = "Previous_Clearing_For_Agriculture",
             "6"  = "Clearing_For_Bare_Soil",
             "7"  = "Previous_Clearing_For_Bare_Soil",
             "8"  = "Herbaceous_Dry_High_Biomass",
             "9"  = "Herbaceous_Dry_Low_Biomass",
             "10" = "Herbaceous_Dry_Post_Fire",
             "11" = "Herbaceous_Wet",
             "12" = "Woodland",
             "13" = "Ecotone",
             "14" = "Palm_Swamps"),
  tiles =  tiles,
  start_date = start_date,
  end_date = end_date,
  version = version,
  data_dir = cube_dirs,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# Step 1.2 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep(".*_all-classes_*.",
                                   samples_validation_list, value = TRUE))

# Step 1.3 -- Calculate accuracy
full_map_acc <- sits_accuracy(cube,
                              validation = samples_validation,
                              memsize = 180,
                              multicores = 28) # adapt to your computer CPU core availability

# Step 1.4 -- Print the area estimated accuracy
full_map_acc

# Step 1.5 -- Show confusion matrix
full_map_acc$error_matrix

# Step 1.6 -- Plotting Full Map Accuracy
plot_accuracy(
  acc       = full_map_acc,
  version   = version,
  tile      = tiles,
  plots_dir = plots_dir,
  prefix    = "full-map-acc"
)

# ============================================================
# 2. Accuracy assessment of PRODES Adjusted Map classified
# ============================================================

# Step 2.1 -- Retrieve local cube of PRODES adjusted map classified
class_cube <- sits_cube(
  source = "BDC",
  collection = "SENTINEL-2-16D",
  bands = "class",
  labels = c("15" = "Deforestation",  # list the grouped classes according to each number they appear in raster 
             "16" = "Water", 
             "17" = "Grassland",
             "18" = "Forest"),
  tiles =  tiles,
  start_date = start_date,
  end_date = end_date,
  version = paste0(version, "-mosaic"),
  data_dir = cube_dirs,
  parse_info = c("satellite", "sensor", "tile", "start_date", "end_date", 
                 "band", "version"))

# Step 2.2 -- Get validation samples points (in geographical coordinates - lat/long)
samples_validation <- st_read(grep(".*_prodes_*.",
                                   samples_validation_list, value = TRUE))

# Step 2.3 -- Calculate accuracy
prodes_acc <- sits_accuracy(class_cube, 
                            validation = samples_validation,
                            memsize = 180,
                            multicores = 28) # adapt to your computer CPU core availability

# Step 2.4 -- Print the area estimated accuracy
prodes_acc

# Step 2.5 -- Show confusion matrix
prodes_acc$error_matrix

# Step 2.6 -- Plotting PRODES Adjusted Map Accuracy
plot_accuracy(
  acc       = prodes_acc,
  version   = paste0(version, "-mosaic"),
  tile      = tiles,
  plots_dir = plots_dir,
  prefix    = "prodes-acc"
)
