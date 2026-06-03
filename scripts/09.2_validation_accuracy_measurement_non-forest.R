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
plots_dir       <- "data/plots/accuracy"
mask_dir         <- "data/raw/auxiliary/masks"
version          <- paste(stringr::str_split_i(model_name, "-", 1),
                          stringr::str_split_i(model_name, "_", 4),
                          stringr::str_split_i(model_name, "_", 7),
                          sep = "-")

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

# =============================================================================
# Análise de Incerteza por Classe — SITS  (inputs em .gpkg)
# =============================================================================
#Readinf Files
pattern_tif <- paste0(".*_", tile, ".*_entropy_", version, "\\.tif$")

uncertainty_raster_path <- list.files(class_dir,
                                      pattern = pattern_tif, 
                                      recursive = TRUE,
                                      full.names = TRUE)

unc_r <- terra::rast(uncertainty_raster_path)
names(unc_r) <- "uncertainty"

#Extracting entropy
message("\n[3/6] Extracting Entropy...")
unc_vals        <- exactextractr::exact_extract(unc_r,
                                                samples_validation,
                                                fun = "mode",
                                                append_cols = "pol_id",
                                                force_df = TRUE,
                                                max_cells_in_memory = 1e+09)
names(unc_vals) <- c("pol_id", "uncertainty")

sf_cls <- dplyr::left_join(samples_validation, unc_vals, by = "pol_id") |>
  sf::st_drop_geometry()


#Statistcs
message("\n[5/6] Calculando estatísticas...")
stats_class <- sf_cls |>
  group_by(class) |>
  summarise(
    n_segs          = n(),
    inc_mean       = mean(uncertainty,     na.rm = TRUE),
    inc_median     = median(uncertainty,   na.rm = TRUE),
    inc_sd          = sd(uncertainty,       na.rm = TRUE),
    inc_q25         = quantile(uncertainty, 0.25, na.rm = TRUE),
    inc_q75         = quantile(uncertainty, 0.75, na.rm = TRUE),
    inc_max         = max(uncertainty,      na.rm = TRUE),
    #pct_alta_incert = mean(uncertainty >= limiar_alto, na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  arrange(desc(inc_mean))
readr::write_excel_csv2(stats_class,
                        file.path(output_dir,
                                  "stats_incerteza_por_classe.csv"))

#Generate plots
message("[6/6] Gerando gráficos...")
n_cls <- nlevels(as.factor(sf_cls$class))
colors <- scales::hue_pal()(n_cls)

# ── 8.1  Violin + Boxplot — distribuição de entropia por classe ──────────────
p_box <- sf_cls |>
  mutate(class = fct_reorder(class, uncertainty, .fun = median, .desc = TRUE)) |>
  ggplot(aes(x = class, y = uncertainty, fill = class)) +
  geom_violin(alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.3, outlier.size = 0.4, outlier.alpha = 0.3) +
  scale_fill_manual(values = colors, guide = "none") +
  labs(title    = "Distribuição de Entropia por Classe",
       #subtitle = "Linha vermelha = limiar de alta incerteza",
       x = NULL,
       y = "Entropia") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title  = element_text(face = "bold"))

ggsave(file.path(plots_dir, "boxplot_incerteza_por_classe.png"),
       p_box, width = 12, height = 6, dpi = 150, bg = "white")
plot(p_box)

message("  Gerando heatmap de probabilidades...")

heat_df <- sf_cls |>
  group_by(class) |>
  summarise(across(5:16, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") |>
  pivot_longer(-class,
               names_to = "class_prob",
               values_to = "mean_prob")

p_heat <- heat_df |>
  ggplot(aes(x = class_prob, y = class, fill = mean_prob)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(mean_prob, 2)), size = 2.5) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#d73027",
    midpoint = 1 / 12,
    name = "Prob.\nmédia"
  ) +
  labs(title    = "Probabilidade Média: Classe Predita × Classe Candidata",
       subtitle = "Diagonal = confiança do modelo  |  Fora da diagonal = confusão entre classes",
       x = "Classe candidata", y = "Classe predita") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title  = element_text(face = "bold"))

w <- max(10, 12 * 1.1)
h <- max(7,  n_cls * 0.65)

ggsave(file.path(plots_dir, "heatmap_probabilidades.png"),
       p_heat, width = w, height = h, dpi = 150, bg = "white")

plot(p_heat)
