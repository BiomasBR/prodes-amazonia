# ============================================================
#  Samples quality assessment, filtering and balancing with SOM
# ============================================================

# ============================================================
# 1. Libraries, paths and some initial parameters
# ============================================================

# Step 1.1 -- Load Required Libraries
library(sits)
library(tibble)
library(dplyr)
library(ggplot2)

# Step 1.2 -- Define the date and time for the start of processing
process_version <- paste0(format(Sys.Date(), "%Y-%m-%d_"), format(Sys.time(), "%Hh%Mm", tz = "America/Sao_Paulo"))

# Step 1.3 -- Define the names and paths for files and folders needed in the processing
rds_path      <- "data/rds"
rds_filename  <- "TS-tiles_012014-012015-013014-013015_1y_2024-08-01_2025-07-31_after-apocalypse-agrupado_2026-06-29_10h48m.rds"
config_dir    <- "../scripts"
plots_dir    <- "data/plots/som"

# Step 1.4 -- Identifier to distinguish this model run from previous versions
no.years <- strsplit(rds_filename, "_")[[1]][3]

# Step 1.5 -- Identifier to distinguish this model run from previous versions
var <- "after-apocalypse-agrupado"

# Step 1.6 -- Create output plots directory per var
plots_path <- file.path(plots_dir, paste0(var,"-",no.years))
dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)

# Step 1.7 -- Extracts the list of tiles
tiles <- gsub(".*_(\\d{6}(-\\d{6})*)_.*", "\\1", rds_filename)


# ============================================================
# 2. Load and Translate Training Samples Dataset
# ============================================================

# Step 2.1 -- Load the samples Time Series from a R file
samples <- readRDS(file.path(rds_path, "time_series", rds_filename))

# Step 2.2 -- Function to read class names and their colors::IMPORTANT
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

# Step 2.3 -- Load color palette and classes translation from external config file
config    <- read_class_config("~/grupos/biomasbr-amazonia/sits-prodes/class_config.txt")
class_translation <- config$class_translation
my_colors <- config$my_colors

# Step 2.4 -- Translate the colors names (English to Portuguese)
nomes_em_pt <- names(class_translation)[match(names(my_colors), class_translation)]
names(my_colors) <- nomes_em_pt

# Step 2.5 -- Filter features
my_colors <- my_colors[names(my_colors) %in% unique(samples$label)]


# ============================================================
# 3. Self-Organizing Map (SOM - 1)
# ============================================================

# Step 3.1 -- Find out the ideal grid size for your original samples data
# Run with a 2x2 grid, then observes the interval indicated by SITS and fill the values in the next step
sits_som_map(samples,
  grid_xdim = 2,
  grid_ydim = 2,
  alpha = 1.0, # starting learning rate (decreases according to number of iterations)
  distance = "dtw", # type of similarity measure (distance)
  rlen = 20 # number of iterations to produce the SOM
)

# Step 3.2 -- Define function SOM-1: clustering SOM, plot/save SOM map, evaluate quality, plot/save (mixed labels)
som_analysis_1 <- function(samples,
                             grid_xdim = 25,
                             grid_ydim = 25,
                             alpha = 1.0,
                             distance = "dtw",
                             rlen = 20,
                             tiles,
                             var,
                             no.years,
                             process_version,
                             plots_path,
                             rds_path,
                             my_colors,
                             plot_width = 3529,
                             plot_height = 1578,
                             plot_dpi = 350) {
  
  # Clustering original Time Series Samples using SOM
  sits_som_map_start <- Sys.time()
  
  som_cluster <- sits_som_map(
    samples,
    grid_xdim = grid_xdim,
    grid_ydim = grid_ydim,
    alpha = alpha,          # starting learning rate (decreases according to number of iterations)
    distance = distance,    # type of similarity measure (distance)
    rlen = rlen             # number of iterations to produce the SOM
  )

  # Plot the SOM map
  plot(som_cluster)
  
  # Save SOM map plot
  ggsave(
    filename = paste0("SOM1-evaluation-cluster_", tiles, "_", var,"-", no.years, "_", process_version, ".png"),
    path = plots_path,
    scale = 1,
    width = plot_width,
    height = plot_height,
    units = "px",
    dpi = plot_dpi
  )
  
  # Save the samples Time Series (SOM object) to a R file
  saveRDS(
    som_cluster,
    paste0(
      rds_path, "/som/",
      "SOM1-", som_cluster$som$grid$xdim, "x", som_cluster$som$grid$ydim, "_",
      tiles, "_", var,"-", no.years, "_", process_version, ".rds"
    )
  )
  
  # SOM evaluation (clustered labels) + plot + save
  som_eval <- sits_som_evaluate_cluster(som_cluster)
  
  p <- plot(som_eval) +
    labs(title = 'Confusion by cluster') +
    xlab("Percentage of Mixture") +
    ylab(NULL) +
    scale_fill_manual(values = my_colors, name = "Legend") +
    theme(legend.position = "right")
  
  print(p)
  
  ggsave(
    filename = paste0("SOM1-confusion-cluster_", tiles, "_", var,"-", no.years, "_", process_version, ".png"),
    plot = p,
    path = plots_path,
    scale = 1,
    width = plot_width,
    height = plot_height,
    units = "px",
    dpi = plot_dpi
  )
  
  sits_som_map_end <- Sys.time()
  sits_som_map_time <- as.numeric(sits_som_map_end - sits_som_map_start, units = "secs")
  
  cat(sprintf(
    "SITS SOM map process duration (HH:MM): %02d:%02d\n",
    as.integer(sits_som_map_time / 3600),
    as.integer((sits_som_map_time %% 3600) / 60)
  ))
  
  return(list(
    som_cluster = som_cluster,
    som_eval = som_eval
  ))
}

# Step 3.3 -- Run function SOM-1
som1 <- som_analysis_1(
  samples          = samples,
  grid_xdim        = 25,
  grid_ydim        = 25,
  alpha            = 1.0,
  distance         = "dtw",
  rlen             = 20,
  tiles            = tiles,
  var              = var,
  no.years         = no.years,
  process_version  = process_version,
  plots_path       = plots_path,
  rds_path         = rds_path,
  my_colors        = my_colors
)

# Step 3.4 -- Returns 'som_cluster' and 'som_eval' objects from SOM analysis 1
som_cluster <- som1$som_cluster
som_eval    <- som1$som_eval


# ============================================================
# 4. Analyse Quality (and Clean) Samples
# ============================================================

# Step 4.1 -- Define function to Evaluates the quality of the samples based on the results of the SOM map 1
analyse_quality <- function(som_cluster, tiles, var, no.years, process_version, plots_path,
                                       prior_threshold = 0.10, 
                                       posterior_threshold = 0.10, 
                                       keep = c("clean", "analyze", "remove")) {
  
  # Evaluates the quality of the samples based on the results of the SOM map 1
  original_samples_quality <- sits_som_clean_samples(
    som_map = som_cluster, 
    prior_threshold = prior_threshold,
    posterior_threshold = posterior_threshold,
    keep = keep
  )
  
  eval_col <- "eval"

  # Ensures that all categories exist before plotting
  categorias_presentes <- unique(original_samples_quality[[eval_col]])
  categorias_faltando  <- setdiff(c("clean", "analyze", "remove"), categorias_presentes)
  
  if (length(categorias_faltando) > 0) {
    message(
      "Notice: the categories [", paste(categorias_faltando, collapse = ", "),
      "] did not occur with the current thresholds (prior=", prior_threshold,
      ", posterior=", posterior_threshold, "). ",
      "Generating graph manually (bypass of plot.som_clean_samples from sits), ",
      "which breaks when a category has zero samples)."
    )
  
  df <- original_samples_quality
  df[[eval_col]] <- factor(df[[eval_col]], levels = c("clean", "analyze", "remove"))
  
  quality_summary <- df |>
    dplyr::count(label, .data[[eval_col]], .drop = FALSE, name = "n") |>
    dplyr::group_by(label) |>
    dplyr::mutate(percentage = 100 * n / sum(n)) |>
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(
    quality_summary,
    ggplot2::aes(x = percentage, y = label, fill = .data[[eval_col]])
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(
      values = c(clean = "#2ca02c", analyze = "#ff7f0e", remove = "#d62728"),
      drop = FALSE   # mantém "analyze" na legenda mesmo com 0%
    ) +
    ggplot2::labs(title = "SOM evaluation - samples quality", x = "Percentage (%)", y = NULL, fill = "Category") +
    ggplot2::theme_minimal()
  
  print(p)
  
  } else {
  p <- plot(original_samples_quality)
  print(p)
  }
  
  ggsave(
    filename = paste0("original-samples-quality_", tiles, "_", var,"-", no.years, "_", process_version, ".png"),
    plot = p,
    path = plots_path,
    scale = 1,
    width = 3529,
    height = 1578,
    units = "px",
    dpi = 350
    )
  
  return(original_samples_quality)
  }

# Step 4.2 -- Run function to Evaluates the quality of the samples
samples_quality <- analyse_quality(
  som_cluster     = som_cluster,
  tiles           = tiles,
  var             = var,
  no.years        = no.years,
  process_version = process_version,
  plots_path      = plots_path
)

# Step 4.3 -- Filter samples according to an evaluation class ("clean", "analyze", "remove") or to a specific class with low samples quantity
cleaned_samples <- original_samples_quality %>% filter(eval == "clean")

# Step 4.4 -- Save the cleaned_samples Time Series to a .rds file
saveRDS(cleaned_samples, 
        paste0(rds_path, "/time_series/", 
               "TS-cleaned-samples_", tiles, "_", var,"-", no.years, "_", process_version, ".rds"))


# ============================================================
# 5. Self-Organizing Map (SOM - 2)
# ============================================================

# Step 5.1 -- Find out the ideal grid size for your cleaned samples data
# Run with a 2x2 grid, then observes the interval indicated by SITS and fill the values in the next step
sits_som_map(cleaned_samples,
             grid_xdim = 2,
             grid_ydim = 2,
             alpha = 1.0, # starting learning rate (decreases according to number of iterations)
             distance = "dtw", # type of similarity measure (distance)
             rlen = 20 # number of iterations to produce the SOM
)

# Step 5.2 -- Define function SOM-2: clustering SOM, plot/save SOM map, evaluate quality, plot/save (mixed labels)
som_analysis_2 <- function(cleaned_samples,
                           grid_xdim = 2,
                           grid_ydim = 2,
                           alpha = 1.0,
                           distance = "dtw",
                           rlen = 20,
                           tiles,
                           var,
                           no.years,
                           process_version,
                           plots_path,
                           rds_path,
                           my_colors,
                           plot_width = 3529,
                           plot_height = 1578,
                           plot_dpi = 350) {
  
  sits_som_map_start2 <- Sys.time()
  
  # Clustering cleaned Time Series Samples
  som_cluster_cleaned <- sits_som_map(
    data       = cleaned_samples,
    grid_xdim  = grid_xdim,
    grid_ydim  = grid_ydim,
    alpha      = alpha,
    distance   = distance,
    rlen       = rlen
    # som_radius = 2,
    # mode = "online" # only for windows' PCs
  )
  
  # Plot the SOM map 2
  plot(som_cluster_cleaned)
  
  # Save the plot
  ggsave(
    filename = paste0("SOM2-evaluation-samples-cleaned_", tiles, "_", var,"-", no.years, "_", process_version, ".png"),
    path = plots_path,
    scale = 1,
    width = plot_width,
    height = plot_height,
    units = "px",
    dpi = plot_dpi
  )
  
  # Save the cleaned Time Series Samples to a R file
  saveRDS(
    som_cluster_cleaned,
    paste0(
      rds_path, "/som/",
      "SOM2-", som_cluster_cleaned$som$grid$xdim, "x", som_cluster_cleaned$som$grid$ydim, "_",
      tiles, "_", var,"-", no.years, "_", process_version, ".rds"
    )
  )
  
  # -- Produce a tibble with a summary of the mixed labels
  som2_eval_clean <- sits_som_evaluate_cluster(som_cluster_cleaned)
  
  # -- Plot the result of summary of the mixed labels
  p <- plot(som2_eval_clean) +
    labs(title = 'Confusion by cluster') +
    xlab("Percentage of Mixture") +
    ylab(NULL) +
    scale_fill_manual(values = my_colors, name = "Legend") +
    theme(legend.position = "right")
  
  print(p)
  
  # -- Save the plot of summary of the mixed labels
  ggsave(
    filename = paste0("SOM2-confusion-cluster-cleaned_", tiles, "_", var,"-", no.years, "_", process_version, ".png"),
    plot = p,
    path = plots_path,
    scale = 1,
    width = plot_width,
    height = plot_height,
    units = "px",
    dpi = plot_dpi
  )
  
  sits_som_map_end2 <- Sys.time()
  sits_som_map_time2 <- as.numeric(sits_som_map_end2 - sits_som_map_start2, units = "secs")
  
  cat(sprintf(
    "SITS SOM map 2 process duration (HH:MM): %02d:%02d\n",
    as.integer(sits_som_map_time2 / 3600),
    as.integer((sits_som_map_time2 %% 3600) / 60)
  ))
  
  return(list(
    som_cluster_cleaned = som_cluster_cleaned,
    som2_eval_clean   = som2_eval_clean,
    duration_secs     = sits_som_map_time2
  ))
}

# Step 5.3 -- Run function SOM-2
som2 <- som_analysis_2(
  cleaned_samples  = cleaned_samples,
  grid_xdim        = 2,
  grid_ydim        = 2,
  alpha            = 1.0,
  distance         = "dtw",
  rlen             = 20,
  tiles            = tiles,
  var              = var,
  no.years         = no.years,
  process_version  = process_version,
  plots_path       = plots_path,
  rds_path         = rds_path,
  my_colors        = my_colors
)

# Step 5.4 -- Returns 'som_cluster_cleaned' and 'som2_eval_clean' objects from SOM analysis 2
som_cluster_cleaned <- som2$som_cluster_clean
som2_eval_clean    <- som2$som2_eval_clean


# ============================================================
# 6. Reduce imbalance among classes
# ============================================================

# Step 6.1 -- Find out the total values by classes
table(samples$label)

# Step 6.2 -- Define function to reduce imbalance among the classes
reduce_imbalance <- function(samples = samples,
                             n_samples_over = 3000,
                             n_samples_under = 3000,
                             tiles,
                             var,
                             no.years,
                             process_version,
                             rds_path) {
  
  sits_reduce_imbalance_start <- Sys.time()
  
  clean_samples_balanced <- sits_reduce_imbalance(
    samples          = samples,
    n_samples_over   = n_samples_over,
    n_samples_under  = n_samples_under
  )
  
  # Removing columns that contain NA values
  clean_samples_balanced <- clean_samples_balanced[, colSums(is.na(clean_samples_balanced)) == 0]
  
  # Save the new Time Series Samples Balanced to a R file
  saveRDS(
    clean_samples_balanced,
    paste0(rds_path, "/time_series/", "TS-balanced", "_", tiles, "_", var, "-", no.years, "_", process_version, ".rds")
  )
  
  sits_reduce_imbalance_end <- Sys.time()
  sits_reduce_imbalance_time <- as.numeric(sits_reduce_imbalance_end - sits_reduce_imbalance_start, units = "secs")
  
  cat(sprintf(
    "SITS reduce imbalance process duration (HH:MM): %02d:%02d\n",
    as.integer(sits_reduce_imbalance_time / 3600),
    as.integer((sits_reduce_imbalance_time %% 3600) / 60)
  ))
  
  return(list(clean_samples_balanced = clean_samples_balanced))
}

# Step 6.3 -- Run the function to reduce imbalance among the classes
balancing <- reduce_imbalance(
  samples          = samples,
  n_samples_over   = 3000,
  n_samples_under  = 3000,
  tiles            = tiles,
  var              = var,
  no.years         = no.years,
  process_version  = process_version,
  rds_path         = rds_path
)

# Step 6.4 -- Returns 'balanced_samples' objects from function "reduce_imbalance"
balanced_samples <- balancing$clean_samples_balanced

# ============================================================
# 7. Self-Organizing Map (SOM - 3)
# ============================================================

# Step 7.1 -- Find out the ideal grid size for your balanced samples data
# Run with a 2x2 grid, then observes the interval indicated by SITS and fill the values in the next step
sits_som_map(balanced_samples,
             grid_xdim = 2,
             grid_ydim = 2,
             alpha = 1.0, # starting learning rate (decreases according to number of iterations)
             distance = "dtw", # type of similarity measure (distance)
             rlen = 20 # number of iterations to produce the SOM
)

# Step 7.2 -- Define function SOM-3: clustering SOM, plot/save SOM map, evaluate quality, plot/save
som_analysis_3 <- function(balanced_samples,
                           grid_xdim = 25,
                           grid_ydim = 25,
                           alpha = 1.0,
                           distance = "dtw",
                           rlen = 20,
                           tiles,
                           var,
                           no.years,
                           process_version,
                           plots_path,
                           my_colors,
                           plot_width = 3529,
                           plot_height = 1578,
                           plot_dpi = 350) {
  
  sits_som_map_start3 <- Sys.time()
  
  # Clustering SOM
  som_cluster_clean_balanced <- sits_som_map(
    balanced_samples,
    grid_xdim = grid_xdim,
    grid_ydim = grid_ydim,
    alpha     = alpha,
    distance  = distance,
    rlen      = rlen
  )
  
  saveRDS(
    som_cluster_clean_balanced,
    paste0(
      rds_path, "/som/",
      "SOM3-", som_cluster_clean_balanced$som$grid$xdim, "x", som_cluster_clean_balanced$som$grid$ydim, "_",
      tiles, "_", var,"-", no.years, "_", process_version, ".rds"
    )
  )
  
  # Plot the SOM map
  plot(som_cluster_clean_balanced)
  
  # -- Save the plot
  ggsave(
    filename = paste0("SOM3-evaluation-cluster-clean-balanced_", tiles, "_", var, "-", no.years, "_", process_version, ".png"),
    path = plots_path,
    scale = 1,
    width = plot_width,
    height = plot_height,
    units = "px",
    dpi = plot_dpi
  )

    # -- Show the summary of the balanced time series sample data
  print(summary(balanced_samples))
    
  # -- Produce a tibble with a summary of the mixed labels
  som_eval_clean_balanced <- sits_som_evaluate_cluster(som_cluster_clean_balanced)

  # -- Plot the result
  p <- plot(som_eval_clean_balanced) +
    labs(title = 'Confusion by cluster') +
    xlab("Percentage of Mixture") +
    ylab(NULL) +
    scale_fill_manual(values = my_colors, name = "Legend") +
    theme(legend.position = "right")
  
  print(p)
  
  # -- Save the plot
  ggsave(
    filename = paste0("SOM3-confusion-cluster-clean-balanced_", tiles, "_", var, "-", no.years, "_", process_version, ".png"),
    plot = p,
    path = plots_path,
    scale = 1,
    width = plot_width,
    height = plot_height,
    units = "px",
    dpi = plot_dpi
  )
  
  sits_som_map_end3 <- Sys.time()
  sits_som_map_time3 <- as.numeric(sits_som_map_end3 - sits_som_map_start3, units = "secs")
  
  cat(sprintf(
    "SITS SOM map 3 process duration (HH:MM): %02d:%02d\n",
    as.integer(sits_som_map_time3 / 3600),
    as.integer((sits_som_map_time3 %% 3600) / 60)
  ))
  
  return(list(
    som_cluster_clean_balanced = som_cluster_clean_balanced,
    som_eval_clean_balanced    = som_eval_clean_balanced
  ))
}

# Step 7.3 -- Run function SOM-3
som3 <- som_analysis_3(
  balanced_samples        = balanced_samples,
  grid_xdim               = 22,
  grid_ydim               = 22,
  alpha                   = 1.0,
  distance                = "dtw",
  rlen                    = 20,
  tiles                   = tiles,
  var                     = var,
  no.years                = no.years,
  process_version         = process_version,
  plots_path              = plots_path,
  my_colors               = my_colors
)

# Step 7.4 -- Returns 'som_cluster_clean_balanced' and 'som3_eval_clean_balanced' objects from SOM analysis 2
som_cluster_clean_balanced <- som3$som_cluster_clean_balanced
som3_eval_clean_balanced    <- som3$som_eval_clean_balanced
