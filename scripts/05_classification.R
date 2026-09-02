# ============================================================
# Vector Data Cube Classification
# ============================================================

# Load required libraries
library(sits)
library(tibble)
library(ggplot2)
library(terra)
library(RColorBrewer)

# Define the parameters: These are user-defined variables
tiles           <- c('016016','018016')
model_name      <- "tcnn-model_2y_2023-08-01_2025-07-28_2026-08-03_eco-3-mt-47d-vsits2_2026-08-18_00h22m.rds"
seg_version     <- "lsmm-snic-spac10-comp03-pad0-rectangular" #SITS recognizes the underscore (_) character as a separator.
label_method    <- "mean"



# Extract the date of the string separated by "_"
start_date     <- "2024-08-01"
end_date       <- "2026-07-31"

# File and folder paths 
models <- c("rf"   = "random_forest",
            "xgb"  = "xgboost",
            "ltae" = "ltae",
            "tcnn" = "temp_cnn",
            "rnet" = "res_net",
            "lstm" = "ltsm")
model_type   <- stringr::str_split_i(model_name, "-", 1)
model_path   <- file.path("data/rds/model", models[model_type], model_name)
vector_path  <- "data/segments"
class_path   <- "data/class"
plots_path   <- "data/plots/accuracy"
log_path     <- "data/logs"
n_cores      <- 16

dir.create(log_path, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_path, paste0("erros_classificacao_",
                                       format(Sys.time(), "%Y-%m-%d_%Hh%Mm"),
                                       ".txt"))

sits_parallel(workers = n_cores)

# Identifier to distinguish this model run from previous runs
var <- stringr::str_split_i(model_name, "_", 6)

# Cube duration in years (independent of the tile)
no.years <- paste0(floor(lubridate::year(end_date) - lubridate::year(start_date)), "y")

# Load the model once (shared across all tiles)
model <- readRDS(model_path)

# ============================================================
# 1. Helper function: logs errors to the log file and the console
# ============================================================

log_erro <- function(tile, etapa, cond) {
  msg <- sprintf(
    "[%s] TILE: %s | ETAPA: %s | ERRO: %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tile, etapa, conditionMessage(cond)
  )
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
}

# ============================================================
# 2. Loop: a complete classification for each tile
# ============================================================

resultados <- list()
for (tile in tiles) {
  
  cat("\n==============================\n")
  cat("Iniciando tile:", tile, "\n")
  cat("==============================\n")
  
  # tryCatch wraps the entire tile pipeline.
  # If any step (cube, segmentation, classification, or labeling)
  # fails, it goes to the "error =" handler, logs the issue, and the for() loop continues to the next tile.
  resultado_tile <- tryCatch({
    
    # --- 2.1 Raw tile cube ---
    cube <- sits_cube(
      source     = "BDC",
      collection = "SENTINEL-2-16D",
      bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A',
                     'B11', 'B12', 'NDVI', 'NBR', 'EVI', 'CLOUD'),
      tiles      = tile,
      start_date = start_date,
      end_date   = end_date,
      progress   = TRUE
    )
    
    # --- 2.2 Local segmented cube (based on previous segmentation) ---
    local_segs_cube <- sits_cube(
      source      = "BDC",
      collection  = "SENTINEL-2-16D",
      raster_cube = cube,
      vector_dir  = vector_path,
      version     = seg_version,
      parse_info  = c("satellite", "sensor", "tile", "start_date",
                      "end_date", "band", "version", "X1")
    )
    
    # --- 2.3 Tile output directory ---
    tile_period_dir <- file.path(class_path, tile, "original_class",
                                 models[model_type],
                                 paste(no.years, var, sep = "-"))
    dir.create(tile_period_dir, recursive = TRUE, showWarnings = FALSE)
    
    # --- 2.4 Probability file version ---
    version <- paste(model_type, no.years, var, sep = "-")
    
    # --- 2.5 Seed for reproducibility ---
    set.seed(88)
    
    cat("Classifing...")
    # --- 2.6 Classification ---
    class_prob <- sits_classify(
      data       = local_segs_cube,
      ml_model   = model,
      multicores = n_cores,
      memsize    = 120,
      output_dir = tile_period_dir,
      version    = version,
      verbose    = TRUE,
      progress   = TRUE
    )
    gc()
              
    cat("Labeling...")
    # --- 2.7 Final labeling ---
    class_map <- sits_label_classification(
      cube         = class_prob,
      output_dir   = tile_period_dir,
      label_method = label_method,
      version      = paste(version, label_method, sep = "-"),
      multicores   = n_cores,
      memsize      = 120,
      max_cells_in_memory = 01e+09,
      progress     = TRUE
    )
    cat("Tile", tile, "finalizada com sucesso!\n")

    # Clean up large objects before returning the status list,
    # so the list() call stays the last (returned) expression.
    rm(class_prob)
    rm(class_map)
    gc()

    list(status = "ok", tile = tile)

  }, error = function(e) {
    log_erro(tile, "pipeline_completa", e)

    # Only remove objects that actually exist at the point of failure,
    # so a partial pipeline doesn't throw a second error inside the handler.
    if (exists("class_prob")) rm(class_prob)
    if (exists("class_map"))  rm(class_map)
    gc()

    list(status = "erro", tile = tile, mensagem = conditionMessage(e))
  })
  
  # Store the result indexed by tile name, so the summary below can
  # look tiles up by name.
  resultados[[tile]] <- resultado_tile
}

status_vec <- vapply(resultados, function(x) x$status, character(1))
cat("\n================ RESUMO FINAL ================\n")
cat("Tiles processadas com sucesso:",
    paste(names(status_vec[status_vec == "ok"]), collapse = ", "), "\n")
cat("Tiles com erro:",
    paste(names(status_vec[status_vec == "erro"]), collapse = ", "), "\n")
if (any(status_vec == "erro")) {
  cat("Detalhes dos erros disponíveis em:", log_file, "\n")
}
cat("================================================\n")
print("Classificação de múltiplos tiles finalizada!")

# ============================================================
# 3. Uncertainty
# ============================================================

# Step 3.1 -- Define function to calculate entropy, rasterize and exclude .gpkg
compute_uncertainty_raster <- function(
    vector_cube,
    tile_period_dir,
    version,
    multicores = n_cores,
    memsize    = 180,
    delete_gpkg = TRUE
) {
  
  # Calculate uncertainty vector cube
  sits_uncertainty(
    vector_cube,
    type       = "entropy",
    multicores = n_cores,
    memsize    = memsize,
    output_dir = tile_period_dir,
    version    = version,
    progress   = TRUE
  )
  
  # List entropy .gpkg files and get the most recent one
  uncertainty_files <- list.files(
    path      = tile_period_dir,
    pattern   = "entropy.*\\.gpkg$",
    full.names = TRUE
  )
  uncertainty_file <- uncertainty_files[which.max(file.info(uncertainty_files)$mtime)]
  
  # Read the segment polygons file with entropy
  uncertainty_polygons <- terra::vect(uncertainty_file)
  
  # Create a raster template based on uncertainty_polygons
  raster_template <- terra::rast(
    terra::ext(uncertainty_polygons),
    res = terra::res(terra::rast(vector_cube$file_info[[1]]$path[1])),
    crs = terra::crs(uncertainty_polygons)
  )
  
  # Rasterize entropy values and scale to UINT16
  uncertainty_raster       <- terra::rasterize(uncertainty_polygons, raster_template, field = "entropy")
  uncertainty_raster_uint16 <- round(uncertainty_raster * 10000)
  
  # Plot
  plot(
    uncertainty_raster_uint16,
    col     = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100),
    maxcell = terra::ncell(uncertainty_raster_uint16),
    main    = "Uncertainty Map - Full Resolution"
  )
  
  # Save as .tif (UINT16, DEFLATE compressed)
  tif_path <- file.path(
    tile_period_dir,
    paste0(tools::file_path_sans_ext(basename(uncertainty_file)), ".tif")
  )
  
  terra::writeRaster(
    uncertainty_raster_uint16,
    filename  = tif_path,
    datatype  = "INT2U",
    overwrite = TRUE,
    gdal      = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9"),
    progress  = TRUE
  )
  
  # Delete the source .gpkg files (optional)
  if (delete_gpkg) {
    removed <- file.remove(uncertainty_files)
    message("Deleted .gpkg files: ", paste(uncertainty_files[removed], collapse = ", "))
  }
  
  message("Uncertainty raster saved: ", tif_path)
  invisible(tif_path)
}

# Step 3.2 -- Run function to calculate entropy, rasterize and exclude .gpkg
compute_uncertainty_raster(
  vector_cube     = vector_cube,
  tile_period_dir = tile_period_dir,
  version         = version,
  multicores      = n_cores, # adapt to your computer CPU core availability
  memsize         = 180, # adapt to your computer CPU core availability
  delete_gpkg     = TRUE  # Keep the .gpkg file if you want to inspect it beforehand
)

# ============================================================
# 4. Plot Uncertainty Blox-Plots
# ============================================================

# Step 4.1 -- Define function to plot uncertainty by class box-plots
plot_uncertainty_boxplot <- function(
    output_dir,
    band_uncertainty  = "entropy",
    class_column      = "class",
    tile,
    version,
    high_uncert_perc  = 0.90,
    scale_factor      = 10000,
    save_plot         = TRUE,
    plots_path,
    width             = 12,
    height            = 6,
    dpi               = 150,
    verbose           = TRUE
)
{
  
  # ---------------------------------------------------------------------------
  # Helper: robust conversion of class_map to an sf object
  # ---------------------------------------------------------------------------
  read_class_map <- function(x) {
    if (inherits(x, "sf")) return(x)
    if (is.character(x) && file.exists(x)) return(sf::st_read(x, quiet = TRUE))
    stop("class_map must be a file path (.gpkg) or an sf object.")
  }
  
  # ---------------------------------------------------------------------------
  # 0. Find entropy .tif and class .gpkg from output_dir
  # ---------------------------------------------------------------------------
  ent_files <- list.files(output_dir, pattern = paste0(tile, ".*entropy.*", version, "\\.tif$"),
                          full.names = TRUE, ignore.case = TRUE)
  if (length(ent_files) == 0) stop("No entropy .tif found for tile '", tile, "' version '", version, "' in: ", output_dir)
  raster_uncertainty <- ent_files[which.max(file.info(ent_files)$mtime)]
  if (verbose) message("  Using entropy raster: ", basename(raster_uncertainty))
  
  cls_files <- list.files(output_dir, pattern = paste0(tile, ".*class.*", version, "\\.gpkg$"),
                          full.names = TRUE, ignore.case = TRUE)
  if (length(cls_files) == 0) stop("No class .gpkg found for tile '", tile, "' version '", version, "' in: ", output_dir)
  class_map <- cls_files[which.max(file.info(cls_files)$mtime)]
  if (verbose) message("  Using class map: ", basename(class_map))
  
  # ---------------------------------------------------------------------------
  # 1. Read entropy raster
  # ---------------------------------------------------------------------------
  if (verbose) message("[1/3] Reading entropy raster...")
  r_ent <- terra::rast(raster_uncertainty)
  idx   <- grep(band_uncertainty, names(r_ent), ignore.case = TRUE)
  if (length(idx) == 0)
    stop("Band '", band_uncertainty, "' not found. Available bands: ",
         paste(names(r_ent), collapse = ", "))
  r_ent        <- r_ent[[idx[1]]]
  names(r_ent) <- "entropy"
  
  # ---------------------------------------------------------------------------
  # 2. Read classification map using the robust helper
  # ---------------------------------------------------------------------------
  if (verbose) message("[2/3] Reading classification map...")
  sf_cls <- read_class_map(class_map)
  
  if (verbose) message("  Class column: '", class_column, "'")
  
  sf_cls <- sf_cls |>
    dplyr::rename(class = dplyr::all_of(class_column)) |>
    dplyr::mutate(class = as.character(class))
  
  if (verbose) {
    unique_classes <- sort(unique(sf_cls$class))
    message("  Classes: ", paste(unique_classes, collapse = ", "))
    message("  Total segments: ", nrow(sf_cls))
  }
  
  # ---------------------------------------------------------------------------
  # 3. Extraction of entropy at segment centroids
  # ---------------------------------------------------------------------------
  if (verbose) message("[3/3] Extracting entropy at segment centroids...")
  centroids <- sf::st_centroid(sf_cls)
  ent_vals  <- terra::extract(
    r_ent,
    terra::vect(centroids)
  )
  sf_cls$entropy_raw <- ent_vals$entropy
  
  sf_cls <- sf_cls[!is.na(sf_cls$entropy_raw), ]
  if (verbose) message("  Segments with valid entropy: ", nrow(sf_cls))
  
  # ---------------------------------------------------------------------------
  # 4. Scale entropy values if necessary
  # ---------------------------------------------------------------------------
  max_val <- max(sf_cls$entropy_raw, na.rm = TRUE)
  if (max_val > 10 && scale_factor != 1) {
    if (verbose) message("  Entropy values appear scaled (max = ", max_val,
                         "), dividing by ", scale_factor, "...")
    sf_cls$entropy <- sf_cls$entropy_raw / scale_factor
  } else {
    sf_cls$entropy <- sf_cls$entropy_raw
  }
  
  # ---------------------------------------------------------------------------
  # 5. Prepare data for plotting
  # ---------------------------------------------------------------------------
  pts <- sf_cls |>
    sf::st_drop_geometry() |>
    dplyr::select(class, entropy) |>
    dplyr::mutate(class = as.factor(class))
  
  threshold <- stats::quantile(pts$entropy, high_uncert_perc, na.rm = TRUE)
  
  n_cls <- nlevels(pts$class)
  colors <- scales::hue_pal()(n_cls)
  
  p <- pts |>
    dplyr::mutate(
      class = forcats::fct_reorder(class, entropy, .fun = stats::median, .desc = TRUE)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = class, y = entropy, fill = class)) +
    ggplot2::geom_violin(alpha = 0.35, color = NA) +
    ggplot2::geom_boxplot(width = 0.3, outlier.size = 0.4, outlier.alpha = 0.3) +
    ggplot2::geom_hline(
      yintercept = threshold, linetype = "dashed",
      color = "red", linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = threshold, hjust = 1.05, vjust = -0.5,
      label = paste0("P", high_uncert_perc * 100),
      color = "red", size = 3.2
    ) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::labs(
      title    = "Entropy Distribution by Class",
      subtitle = "Red dashed line = high uncertainty threshold",
      x        = NULL,
      y        = "Entropy distribution"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      plot.title  = ggplot2::element_text(face = "bold")
    )
  
  # ---------------------------------------------------------------------------
  # 5. Save plot
  # ---------------------------------------------------------------------------
  if (save_plot) {
    dir.create(plots_path, showWarnings = FALSE, recursive = TRUE)
    file_name <- paste0("bxplt-uncertainty-by-class_", tile, "_", version, ".png")
    file_path <- file.path(plots_path, file_name)
    ggplot2::ggsave(file_path, p, width = width, height = height, dpi = dpi, bg = "white")
    if (verbose) message("Plot saved to: ", normalizePath(file_path))
  }
  print(p)
  invisible(p)
}

# Step 5.1 -- Run function
plot_uncertainty_boxplot(
  output_dir        = tile_period_dir,
  tile              = tile,
  version           = version,
  high_uncert_perc  = 0.90,
  plots_path        = plots_path,
  verbose           = TRUE
)
