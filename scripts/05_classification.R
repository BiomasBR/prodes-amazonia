# ============================================================
#  Classification of Vector Data Cube
# ============================================================

# Load required libraries
library(sits)
library(tibble)
library(ggplot2)
library(terra)
library(RColorBrewer)

# Define the parameters: These are user-defined variables
tiles           <- c('012014', '012015', '013014', '013015')
model_name      <- "tcnn-model_2y_2023-08-01_2025-07-13_2026-08-03_eco-3-mt-46d_2026-08-03_16h02m.rds"
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
n_cores      <- 28

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
      vector_band = "segments",
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
      memsize    = 180,
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
      memsize      = 180,
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
