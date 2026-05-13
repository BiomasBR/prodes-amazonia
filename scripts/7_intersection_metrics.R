# Load required libraries
library(sf)
library(segmetric)
library(furrr)
library(future)
library(parallelly)
library(openxlsx)

# ============================================================
# Define parameters and paths
# ============================================================

tile <- "012014"
version <- "camara-agrupado-reclassificado"

prodes_dir  <- "data/raw/samples/prodes-2025"
class_dir   <- paste0("data/class/",tile,"/post_processed")

# List the classified polygons 
class <- list.files(
  class_dir,
  pattern    = paste0("^class-post-processed_", tile, "_[^_]+_[^_]+_", version, "\\.gpkg$"),
  full.names = TRUE,
  recursive  = TRUE
)

# List the PRODES reference polygons 
ref_prodes <- list.files(
  prodes_dir,
  pattern    = paste0(".*", tile, ".*\\.gpkg$"),
  full.names = TRUE,
  recursive  = TRUE
)

# ============================================================
# Define intersection metrics function
# ============================================================

compute_metrics <- function(ref_prodes,
                                class,
                                n_workers = max(1L, floor(parallelly::availableCores() / 1.5))) {
  
  # --- Configure parallelization ----------------------------
  plan(multisession, workers = n_workers)
  on.exit(plan(sequential), add = TRUE)
  message(sprintf("[PARALLEL] Using %d workers.", n_workers))
  
  # --- Helper: secure reading with correction via GEOS --------
  safe_read_sf <- function(path, layer_name) {
    sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(TRUE), add = TRUE)
    
    layer <- sf::read_sf(path)
    
    # Corrige topologia via GEOS (funciona em lon/lat com s2 desligado)
    layer <- sf::st_make_valid(layer)
    
    # Remove geometrias ainda inválidas ou vazias
    valid_mask <- sf::st_is_valid(layer) & !sf::st_is_empty(layer)
    n_dropped  <- sum(!valid_mask)
    if (n_dropped > 0)
      warning(sprintf("[%s] %d geometria(s) inválida(s) removida(s) após correção.",
                      layer_name, n_dropped))
    
    layer[valid_mask, ]
  }
 
  # --- Load PRODES -------------------------------------
  if (length(ref_prodes) > 1)
    warning(sprintf("[PRODES] %d files found; using the first one.", length(ref_prodes)))
  
  pol_ref <- safe_read_sf(ref_prodes[1], "PRODES")
  
  # --- Load SITS classification -------------------------
  if (length(class) > 1)
    warning(sprintf("[SITS] %d files found; using the first one.", length(class)))
  
  pol_class <- safe_read_sf(class[1], "SITS")
  
  # --- Align CRS -----------------------------------------
  if (!sf::st_crs(pol_ref) == sf::st_crs(pol_class)) {
    message("[CRS] Reprojecting pol_class for the CRS of pol_ref.")
    pol_class <- sf::st_transform(pol_class, sf::st_crs(pol_ref))
  }
  
  # --- Segmetric Object -----------------------
  seg_obj <- segmetric::sm_read(ref_sf = pol_ref, seg_sf = pol_class)
  
  # --- Calculate metrics in parallel (segmetric) -----------
  metricas_ids <- c(
    "OS2", "US2", "AFI", "D_index",
    "precision", "recall", "M", "IoU", "Dice"
  )
  
  metricas_list <- furrr::future_map(
    metricas_ids,
    .f        = \(m) list(name = m, result = segmetric::sm_compute(seg_obj, m)),
    .options  = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  )
  
  # --- Summarize metrics in data.frame (structure: 1 line per metric) ----
  r_df <- do.call(rbind, lapply(metricas_list, function(x) {
    
    # Extracts only the numeric values from the segmetric object 
    vals <- unlist(x$result[sapply(x$result, is.numeric)])
    vals <- vals[is.finite(vals)]   # remove NA/Inf
    
    data.frame(
      Minimo        = min(vals,            na.rm = TRUE),
      Maximo        = max(vals,            na.rm = TRUE),
      Media         = mean(vals,           na.rm = TRUE),
      Mediana       = median(vals,         na.rm = TRUE),
      Desvio_Padrao = sd(vals,             na.rm = TRUE),
      row.names     = x$name
    )
  }))
  
  # --- Calculate the area of intersection -------------------------
  
  # Total area of PRODES e total area of classification
  area_ref_ha   <- sum(as.numeric(sf::st_area(pol_ref)),   na.rm = TRUE) / 10000
  area_class_ha <- sum(as.numeric(sf::st_area(pol_class)), na.rm = TRUE) / 10000
  
  intersec_sf  <- sf::st_intersection(pol_ref, pol_class)
  areas_m2     <- as.numeric(sf::st_area(intersec_sf))
  areas_ha     <- areas_m2 / 10000
  total_area_ha <- sum(areas_ha, na.rm = TRUE)
  
  # % da interseção em relação ao PRODES  → recall geométrico
  pct_ref   <- (total_area_ha / area_ref_ha)   * 100
  
  # % da interseção em relação à classificação → precisão geométrica
  pct_class <- (total_area_ha / area_class_ha) * 100
  
  area_df <- data.frame(
    polygon_id      = seq_along(areas_ha),
    area_intersecao_ha = areas_ha
  )
  
  summary_area_df <- data.frame(
    descricao = c(
      "Area PRODES (ha)",
      "Area classificacao (ha)",
      "Area intersecao total (ha)",
      "Intersecao / PRODES (%)",
      "Intersecao / Classificacao (%)"
    ),
    valor = c(
      area_ref_ha,
      area_class_ha,
      total_area_ha,
      pct_ref,
      pct_class
    )
  )
  
  message(sprintf("Total area of PRODES:  %.4f ha", area_ref_ha))
  message(sprintf("Total area of classification: %.4f ha", area_class_ha))
  message(sprintf("Total area of intersection: %.4f ha", total_area_ha))
  message(sprintf("intersection / PRODES:  %.2f %%", pct_ref))
  message(sprintf("intersection / Class.:  %.2f %%", pct_class))

  # --- Export to Excel ---------------------------------
  class_filename <- tools::file_path_sans_ext(basename(class[1]))
  metrics_filename <- sub("^class", "metrics", class_filename)
  
  wb <- openxlsx::createWorkbook()
  
  sheet_name <- tools::file_path_sans_ext(basename(ref_prodes[1]))
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, r_df, rowNames = TRUE)
  
  openxlsx::addWorksheet(wb, "Polygons_Intersection")
  openxlsx::writeData(wb, "Polygons_Intersection", area_df)
  
  openxlsx::addWorksheet(wb, "Summary_Intersection")
  openxlsx::writeData(wb, "Summary_Intersection", summary_area_df)
  
  out_path <- file.path(
    dirname(class[1]),
    paste0(metrics_filename, ".xlsx")
  )
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[OK] File saved in: %s", out_path))
  
  invisible(list(
    metricas        = r_df,
    area_intersecao = area_df,
    total_area_ha      = total_area_ha,
    out_path        = out_path
  ))
}

# ============================================================
# Use intersection metrics function
# ============================================================

metrics_results <- compute_metrics(
  ref_prodes = ref_prodes,
  class   = class
)
