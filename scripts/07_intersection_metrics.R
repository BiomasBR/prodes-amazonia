# Load required libraries
library(sf)
library(segmetric)
library(openxlsx)
library(dplyr)

# ============================================================
# Define parameters and paths
# ============================================================

tile <- "012014"
version <- "rf-1y-all-samples-new-pol-avg-false-pixel-mode"

prodes_dir  <- "data/raw/prodes-2025"
class_dir   <- file.path("data/class", tile, "post_processed")

# List the classified polygons 
class_files <- list.files(
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

# Verify that at least one file was found
if (length(ref_prodes) == 0) stop("No PRODES reference file found.")
if (length(class_files) == 0) stop("No classification file found.")

# ============================================================
# Define intersection metrics function (with polygon disaggregation)
# ============================================================

compute_metrics <- function(ref_prodes, class_files, overwrite = TRUE) {
  
  old_s2_state <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE) # Disable spherical geometry (S2), use planar (GEOS)
  on.exit(sf::sf_use_s2(old_s2_state), add = TRUE)
  
  # --- Helper: secure reading with correction (only for classification) --------
  safe_read_sf <- function(path, layer_name, extract_polygons = FALSE) {
    layer <- sf::read_sf(path)
    layer <- sf::st_make_valid(layer)
    
    if (extract_polygons) {
      # Only extract if necessary, avoiding warning "x is already of type POLYGON"
      geom_types <- unique(as.character(sf::st_geometry_type(layer)))
      if (any(grepl("MULTI|COLLECTION", geom_types))) {
        layer <- sf::st_collection_extract(layer, "POLYGON")
      }
    }
    
    valid_mask <- sf::st_is_valid(layer) & !sf::st_is_empty(layer)
    n_dropped  <- sum(!valid_mask)
    if (n_dropped > 0) {
      warning(sprintf("[%s] %d invalid geometry(ies) removed after correction.",
                      layer_name, n_dropped))
      layer <- layer[valid_mask, ]
    }
    return(layer)
  }
  
  # --- Load PRODES (Reference) - keep original geometry types initially ----
  if (length(ref_prodes) > 1)
    warning(sprintf("[PRODES] %d files found; using the first one.", length(ref_prodes)))
  pol_ref <- safe_read_sf(ref_prodes[1], "PRODES", extract_polygons = FALSE)
  
  # --- NEW STEP: Disaggregate MULTIPOLYGON into individual POLYGONs ----------
  # This ensures that each reference object is a single polygon,
  # allowing meaningful per-polygon metrics.
  geom_types_ref <- unique(as.character(sf::st_geometry_type(pol_ref)))
  if (any(grepl("MULTI|COLLECTION", geom_types_ref))) {
    pol_ref <- sf::st_cast(pol_ref, "POLYGON")
    pol_ref <- sf::st_make_valid(pol_ref)
    message(sprintf("PRODES disaggregated into %d individual polygons.", nrow(pol_ref)))
  } else {
    message(sprintf("PRODES already contains %d simple polygons.", nrow(pol_ref)))
  }
  
  # --- Load SITS classification (extract polygons) ---
  if (length(class_files) > 1)
    warning(sprintf("[SITS] %d files found; using the first one.", length(class_files)))
  pol_class <- safe_read_sf(class_files[1], "SITS", extract_polygons = TRUE)
  
  # --- Align CRS -----------------------------------------
  if (!sf::st_crs(pol_ref) == sf::st_crs(pol_class)) {
    message("[CRS] Reprojecting pol_class to match CRS of pol_ref.")
    pol_class <- sf::st_transform(pol_class, sf::st_crs(pol_ref))
  }
  
  # --- Project to a equal-area coordinate system for correct area calculation ----
  equal_area_crs <- "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs"
  
  pol_ref_proj <- sf::st_transform(pol_ref, equal_area_crs)
  pol_class_proj <- sf::st_transform(pol_class, equal_area_crs)
  
  # --- Segmetric Object (uses original CRS) ---
  seg_obj <- segmetric::sm_read(ref_sf = pol_ref, seg_sf = pol_class)
  
  # --- Calculate metrics ---
  metricas <- seg_obj %>%
    segmetric::sm_compute("OS2") %>%
    segmetric::sm_compute("US2") %>%
    segmetric::sm_compute("AFI") %>%
    segmetric::sm_compute("D_index") %>%
    segmetric::sm_compute("precision") %>%
    segmetric::sm_compute("recall") %>%
    segmetric::sm_compute("M") %>%
    segmetric::sm_compute("IoU") %>%
    segmetric::sm_compute("Dice")
  
  # --- Summarize only the intended metrics ---
  wanted_metrics <- c("OS2", "US2", "AFI", "D_index", "precision", 
                      "recall", "M", "IoU", "Dice")
  r_df <- data.frame()
  
  metrics_evaluation <- function(x) {
    data.frame(
      Minimo        = round(min(x,    na.rm = TRUE), 4),
      Maximo        = round(max(x,    na.rm = TRUE), 4),
      Media         = round(mean(x,   na.rm = TRUE), 4),
      Mediana       = round(median(x, na.rm = TRUE), 4),
      Desvio_Padrao = round(ifelse(length(x) > 1, sd(x, na.rm = TRUE), NA), 4)
    )
  }
  
  for (m in wanted_metrics) {
    if (!m %in% names(metricas)) next
    vals <- unlist(metricas[[m]])
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) next
    r <- metrics_evaluation(vals)
    row.names(r) <- m
    r_df <- rbind(r_df, r)
  }
  
  # --- Calculate areas (in hectares) using projected data ---
  area_ref_ha   <- sum(as.numeric(sf::st_area(pol_ref_proj)),   na.rm = TRUE) / 10000
  area_class_ha <- sum(as.numeric(sf::st_area(pol_class_proj)), na.rm = TRUE) / 10000
  
  intersec_sf  <- sf::st_intersection(pol_ref_proj, pol_class_proj)
  areas_m2     <- as.numeric(sf::st_area(intersec_sf))
  areas_ha     <- areas_m2 / 10000
  total_area_ha <- sum(areas_ha, na.rm = TRUE)
  
  pct_ref   <- (total_area_ha / area_ref_ha)   * 100
  pct_class <- (total_area_ha / area_class_ha) * 100
  
  area_df <- data.frame(
    polygon_id      = seq_along(areas_ha),
    area_intersecao_ha = areas_ha
  )
  
  summary_area_df <- data.frame(
    descricao = c(
      "PRODES area (ha)",
      "Total Classification Area (ha)",
      "Total Intersection Area (ha)",
      "Intersection / PRODES (%)",
      "Intersection / Classification (%)"
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
  message(sprintf("Intersection / PRODES:  %.2f %%", pct_ref))
  message(sprintf("Intersection / Class.:  %.2f %%", pct_class))
  
  # --- Export to Excel ---------------------------------
  class_filename <- tools::file_path_sans_ext(basename(class_files[1]))
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
    dirname(class_files[1]),
    paste0(metrics_filename, ".xlsx")
  )
  
  if (!overwrite && file.exists(out_path)) {
    stop(sprintf("File already exists and overwrite = FALSE: %s", out_path))
  }
  
  openxlsx::saveWorkbook(wb, out_path, overwrite = overwrite)
  message(sprintf("[OK] File saved in: %s", out_path))
  
  invisible(list(
    metricas        = r_df,
    area_intersecao = area_df,
    total_area_ha   = total_area_ha,
    out_path        = out_path
  ))
}

# ============================================================
# Use intersection metrics function
# ============================================================
metrics_results <- compute_metrics(
  ref_prodes  = ref_prodes,
  class_files = class_files,
  overwrite   = TRUE)
