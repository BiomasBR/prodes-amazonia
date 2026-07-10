# ============================================================
#  Stratification of validation samples
# ============================================================

# Load required libraries
library(tibble)
library(sits)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(fs)
library(stringr)
library(purrr)

# Define the parameters: These are user-defined variables
model_name      <- "rf-model_4t_012014-012015-013014-013015_2y_2023-08-01_2025-07-31_after-apocalypse-agrupado_2026-06-29_11h47m.rds"
tiles           <- c("012014")
version         <- "rf-2y-after-apocalypse-agrupado-mean"

# define and load model path
models <- c("rf"   = "random_forest",
            "xgb"  = "xgboost",
            "ltae" = "ltae",
            "tcnn" = "temp_cnn",
            "rnet" = "res_net",
            "lstm" = "ltsm")
model_type       <- stringr::str_split_i(model_name, "-", 1)
model_path       <- file.path("data/rds/model", models[model_type], model_name)
model            <- readRDS(model_path)

# define classification, mask and output paths
class_dir       <- "data/class"
output_dir      <- "data/raw/samples/validation_samples"
mask_dir        <- "data/raw/auxiliary/masks"
prodes_dir      <- "data/raw/prodes-2025"

# define segments path
pattern          <- sprintf(".*_(%s)_", paste(tiles, collapse = "|"))
seg_path         <- list.files("data/segments",
                               pattern = pattern,
                               full.names = TRUE)

# List the PRODES reference polygons 
ref_prodes <- list.files(
  prodes_dir,
  pattern    = paste0(".*", tiles, ".*\\.gpkg$"),
  full.names = TRUE,
  recursive  = TRUE
)

# Define function to create validation
sits_validation_sampling <- function(
    cube            = NULL,
    sampling_design = NULL,
    validation_type = NULL,
    alloc         = NULL,
    overhead      = 1.2,
    progress      = TRUE,
    multicores    = 8,
    polygons      = polygons,
    prodes        = read_sf(ref_prodes),
    output_dir,
    version,
    date_process  = format(Sys.Date(), "%Y-%m-%d")
) {
  
  # ── 0. Basic validations ────────────────────────────────────────────────────
  stopifnot(is.character(alloc))
  stopifnot(is.numeric(overhead) && overhead >= 1)
  fs::dir_create(output_dir, recurse = TRUE)
  
  tile_id           <- paste(cube$tile, collapse = "-")
  overhead_col_name <- paste0("overhead_", overhead)
  file_suffix       <- paste(validation_type, tile_id, version, date_process,
                             sep = "_")
  
  # ── 1. Generate stratified samples ────────────────────────────────────────
  cli::cli_inform("Generating stratified samples ({alloc}, overhead = {overhead})...")
  
  samples_sf <- sits::sits_stratified_sampling(
    cube            = cube,
    sampling_design = sampling_design,
    alloc           = alloc,
    overhead        = overhead,
    progress        = progress,
    multicores      = multicores
  ) |>
    dplyr::rename(sits_label = label) |>
    dplyr::mutate(validation_label = NA_character_)
  
  # ── 2. Save points (gpkg) ──────────────────────────────────────────────────
  points_path <- file.path(
    output_dir,
    paste0("validation-samples_points_", file_suffix, ".gpkg")
  )
  
  cli::cli_inform("Saving points: {fs::path_file(points_path)}")
  sf::st_write(samples_sf, points_path, delete_dsn = TRUE, append = FALSE)
  
  # ── 3. Save polygons (gpkg) ───────────────────────────────────────────────
  polygons_path <- NULL
  
  if (!is.null(polygons)) {
    cli::cli_inform("Filtering and saving polygons...")
    
    validation_polygons <- polygons |>
      sf::st_join(
        samples_sf |>
          sf::st_transform(sf::st_crs(polygons)) |>
          dplyr::select(sits_label),
        join = sf::st_intersects,
        left = FALSE
      )
    
    # ── Análise de Interseção com o PRODES ────────────────────────────────────
    if (!is.null(prodes)) {
      cli::cli_inform("Calculating PRODES intersection percentages...")
      
      # Garantir que o PRODES use a mesma projeção dos polígonos
      prodes_crs <- sf::st_transform(prodes, sf::st_crs(validation_polygons))
      
      if (nrow(prodes_crs) > 0) {
        # Corrigir geometrias inválidas e converter em partes simples (POLYGON)
        prodes_single <- prodes_crs |> 
          sf::st_make_valid() |> 
          sf::st_cast("POLYGON")
        
        # Criar IDs temporários e calcular a área original de cada polígono de validação
        validation_polygons$tmp_id    <- 1:nrow(validation_polygons)
        validation_polygons$poly_area <- as.numeric(sf::st_area(validation_polygons))
        
        # Realizar a interseção geométrica detalhada
        suppressWarnings({
          intersections <- sf::st_intersection(
            validation_polygons[, c("tmp_id", "poly_area")], 
            prodes_single
          )
        })
        
        if (nrow(intersections) > 0) {
          # Calcular a área de cada fragmento intersectado
          intersections$int_area <- as.numeric(sf::st_area(intersections))
          
          # Somar as áreas de interseção caso um polígono toque em mais de uma parte simples do PRODES
          prodes_cover <- intersections |>
            sf::st_drop_geometry() |>
            dplyr::group_by(tmp_id) |>
            dplyr::summarise(total_int_area = sum(int_area), .groups = "drop")
          
          # Cruzar de volta os dados e calcular o percentual final na coluna "prodes"
          validation_polygons <- validation_polygons |>
            dplyr::left_join(prodes_cover, by = "tmp_id") |>
            dplyr::mutate(
              prodes = (tidyr::replace_na(total_int_area, 0) / poly_area) * 100,
              prodes = round(prodes, 2) # Limita a duas casas decimais
            ) |>
            dplyr::select(-tmp_id, -poly_area, -total_int_area)
        } else {
          # Se não houve interseção real de área
          validation_polygons$prodes <- 0
          validation_polygons <- validation_polygons |> dplyr::select(-tmp_id, -poly_area)
        }
      } else {
        # Se o objeto PRODES fornecido estiver vazio
        validation_polygons$prodes <- 0
      }
    }
    
    polygons_path <- file.path(
      output_dir,
      paste0("validation-samples_polygons_", file_suffix, ".gpkg")
    )
    
    sf::st_write(validation_polygons, polygons_path, delete_dsn = TRUE, append = FALSE)
  }
  
  # ── 4. Save Excel spreadsheet with sampling design summary ────────────────────────────
  cli::cli_inform("Generating a summary of the sampling design...")
  
  # Actual count of samples generated by class
  samples_count <- samples_sf |>
    sf::st_drop_geometry() |>
    dplyr::group_by(sits_label) |>
    dplyr::summarise(total_generated = dplyr::n(), .groups = "drop")
  
  design_df <- as.data.frame(sampling_design) |>
    tibble::rownames_to_column("class") |>
    dplyr::select(class, prop, expected_ua, std_dev,
                  alloc_chosen = dplyr::all_of(alloc)) |>
    dplyr::rename(!!alloc := alloc_chosen) |>
    dplyr::mutate(
      alloc_val            = as.integer(.data[[alloc]]),
      !!overhead_col_name := ceiling(alloc_val * overhead) - alloc_val
    ) |>
    dplyr::select(-alloc_val) |>
    dplyr::left_join(
      samples_count,
      by = dplyr::join_by(class == sits_label)
    )
  
  xlsx_path <- file.path(
    output_dir,
    paste0("sampling-design_", file_suffix, ".xlsx")
  )
  
  # Style: highlighted header
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "sampling_design")
  
  header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF",
    fgFill     = "#2E75B6",
    halign     = "CENTER",
    textDecoration = "bold",
    border     = "Bottom"
  )
  
  openxlsx::writeData(wb, sheet = 1, x = design_df, headerStyle = header_style)
  openxlsx::setColWidths(wb, sheet = 1, cols = 1:ncol(design_df), widths = "auto")
  openxlsx::freezePane(wb, sheet = 1, firstRow = TRUE)
  openxlsx::saveWorkbook(wb, file = xlsx_path, overwrite = TRUE)
  
  cli::cli_inform("Excel file saved in: {fs::path_file(xlsx_path)}")
  
  # ── 5. Returns paths and objects generated  ───────────────────────────────────
  invisible(list(
    samples_sf    = samples_sf,
    design_df     = design_df,
    points_path   = points_path,
    polygons_path = polygons_path,
    xlsx_path     = xlsx_path
  ))
}

# ============================================================
# 1. Create raster file from classified vector map
# ============================================================

# Step 1.1 -- Function to rasterize
sits_rasterize_segments <- function(file, res, class_raster_dir, style = NULL) {
  
  stopifnot(!is.null(res))
  stopifnot(!is.null(file))
  stopifnot(!is.null(class_raster_dir))
  
  fs::dir_create(class_raster_dir, recurse = TRUE)
  
  file <- fs::path_expand(file)
  class_raster_dir <- fs::path_expand(class_raster_dir)
  
  output_file_base <- fs::path(class_raster_dir) / fs::path_file(file) |>
    fs::path_ext_remove()
  
  output_file <- stringr::str_c(output_file_base, ".tif")
  output_style <- stringr::str_c(output_file_base, ".qml")
  
  if (fs::file_exists(output_file)) {
    return(output_file)
  }
  
  file_sf <- sf::st_read(file, quiet = TRUE)
  
  if (is.null(style)) {
    file_sf <- file_sf |>
      dplyr::mutate(
        class_num = .data[["class"]] |>
          as.factor() |>
          as.numeric()
      )
    
    style <- file_sf |>
      tibble::as_tibble() |>
      dplyr::select(dplyr::all_of(c("class", "class_num"))) |>
      dplyr::distinct(.data[["class"]], .data[["class_num"]]) |>
      dplyr::mutate(color = RColorBrewer::brewer.pal(n = dplyr::n(), name = "Set3")) |>
      dplyr::rename(
        index = class_num,
        name = class
      ) |>
      dplyr::arrange(.data[["index"]])
    
  } else {
    
    file_sf <- file_sf |>
      dplyr::rename(name = class) |>
      dplyr::left_join(
        style |> dplyr::select(dplyr::all_of(c("name", "index")))
      ) |>
      dplyr::mutate(
        class_num = .data[["index"]]
      )
  }
  
  file_bbox <- sf::st_bbox(file_sf)
  
  file_gpkg <- fs::file_temp(ext = ".gpkg")
  
  sf::st_write(obj = file_sf, dsn = file_gpkg, quiet = TRUE)
  
  a_srs <- readRDS(
    url("https://github.com/restore-plus/restore-utils/raw/refs/heads/main/inst/extdata/crs/bdc.rds")
  )
  
  gdalUtilities::gdal_rasterize(
    src_datasource = file_gpkg,
    dst_filename = output_file,
    a = "class_num",
    at = TRUE,
    tr = c(res, res),
    te = file_bbox,
    ot = "Int16",
    init = 255,
    a_nodata = 255,
    co = c(
      "COMPRESS=ZSTD",
      "PREDICTOR=2",
      "ZSTD_LEVEL=1",
      "BIGTIFF=YES"
    ),
    a_srs = a_srs
  )
  
  sits:::.colors_qml(
    color_table = style,
    file = output_style
  )
  return(output_file)
}

# Step 1.2 -- Style from ML model
style <- tibble::tibble(
  name = sits_labels(model),
  index = 1:length(sits_labels(model)),
  color = pals::cols25(length(sits_labels(model)))
)

# Step 1.3 -- Rasterize classified vectors
to_raster <- paste0(".*_class_", version, ".*\\.gpkg$")
class_files <- list.files(
  path = class_dir,
  pattern = to_raster,
  full.names = TRUE,
  recursive = TRUE
)
raster_files <- purrr::map(class_files, function(file) {
  file_name <- fs::path_file(file)
  cli::cli_inform("Processing: {file_name}")
  tile_id <- stringr::str_extract(file_name, "\\d{6}")
  tile_period_dir <- file.path(
    class_dir,
    tile_id,
    "accuracy"
  )
  
  fs::dir_create(tile_period_dir, recurse = TRUE)
  sits_rasterize_segments(
    file = file,
    res = 10,
    style = style,
    class_raster_dir = tile_period_dir
  )
})

# ============================================================
# 2. SITS Cube and Segments
# ============================================================

# Step 2.1 -- Get labels associated to the trained model data set (Enumerate them in the order they appear according to "sits_labels(model)")
pattern <- paste0(".*", tiles, ".*", version, ".*\\.tif$")
cube_dirs <- grep("accuracy",
                  list.dirs(class_dir, recursive = TRUE) |> 
                    purrr::keep(~ length(list.files(.x, pattern = pattern)) > 0),
                  value = TRUE)

# Step 2.2 -- Create an array with the labels for the model classes
labels <- c(
  x = sits_labels(model)
)
names(labels) <- 1:length(labels)

# Step 2.3 -- Load the original cube with classified raster file
cube_list <- map(cube_dirs, function(dir) {
  sits_cube(
    source      = "BDC",
    collection  = "SENTINEL-2-16D",
    bands       = "class",
    labels      = labels,
    data_dir    = dir, # Takes one path from 'cube_dirs' at a time
    version     = version,
    parse_info  = c("satellite", "sensor", "tile", "start_date", "end_date", 
                    "band", "version")
  )
})

# Step 2.4 -- Combine the list of tibbles into a single multi-row sits cube
cube <- do.call(rbind, cube_list)

# Step 2.5 -- Iterates over the paths of the segments files, then stacks them all into a single spatial data frame
polygons <- seg_path |> 
  map(read_sf) |> 
  bind_rows()

# ============================================================
# 3. Full Map Stratified Random Sampling
# ============================================================

# Step 3.1 -- Full Map Sampling design
sampling_design <- sits_sampling_design(
  cube = cube,
  expected_ua = c(
    "Corpo_Dagua"                           = 0.95,
    "Corte_Raso_Com_Herbacea"               = 0.70,
    "Corte_Raso_Com_Solo_Exposto"           = 0.70,
    "Degradacao"                            = 0.70,
    "Floresta"                              = 0.95,
    "Vegetacao_Natural_Nao_Florestal"       = 0.70
  ),
  alloc_options = c(120, 100, 75, 50, 30),
  std_err = 0.01,
  rare_class_prop = 0.015
)

# Step 3.2 -- Show Full Map sampling design
sampling_design

# Step 3.3 -- Run function to create validation samples all classes
result_all_classes <- sits_validation_sampling(
  cube            = cube,
  sampling_design = sampling_design,
  validation_type = "all-classes",
  alloc           = "alloc_100",
  overhead        = 1.2,
  progress        = TRUE,
  multicores      = 12,
  polygons        = polygons,
  prodes          = sf::read_sf(ref_prodes),
  output_dir     = output_dir,
  version         = version,
  date_process    = format(Sys.Date(), "%Y-%m-%d")
)

# ============================================================
# 4. Grouped Adjusted Map Accuracy
# ============================================================

# Step 4.1 -- Create a data cube of type mask
counter_mask <- c("1" = "Natural Vegetation", "0" = "Deforestation Mask")
prodes_mask <- sits_cube(source = "BDC",
                         collection = "SENTINEL-2-16D",
                         tiles = cube$tile,
                         data_dir = mask_dir,
                         parse_info = c("X1", "tile", "start_date",
                                        "end_date", "band", "version"),
                         bands = "class",
                         version = "contra-mask-geral-amz",
                         labels = counter_mask)

# Step 4.2 -- Define and create a directory to store the regrouped (reclassified) cube
dir_path <- file.path(
  class_dir,
  cube$tile,
  "grouped"
)
fs::dir_create(dir_path, recurse = TRUE)

# Step 4.3 -- Reclassifies the original class cube into groups
cube_reclass <- sits_reclassify(
        cube = cube,
        mask = prodes_mask,
        rules = list(
          "Deforestation" =
            cube %in% c(
              "Corte_Raso_Com_Arvores_Remanescentes",
              "Corte_Raso",
              "Corte_Raso_Com_Vegetacao"
            ),
          "Degradation" =
            cube %in% c(
              "Degradacao",
              "Degradacao_Por_Fogo"
            ),
          "Other_Classes" =
            cube %in% c(
              "Corpo_Dagua",
              "Corte_Raso_Antigo",
              "Corte_Raso_Antigo_Com_Vegetacao",
              "Floresta",
              "Floresta_Transicional",
              "Vegetacao_Natural_Nao_Florestal",
              "Area_Inundavel"
            )
        ),
        multicores = 24,
        memsize = 180,
        version = paste("grouped", version, sep = "-"),
        output_dir = dir_path,
        progress = TRUE
  )

# 4.4 -- Sampling design degradation
sampling_design_grouped <- sits_sampling_design(
  cube = cube_reclass,
  expected_ua = c(
    "Deforestation" = 0.70,
    "Degradation" = 0.60, 
    "Other_Classes" = 0.95
  ),
  alloc_options = c(120, 100),
  std_err = 0.01,
  rare_class_prop = 0.05
)

# 4.5 -- Show sampling design
sampling_design_grouped

# 4.6 -- Run function to create validation samples grouped
result_grouped <- sits_validation_sampling(
  cube            = cube_reclass,
  sampling_design = sampling_design_grouped,
  validation_type = "grouped",
  alloc           = "alloc_100",
  overhead        = 1.2,
  progress        = TRUE,
  multicores      = 12,
  polygons        = polygons,
  prodes          = sf::read_sf(ref_prodes),
  output_dir      = output_dir,
  version         = version,
  date_process    = format(Sys.Date(), "%Y-%m-%d")
)
