# ============================================================
# Forest Post-processing
# ============================================================

# Load required libraries
library(sf)
library(dplyr)
library(sits)
library(terra)
library(units)
library(smoothr)
library(purrr)
library(stringr)
library(RANN)

# Define the parameters: These are user-defined variables
model_name  <- "tcnn-model_2y_2023-08-01_2025-07-28_2026-08-03_eco-3-mt-47d-vsits2_2026-08-18_00h22m.rds"
version     <- "tcnn-2y-eco-3-mt-47d-vsits2-mean"
tiles       <- c('020016')

# File and folder paths
seg_version <- "lsmm-snic-spac10-comp03-pad0-rectangular"
class_path  <- "data/class"
mask_path   <- "data/raw/auxiliary/mascara_geral_amz_v2025_postgis_nb.gpkg"
config_dir  <- ".."

# Brazil Albers Equal Area (EPSG 10857)
crs_proc <- "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"Unknown based on GRS80 ellipsoid\",\n            ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,\n                LENGTHUNIT[\"metre\",1],\n                ID[\"EPSG\",7019]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Albers Equal Area\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",-12,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-54,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",-2,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",-22,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",5000000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",10000000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]"
# SIRGAS 2000
crs_final <- 4674

# Precision in metric CRS
precision <- units::set_units(1, "mm")

models <- c("rf"   = "random_forest",
            "xgb"  = "xgboost",
            "ltae" = "ltae",
            "tcnn" = "temp_cnn",
            "rnet" = "res_net",
            "lstm" = "ltsm")
model_type <- stringr::str_split_i(model_name, "-", 1)
model_path <- file.path("data/rds/model", models[model_type], model_name)
model      <- readRDS(model_path)
years <- regmatches(version, regexpr("\\d+y", version))

# Biome boundary (shared by all tiles, loaded only once)
biome <- read_sf("data/raw/auxiliary/amazon-biome-border-epsg10857.gpkg") |>
  st_make_valid() |>
  st_transform(crs_proc)

edge_tiles <- c(
  "001014", "002011", "002012", "002013", "002014", "002015", "002016",
  "003011", "003015", "003016", "004010", "004011", "004016", "005004",
  "005005", "005006", "005007", "005009", "005010", "005016", "005017",
  "006004", "006005", "006006", "006007", "006008", "006009", "006017",
  "007003", "007004", "007017", "008003", "008004", "008005", "008017",
  "009005", "009016", "009017", "010001", "010004", "010005", "010016",
  "010017", "010018", "011001", "011002", "011003", "011004", "011017",
  "011018", "011019", "012001", "012002", "012003", "012019", "013001",
  "013002", "013019", "014001", "014019", "014020", "015000", "015001",
  "015019", "015020", "015021", "016000", "016001", "016002", "016003",
  "016004", "016018", "016019", "016020", "016021", "016022", "016023",
  "017004", "017018", "017021", "017022", "017023", "018003", "018004",
  "018018", "018019", "018020", "018021", "018022", "018023", "019003",
  "019019", "019020", "019021", "019022", "020003", "020018", "020019",
  "020020", "021003", "021019", "021020", "022003", "022020", "023003",
  "023019", "023020", "024001", "024002", "024003", "024017", "024018",
  "024019", "024020", "025001", "025002", "025003", "025016", "025017",
  "025018", "025019", "026003", "026004", "026005", "026013", "026014",
  "026015", "026016", "027005", "027006", "027013", "027014", "027015",
  "028006", "028011", "028012", "028013", "028014", "028015", "029006",
  "029011", "030006", "030007", "030011", "031007", "031010", "031011",
  "032007", "032008", "032009", "032010", "033008", "033009"
)

# ============================================================
# 3. Helper functions
# ============================================================

# Extracting the cloud mask
extract_cloud_mask <- function(
    sits_classification_path,
    sits_reclassification,
    cloud_values = c(3, 8, 9, 10),
    output_dir = NULL,
    collection = "SENTINEL-2-16D",
    date_window_days = 1
) {
  
  # ----------------------------------------------------------
  # 1. Extract metadata from the filename
  # ----------------------------------------------------------
  filename_base <- basename(sits_classification_path)
  
  last_date_str <- regmatches(
    filename_base,
    gregexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", filename_base)
  )[[1]]
  
  if (length(last_date_str) == 0) {
    stop("No date in YYYY-MM-DD format found in file name: ", filename_base)
  }
  
  last_date <- as.Date(tail(last_date_str, 1))
  
  start_date_scl <- format(last_date - date_window_days, "%Y-%m-%d")
  end_date_scl   <- format(last_date, "%Y-%m-%d")
  
  tile_id <- regmatches(
    filename_base,
    regexpr("(?<=SENTINEL-2_MSI_)[0-9]+", filename_base, perl = TRUE)
  )
  
  if (length(tile_id) == 0 || tile_id == "") {
    stop("Could not extract tile_id from file name: ", filename_base)
  }
  
  message(" -> SCL time window: ", start_date_scl, " to ", end_date_scl)
  
  # ----------------------------------------------------------
  # 1.1 Check if the output file already exists
  # ----------------------------------------------------------
  if (!is.null(output_dir)) {
    output_filename <- paste0("cloud-vec_", tile_id, "_", end_date_scl, ".gpkg")
    output_path <- file.path(output_dir, output_filename)
    
    if (file.exists(output_path)) {
      message(" -> File already exists: ", output_filename, ". Skipping processing.")
      
      cloud_vec <- sf::st_read(output_path, quiet = TRUE)
      
      return(invisible(list(
        cloud_vec    = cloud_vec,
        tile_id      = tile_id,
        end_date_scl = end_date_scl
      )))
    }
  }
  
  # ----------------------------------------------------------
  # 2. Build the BDC cube
  # ----------------------------------------------------------
  scl_cube <- sits::sits_cube(
    source = "BDC",
    collection = collection,
    tiles = tile_id,
    bands = "SCL",
    start_date = start_date_scl,
    end_date = end_date_scl
  )
  
  # ----------------------------------------------------------
  # 3. Load the SCL raster
  # ----------------------------------------------------------
  scl_files <- scl_cube$file_info[[1]] |>
    dplyr::filter(
      band == "CLOUD",
      date == as.Date(end_date_scl)
    ) |>
    dplyr::pull(path)
  
  if (length(scl_files) == 0) {
    stop("No SCL file found for date: ", end_date_scl)
  }
  
  scl_raster <- terra::rast(scl_files[1])
  
  message(" -> SCL file loaded: ", scl_files[1])
  
  # ----------------------------------------------------------
  # 4. Create the cloud mask
  # ----------------------------------------------------------
  scl_mask <- terra::classify(
    scl_raster,
    rcl = cbind(cloud_values, rep(1, length(cloud_values))),
    others = NA
  )
  
  # ---------------------------------------------------------------------------
  # 4.1 Check if any cloud value was found
  # ---------------------------------------------------------------------------
  if (all(is.na(terra::values(scl_mask)))) {
    message("  -> No clouds identified")
    return(invisible(list(
      cloud_vec    = NULL,
      tile_id      = tile_id,
      end_date_scl = end_date_scl
    )))
  }
  
  # ----------------------------------------------------------
  # 5. Clip to the classification extent
  # ----------------------------------------------------------
  class_bbox <- sits_reclassification |>
    sf::st_transform(terra::crs(scl_mask)) |>
    terra::vect() |>
    terra::ext()
  
  scl_raster_crop <- terra::crop(scl_mask, class_bbox)
  
  # ----------------------------------------------------------
  # 6. Vectorize the cloud mask
  # ----------------------------------------------------------
  cloud_vec <- terra::as.polygons(scl_raster_crop, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_transform(sf::st_crs(sits_reclassification)) |>
    smoothr::fill_holes(threshold = Inf)
  
  # ----------------------------------------------------------
  # 7. Save 
  # ----------------------------------------------------------
  if (!is.null(output_dir)) {
    cloud_vec |>
      sf::st_transform(4674) |>
      sf::st_collection_extract("POLYGON") |>
      sf::st_cast("POLYGON") |>
      sf::st_write(output_path, append = FALSE)
    
    message(" -> Cloud vector saved (EPSG:4674): ", output_path)
  }
  
  # ----------------------------------------------------------
  # 8. Return
  # ----------------------------------------------------------
  return(
    invisible(
      list(
        cloud_vec   = cloud_vec,
        tile_id     = tile_id,
        end_date_scl = end_date_scl
      )
    )
  )
}

# Cloud/shadow difference
remove_cloud_areas <- function(
    sits_reclassification,
    cloud_vec,
    buffer_dist = 100
) {
  # Check if cloud_vec exists and contains features
  if (is.null(cloud_vec) || nrow(cloud_vec) == 0) {
    message("  -> No cloud vectors were found")
    return(invisible(sits_reclassification))
  }
  
  # Ensure the same CRS before any geometric operation
  if (st_crs(cloud_vec) != st_crs(sits_reclassification)) {
    cloud_vec <- st_transform(cloud_vec, st_crs(sits_reclassification))
  }
  
  # Dissolve
  cloud_union <- sf::st_union(sf::st_make_valid(cloud_vec))
  
  # Buffer
  cloud_vec_buffer <- sf::st_buffer(cloud_union, dist = buffer_dist)
  
  # Remove cloud/shadow areas from the classification
  sits_classification_cloud_cleaned <- sf::st_difference(
    sf::st_make_valid(sits_reclassification),
    cloud_vec_buffer
  ) |>
    sf::st_cast("MULTIPOLYGON")
  
  return(invisible(sits_classification_cloud_cleaned))
}

# Calculate area, perimeter, shared boundaries and equivalent radius
calculate_edge_metrics <- function(class, prodes_mask, crs_planar) {
  
  # Preserves the original state of S2 and ensures restoration upon completion of execution
  s2_state <- sf_use_s2()
  on.exit(sf_use_s2(s2_state), add = TRUE)
  
  # Assigns a unique temporary ID for control purposes
  class$id_feicao <- seq_len(nrow(class))
  
  # Geometric Metrics
  area_vec <- as.numeric(st_area(class))
  perim_vec <- as.numeric(st_length(st_boundary(class)))
  
  class$area <- area_vec
  class$perimetro_total <- perim_vec
  class$raio_equivalente <- 2 * (area_vec / perim_vec)
  class$raio_equivalente[!is.finite(class$raio_equivalente)] <- 0
  
  # Shared Edge Calculation
  class_linhas <- st_cast(class, "MULTILINESTRING")
  
  sf_use_s2(FALSE)
  borda_compartilhada <- st_intersection(class_linhas, prodes_mask)
  borda_compartilhada$comp_compartilhado <- as.numeric(st_length(borda_compartilhada))
  
  # Grouping of segments by feature
  borda_resumo <- borda_compartilhada |>
    st_drop_geometry() |>
    group_by(id_feicao) |>
    summarise(comp_compartilhado = sum(comp_compartilhado), .groups = "drop")
  
  # Combines the results and calculates the final proportion
  class <- class |>
    left_join(borda_resumo, by = "id_feicao") |>
    mutate(
      comp_compartilhado = coalesce(comp_compartilhado, 0),
      prop_comp = comp_compartilhado / perimetro_total,
      prop_comp = ifelse(!is.finite(prop_comp), 0, prop_comp)
    )
  
  # Removes the temporary ID column
  class$id_feicao <- NULL
  
  return(class)
}

# Chop polygons
chop_polygons <- function(pol, class, mask, dist){
  
  buf_neg <- st_buffer(
    pol,
    dist = dist,
    joinStyle = "MITRE",
    mitreLimit = 2
  )
  
  # Remove null or invalid geometries that may appear
  buf_neg <- buf_neg[!st_is_empty(buf_neg), ]
  buf_neg <- st_make_valid(buf_neg) |>
    st_collection_extract("POLYGON") |>
    st_cast("POLYGON")
  # ------------------------------------------------------------
  # 2. Distance-based allocation (competitive growth)
  # ------------------------------------------------------------
  # Convert to SpatVector
  orig_v <- vect(pol)
  buf_v  <- vect(buf_neg)
  
  # Unique field in the seeds
  buf_v$id_seed <- 1:nrow(buf_v)
  
  # Empty raster covering the original extent
  r_template <- rast(class)
  
  # Rasterize the seeds
  sementes <- rasterize(
    buf_v,
    r_template,
    field = "id_seed",
    background = NA
  )
  
  seed_cells <- which(!is.na(values(sementes)))
  
  # Coordinates of these cells
  seed_centroid <- xyFromCell(
    sementes,
    seed_cells
  )
  
  # ID corresponding to each seed cell
  seed_ids <- values(sementes)[seed_cells]
  
  # Generate the centroids of the template raster
  xy_pontos <- xyFromCell(
    r_template,
    1:ncell(r_template)
  )
  
  # Nearest-neighbor interpolation
  nn <- RANN::nn2(
    data = seed_centroid,
    query = xy_pontos,
    k = 1
  )
  # ---------------------------------------------------------
  # 6. Raster
  # ---------------------------------------------------------
  valores <- seed_ids[nn$nn.idx[, 1]]
  values(r_template) <- valores
  
  aloc_final <- mask(r_template, orig_v)
  
  # Convert raster to vector polygons
  allocated_polygons <- disagg(as.polygons(aloc_final, aggregate=TRUE))
  
  allocated_polygons$area_ha <- expanse(allocated_polygons, unit = "ha")
  
  large  <- allocated_polygons[allocated_polygons$area_ha >= 1, ]
  small <- disagg(
    aggregate(
      allocated_polygons[allocated_polygons$area_ha <  1, ]
    )
  )
  
  combined1 <- combineGeoms(
    x        = large,
    y        = small,
    overlap  = FALSE,
    boundary = TRUE,
    distance = FALSE,
    dissolve = TRUE,
    erase    = TRUE,
    append   = TRUE       # includes geometries of y that do not match
  )
  
  combined <- disagg(combined1)
  
  combined$area_ha2 <- expanse(combined, unit = "ha")
  
  precision <- units::set_units(1, "mm")
  
  combined <-  combined |>
    st_as_sf() |> 
    st_cast("MULTIPOLYGON") |> 
    st_cast("POLYGON") |>
    sf::st_set_precision(precision) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON")
  
  combined$touches_mask <- lengths(
    st_intersects(combined, mask)
  ) > 0
  
  combined <- combined |>
    dplyr::filter(area_ha2 >= 1 | touches_mask == TRUE)
  
  return(combined)
}

# Calculate and display the elapsed time
log_step_time <- function(step_name, start_time) {
  elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
  message("--> [Processing Time ", step_name, "]: ", elapsed, " seconds")
}

# ============================================================
# 4. Main function: process ONE tile
# ============================================================

process_tile <- function(tile) {
  
  t_total_start <- Sys.time()
  
  message("\n======================================")
  message("Starting tile post-processing: ", tile)
  message("========================================")
  
  # ----------------------------------------------------------
  # 1. Reading Classification File
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 1 of 10 -> Reading classification file.")
  
  raw_class_path <- list.files(
    class_path,
    pattern = paste0(".*_", tile, "_.*_class_", version, "\\.tif$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(raw_class_path) == 0) {
    stop("No classification raster found for the tile ", tile)
  }
  
  if (length(raw_class_path) > 1) {
    stop(
      "More than one classification raster found for the tile ", tile, ":\n",
      paste(" -", raw_class_path, collapse = "\n"),
      "\nAdjust the search pattern (or remove duplicate files) so that only 1 remains."
    )
  }
  
  post_class_path <- file.path(class_path, tile, "post_processed", version)
  dir.create(post_class_path, showWarnings = FALSE, recursive = TRUE)
  
  raw_class <- rast(raw_class_path)
  levels(raw_class) <- data.frame(
    ID = seq_along(sits_labels(model)),
    classe = sits_labels(model)
  )
  
  crs_proc <- crs(raw_class)
  log_step_time("Step 1", t_step)
  
  # ----------------------------------------------------------
  # 2. Classification Classes and Vectorization
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 2 of 10 -> Vectorizing deforestation classes.")
  
  labels <- c('Clear_Cut', 'Clear_Cut_Herbaceous', 'Mininig')
  labels_ids <- match(labels, sits_labels(model))
  
  if (anyNA(labels_ids)) {
    stop(
      "The following labels were not found in sits_labels(model): ",
      paste(labels[is.na(labels_ids)], collapse = ", "),
      ". Labels available in the model: ",
      paste(sits_labels(model), collapse = ", ")
    )
  }
  
  deforest_class <- ifel(
    raw_class %in% labels_ids,
    raw_class,
    NA
  ) |>
    categories(value = data.frame(
      ID = labels_ids,
      classe = labels
    ))
  
  vector_class <- as.polygons(deforest_class, aggregate = TRUE)
  names(vector_class) <- "class"
  vector_multipolygons <- aggregate(vector_class, by = "class")
  vector_multipolygons <- sf::st_as_sf(vector_multipolygons) |>
    sf::st_make_valid()
  
  rm(deforest_class, vector_class)
  gc()
  log_step_time("Step 2", t_step)
  
  # ----------------------------------------------------------
  # 3. Remove polygons outside the biome border
  # ----------------------------------------------------------
  t_step <- Sys.time()
  if (tile %in% edge_tiles) {
    message("Step 3 of 10 -> The tile ", tile, " is an edge tile. Running intersection.")
    class_biome <- st_intersection(vector_multipolygons, biome)
  } else {
    message("Step 3 of 10 -> The tile ", tile, " is not an edge tile. Intersection ignored.")
    class_biome <- vector_multipolygons
  }
  rm(vector_multipolygons)
  gc()
  
  # ----------------------------------------------------------
  # 4. Extraction of cloud features
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 4 of 10 -> Analyzing cloud cover.")
  
  result <- extract_cloud_mask(
    sits_classification_path = raw_class_path,
    sits_reclassification    = class_biome,
    cloud_values               = c(3, 8, 9, 10),
    output_dir                = post_class_path
  )
  
  cloud_vec    <- result$cloud_vec
  end_date_scl <- result$end_date_scl
  log_step_time("Step 4", t_step)
  
  # ----------------------------------------------------------
  # 5. Cloud/shadow difference
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 5 of 10 -> Removing classification in cloud areas (if exist).")
  
  sits_classification_cloud_cleaned <- remove_cloud_areas(
    sits_reclassification = class_biome,
    cloud_vec             = cloud_vec,
    buffer_dist           = 100
  )
  
  rm(result, cloud_vec, class_biome)
  gc()
  log_step_time("Step 5", t_step)
  
  # ----------------------------------------------------------
  # 6. Remove polygons < 1.3 hectare
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 6 of 10 -> Removing polygons < 1 hectare - keeping those that intersect the PRODES cumulative mask.")
  
  query <- sprintf("SELECT * FROM mascara_geral_amz_v2025_nb WHERE tile = '%s'", tile)
  
  mask_union <- read_sf(mask_path, query = query) |>
    sf::st_transform(crs_proc) |>
    sf::st_set_precision(precision) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_union() |>
    sf::st_sf() |>
    sf::st_make_valid()
  
  sits_classification_cloud_cleaned$area_m2 <- as.numeric(sf::st_area(sits_classification_cloud_cleaned))
  sits_classification_cloud_cleaned$area_ha <- sits_classification_cloud_cleaned$area_m2 / 10000
  
  sits_classification_cloud_cleaned <- sf::st_set_precision(sits_classification_cloud_cleaned, precision)|>
    sf::st_make_valid()|>
    sf::st_collection_extract("POLYGON")
  
  sits_classification_cloud_cleaned$touches_mask <- lengths(
    sf::st_intersects(sits_classification_cloud_cleaned, mask_union)
  ) > 0
  
  class_filtered <- sits_classification_cloud_cleaned |>
    dplyr::filter(area_ha >= 1 | touches_mask == TRUE)|>
    sf::st_cast("POLYGON") |>
    sf::st_make_valid()|>
    sf::st_collection_extract("POLYGON")
  
  rm(sits_classification_cloud_cleaned)
  gc()
  log_step_time("Step 6", t_step)
  
  # ----------------------------------------------------------
  # 7. Fill holes < 1.3 hectare
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 7 of 10 -> Merging polygons with the PRODES cumulative mask.")
  
  merged <- list(class_filtered, mask_union) |>
    purrr::map(\(x) {
      sf::st_geometry(x) <- "geom"
      x
    }) |>
    purrr::map(\(x) sf::st_cast(x, "MULTIPOLYGON")) |>
    dplyr::bind_rows() |>
    sf::st_union() |>
    sf::st_set_precision(precision) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON")
  
  message("Step 7 of 10 -> Filling holes < 1.3 ha")
  
  smoothed <- smoothr::fill_holes(
    merged,
    threshold = units::set_units(13000, "m^2")) |>
    sf::st_set_precision(precision) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON")
  
  rm(class_filtered, merged)
  gc()
  
  message("Step 7 of 10 -> Taking off the PRODES cumulative mask.")
  
  class_diff_mask <- sf::st_difference(smoothed, mask_union) |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_cast("POLYGON") |>
    sf::st_sf()|>
    sf::st_set_precision(precision) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON")
  
  log_step_time("Step 7", t_step)
  
  # ----------------------------------------------------------
  # 8. Remove old boundaries polygons
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 8 of 10 -> Calculating shape metrics of polygons.")
  
  supression_polygons <- calculate_edge_metrics(
    class = class_diff_mask,
    prodes_mask = mask_union,
    crs_planar = crs_proc
  )
  
  message("Step 8 of 10 -> Removing old boundaries polygons.")
  
  supression_polygons <- supression_polygons |>
    dplyr::filter(!(prop_comp > 0.1 & prop_comp < 0.9 & raio_equivalente < 35)) |>
    sf::st_set_precision(precision) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON")
  
  rm(smoothed, class_diff_mask)
  gc()
  log_step_time("Step 8", t_step)
  
  # ----------------------------------------------------------
  # Chopping polygons
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 9 of 10-> Chopping polygons")
  
  chopped_polygons <- chop_polygons(supression_polygons, raw_class, mask_union, -51)
  chopped_polygons <- chop_polygons(chopped_polygons, raw_class, mask_union, -16)
  
  log_step_time("Step 9", t_step)
  
  rm(supression_polygons, raw_class, mask_union)
  gc()
  
  # ----------------------------------------------------------
  # 10. Save final result
  # ----------------------------------------------------------
  t_step <- Sys.time()
  message("Step 10 of 10 -> Saving geopackage file.")
  
  final <- chopped_polygons |>
    st_as_sf(sf_column_name = "geom") |>
    dplyr::select(
      any_of(c(
        "fid"))) |>
    sf::st_collection_extract("POLYGON") |>
    sf::st_cast("POLYGON") |>
    sf::st_transform(crs_final)
  
  output_file <- file.path(
    post_class_path,
    paste0("rascunho-sits_t",
           tile, "_",
           end_date_scl, 
           ".gpkg")
  )
  
  sf::st_write(final, dsn = output_file, delete_dsn = TRUE)
  
  log_step_time("Step 10", t_step)
  
  message("Tile ", tile, " successfully processed -> ", output_file)
  
  rm(chopped_polygons, final)
  gc()
  
  return(invisible(output_file))
}

# ============================================================
# 5. # Loop over tiles 
# ============================================================
results <- vector("list", length(tiles))
names(results) <- tiles

for (tile in tiles) {
  results[tile] <- list(
    tryCatch(
      {
        process_tile(tile)
      },
      error = function(e) {
        message("ERROR in tile ", tile, ": ", conditionMessage(e))
        NULL
      }
    )
  )
}

# ============================================================
# 6. Final summary
# ============================================================

success <- names(results)[!vapply(results, is.null, logical(1))]
failure   <- names(results)[vapply(results, is.null, logical(1))]

message("\n========== PROCESSING SUMMARY ==========")
message("Total tiles: ", length(tiles))
message("Success (", length(success), "): ", paste(success, collapse = ", "))
if (length(failure) > 0) {
  message("Failure (", length(failure), "): ", paste(failure, collapse = ", "))
} else {
  message("No faults recorded.")
}
