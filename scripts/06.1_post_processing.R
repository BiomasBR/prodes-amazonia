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

# Define the parameters: These are user-defined variables
model_name  <- "tcnn-model_2y_2023-08-01_2025-07-13_2026-08-03_eco-3-mt-46d_2026-08-03_16h02m.rds"
version     <- "tcnn-2y-eco-3-mt-46d-mean"
tiles       <- c('024018','023019','024017','024013','025014','025015','024015','023014',
                 '019015','025016','023015','025017','020019','010017','025019','023020',
                 '022014','018020','017023','016023','022013','025018','022020','021020',
                 '017017','024019','026016','016017','026014','010016','018014','027013',
                 '018019','027014','012015','015021','020018','019016','015018','019020',
                 '026015','027012')

# File and folder paths
seg_version <- "lsmm-snic-spac10-comp03-pad0-rectangular"
class_path  <- "data/class"
mask_path   <- "data/raw/auxiliary/mask_geral_amz_v2025.gpkg"
config_dir  <- ".."

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
  st_make_valid()

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
  
  message(" -> Tile: ", tile_id)
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

# 3.2 Cloud/shadow difference
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

# 3.3 Calculate area, perimeter, shared boundaries and equivalent radius
calculate_edge_metrics <- function(class, prodes_mask, crs_planar = 5880) {

  # Preserves the original state of S2 and ensures restoration upon completion of execution
  s2_state <- sf_use_s2()
  on.exit(sf_use_s2(s2_state), add = TRUE)
  
  # Redesigns and validates geometries
  prodes_mask <- st_transform(prodes_mask, crs_planar) |>
    st_buffer(1) |>
    st_make_valid()
  
  class <- st_transform(class, crs_planar) |>
    st_make_valid()
  
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

# ============================================================
# 4. Main function: process ONE tile
# ============================================================

process_tile <- function(tile) {
  
  message("\n==============================")
  message("Starting tile post-processing: ", tile)
  message("==============================")
  
  # ---- Step 1.3 -- defines the path for the classification raster ----
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
  
  # ---- Step 1.5 -- define and create the post-classification path ----
  post_class_path <- file.path(class_path, tile, "post_processed", version)
  dir.create(post_class_path, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------
  # 2. Classification Classes
  # ----------------------------------------------------------
  raw_class <- rast(raw_class_path)
  levels(raw_class) <- data.frame(
    ID = seq_along(sits_labels(model)),
    classe = sits_labels(model)
  )
  
  labels <- c('Clear_Cut', 'Clear_Cut_Herbaceous','Mininig')
  
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
  vector_multipolygons <- sf::st_as_sf(vector_multipolygons)
  
  # ----------------------------------------------------------
  # 3. Extraction of cloud features
  # ----------------------------------------------------------
  result <- extract_cloud_mask(
    sits_classification_path = raw_class_path,
    sits_reclassification    = vector_multipolygons,
    cloud_values             = c(3, 8, 9, 10),
    output_dir               = post_class_path
  )
  
  cloud_vec    <- result$cloud_vec
  end_date_scl <- result$end_date_scl
  
  # ----------------------------------------------------------
  # 4. Cloud/shadow difference
  # ----------------------------------------------------------
  sits_classification_cloud_cleaned <- remove_cloud_areas(
    sits_reclassification = vector_multipolygons,
    cloud_vec             = cloud_vec,  # NULL se nao houver nuvens
    buffer_dist           = 100
  )
  
  # ----------------------------------------------------------
  # 5. Fill holes < 1 hectare
  # ----------------------------------------------------------
  query <- sprintf("SELECT * FROM mascara_geral_amz_v2025 WHERE tile = '%s'", tile)
  prodes_mask <- read_sf(mask_path, query = query)
  
  prodes_mask <- sf::st_transform(
    prodes_mask,
    sf::st_crs(sits_classification_cloud_cleaned)
  )
  
  merged <- list(sits_classification_cloud_cleaned, prodes_mask) |>
    purrr::map(sf::st_make_valid) |>
    purrr::map(\(x) sf::st_transform(x, sf::st_crs(sits_classification_cloud_cleaned))) |>
    purrr::map(\(x) {
      sf::st_geometry(x) <- "geom"
      x
    }) |>
    purrr::map(\(x) sf::st_cast(x, "MULTIPOLYGON")) |>
    dplyr::bind_rows() |>
    sf::st_union()
  
  smoothed <- smoothr::fill_holes(
    merged,
    threshold = units::set_units(10000, "m^2")
  )
  
  # ----------------------------------------------------------
  # 6. Difference with deforestation mask
  # ----------------------------------------------------------
  smoothed <- sf::st_transform(smoothed, sf::st_crs(prodes_mask)) |>
    st_make_valid()
  
  prodes_mask <- prodes_mask |>
    st_make_valid()
  
  mask_union <- prodes_mask |>
    st_union() |>
    st_make_valid()
  
  class_diff_mask <- sf::st_difference(
    smoothed,
    mask_union
  ) |> sf::st_collection_extract("POLYGON") |>
    sf::st_cast("POLYGON") |>
    sf::st_sf()
  
  # ----------------------------------------------------------
  # 7. Remove polygons outside the biome border
  # ----------------------------------------------------------
  biome_tile <- st_transform(biome, st_crs(class_diff_mask))
  
  if (tile %in% edge_tiles) {
    message("The tile ", tile, " is an edge tile. Running intersection.")
    class_biome <- st_intersection(class_diff_mask, biome_tile)
  } else {
    message("The tile ", tile, " is not an edge tile. Intersection ignored.")
    class_biome <- class_diff_mask
  }
  
  # ----------------------------------------------------------
  # 8. Remove polygons < 1 hectare
  # ----------------------------------------------------------
  class_biome$area_m2 <- as.numeric(sf::st_area(class_biome))
  class_biome$area_ha <- class_biome$area_m2 / 10000
  
  class_biome_bigger_than_1ha <- class_biome |>
    dplyr::filter(area_ha >= 1)
  
  supression_polygons <- st_transform(
    class_biome_bigger_than_1ha,
    crs = 4674
  ) |>
    sf::st_cast("POLYGON") |>
    sf::st_make_valid()
  
  # # ----------------------------------------------------------
  # # 9. Calculate old boundaries polygons
  # # ----------------------------------------------------------

  supression_polygons <- calculate_edge_metrics(
    class = supression_polygons,
    prodes_mask = prodes_mask,
    crs_planar = 5880
  )

  # # ----------------------------------------------------------
  # # 10. Remove old boundaries polygons
  # # ----------------------------------------------------------
  # 
   supression_polygons <- supression_polygons |>
      dplyr::filter(!(prop_comp > 0.1 & prop_comp < 0.9 & raio_equivalente < 35)) |>
      sf::st_transform(4674)
  
  # ----------------------------------------------------------
  # 11. Assigns names of the classes with the greatest spatial intersection
  # ----------------------------------------------------------
  sf_use_s2(FALSE)
  
  vector_multipolygons_valid <- vector_multipolygons |>
    st_transform(st_crs(supression_polygons)) |>
    st_make_valid()
  
  supression_polygons <- supression_polygons |>
    mutate(.id_temp = row_number())
  
  intersecao <- st_intersection(supression_polygons, vector_multipolygons_valid)
  intersecao$area_intersec <- st_area(intersecao)
  maior_classe <- intersecao |>
    st_drop_geometry() |>
    group_by(.id_temp) |>
    slice_max(order_by = area_intersec, n = 1, with_ties = FALSE) |>
    select(.id_temp, class)
  
  supression_polygons <- supression_polygons |>
    left_join(maior_classe, by = ".id_temp") |>
    select(-.id_temp)
  
  sf_use_s2(TRUE)
  
  # ----------------------------------------------------------
  # 12. Select Boundaries Segments
  # ----------------------------------------------------------
  vector_path <- list.files(
    "data/segments",
    pattern = paste0("SENTINEL-2_MSI_", tile, "_.*_segments_", seg_version, "(_\\d{4}-\\d{2}-\\d{2})?\\.gpkg$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(vector_path) == 0) {
    stop("No segment file found for the tile ", tile)
  }
  
  if (length(vector_path) > 1) {
    stop(
      "More than one segment file found for the tile ", tile, ":\n",
      paste(" -", vector_path, collapse = "\n"),
      "\nAdjust the search pattern (or remove duplicate files) so that only 1 remains."
    )
  }
  
  sf_use_s2(FALSE)
  
  segments <- read_sf(
    vector_path,
    wkt_filter = st_as_text(st_combine(st_transform(st_make_valid(supression_polygons), st_layers(vector_path)$crs[[1]])))
  ) |>
    st_transform(st_crs(supression_polygons)) |>
    st_make_valid() |>
    st_difference(st_union(st_union(st_make_valid(supression_polygons)), st_make_valid(st_transform(mask_union, st_crs(supression_polygons))))) |>
    dplyr::filter(!st_is_empty(geom))
  
  sf_use_s2(TRUE)
  
  # ----------------------------------------------------------
  # 13. Save final result
  # ----------------------------------------------------------
  st_geometry(segments) <- "geom"
  st_geometry(supression_polygons) <- "geom"
  
  merged_polygons <- bind_rows(supression_polygons, segments) |>
    st_as_sf(sf_column_name = "geom") |>
    dplyr::select(
      any_of(c(
        "fid", 
        "class"
      ))
    )
  
  output_file <- file.path(
    post_class_path,
    paste0("class-post-processed_",
           tile, "_",
           years, "_",
           end_date_scl, "_",
           version, ".gpkg")
  )
  
  sf::st_write(merged_polygons, dsn = output_file, delete_dsn = TRUE)
  
  message("Tile ", tile, " successfully processed -> ", output_file)
  
  return(invisible(output_file))
}

# ============================================================
# 5. # Loop over tiles 
# ============================================================

resultados <- vector("list", length(tiles))
names(resultados) <- tiles

for (tile in tiles) {
  
  resultados[tile] <- list(
    tryCatch(
      {
        withCallingHandlers(
          {
            process_tile(tile)
          },
          warning = function(w) {
            message("AVISO no tile ", tile, ": ", conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) {
        message("ERRO no tile ", tile, ": ", conditionMessage(e))
        NULL  # marca falha e permite que o loop continue para o proximo tile
      }
    )
  )
}

# ============================================================
# 6. Final summary
# ============================================================

sucesso <- names(resultados)[!vapply(resultados, is.null, logical(1))]
falha   <- names(resultados)[vapply(resultados, is.null, logical(1))]

message("\n========== PROCESSING SUMMARY ==========")
message("Total tiles: ", length(tiles))
message("Success (", length(sucesso), "): ", paste(sucesso, collapse = ", "))
if (length(falha) > 0) {
  message("Failure (", length(falha), "): ", paste(falha, collapse = ", "))
} else {
  message("No faults recorded.")
}
