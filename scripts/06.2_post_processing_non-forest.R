# ============================================================
# Post-Processing in Non-Forest areas 
# ============================================================

# =============================================================================
# 1. Libraries, paths and some initial parameters
# =============================================================================

# Step 1.1 -- Load Required Libraries
library(sf)
library(dplyr)
library(sits)
library(terra)
library(units)
library(smoothr)

# Step 1.2 -- Define paths for files and folders
tile            <- "012014"
version         <- "rf-3y-all-samples-new-pol-avg-false"
class_path      <- "data/class"
mask_path       <- "data/raw/auxiliary/mask_geral_amz_v2024.gpkg" #nome da máscara em gpkg geral

# Step 1.3 -- define raw classification path and load the file
raw_class_path <- list.files(class_path,
                             pattern = paste0(".*", tile, ".*", "class",
                                              "*.*", version, ".*\\.gpkg$"),
                             full.names = TRUE,
                             recursive = TRUE)
raw_class <- read_sf(raw_class_path)

# Step 1.4 -- extract the number of 'years'
years <- regmatches(version, regexpr("\\d+y", version))

# Step 1.5 -- define and create post classification path
post_class_path <- file.path(class_path, tile, "post_processed")
dir.create(post_class_path,
           showWarnings = FALSE,
           recursive = TRUE)

# ============================================================
# 2. Probabilistic reclassification
# ============================================================

# Step 2.1 -- Create column "Soma_Supressao"
sits_classification_raw <- sits_classification_raw |>
  mutate(
    Soma_Supressao = rowSums(across(all_of(c(
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto",
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
      "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo"
    ))), na.rm = TRUE)
  )

# Step 2.2 -- Create column "Soma_NF"
sits_classification_raw <- sits_classification_raw |>
  mutate(
    Soma_NF = rowSums(across(all_of(c(
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Mais_Biomassa",
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Menos_Biomassa",
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Umida",
      "Vegetacao_Natural_Nao_Florestal_Vereda",
      "Vegetacao_Natural_Nao_Florestal_Mata",
      "Vegetacao_Natural_Nao_Florestal_Transicao_Florestal",
      "Vegetacao_Natural_Nao_Florestal_Herbacea_Seca_Pos_Fogo"
    )
    )), na.rm = TRUE)
  )

# Step 2.3 -- Create column "Soma_Agua"
sits_classification_raw <- sits_classification_raw |>
  mutate(
    Soma_Agua = rowSums(across(all_of(c(
      "Hidrografia_Lago",
      "Hidrografia_Rio"
    ))), na.rm = TRUE)
  )

# Step 2.4 -- Reclassify all segments where Sum_Suppression > Sum_NF & Sum_Suppression >= Sum_Water
sits_classification_raw <- sits_classification_raw |>
  mutate(
    class = if_else(Soma_Supressao >= Soma_NF & Soma_Supressao >= Soma_Agua, "Soma_Supressao", class)
  )

# Step 2.5 -- Filter the dataset to keep only the classes related to suppression
sits_classification <- sits_classification_raw |>
  dplyr::filter(class %in% c(
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura",
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto",
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Agricultura_Antigo",
    "Supressao_de_Vegetacao_Natural_Nao_Florestal_Com_Solo_Exposto_Antigo",
    "Soma_Supressao" # added in the previous step
  ))

# Step 2.6 -- Standardize the nomenclature of the classes
sits_reclassification <- sits_classification |>
  dplyr::mutate(class = "supressao")

# Step 2.7 -- Dissolving and aggregating geometries
sits_reclassification <- sits_reclassification |>
  dplyr::group_by(class) |>
  dplyr::summarise(geometry = sf::st_union(geom), .groups = "drop") |>
  sf::st_cast("MULTIPOLYGON")

# ============================================================
# 3. Extraction of cloud features
# ============================================================
# Extracts vector cloud/shadow mask from the SCL band at the latest date
# from the BDC cube that generated the sits classification

# Step 3.1 -- Defining the 'extract_cloud_mask' function
extract_cloud_mask <- function(
    sits_classification_path,
    sits_reclassification,
    cloud_values    = c(3, 8, 9, 10),
    output_dir      = NULL,
    collection      = "SENTINEL-2-16D",
    date_window_days = 1
) {
  # ---------------------------------------------------------------------------
  # 1. Extract metadata from the file name
  # ---------------------------------------------------------------------------
  filename_base <- basename(sits_classification_path)
  
  last_date_str <- regmatches(
    filename_base,
    gregexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", filename_base)
  )[[1]]
  
  if (length(last_date_str) == 0) {
    stop("No date in YYYY-MM-DD format found in the file name: ", filename_base)
  }
  
  last_date      <- as.Date(tail(last_date_str, 1))
  start_date_scl <- format(last_date - date_window_days, "%Y-%m-%d")
  end_date_scl   <- format(last_date, "%Y-%m-%d")
  
  tile_id <- regmatches(
    filename_base,
    regexpr("(?<=SENTINEL-2_MSI_)[0-9]+", filename_base, perl = TRUE)
  )
  
  if (length(tile_id) == 0 || tile_id == "") {
    stop("Could not extract tile_id from the file name: ", filename_base)
  }
  
  message("  -> Tile: ", tile_id)
  message("  -> SCL time window: ", start_date_scl, " to ", end_date_scl)
  
  # ---------------------------------------------------------------------------
  # 2. Build BDC cube with SCL band
  # ---------------------------------------------------------------------------
  scl_cube <- sits::sits_cube(
    source     = "BDC",
    collection = collection,
    tiles      = tile_id,
    bands      = "SCL",
    start_date = start_date_scl,
    end_date   = end_date_scl
  )
  
  # ---------------------------------------------------------------------------
  # 3. Load SCL raster
  # ---------------------------------------------------------------------------
  scl_files <- scl_cube$file_info[[1]] |>
    dplyr::filter(band == "CLOUD", date == as.Date(end_date_scl)) |>
    dplyr::pull(path)
  
  if (length(scl_files) == 0) {
    stop("No SCL file found for the date: ", end_date_scl)
  }
  
  scl_raster <- terra::rast(scl_files[1])
  message("  -> SCL file loaded: ", scl_files[1])
  
  # ---------------------------------------------------------------------------
  # 4. Create binary cloud/shadow mask
  # ---------------------------------------------------------------------------
  scl_mask <- terra::classify(
    scl_raster,
    rcl    = cbind(cloud_values, rep(1, length(cloud_values))),
    others = NA
  )

  # ---------------------------------------------------------------------------
  # 4.1 Check if any cloud values were found
  # ---------------------------------------------------------------------------
  if (all(is.na(terra::values(scl_mask)))) {
    message("  -> No clouds identified")
    return(invisible(list(
      cloud_vec    = NULL,
      tile_id      = tile_id,
      end_date_scl = end_date_scl
    )))
  }
  
  # ---------------------------------------------------------------------------
  # 4.1 Check if any cloud values were found
  # ---------------------------------------------------------------------------
  if (all(is.na(terra::values(scl_mask)))) {
    message("  -> No clouds identified")
    return(invisible(list(
      cloud_vec    = NULL,
      tile_id      = tile_id,
      end_date_scl = end_date_scl
    )))
  }
  
  # ---------------------------------------------------------------------------
  # 5. Crop by the classification extent
  # ---------------------------------------------------------------------------
  class_bbox <- sits_reclassification |>
    sf::st_transform(terra::crs(scl_mask)) |>
    terra::vect() |>
    terra::ext()
  
  scl_raster_crop <- terra::crop(scl_mask, class_bbox)
  
  # ---------------------------------------------------------------------------
  # 6. Vectorize and fill holes in the cloud mask
  # ---------------------------------------------------------------------------
  cloud_vec <- terra::as.polygons(scl_raster_crop, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_transform(sf::st_crs(sits_reclassification)) |>
    smoothr::fill_holes(threshold = Inf)  # fills holes of all sizes
  
  # ---------------------------------------------------------------------------
  # 7. Save result (optional)
  # ---------------------------------------------------------------------------
  if (!is.null(output_dir)) {
    output_filename <- paste0(
      "Nuvem_",
      tile_id, "_",
      end_date_scl,
      ".gpkg"
    )
    output_path <- file.path(output_dir, output_filename)
    
    cloud_vec |>
      sf::st_transform(4674) |>   # reproject to EPSG:4674 before saving
      sf::st_write(output_path, append = FALSE)
    
    message("  -> Cloud vector saved in EPSG:4674: ", output_path)
  }
  
  # ---------------------------------------------------------------------------
  # 8. Return list with cloud vector and extracted metadata
  # ---------------------------------------------------------------------------
  return(invisible(list(
    cloud_vec    = cloud_vec,
    tile_id      = tile_id,
    end_date_scl = end_date_scl
  )))
}

# Step 3.2 -- Applying the 'extract_cloud_mask' function
result <- extract_cloud_mask(
  sits_classification_path = sits_classification_path,
  sits_reclassification      = sits_reclassification,
  cloud_values             = c(3, 8, 9, 10),
  output_dir               = output_dir
)

# Step 3.3 -- Access each object in the result individually
cloud_vec    <- result$cloud_vec
tile_id      <- result$tile_id
end_date_scl <- result$end_date_scl

# ============================================================
# 4. Difference with cloud/shadow
# ============================================================

# Step 4.1 -- Defining the 'remove_cloud_areas' function
remove_cloud_areas <- function(
    sits_reclassification,
    cloud_vec,
    buffer_dist = buffer_dist
) {
  # Check if cloud_vec exists and has features
  if (is.null(cloud_vec) || nrow(cloud_vec) == 0) {
    message("  -> No cloud vectors were found")
    return(invisible(sits_reclassification))
  }
  
  # Dissolve
  cloud_union <- sf::st_union(cloud_vec) 
 
   # Buffer
  cloud_vec_buffer <- sf::st_buffer(cloud_union, dist = buffer_dist) 
  
  # Remove cloud/shadow areas from classification
  sits_classification_cloud_cleaned <- sf::st_difference(
    sits_reclassification,
    cloud_vec_buffer
  ) |>
    sf::st_cast("MULTIPOLYGON")
  
  return(invisible(sits_classification_cloud_cleaned))
}

# Step 4.2 -- Run 'remove_cloud_areas' function
sits_classification_cloud_cleaned <- remove_cloud_areas(
  sits_reclassification = sits_reclassification,
  cloud_vec             = cloud_vec,  # NULL if there are no clouds
  buffer_dist           = 100
)

# ============================================================
# 5. Fill holes < 1 hectares / First round
# ============================================================

# Step 5.1 -- Read accumulated deforestation data, filtering only for the tile under analysis
mask <- sf::st_read(mask_path, quiet = TRUE) |>
  dplyr::filter(tile == tile_id)

# Step 5.2 -- Ensure that the mask is in the same spatial reference system
mask <- sf::st_transform(
  mask,
  sf::st_crs(sits_classification_cloud_cleaned)
)

# Step 5.3 -- Merge reclassification with the accumulated mask,
# dissolving all geometries into a single continuous set
merged <- sf::st_union(
  c(
    sf::st_geometry(sf::st_make_valid(sits_classification_cloud_cleaned)),
    sf::st_geometry(sf::st_make_valid(mask))
  )
)

# Step 5.4 -- Fills internal holes smaller than 1 hectare (10000 m²)
smoothed <- merged |>
  sf::st_collection_extract("POLYGON") |>
  sf::st_union() |>
  smoothr::fill_holes(threshold = units::set_units(10000, "m^2"))

# ============================================================
# 6. Difference with the deforestation mask / First round
# ============================================================

# Step 6.1 -- Ensure consistent CRS prior to space operation
smoothed <- sf::st_transform(smoothed, sf::st_crs(mask))

# Step 6.2 -- Fix geometries
smoothed <- smoothed |>
  st_make_valid()
mask <- mask |>
  st_make_valid()

# Step 6.3 -- Make a union
mask_union <- mask |>
  st_union() |>
  st_make_valid()

# Step 6.4 -- Remove areas already present in the deforestation mask
class_diff_mask <- sf::st_difference(
  smoothed,
  mask_union
) |>
  sf::st_cast("POLYGON") |>
  sf::st_sf()

# ============================================================
# 7. Fill in "bays" and smooth edges
# ============================================================

# Step 7.1 -- Mathematical morphology / Positive buffer then negative buffer
class_diff_mask_filled_bays <- class_diff_mask |>
  sf::st_buffer(dist =  50)                    |>
  sf::st_buffer(dist = -50)                    |>
  sf::st_cast("POLYGON")                       |>
  tibble::rowid_to_column("id")

# ============================================================
# 8. Fill holes < 1 hectares / Second round
# ============================================================

# Step 8.1 -- Merge with the accumulated mask, dissolving all geometries into a single continuous set
merged_2 <- sf::st_union(
  c(
    sf::st_geometry(sf::st_make_valid(class_diff_mask_filled_bays)),
    sf::st_geometry(sf::st_make_valid(mask))
  )
)

# Step 8.2 -- Fills internal holes smaller than 1 hectare (10000 m²)
smoothed_2 <- merged_2 |>
  sf::st_collection_extract("POLYGON") |>
  sf::st_union() |>
  smoothr::fill_holes(threshold = units::set_units(10000, "m^2"))

# ============================================================
# 9. Difference with the deforestation mask / Second round
# ============================================================

# Step 9.1 -- Remove areas already present in the deforestation mask
class_diff_mask_2 <- sf::st_difference(
  smoothed_2,
  mask_union) |>
  sf::st_cast("POLYGON") |>
  sf::st_sf()

# ============================================================
# 10. Remove polygons outside the biome border
# ============================================================

biome <- read_sf("data/raw/amazon-nf-biome-border-epsg10857.gpkg") |>
  st_make_valid()

class_biome <- st_intersection(class_diff_mask_2, biome)

# ============================================================
# 11. Remove polygons < 1 hectare
# ============================================================

# Step 11.1 -- Calculate the area of each polygon (in m²)
class_biome$area_m2 <- as.numeric(sf::st_area(class_biome))

# Step 11.2 -- Convert square meters to hectares
class_biome$area_ha <- class_biome$area_m2 / 10000

# Step 11.3 -- Keep only polygons with an area greater than or equal to 1 hectare
class_biome_bigger_than_1ha <- class_diff_mask_2 |>
  dplyr::filter(area_ha >= 1)

# ============================================================
# 12. Save final result
# ============================================================

# Step 12.1 -- Reproject to EPSG:4674 (SIRGAS 2000)
supression_polygons <- st_transform(
  class_biome_bigger_than_1ha,
  crs = 4674
) |>
  sf::st_cast("POLYGON") |>
  sf::st_make_valid()

# Step 12.2 -- Save the final result
sf::st_write(
  supression_polygons,
  file.path(
    post_class_path,
    paste0("class-post-processed_",
           tile,"_",years,"_",
           end_date_scl,"_",
           version,".gpkg")
  )
)
