# ============================================================
# Post-Processing in Forest Areas
# ============================================================

# ============================================================
# 1. Libraries, paths, and initial parameters
# ============================================================

# Step 1.1 -- Load required libraries
library(sf)
library(dplyr)
library(sits)
library(terra)
library(units)
library(smoothr)

# Step 1.2 -- Define paths for files and folders
tile            <- "012014"
version         <- "rf-1y-all-samples-new-pol-avg-false-final-seg"
class_path      <- "data/class"
mask_path       <- "data/raw/auxiliary/mask_geral_amz_v2024.gpkg" #nome da máscara em gpkg geral
config_dir      <- ".."

# Step 1.3 -- define raw classification path and load the file
raw_class_path <- list.files(class_path,
                             pattern <- paste0(".*_", tile, "_.*_class_", version, "\\.gpkg$"),
                             full.names = TRUE,
                             recursive = TRUE)
raw_class <- read_sf(raw_class_path)

# Step 1.4 -- extract the number of 'years'
years <- regmatches(version, regexpr("\\d+y", version))

# Step 1.5 -- define and create post classification path
post_class_path <- file.path(class_path, tile, "post_processed", version)
dir.create(post_class_path,
           showWarnings = FALSE,
           recursive = TRUE)

# ============================================================
# . Translation Stage
# ============================================================

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

config     <- read_class_config(file.path(config_dir, "class_config.txt"))
class_translation <- config$class_translation

raw_class$class <- ifelse(
       raw_class$class %in% names(class_translation),
       class_translation[raw_class$class],
       raw_class$class
   )

raw_class <- raw_class |>
  rename_with(~ class_translation[.x],
              .cols = any_of(names(class_translation)))

# ============================================================
# 2. Probabilistic reclassification
# ============================================================

# Step 2.1 -- Create sum columns
raw_class <- raw_class |>
  mutate(
    Suppression_sum = rowSums(
      across(
        all_of(c(
          "Clear_Cut_Bare_Soil",
          "Clear_Cut_Vegetation",
          "Clear_Cut_Trees"
        ))
      ),
      na.rm = TRUE
    ),
    Degrad_sum = rowSums(
      across(
        all_of(c(
          "Degradation",
          "Degradation_Fire"
        ))
      ),
      na.rm = TRUE
    ),
    Natural = rowSums(
     across(
        all_of(c(
          "Forest",
          "Non_Forest_Natural_Vegetation",
          "Transition_Forest",
          "Wetland"
        ))
      ),
      na.rm = TRUE
    )
  )

# Step 2.4 -- Reclassification
raw_class <- raw_class |>
  mutate(
    class = if_else(
      Suppression_sum >= Degrad_sum &
        Suppression_sum >= Natural &
        Suppression_sum >= Water,
      "Suppression_sum",
      class
    )
  )

# Step 2.5 -- Filter only suppression-related classes
post_class <- raw_class |>
  filter(
    class %in% c(
      "Clear_Cut_Bare_Soil",
      "Clear_Cut_Vegetation",
      "Clear_Cut_Burned_Area",
      "Clear_Cut_Trees",
      "Suppression_sum"
    )
  ) |>
  mutate(class = "supression")

# Step 2.7 -- Dissolve and aggregate geometries
post_class <- post_class |>
  group_by(class) |>
  summarise(
    geometry = st_union(geom),
    .groups = "drop"
  ) |>
  st_cast("MULTIPOLYGON")

# ============================================================
# 3. Extraction of cloud features
# ============================================================

# Step 3.1 -- Function definition
extract_cloud_mask <- function(
    sits_classification_path,
    sits_reclassification,
    cloud_values = c(3, 8, 9, 10),
    output_dir = NULL,
    collection = "SENTINEL-2-16D",
    date_window_days = 1
) {
  
  # ----------------------------------------------------------
  # 1. Extract metadata from filename
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
  # 2. Build BDC cube
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
  # 3. Load SCL raster
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
  # 4. Create cloud mask
  # ----------------------------------------------------------
  scl_mask <- terra::classify(
    scl_raster,
    rcl = cbind(cloud_values, rep(1, length(cloud_values))),
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
  
  # ----------------------------------------------------------
  # 5. Crop to classification extent
  # ----------------------------------------------------------
  class_bbox <- sits_reclassification |>
    sf::st_transform(terra::crs(scl_mask)) |>
    terra::vect() |>
    terra::ext()
  
  scl_raster_crop <- terra::crop(scl_mask, class_bbox)
  
  # ----------------------------------------------------------
  # 6. Vectorize cloud mask
  # ----------------------------------------------------------
  cloud_vec <- terra::as.polygons(scl_raster_crop, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_transform(sf::st_crs(sits_reclassification)) |>
    smoothr::fill_holes(threshold = Inf)
  
  # ----------------------------------------------------------
  # 7. Save (optional)
  # ----------------------------------------------------------
  if (!is.null(output_dir)) {
    output_filename <- paste0(
      "cloud-vec_", tile_id, "_", end_date_scl, ".gpkg"
    )
    
    output_path <- file.path(output_dir, output_filename)
    
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

# Step 3.2 -- Apply function
result <- extract_cloud_mask(
  sits_classification_path = raw_class_path,
  sits_reclassification    = post_class,
  cloud_values             = c(3, 8, 9, 10),
  output_dir               = post_class_path
)

# Step 3.3 -- Extract outputs
cloud_vec   <- result$cloud_vec
end_date_scl <- result$end_date_scl

# ============================================================
# 4. Cloud/shadow difference
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
  sits_reclassification = post_class,
  cloud_vec             = cloud_vec,  # NULL if there are no clouds
  buffer_dist           = 100
)

# ============================================================
# 5. Fill holes < 1 hectare (first round)
# ============================================================

query <- sprintf("SELECT * FROM mask_geral_amz_v2024 WHERE tile = '%s'", tile) #nome da máscara em gpkg geral
prodes_mask <- read_sf(mask_path,
                       query = query) 

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

# ============================================================
# 6. Difference with deforestation mask (first round)
# ============================================================

smoothed <- sf::st_transform(smoothed, sf::st_crs(prodes_mask))

smoothed <- smoothed |> 
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

# ============================================================
# 7. Fill bays and smooth edges 
# ============================================================

class_diff_mask_filled_bays <- class_diff_mask |>
  sf::st_buffer(dist = 50) |>
  sf::st_buffer(dist = -50) |>
  sf::st_cast("POLYGON") |>
  tibble::rowid_to_column("id")

# ============================================================
# 8. Fill holes < 1 hectare (second round)
# ============================================================

merged_2 <- sf::st_union(
  c(
    sf::st_geometry(sf::st_make_valid(class_diff_mask_filled_bays)),
    sf::st_geometry(mask_union)
  )
) |>
  sf::st_collection_extract("POLYGON") |>  # pull only polygon parts
  sf::st_union() |>                         # re-union into one geometry
  sf::st_make_valid()

smoothed_2 <- smoothr::fill_holes(
  merged_2,
  threshold = units::set_units(10000, "m^2")
) |>
  sf::st_make_valid()

# ============================================================
# 9. Difference with deforestation mask (second round)
# ============================================================

class_diff_mask_2 <- sf::st_difference(
  smoothed_2,
  mask_union
) |>
  sf::st_make_valid() |>
  sf::st_collection_extract(type = "POLYGON") |>
  sf::st_cast("POLYGON") |>
  sf::st_sf()|>
  st_make_valid()

# ============================================================
# 10. Remove polygons outside the biome border
# ============================================================

biome <- read_sf("data/raw/auxiliary/borders/amazon-biome-border-epsg10857.gpkg") |>
  st_make_valid() |>
  st_transform(st_crs(class_diff_mask_2))

class_biome <- st_intersection(class_diff_mask_2, biome)

# ============================================================
# 11. Remove polygons < 1 hectare
# ============================================================

class_biome$area_m2 <- as.numeric(sf::st_area(class_biome))
class_biome$area_ha <- class_biome$area_m2 / 10000

class_biome_bigger_than_1ha <- class_biome |>
  dplyr::filter(area_ha >= 1)

# ============================================================
# 12. Save final result
# ============================================================

supression_polygons <- st_transform(
  class_biome_bigger_than_1ha,
  crs = 4674
) |>
  sf::st_cast("POLYGON") |>
  sf::st_make_valid()

sf::st_write(
  supression_polygons,
  file.path(
    post_class_path,
    paste0("class-post-processed_",
           tile,"_",years,"_",
           end_date_scl,"_",
           version, ".gpkg")
  )
)
