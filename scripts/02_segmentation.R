# ============================================================
# Segmentation using Simple Non-Iterative Clustering (SNIC) algorithm
# ============================================================

# Load required libraries
library(sf, lib.loc = "/opt/r/R/x86_64-pc-linux-gnu-library/4.4")
library(lubridate)
library(sits)

# Define the parameters: These are user-defined variables
start_date <- '2025-07-01'
end_date   <- '2025-07-31'
tiles      <- '014002'

# Parameters for SNIC segmentation
grid_seeding <- "rectangular"
spacing      <- 10
compactness  <- 0.3
padding      <- 0

# File and folder paths
mixture_path  <- "data/raw/mixture_model"
segments_path <- "data/segments"

# Store processing date
date_process  <- format(Sys.Date(), "%Y-%m-%d")

# ============================================================
# 1. Load Data Cubes (Spectral Bands + Mixture Model Fractions)
# ============================================================

# Step 1.1 -- Create a training cube from a collection
cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07',
                 'B08', 'B8A', 'B11', 'B12'),
  tiles      = tiles,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE
)

# Step 1.2 -- Retrieve Mixture Model Cube from a predefined repository
mm_cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c("SOIL", "VEG", "DARK"),
  tiles      = tiles,
  data_dir   = mixture_path,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE
)

# Step 1.3 -- Merge the Training Cube with Mixture Model Cube
mm_cube_fraction <- sits_merge(mm_cube, cube)

# ============================================================
# 2. Image Segmentation
# ============================================================

# Step 2.1 -- Define the version of your segmentation file
version <- paste0(
  "LSMM-SNIC-spac", spacing,
  "-comp", gsub("[.]", "", as.character(compactness)),
  "-pad", padding,
  "-", grid_seeding,
  "_", date_process
)

# Step 2.2 -- Segment sits cube using SNIC algorithm
mm_cube_segments <- sits_segment(
  cube      = mm_cube_fraction,
  seg_fn    = sits_snic(
    grid_seeding = grid_seeding,
    spacing      = spacing,
    compactness  = compactness,
    padding      = padding
  ),
  memsize    = 180,
  multicores = 28,
  output_dir = segments_path,
  version    = version
)

# ============================================================
# 3. Tile Area and Segmentation Coverage Utilities
# ============================================================

# Step 3.1 -- Compute the bounding box area of the tile
calculate_bbox_area <- function(cube) {
  width    <- as.numeric(cube$xmax[[1]] - cube$xmin[[1]])
  height   <- as.numeric(cube$ymax[[1]] - cube$ymin[[1]])
  area_km2 <- (width * height) / 1e6
  return(round(area_km2))
}

# Step 3.2 -- Retrieve the most recent segmentation file for a tile
get_latest_segment_file <- function(segments_path, tile_number) {
  seg_files <- list.files(
    segments_path,
    pattern    = paste0("_", tile_number, "_.*\\.gpkg$"),
    full.names = TRUE,
    recursive  = FALSE
  )
  
  if (length(seg_files) == 0) return(NULL)
  
  seg_files[order(file.mtime(seg_files), decreasing = TRUE)][1]
}

# Step 3.3 -- Calculate total segmented area from a file
calculate_segments_area <- function(seg_file) {
  if (is.null(seg_file) || !file.exists(seg_file)) return(0)
  
  segments <- tryCatch({
    s <- sf::st_read(seg_file, quiet = TRUE)
    if (!is.null(s) && nrow(s) > 0) s else NULL
  }, error = function(e) {
    message(sprintf("Error reading file: %s\n  -> %s", seg_file, e$message))
    NULL
  })
  
  if (is.null(segments)) return(0)
  
  areas_m2  <- as.numeric(sf::st_area(segments))
  total_km2 <- sum(areas_m2, na.rm = TRUE) / 1e6
  return(round(total_km2))
}

# ============================================================
# 4. Iterative Temporal Expansion for Cloud-Free Segmentation
# ============================================================

# Step 4.1 -- Define initial temporal boundaries
end_date_fixed   <- as.Date(end_date)
start_date_month <- as.Date(start_date, "%Y-%m-%d")

# Step 4.2 -- Remove any existing segmentation files for the tile
seg_files_init <- list.files(
  segments_path,
  pattern = paste0("_", tiles, "_.*\\.gpkg$"),
  full.names = TRUE,
  recursive = FALSE
)

if (length(seg_files_init) > 0) {
  file.remove(seg_files_init)
  cat("Initial segmentation files removed before loop.\n")
}

# Step 4.3 -- Compute fixed bounding box area (reference for coverage)
bbox_area <- calculate_bbox_area(cube)
cat(sprintf("BBOX area (fixed): %d km2\n", bbox_area))
prev_seg_file <- NULL

# Step 4.4 -- Iterative loop (maximum of 5 iterations)
max_iter <- 5
iter <- 0
while (iter < max_iter) {
  iter <- iter + 1
  cat(sprintf("\nIteration %d / %d\n", iter, max_iter))
  cat(sprintf("Period    : %s to %s\n", start_date_month, end_date_fixed))
  cat(sprintf("BBOX area : %d km2\n", bbox_area))
  # Expand temporal window one month backward
  start_date_month <- floor_date(start_date_month - days(1), "month")
  cat(sprintf("Processing period: %s to %s\n",
              start_date_month, end_date_fixed))
  # Load spectral cube
  new_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c('B02','B03','B04','B05','B06','B07',
                   'B08','B8A','B11','B12','CLOUD'),
    tiles      = tiles,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed, "%Y-%m-%d"),
    progress   = TRUE
  )
  # Load fraction (LSMM) cube
  mm_cube <- sits_cube(
    source     = "BDC",
    collection = "SENTINEL-2-16D",
    bands      = c("SOIL", "VEG", "WATER"),
    tiles      = tiles,
    data_dir   = mixture_path,
    start_date = format(start_date_month, "%Y-%m-%d"),
    end_date   = format(end_date_fixed, "%Y-%m-%d"),
    progress   = TRUE
  )
  # Merge cubes
  mm_cube_fraction <- sits_merge(mm_cube, new_cube)
  # Define version identifier
  version <- paste0(
    "LSMM-SNIC-spac", spacing,
    "-comp", gsub("[.]", "", as.character(compactness)),
    "-pad", padding,
    "-", grid_seeding,
    "_", date_process
  )
  # Run segmentation
  sits_segment(
    cube      = mm_cube_fraction,
    seg_fn    = sits_snic(
      grid_seeding = grid_seeding,
      spacing      = spacing,
      compactness  = compactness,
      padding      = padding
    ),
    memsize    = 60,
    multicores = 14,
    output_dir = segments_path,
    version    = version
  )
  # Retrieve latest segmentation result
  new_seg_file <- get_latest_segment_file(segments_path, tiles)
  if (is.null(new_seg_file)) {
    cat("Warning: no segment file found. Stopping.\n")
    break
  }
  cat(sprintf("Generated file: %s\n", basename(new_seg_file)))
  # Remove previous incomplete segmentation
  if (!is.null(prev_seg_file) &&
      file.exists(prev_seg_file) &&
      normalizePath(prev_seg_file) != normalizePath(new_seg_file)) {
    file.remove(prev_seg_file)
    cat(sprintf("Removed previous file: %s\n", basename(prev_seg_file)))
  }
  # Compute segmented area
  total_segment_area <- calculate_segments_area(new_seg_file)
  cat(sprintf("Segments area: %d km2\n", total_segment_area))
  # Check if full spatial coverage is reached
  if (bbox_area <= total_segment_area) {
    if (!is.null(prev_seg_file) && file.exists(prev_seg_file)) {
      file.remove(prev_seg_file)
      cat(sprintf("Removed incomplete file: %s\n",
                  basename(prev_seg_file)))
    }
    cat("Full coverage reached. Stopping loop.\n")
    break
  }
  prev_seg_file <- new_seg_file
}
# Step 4.5 -- Final status
if (iter >= max_iter && bbox_area > total_segment_area) {
  cat("Max iterations reached (5) without full coverage. Processing stopped.\n")
}
print("Segmentation complete.")
