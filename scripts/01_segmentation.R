# ============================================================
# Linear Spectral Mixture Model (LSMM) + Simple Non-Iterative Clustering (SNIC) Algorithm
# ============================================================

# Load required libraries
library(sits)
library(tibble)
library(sf, lib.loc = "/opt/r/R/x86_64-pc-linux-gnu-library/4.4")

# Define the parameters: These are user-defined variables
tiles           <- c('004010','026004')
start_date      <- '2026-07-12'
end_date        <- '2026-07-28'
n_cores         <- 28    # adapt to your computer CPU core availability
memsize         <- 180   # adapt to your computer memory availability

sits_parallel(workers = n_cores)

# Parameters for SNIC segmentation
grid_seeding    <- "rectangular"
spacing         <- 10
compactness     <- 0.3
padding         <- 0

# Tolerance (in pixels) between segments BBOX and tile feature BBOX
bbox_tol_px     <- 1

# File and folder paths
mixture_path    <- "data/raw/mixture_model"
segments_path   <- "data/segments"
tiles_bdc_path  <- "~/grupos/biomasbr-amazonia/dados-auxiliares/tiles-bdc-amz.gpkg"
tiles_id_col    <- "tile"   # column with the tile identifier in the gpkg

# Store processing date
date_process    <- format(Sys.Date(), "%Y-%m-%d")

# Spectral bands
bands           <- c('B02', 'B03', 'B04', 'B05', 'B06', 'B07',
                     'B08', 'B8A', 'B11', 'B12','CLOUD')

# Endmember values (class names cannot contain numbers)
endmembers <- tibble::tribble(
  ~class,    ~B02, ~B03, ~B04, ~B05, ~B06, ~B07, ~B08, ~B8A, ~B11, ~B12,
  "soil",    1799, 2154, 3028, 3303, 3472, 3656, 3566, 3686, 5097, 4736,
  "veg",      827,  892,  410, 1070, 4206, 5646, 5495, 6236, 2101,  775,
  "dark",     946,  739,  280,  208,  180,  167,  135,  129,   26,   14,
)
fraction_bands  <- toupper(endmembers$class)

# Version of the segmentation files
version <- paste0(
  "LSMM-SNIC-spac", spacing,
  "-comp", gsub("[.]", "", as.character(compactness)),
  "-pad", padding,
  "-", grid_seeding,
  "_", date_process
)

# ============================================================
# 1. BDC Tiles Grid and Utilities
# ============================================================

# Step 1.1 -- Load the BDC tiles grid (tile identifiers as 6 digits)
tiles_bdc <- sf::st_read(tiles_bdc_path, quiet = TRUE)

if (!tiles_id_col %in% names(tiles_bdc)) {
  stop(sprintf("Column '%s' not found in %s. Available columns: %s",
               tiles_id_col, tiles_bdc_path, paste(names(tiles_bdc), collapse = ", ")))
}
tiles_bdc$tile_id <- formatC(as.integer(gsub("\\D", "", tiles_bdc[[tiles_id_col]])),
                             width = 6, flag = "0")

missing_tiles <- setdiff(tiles, tiles_bdc$tile_id)
if (length(missing_tiles) > 0) {
  warning(sprintf("Tiles not found in the BDC grid: %s", paste(missing_tiles, collapse = ", ")))
}

# Step 1.2 -- Retrieve the most recent segmentation file for a tile
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

# Step 1.3 -- Find the tile of the BDC layer where the segments BBOX is located
#             and compare it with the expected tile
check_bbox_location <- function(seg_file, tile_number, tolerance) {
  segments <- sf::st_read(seg_file, quiet = TRUE)
  seg_bbox <- sf::st_bbox(segments)
  seg_poly <- sf::st_as_sfc(seg_bbox)
  grid     <- sf::st_transform(tiles_bdc, sf::st_crs(segments))
  
  # Tile with the largest overlap with the segments BBOX
  idx <- sf::st_intersects(seg_poly, grid)[[1]]
  if (length(idx) == 0) {
    return(list(status  = "mismatch",
                message = "segments BBOX does not intersect any tile in the BDC tiles layer"))
  }
  overlap  <- vapply(idx, function(i) {
    as.numeric(sf::st_area(sf::st_intersection(sf::st_geometry(grid)[i], seg_poly)))
  }, numeric(1))
  found    <- idx[which.max(overlap)]
  found_id <- grid$tile_id[found]
  
  # Largest offset between the BBOX sides (xmin, ymin, xmax, ymax)
  offset <- max(abs(as.numeric(seg_bbox) - as.numeric(sf::st_bbox(grid[found, ]))))
  
  if (found_id != tile_number) {
    list(status  = "mismatch",
         message = sprintf("segments BBOX is located on tile %s in the BDC tiles layer, not on tile %s",
                           found_id, tile_number))
  } else if (offset > tolerance) {
    list(status  = "partial",
         message = sprintf("segments BBOX is located on tile %s in the BDC tiles layer, but does not cover its full extent (max offset: %.1f m)",
                           found_id, offset))
  } else {
    list(status  = "ok",
         message = sprintf("segments BBOX matches the location of tile %s in the BDC tiles layer",
                           found_id))
  }
}

# ============================================================
# 2. Loop: LSMM + segmentation + BBOX check for each tile
# ============================================================

tile_status  <- setNames(character(length(tiles)), tiles)
tile_message <- setNames(character(length(tiles)), tiles)

for (tile in tiles) {
  
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("Processing tile: %s\n", tile))
  cat(sprintf("============================================================\n"))
  
  tile_result <- tryCatch({
    
    # ------------------------------------------------------------
    # 2.1 Load Data Cube and Create Fraction Images (LSMM)
    # ------------------------------------------------------------
    
    # Step 2.1.1 -- Create a cube from a collection
    cube <- sits_cube(
      source     = "BDC",
      collection = "SENTINEL-2-16D",
      bands      = bands,
      tiles      = tile,
      start_date = start_date,
      end_date   = end_date,
      progress   = TRUE
    )
    
    # Step 2.1.2 -- Generate the mixture model cube and calculate the process duration
    sits_mixture_model_start <- Sys.time()
    sits_mixture_model(
      data       = cube,
      endmembers = endmembers,
      multicores = n_cores,
      memsize    = memsize,
      output_dir = mixture_path
    )
    sits_mixture_model_time <- as.numeric(Sys.time() - sits_mixture_model_start, units = "secs")
    cat(sprintf("SITS LSMM process duration (HH:MM): %02d:%02d\n",
                as.integer(sits_mixture_model_time / 3600),
                as.integer((sits_mixture_model_time %% 3600) / 60)))
    cat("Linear Spectral Mixture Model created successfully.\n")
    
    # Step 2.1.3 -- Retrieve the Mixture Model Cube and merge it with the spectral cube
    mm_cube <- sits_cube(
      source     = "BDC",
      collection = "SENTINEL-2-16D",
      bands      = fraction_bands,
      tiles      = tile,
      data_dir   = mixture_path,
      start_date = start_date,
      end_date   = end_date,
      progress   = TRUE
    )
    mm_cube_fraction <- sits_merge(mm_cube, cube)
    
    # ------------------------------------------------------------
    # 2.2 Image Segmentation
    # ------------------------------------------------------------
    
    sits_segment(
      cube      = mm_cube_fraction,
      seg_fn    = sits_snic(
        grid_seeding = grid_seeding,
        spacing      = spacing,
        compactness  = compactness,
        padding      = padding
      ),
      memsize    = memsize,
      multicores = n_cores,
      output_dir = segments_path,
      version    = version
    )
    
    seg_file <- get_latest_segment_file(segments_path, tile)
    if (is.null(seg_file)) stop("no segment file found")
    cat(sprintf("Generated file: %s\n", basename(seg_file)))
    
    # ------------------------------------------------------------
    # 2.3 Segments BBOX vs Tile Feature (BDC grid)
    # ------------------------------------------------------------
    
    pixel_res <- suppressWarnings(max(cube$file_info[[1]]$xres, na.rm = TRUE))
    if (!is.finite(pixel_res)) pixel_res <- 10
    
    check_bbox_location(seg_file, tile, bbox_tol_px * pixel_res)
    
  }, error = function(e) {
    list(status = "error", message = conditionMessage(e))
  })
  
  tile_status[tile]  <- tile_result$status
  tile_message[tile] <- tile_result$message
  
  label <- c(ok = "OK", partial = "WARNING", mismatch = "WARNING", error = "ERROR")[[tile_result$status]]
  cat(sprintf("[%s] Tile %s: %s.\n", label, tile, tile_result$message))
}

# ============================================================
# 3. File Naming Standardization and Final Report
# ============================================================

# Step 3.1 -- Rename segmentation files to standardize the file naming convention
version_date_regex <- paste0(grid_seeding, "-(\\d{4}-\\d{2}-\\d{2})")
arqs <- list.files(segments_path, pattern = version_date_regex, full.names = TRUE)

if (length(arqs) > 0) {
  renamed <- file.rename(arqs, sub(version_date_regex, paste0(grid_seeding, "_\\1"), arqs))
  if (!all(renamed)) {
    warning(sprintf("Could not rename: %s", paste(basename(arqs[!renamed]), collapse = ", ")))
  }
}

# Step 3.2 -- Report only tiles with errors or problems
problems <- tile_status != "ok"
if (any(problems)) {
  cat("\nTiles with problems:\n")
  cat(sprintf("  - Tile %s [%s]: %s\n", tiles[problems], tile_status[problems], tile_message[problems]),
      sep = "")
}

print("Segmentation complete.")
