# ============================================================
# Linear Spectral Mixture Model (LSMM) for Feature Extraction
# ============================================================

# Load Required Libraries
library(sits)
library(tibble)

# Paths for files and folders
mixture_path <- "data/raw/mixture_model"
tiles      = c('012014', '012015', '013014', '013015')
start_date = '2024-08-01'
end_date   = '2025-07-31'

# ============================================================
# 1. Define and Load Raster Data Cubes from a collection
# ============================================================

cube <- sits_cube(
  source     = "BDC",
  collection = "SENTINEL-2-16D",
  bands      = c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12', 'CLOUD'),
  tiles      = tiles,
  start_date = start_date,
  end_date   = end_date,
  progress   = TRUE)

# ============================================================
# 2. Create Fraction Image Features from the Mixture Model Cube
# ============================================================

# Step 2.1 -- Define endmember values (class names cannot contain numbers)
endmembers <- tibble::tribble(
  ~class,    ~B02, ~B03, ~B04, ~B05, ~B06, ~B07, ~B08, ~B8A, ~B11, ~B12,
  "soil",    1799, 2154, 3028, 3303, 3472, 3656, 3566, 3686, 5097, 4736,
  "veg",      827,  892,  410, 1070, 4206, 5646, 5495, 6236, 2101,  775,
  "water",    946,  739,  280,  208,  180,  167,  135,  129,   26,   14,
)

# Step 2.2 -- Generate the mixture model cube and calculate the process duration
sits_mixture_model_start <- Sys.time()
mm_cube <- sits_mixture_model(
  data       = cube,
  endmembers = endmembers,
  multicores = 28,  # adapt to your computer CPU core availability
  memsize    = 180, # adapt to your computer memory availability
  output_dir = mixture_path
)
sits_mixture_model_end  <- Sys.time()
sits_mixture_model_time <- as.numeric(sits_mixture_model_end - sits_mixture_model_start, units = "secs")
sprintf("SITS LSMM process duration (HH:MM): %02d:%02d", as.integer(sits_mixture_model_time / 3600), as.integer((sits_mixture_model_time %% 3600) / 60))
print("Linear Spectral Mixture Model created successfully.")
