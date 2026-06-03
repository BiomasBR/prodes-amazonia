# ============================================================
#  Creation of uncertainty samples
# ============================================================

# Load required libraries
library(sits)
library(sf)
library(purrr)

# Define the parameters: These are user-defined variables
n_samples   <- 100
min_uncert  <- 0.7
tiles       <- c("012014", "012015", "013014", "013015")
version     <- "rf-1y-all-classes"

# Date and time of the start of processing
date_process    <- format(Sys.Date(), "%Y-%m-%d")

# File and folder paths
samples_dir       <- "data/raw/samples/uncertainty_samples"
uncertainty_dir   <- "data/class"
pattern           <- sprintf(".*_(%s)_", paste(tiles, collapse = "|"))
seg_path          <- list.files("data/segments",
                                 pattern = pattern,
                                 full.names = TRUE)

# ============================================================
# 1. SITS Cube
# ============================================================
cube_dirs <- list.dirs(uncertainty_dir, recursive = TRUE)

cube_dirs <- cube_dirs[
  sapply(cube_dirs, function(x) {
    files <- list.files(x, pattern = paste0(pattern,
                                            ".*_entropy",
                                            ".*\\.tif$"),
                        full.names = TRUE)
    return(length(files) > 0)
  })
]

cube_list <- map(cube_dirs, function(dir) {
  sits_cube(
    source      = "BDC",
    collection  = "SENTINEL-2-16D",
    bands       = "entropy",
    data_dir    = dir, # Takes one path from 'cube_dirs' at a time
    version     = version,
    parse_info  = c("satellite", "sensor", "tile", "start_date", "end_date", 
                    "band", "version")
  )
})

# 1.2 -- Combine the list of tibbles into a single multi-row sits cube
uncert_cube <- do.call(rbind, cube_list)
                
# ============================================================
# 2. Generate Uncertainty points samples
# ============================================================
uncert_sf <- sits_uncertainty_sampling(
  uncert_cube,
  n = n_samples,
  min_uncert = min_uncert,
  sampling_window = 10L,
  multicores = 28,
  memsize = 180
)

uncert_sf <- st_as_sf(uncert_sf,
                       coords = c("longitude", "latitude"),
                       crs = "EPSG:4674") 

uncert_sf_file_path <- file.path(samples_dir,
                                   paste0("uncertainty-samples_all-classes_",
                                          paste(uncert_cube$tile, collapse = "-"),
                                          "_", version, "_",
                                          date_process, ".gpkg"))

sf::st_write(uncert_sf,
             uncert_sf_file_path,
             delete_dsn = TRUE, append = FALSE)

# 2 -- Polygons Samples
polygons <- seg_path |> 
  map(read_sf) |> 
  bind_rows()

uncertainty_polygons <- st_filter(polygons,
                                  st_transform(uncert_sf,
                                               crs(polygons))
                                  )

polygons_sf_file_path <- file.path(samples_dir,
                                   paste0("uncertainty-samples-polygons_all-classes_",
                                          paste(uncert_cube$tile, collapse = "-")
                                          ,"_", version, "_",
                                          date_process, ".gpkg"))

sf::st_write(uncertainty_polygons,
             polygons_sf_file_path,
             delete_dsn = TRUE, append = FALSE)