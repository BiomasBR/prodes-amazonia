library(rstac)
library(terra)
library(lubridate)

tiles <- c('024005','024006','025001','025002','025003','025004','025005','025006','026005','026006','027006','028006','029006','029007','030006','031007','031008','032007','032008')

ano <- '2026' #formato YYY-MM-DD

upper <- -9
lower <- -20
method <- "near"

ref_date <- as.Date(paste0(ano, "-07-01"))

data_inicio_4m <- as.character(ref_date %m-% months(3))
data_inicio_8m <- as.character(ref_date %m-% months(6))



crs <- "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"Unknown based on GRS80 ellipsoid\",\n            ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,\n                LENGTHUNIT[\"metre\",1],\n                ID[\"EPSG\",7019]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Albers Equal Area\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",-12,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-54,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",-2,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",-22,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",5000000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",10000000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]"

#===============================================================================
# Download Function
#===============================================================================

build_vh_mosaic <- function(
    target_date,
    bbox,
    search_window = 30){
  
  message("----------------------------------------")
  message("Searching images from ", target_date)
  
  date_search <- paste(
    target_date,
    as.character(as.Date(target_date) %m+%
                   lubridate::days(search_window)),
    sep = "/"
  )
  
  items <- rstac::stac("https://data.inpe.br/bdc/stac/v1/") |>
    rstac::stac_search(
      collections = "sentinel-1-rtc-1",
      bbox = as.numeric(bbox),
      datetime = date_search) |>
    rstac::ext_query(orbit_direction == "DESCENDING") |>
    rstac::post_request() |>
    rstac::items_fetch()
  
  if(length(items$features) == 0){
    stop("No Sentinel-1 images found.")
  }
  
  raw_datetimes <- sapply(items$features, function(x){
    
    if(!is.null(x$properties$datetime)){
      return(x$properties$datetime)
    }
    x$datetime
  })
  
  calendar_dates <- as.Date(raw_datetimes)
  
  unique_dates <- sort(unique(calendar_dates))
  
  message(length(unique_dates),
          " acquisition dates found.")
  
  message(length(items$features),
          " images found.")
  
  tile_ext <- vect(sf::st_as_sfc(bbox))
  
  ## ---------------------------------------------------------------------
  ## 1) Crop de TODAS as cenas (todas as datas),
  ##    junto com a órbita relativa de cada cena ("Relative Orbit" no STAC).
  ##    Isso precisa ser feito para todo o conjunto ANTES do filtro por
  ##    órbita, já que o "outcore" de cada órbita é uma média temporal
  ##    calculada sobre todas as cenas daquela órbita, não só as de um dia.
  ## ---------------------------------------------------------------------
  
  n_feat <- length(items$features)
  
  scene_names   <- character(n_feat)
  scene_orbits  <- character(n_feat)
  ratio_rasters <- vector("list", n_feat)
  
  for(i in seq_len(n_feat)){
    
    feature <- items$features[[i]]
    
    url_vh <- feature$assets$Gamma0_VH$href
    
    scene_names[i]  <- tools::file_path_sans_ext(basename(url_vh))
    scene_orbits[i] <- as.character(feature$properties$relative_orbit)
    
    r <- terra::rast(paste0("/vsicurl/", url_vh))
    
    message("Cropping image: ", scene_names[i])
    
    r_cropped <- terra::crop(r, tile_ext)
    
    if(!exists("grid_template")){
      grid_template <- terra::rast(
        tile_ext,
        resolution = terra::res(r_cropped),
        crs = terra::crs(r_cropped)
      )
    }
    
    if(!isTRUE(terra::compareGeom(r_cropped, grid_template, stopOnError = FALSE))){
      r_cropped <- terra::resample(r_cropped, grid_template, method = "bilinear")
    }
    
    names(r_cropped) <- ("VH")
    
    ratio_rasters[[i]] <- r_cropped
    
    rm(r, r_cropped)
    gc()
  }
  
  ## ---------------------------------------------------------------------
  ## 2) Filtro por órbita (equivalente a getFilteredColbyOrbit do GEE)
  ## ---------------------------------------------------------------------
  message("Applying orbit filter...")
  
  filtered_rasters <- apply_filter_by_orbit(ratio_rasters, scene_orbits)
  gc()
  
  ## ---------------------------------------------------------------------
  ## 3) Para CADA data única, faz o merge das cenas daquele dia -> 1 tile
  ##    por data (mesma lógica de antes: "imagem de data 1", "data 2" ...
  ##    "data n"). Nenhuma data é descartada e não há mais lógica de
  ##    cobertura/early-stop: todas as datas encontradas na janela de busca
  ##    entram no cálculo.
  ## ---------------------------------------------------------------------
  
  date_rasters <- vector("list", length(unique_dates))
  
  for(k in seq_along(unique_dates)){
    
    current_date <- unique_dates[k]
    
    message("----------------------------------------")
    message("Processing date: ", as.character(current_date))
    
    idx <- which(calendar_dates == current_date)
    
    tile_list <- filtered_rasters[idx]
    
    ## Apenas uma cena neste dia ##
    
    if(length(tile_list) == 1){
      day_raster <- tile_list[[1]]
    } else{
      
      ## Merge apenas das cenas do mesmo dia -> forma o tile daquela data ##
      
      message("Merging ", length(tile_list), " same-day scenes.")
      day_raster <- terra::merge(terra::sprc(tile_list))
    }
    
    date_rasters[[k]] <- day_raster
    gc()
  }
  
  ## ---------------------------------------------------------------------
  ## 4) Raster final: MÉDIA de todas as imagens de data (date 1, date 2,
  ##    ..., date n), equivalente ao .mean() do GEE sobre a ImageCollection
  ##    do período.
  ## ---------------------------------------------------------------------
  
  message("Averaging ", length(date_rasters), " date images...")
  
  final_raster <- terra::mosaic(
    terra::sprc(date_rasters),
    fun = "mean"
  )
  gc()
  return(
    list(
      raster = final_raster,
      dates = as.character(unique_dates),
      scenes = scene_names,
      target_date = target_date,
      n_dates = length(unique_dates),
      n_scenes = length(scene_names)
    )
  )
}

#===============================================================================
# Auxiliary Functions
#===============================================================================

apply_boxcar <- function(r){
  boxcar_window <- matrix(1, nrow = 3, ncol = 3)
  
  terra::focal(r, w = boxcar_window,
               fun = mean,
               na.rm = TRUE)
}

# Equivalente a getFilteredColbyOrbit do GEE:
# recebe uma lista de rasters (uma por cena) e um vetor de mesma ordem/tamanho
# com o número da órbita relativa de cada raster ("Relative Orbit" no STAC),
# e devolve a lista de rasters filtradas, na MESMA ordem de entrada.
apply_filter_by_orbit <- function(rasters, orbits){
  
  filtered <- vector("list", length(rasters))
  unique_orbits <- unique(orbits)
  
  for(orb in unique_orbits){
    
    idx <- which(orbits == orb)
    input_col <- rasters[idx]
    
    ref <- input_col[[1]]
    input_col <- lapply(input_col, function(img){
      if(!isTRUE(terra::compareGeom(img, ref, stopOnError = FALSE))){
        img <- terra::resample(img, ref, method = "bilinear")
      }
      img
    })
    
    # versão suavizada (boxcar) de cada imagem da órbita
    smooth_col <- lapply(input_col, apply_boxcar)
    
    # razão imagem_original / imagem_suavizada, para cada imagem da órbita
    ratio_col <- mapply(
      function(img, sm) img / sm,
      input_col, smooth_col,
      SIMPLIFY = FALSE
    )
    
    # média pixel-a-pixel dessas razões (equivalente a sum()/count() no GEE)
    outcore <- terra::app(terra::rast(ratio_col), fun = mean, na.rm = TRUE)
    
    # multiplica cada imagem já suavizada pelo "outcore"
    filtered[idx] <- lapply(smooth_col, function(sm) sm * outcore)
    rm(input_col, smooth_col, ratio_col, outcore)
    gc()
  }
  filtered
}

to_db <- function(r){
  10 * log10(r)
}

stretch_vh <- function(r, upper, lower){

  r <- terra::clamp(r,
                    lower = lower,
                    upper = upper
                   )
  
  ((r - lower) / (upper - lower)) * 254
}

# Função auxiliar: gera o bloco de texto de um período
make_period_report <- function(title, mosaic) {
  c(
    sprintf("---------------- %s ----------------", title),
    paste("Reference date:", mosaic$target_date),
    paste("Dates used:", paste(mosaic$dates, collapse = ", ")),
    paste("Number of dates:", mosaic$n_dates),
    paste("Number of scenes:", mosaic$n_scenes),
    "",
    "Scenes:",
    paste("  -", mosaic$scenes),
    ""
  )
}


#===============================================================================
# BUILD THE THREE TEMPORAL MOSAICS
#===============================================================================

results <- vector("list", length(tiles))
names(results) <- tiles

ind <- 0
for (tile in tiles) {
  
  ind <- ind + 1
  
  results[tile] <- list(
    tryCatch(
      {
        cat("Starting ", tile, "\n",
            "Scene ", ind, "/", length(tiles), "\n")
        
        class_dir <- "data/class"
        
        radar_dir <- file.path(class_dir, tile, "radar_image")
        dir.create(radar_dir, showWarnings = FALSE, recursive = TRUE)
        
        bbox <- sits::sits_tiles_to_roi(tiles = tile,
                                        grid_system = "BDC_SM_V2")
        
        # ---------------- CURRENT ----------------
        message("Current image")
        vh_atual <- build_vh_mosaic(
          target_date = as.character(ref_date),
          bbox = bbox
        )
        
        if (is.null(vh_atual))
          stop("Current mosaic could not be created.")
        
        vh_atual_img <- to_db(vh_atual$raster)
        report_atual <- make_period_report("CURRENT", vh_atual)
        data_ref <- base::sort(vh_atual$dates, decreasing = TRUE)[1]
        
        rm(vh_atual)
        
        vh_atual_img <- stretch_vh(vh_atual_img, upper, lower)
        vh_atual_10m <- project(vh_atual_img,
                                crs,
                                res = 10,
                                method = method,
                                threads = TRUE)
        rm(vh_atual_img)
        gc()
        
        # ---------------- 4 MONTHS ----------------
        message("4 months")
        vh_4m <- build_vh_mosaic(
          target_date = data_inicio_4m,
          bbox = bbox
        )
        
        if (is.null(vh_4m))
          stop("4-month mosaic could not be created.")
        
        vh_4m_img <- to_db(vh_4m$raster)
        report_4m <- make_period_report("4 MONTHS", vh_4m)
        
        rm(vh_4m)
        
        vh_4m_img <- stretch_vh(vh_4m_img, upper, lower)
        vh_4m_10m <- project(vh_4m_img,
                             vh_atual_10m,
                             method = method,
                             threads = TRUE)
        rm(vh_4m_img)
        gc()
        
        # ---------------- 8 MONTHS ----------------
        message("8 months")
        vh_8m <- build_vh_mosaic(
          target_date = data_inicio_8m,
          bbox = bbox
        )
        
        if (is.null(vh_8m))
          stop("8-month mosaic could not be created.")
        
        vh_8m_img <- to_db(vh_8m$raster)
        report_8m <- make_period_report("8 MONTHS", vh_8m)
        
        rm(vh_8m)
        gc()
        
        vh_8m_img <- stretch_vh(vh_8m_img, upper, lower)
        vh_8m_10m <- project(vh_8m_img,
                             vh_atual_10m,
                             method = method,
                             threads = TRUE)
        rm(vh_8m_img)
        gc()
        
        #=======================================================================
        # STACK RGB
        #=======================================================================
        
        temporal_rgb_10m <- c(vh_atual_10m,
                              vh_4m_10m,
                              vh_8m_10m)
        
        rm(vh_atual_10m,
           vh_4m_10m,
           vh_8m_10m)
        gc()
        
        temporal_rgb_10m <- project(temporal_rgb_10m,
                                    "EPSG:4674",
                                    method = method,
                                    threads = TRUE)
        
        names(temporal_rgb_10m) <- c(
          "VH_Current",
          "VH_4Months",
          "VH_8Months"
        )
        
        #=======================================================================
        # EXPORT
        #=======================================================================
        
        name_export <- sprintf("%s/S1_%s_VH_R07G04B01_mean-temporal_%s.tif",
                               radar_dir, tile, data_ref)
        
        terra::writeRaster(temporal_rgb_10m,
                           filename = name_export,
                           datatype = "INT1U",
                           overwrite = TRUE,
                           gdal = c("COMPRESS=DEFLATE",
                                    "PREDICTOR=2",
                                    "ZLEVEL=9"),
                           progress = TRUE)
        
        rm(temporal_rgb_10m)
        gc()
        
        #=======================================================================
        # Sentinel-1 Temporal RGB Report
        #=======================================================================
        
        report_header <- c(
          "===================================================",
          "Sentinel-1 Temporal RGB Report",
          "===================================================",
          "",
          paste("Tile:", tile),
          paste("Processing date:", Sys.time()),
          "",
          ""
        )
        
        report_final <- c(
          report_header,
          report_atual,
          report_4m,
          report_8m
        )
        
        # somente o report final é salvo
        writeLines(
          report_final,
          sprintf("%s/S1_%s_VH_R07G04B01_mean-temporal_%s_report.txt",
                  radar_dir, tile, ano)
        )
        
        rm(report_atual, report_4m, report_8m, report_header, report_final)
        gc()
        
        cat(tile, " Finished\n")
        
        # valor de sucesso (não-NULL) para o resumo final
        name_export
      },
      error = function(e) {
        message("ERROR in tile ", tile, ": ", conditionMessage(e))
        NULL
      },
      finally = {
        gc()
      }
    )
  )
}

# ============================================================
# 6. Final summary
# ============================================================

success <- names(results)[!vapply(results, is.null, logical(1))]
failure <- names(results)[vapply(results, is.null, logical(1))]

message("\n========== PROCESSING SUMMARY ==========")
message("Total tiles: ", length(tiles))
message("Success (", length(success), "): ", paste(success, collapse = ", "))
if (length(failure) > 0) {
  message("Failure (", length(failure), "): ", paste(failure, collapse = ", "))
} else {
  message("No faults recorded.")
}
