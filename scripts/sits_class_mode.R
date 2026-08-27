library(terra)
library(sf)
library(exactextractr)

class_sits_path <- "~/grupos/biomasbr-amazonia/sits-prodes/prodes.amz/data/class/012014/original_class/random_forest/2y-after-apocalypse-agrupado/SENTINEL-2_MSI_012014_2023-07-28_2025-07-28_probs_rf-2y-after-apocalypse-agrupado.tif"
segs_path <- "~/grupos/biomasbr-amazonia/sits-prodes/prodes.amz/data/segments/SENTINEL-2_MSI_012014_2025-07-12_2025-07-28_segments_lsmm-snic-spac10-comp03-pad0-rectangular_2026-06-25.gpkg"

class_sits <- rast(class_sits_path)
segs <- read_sf(segs_path)

segs <- st_transform(segs, crs(class_sits))

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

nomes <- c("pol_id",
           "Corpo_Dagua",
           "Corte_Raso_Com_Vegetacao",
           "Corte_Raso",
           "Degradacao",
           "Floresta",
           "Vegetacao_Natural_Nao_Florestal")

# if (grepl("class", basename(class_sits_path))){
# class_seg <- exact_extract(class_sits,
#                            segs,
#                            fun = "mode",
#                            append_cols = "pol_id",
#                            force_df = TRUE)
# } else {
class_seg <- exact_extract(class_sits,
                           segs,
                           fun = "mean",
                           stack_apply = TRUE,
                           append_cols = "pol_id",
                           force_df = TRUE)

names(class_seg) <- nomes
class_seg$class <- colnames(class_seg)[max.col(class_seg[, -1])+1]
#head(class_seg)

# config     <- read_class_config(file.path(config_dir, "class_config.txt"))
# class_translation <- config$class_translation

# class_seg <- class_seg |> 
#   mutate(mode = nomes[mode]) |>
#   mutate(mode = class_translation[mode]) |> 
#   rename(class = mode)

class_seg <- dplyr::left_join(segs, class_seg, by = "pol_id")

sf::st_write(
  class_seg,
  paste0(tools::file_path_sans_ext(class_sits_path),
         "_mean",
         ".gpkg"))