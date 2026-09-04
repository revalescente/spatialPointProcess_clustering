library(anndataR)
library(purrr)
library(dplyr)


data_path <- "/mnt/ganimede/shared/imageFeatureTCGA_OV/hovernet/h5ad"
img_path <- "/mnt/ganimede/shared/imageFeatureTCGA_OV/hovernet/thumb"

# Convert h5ad to ppp object

#' Barebone function to read AnnData and create a spatstat ppp object
#'
#' @param data_path Path to the .h5ad file.
#' @param spat_coords_name Names of the spatial coordinates in obsm$spatial (default: c("x", "y")).
#' @param marks_col The column name in `obs` to use as the point marks (default: "type").
#' @return A basic spatstat `ppp` object.
ann2ppp_base <- function(data_path, spat_coords_name = c("x", "y"), marks_col = "type") {
  
  # 1. Read AnnData
  ann <- read_h5ad(data_path)
  
  # 2. Extract spatial coordinates
  spat <- ann$obsm$spatial
  x_coords <- spat[, 1]
  y_coords <- spat[, 2]
  
  # 3. Extract marks (e.g., cell types)
  marks_data <- NULL
  if (marks_col %in% colnames(ann$obs)) {
    marks_data <- as.factor(ann$obs[[marks_col]])
  }
  
  # 4. Create a basic bounding box window 
  # (You can replace this with your custom mask/owin logic in another function)
  basic_win <- spatstat.geom::owin(
    xrange = range(x_coords, na.rm = TRUE), 
    yrange = range(y_coords, na.rm = TRUE)
  )
  
  # 5. Build the PPP object
  p <- spatstat.geom::ppp(
    x = x_coords, 
    y = y_coords, 
    window = basic_win,
    marks = marks_data
  )
  
  return(p)
}

save_all_h5ad_to_ppp <- function(h5ad_files, out_dir = "dataset/samples_ppp") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (sample_name in names(h5ad_files)) {
    file_path <- h5ad_files[[sample_name]]
    mask <- mask_hovernet[[sample_name]]
    
    # Extract metadata
    meta_row <- metadata[metadata$`Image Name` == sample_name, ]
    app_mag <- ifelse(nrow(meta_row) > 0 && !is.na(meta_row$AppMag[1]), as.numeric(meta_row$AppMag[1]), 40)
    mpp     <- ifelse(nrow(meta_row) > 0 && !is.na(meta_row$MPP[1]), as.numeric(meta_row$MPP[1]), 1)
    
    # Convert and save
    my_ppp <- ann2ppp(data_path = file_path, mask = mask, app_mag = app_mag, mpp = mpp)
    saveRDS(my_ppp, file = file.path(out_dir, paste0(sample_name, ".rds")))
  }
}