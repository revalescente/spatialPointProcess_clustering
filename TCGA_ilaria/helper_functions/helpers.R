library(anndataR)
library(SpatialExperiment)

# Function to transform from anndata to spatialExperiment
ann2spe <- \(data_path, 
             img_path = NULL, 
             spat_coords_name = c("x", "y"), 
             mpp_value = NULL, ...) {
  
  h5_file <- file.path(data_path, sub("\\.svs$", ".h5ad", sample))
  img_file <- file.path(img_path, sub("\\.svs$", ".png", sample)) 
  # anndata conversion
  ann <- read_h5ad(h5_file)
  # object to build SpatialExperiment
  counts <- t(ann$X) # transpose it
  cd <- ann$obs
  meta <- ann$uns
  spat <- ann$obsm$spatial
  colnames(spat) <- spat_coords_name
  cd <- cbind(cd, spat)
  spe <- SpatialExperiment(assays = base::list("counts" = counts), 
                                              colData = cd,
                                              metadata = meta,
                                              spatialCoordsNames = spat_coords_name,
                                              imageSources = img_file, ...)
  metadata(spe)$MPP <- mpp_value
  spe
}

spe_list <- map(names(mpp_values), \(name){
  current_mpp <- mpp_values[name]
  ann2spe()
})