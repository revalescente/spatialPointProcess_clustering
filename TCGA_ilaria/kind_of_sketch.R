sketch_spatial_R <- function(ppp_obj, idx, n_target) {
  if (length(idx) <= n_target) return(idx)
  
  # Extract X and Y coordinates for these specific points
  coords <- cbind(x = ppp_obj$x[idx], y = ppp_obj$y[idx])
  
  # Cluster the physical space into 'n_target' even regions
  # Using MacQueen algorithm because it is blazingly fast for 2D coordinates
  km <- kmeans(coords, centers = n_target, algorithm = "MacQueen", iter.max = 15)
  
  # Pick one random cell from each spatial cluster
  # This forces the downsampled cells to perfectly trace the tissue structure
  sketched_rel_idx <- tapply(seq_along(idx), km$cluster, function(x) sample(x, 1))
  
  return(idx[sketched_rel_idx])
}

prova <- ppp_list$`TCGA-13-A5FU-01Z-00-DX1.9AD9E4B9-3F87-4879-BC0F-148B12C09036`

ind_i <- which(marks(prova)$merged_type == "neoplastic")
ind_j <- which(marks(prova)$merged_type == "stromal_benign")
# sketching
keep_i <- sketch_spatial_R(prova, ind_i, 10000)
keep_j <- sketch_spatial_R(prova, ind_j, 10000)
prova_sub <- prova[c(keep_i, keep_j)]

# downsample
dkep_i <- sample(ind_i, min(length(ind_i), 10000))
dkep_j <- sample(ind_j, min(length(ind_j), 10000))
prova_sub2 <- prova[c(dkep_i, dkep_j)]

p1 <- poinside_plot(prova_sub, title_p = "sketched sample")

p3 <- poinside_plot(prova_sub2, title_p = "downsample")

p1 | p2
p1 | p3
  

sketch_spatial_fast <- function(ppp_obj, idx, n_target) {
  if (length(idx) <= n_target) return(idx)
  
  # Determine a high-resolution grid size (creating ~5x more boxes than n_target)
  # This ensures we capture all spatial nuances before thinning down to exactly n_target
  grid_res <- ceiling((n_target * 5)^2)
  
  # Get coordinates
  x_coords <- ppp_obj$x[idx]
  y_coords <- ppp_obj$y[idx]
  
  # Assign each point to an X and Y grid bin
  x_bin <- findInterval(x_coords, seq(min(x_coords), max(x_coords), length.out = grid_res))
  y_bin <- findInterval(y_coords, seq(min(y_coords), max(y_coords), length.out = grid_res))
  
  # Create a unique ID for each grid square
  grid_id <- paste(x_bin, y_bin, sep = "_")
  
  # Shuffle the data randomly
  shuffled_order <- sample(length(idx))
  idx_shuffled <- idx[shuffled_order]
  grid_id_shuffled <- grid_id[shuffled_order]
  
  # By taking the first occurrence of each grid_id, we get exactly 1 random cell per grid square!
  kept_idx <- idx_shuffled[!duplicated(grid_id_shuffled)]
  
  # If the grid gave us more points than n_target, randomly downsample to exactly n_target
  if (length(kept_idx) > n_target) {
    kept_idx <- sample(kept_idx, n_target)
  }
  
  return(kept_idx)
}

keep2_i <- sketch_spatial_fast(prova, ind_i, 10000)
keep2_j <- sketch_spatial_fast(prova, ind_j, 10000)
prova_sub_fast <- prova[c(keep2_i, keep2_j)]

p4 <- poinside_plot(prova_sub_fast, title_p = "fast sketch")
p1 | p4
p1 | p3

# compare the 3 sampling methods
marks(prova_sub)[,1] <- NULL
marks(prova_sub2)[,1] <- NULL
marks(prova_sub_fast)[,1] <- NULL
L_sub <- Lcross(prova_sub, i = "neoplastic", j = "stromal_benign", correction = "border", r = r_values)
L_sub2 <- Lcross(prova_sub2, i = "neoplastic", j = "stromal_benign", correction = "border", r = r_values)
L_sub_fast <- Lcross(prova_sub_fast, i = "neoplastic", j = "stromal_benign", correction = "border", r = r_values)

lplot <- \(L_obj, titled = "L-cross Function "){
  as.data.frame(L_obj) |> 
    mutate(
      # Use .data[[correction]] to dynamically grab the 'border' column (or 'iso', etc.)
      L_r_minus_r = .data[["border"]] - r,  
      
      # Extract the marks names exactly as you designed
      mark_i = trimws(strsplit(gsub("list\\(|\\)", "", attr(L_obj, "fname")[2]), ",")[[1]][1]),
      mark_j = trimws(strsplit(gsub("list\\(|\\)", "", attr(L_obj, "fname")[2]), ",")[[1]][2])
    ) |> 
    ggplot(aes(x = r, y = L_r_minus_r)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(alpha = 0.7, linewidth = 0.6) + 
    theme_minimal() +
    labs(
      title = titled,
      x = "Radius (r)",
      y = "L(r) - r"
    )
}

lplot(L_sub, "sketching") | lplot(L_sub2, "subsample") | lplot(Lprova, "original?") | lplot(L_sub_fast, "fast sketching")
