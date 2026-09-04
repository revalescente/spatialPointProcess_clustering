library(ggplot2)
library(scales)
library(paletteer)
library(spatstat)
library(SpatialExperiment)

# Use in a ggplot2 chart:
scale_colour_paletteer_d("colorblindr::OkabeIto")
scale_fill_paletteer_d("colorblindr::OkabeIto")
c("#E69F00FF", "#56B4E9FF", "#009E73FF", "#F0E442FF", 
  
  "#0072B2FF", "#D55E00FF", "#CC79A7FF", "#999999FF")

levels(ppp_list$`062921_D0_m3a_0_1`$marks)
levels(ppp_list$`082421_D21_m2_18_1`$marks)

# scale colors manual
my_colors <- c("#E69F00FF", "#56B4E9FF", "#009E73FF", "#F0E442FF", 
               "#0072B2FF", "#D55E00FF", "#CC79A7FF", "#999999FF")

cell_types <- c("Endothelial", "EntericNervous", "Epithelial", 
                "Fibroblast", "ICC", "Immune", 
                "Smooth Muscle Cells", "Adipose") # Added "Unknown" to match 9 colors

# 2. Name the vector
names(my_colors) <- cell_types

# spatial plot of points
points_plot <- \(ppp) {
  ggplot() +
    geom_raster(
      data = as.data.frame(Window(ppp)),
      aes(x = x, y = y),
      fill = "grey92",
      alpha = 1
    ) +
    geom_point(
      data = as.data.frame(ppp),
      aes(x = x, y = y, color = marks),
      size = 0.15,
      alpha = 0.75
    ) +
    coord_fixed(expand = FALSE) +
    scale_x_reverse() +
    scale_y_reverse() +
    #scale_colour_paletteer_d("colorblindr::OkabeIto", name = "Cell types") +
    scale_color_manual(values = my_colors,
                       name = "Cell types") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      #plot.title = element_text(size = 22, face = "bold"),
      #plot.subtitle = element_text(size = 14),
      #axis.title = element_text(size = 16),
      #axis.text  = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 36, face = "bold"),
      legend.text  = element_text(size = 26),
      legend.key.height = unit(10, "mm"),
      legend.key.width  = unit(10, "mm"),
      plot.margin = margin(8, 12, 8, 8, "mm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    labs(
      title = "",
      subtitle = "",
      x = "X coordinate",
      y = "Y coordinate"
    )
}
p1 <- points_plot(ppp_list$`062921_D0_m3a_0_1`)
p3 <- points_plot(ppp_list$`082421_D21_m2_18_1`)


# highlight 2 cell thypes
points_plot_highlight <- \(ppp) {
  df_mask <- as.data.frame(Window(ppp))
  df_pts  <- as.data.frame(ppp)
  
  # keep original labels, but collapse everything else to "Other"
  df_pts$highlight <- as.character(df_pts$marks)
  df_pts$highlight[!df_pts$highlight %in% c("Epithelial", "Fibroblast")] <- "Other"
  df_pts$highlight <- factor(df_pts$highlight, levels = c("Epithelial", "Fibroblast", "Other"))
  
  ggplot() +
    geom_raster(
      data = df_mask,
      aes(x = x, y = y),
      fill = "grey92",
      alpha = 1
    ) +
    geom_point(
      data = df_pts,
      aes(x = x, y = y, color = highlight),
      size = 0.15,
      alpha = 0.75
    ) +
    coord_fixed(expand = FALSE) +
    scale_x_reverse() +
    scale_y_reverse() +
    scale_color_manual(
      values = c(
        "Epithelial" = "#009E73FF", # blue
        "Fibroblast" = "#F0E442FF", # orange
        "Other"      = "lightgrey"
      ),
      breaks = c("Epithelial", "Fibroblast", "Other"),
      name = "Cell types"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      #plot.title = element_text(size = 22, face = "bold"),
      #plot.subtitle = element_text(size = 14),
      #axis.title = element_text(size = 16),
      #axis.text  = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 36, face = "bold"),
      legend.text  = element_text(size = 26),
      legend.key.height = unit(10, "mm"),
      legend.key.width  = unit(10, "mm"),
      plot.margin = margin(8, 12, 8, 8, "mm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    labs(
      title = "",
      x = "X coordinate",
      y = "Y coordinate"
    )
}
p2=points_plot_highlight(ppp_list$`062921_D0_m3a_0_1`)
p4=points_plot_highlight(ppp_list$`082421_D21_m2_18_1`)

# Lfuns of selected samples
plot_fv_list <- function(list_of_funcs, spe, correction = "border") {
  
  # 1. Extract unique mapping of unique_id to sample_type from the SPE object
  meta_map <- as.data.frame(colData(spe)) |> 
    dplyr::select(unique_id, sample_type) |> 
    dplyr::distinct()
  
  # 2. Run the map_dfr block dynamically on whatever list you provide
  df_lcross <- map_dfr(names(list_of_funcs), function(slide_name) {
    
    # Extract the specific L-function object from the provided list
    L_obj <- list_of_funcs[[slide_name]]
    
    as.data.frame(L_obj) |> 
      mutate(
        Slide = slide_name,
        # Use .data[[correction]] to dynamically grab the 'border' column (or 'iso', etc.)
        L_r_minus_r = .data[[correction]] - r,  
        
        # Extract the marks names exactly as you designed
        mark_i = trimws(strsplit(gsub("list\\(|\\)", "", attr(L_obj, "fname")[2]), ",")[[1]][1]),
        mark_j = trimws(strsplit(gsub("list\\(|\\)", "", attr(L_obj, "fname")[2]), ",")[[1]][2])
      )
  })
  
  # 3. Join the sample_type metadata to our plot dataframe
  df_lcross <- df_lcross |> 
    left_join(meta_map, by = c("Slide" = "unique_id"))
  
  # 4. Plot with color mapped to sample_type
  p <- ggplot(df_lcross, aes(x = r, y = L_r_minus_r, group = Slide, color = sample_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(alpha = 0.7, linewidth = 1.5) + 
    theme_minimal() +
    labs(
      title = "",
      subtitle = "",
      x = "Radius (r)",
      y = "L(r) - r",
      color = "Sample Type"
    ) +
    theme(legend.position = "right",
          legend.title = element_text(size = 36, face = "bold"),
          legend.text  = element_text(size = 26),
          axis.title = element_text(size = 26),
          axis.text  = element_text(size = 18),
          ) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1, linewidth = 4), ncol = 1))
  
  return(p)
}
(p5 <- plot_fv_list(list_of_funcs = Lcross_FibEp, spe = spe, correction = "border"))


# fPC scores
plot_fpca_scores <- function(fpca_obj, spe, ns, sz) {
  
  # 1. Extract metadata from SPE
  meta_map <- as.data.frame(colData(spe)) |> 
    dplyr::select(unique_id, sample_type)
  
  # 2. Build dataframe for ggplot
  # Assuming fpca_obj$scores is a matrix where rownames match unique_id
  df_scores <- data.frame(
    PC1 = fpca_obj$scores[, 1],
    PC2 = fpca_obj$scores[, 2],
    unique_id = rownames(fpca_obj$scores)
  ) |> 
    left_join(meta_map, by = "unique_id")
  
  # 4. Generate the plot
  p <- ggplot(df_scores, aes(x = PC1, y = PC2, color = sample_type)) +
    geom_point(shape = ns, size = sz) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    theme_minimal() +
    labs(
      title = "",
      x = "PC1",
      y = "PC2",
      color = "Sample Type"
    ) +
    theme(
      legend.position = "none",
      legend.title = element_text(size = 36, face = "bold"),
      legend.text  = element_text(size = 26),
      #axis.line = element_line(color = "black", linewidth = 1)
    ) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1), ncol = 1))
  
  return(p)
}

p6 <- plot_fpca_scores(fpca, spe, ns = 19, sz = 6)
p6
ggsave(plot = p6, filename = "fpca_scores.png", dpi = 600)
saveRDS(sub_ppp, file = "~/repositories/spatialPointProcess_clustering/dataset/Merfish_dataset/ppp_davide_subset.rds")

# TCGA-OV plots

# spatial pp
points_plot_highlight_tcga <- \(ppp) {
  # select the merged column
  marks(ppp) <- marks(ppp)[["merged_type"]]
  
  df_mask <- as.data.frame(Window(ppp))
  df_pts  <- as.data.frame(ppp)

  # keep original labels, but collapse everything else to "Other"
  df_pts$highlight <- as.character(df_pts$marks)
  df_pts$highlight[!df_pts$highlight %in% c("neoplastic", "stromal_benign")] <- "Other"
  df_pts$highlight <- factor(df_pts$highlight, levels = c("neoplastic", "stromal_benign", "Other"))
  
  ggplot() +
    geom_raster(
      data = df_mask,
      aes(x = x, y = y),
      fill = "grey92",
      alpha = 1
    ) +
    geom_point(
      data = df_pts,
      aes(x = x, y = y, color = highlight),
      size = 0.15,
      alpha = 0.25
    ) +
    coord_fixed(expand = FALSE) +
    scale_x_reverse() +
    scale_y_reverse() +
    scale_color_manual(
      values = c(
        "neoplastic" = "#009E73FF", # blue
        "stromal_benign" = "#F0E442FF", # orange
        "Other"      = "lightgrey"
      ),
      breaks = c("neoplastic", "stromal_benign", "Other"),
      name = "Cell types"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      #plot.title = element_text(size = 22, face = "bold"),
      #plot.subtitle = element_text(size = 14),
      #axis.title = element_text(size = 16),
      #axis.text  = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 36, face = "bold"),
      legend.text  = element_text(size = 26),
      legend.key.height = unit(10, "mm"),
      legend.key.width  = unit(10, "mm"),
      plot.margin = margin(8, 12, 8, 8, "mm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    labs(
      title = "",
      x = "X coordinate",
      y = "Y coordinate"
    )
}

p7 <- points_plot_highlight_tcga(ppp_list[[10]])
p8 <- points_plot_highlight_tcga(ppp_list[[83]])
p7
p8

