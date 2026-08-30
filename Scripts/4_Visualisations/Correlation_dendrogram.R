###############################################
# Load libraries -----------------------------------------------------------
################################################
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ComplexHeatmap"))

library(terra)
library(here)
library(ComplexHeatmap)
library(circlize)
library(grid)

###############################################
# Load variables -----------------------------------------------------------
################################################
inpath <- here("Outputs/NonClimate_Vars")

# list all raster files
raster_files <- list.files(inpath, full.names = TRUE)

# split into two groups based on filename
icefree_files   <- raster_files[grepl("ICEFREE",   raster_files, ignore.case = TRUE)]
coastline_files <- raster_files[grepl("COASTLINE", raster_files, ignore.case = TRUE)]


###############################################
# Helper function: raster group -> correlation matrix ----------------------
################################################
# n_sample = number of random locations to extract values from

compute_cor_matrix <- function(file_list, n_sample = 50000, seed = 123) {
  
  r_stack <- rast(file_list)
  
  # clean names: strip path/extension, keep something readable for plotting
  names(r_stack) <- tools::file_path_sans_ext(basename(file_list))
  
  set.seed(seed)
  
  # sample random cells (non-NA) across the stack
  sampled_vals <- spatSample(r_stack,
                             size = n_sample,
                             method = "random",
                             na.rm = TRUE,
                             as.df = TRUE)
  
  # Pearson correlation, pairwise complete in case some layers have
  # differing NA footprints not fully removed by na.rm above
  cor_mat <- cor(sampled_vals, use = "pairwise.complete.obs", method = "pearson")
  
  return(cor_mat)
}

###############################################
# Compute correlation matrices ----------------------------------------------
################################################
cor_icefree   <- compute_cor_matrix(icefree_files)
cor_coastline <- compute_cor_matrix(coastline_files)

###############################################
# Helper function: replicate superheat-style circle correlation plot -------
################################################
# Clustering is done on Euclidean distance between correlation *profiles*
# (i.e. dist(cor_mat)), matching "hierarchical clustering based on their
# correlations" as described in the source figure caption.

plot_cor_heatmap <- function(cor_mat, title = "") {
  
  col_fun <- colorRamp2(c(-1, 0, 1), c("#67001f", "white", "#053061"))
  
  # Euclidean distance between correlation profiles (rows of cor_mat)
  d <- dist(cor_mat, method = "euclidean")
  hc <- hclust(d, method = "complete")
  
  ht <- Heatmap(
    cor_mat,
    name = "Pearson r",
    col = col_fun,
    
    cluster_rows = hc,
    cluster_columns = hc,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    
    rect_gp = gpar(type = "none"),  # suppress default cell background
    cell_fun = function(j, i, x, y, width, height, fill) {
      r_val <- cor_mat[i, j]
      radius <- abs(r_val) / 2 * min(unit.c(width, height))
      grid.circle(
        x = x, y = y, r = radius,
        gp = gpar(fill = col_fun(r_val), col = NA)
      )
    },
    
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90,
    
    heatmap_legend_param = list(
      at = c(-1, -0.5, 0, 0.5, 1),
      title = "r"
    ),
    
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  return(ht)
}

###############################################
# Build the two plots --------------------------------------------------------
################################################
ht_icefree   <- plot_cor_heatmap(cor_icefree,   title = "ICEFREE variables")
ht_coastline <- plot_cor_heatmap(cor_coastline, title = "COASTLINE variables")

###############################################
# Draw / save to file --------------------------------------------------------
################################################
outpath <- here("Plots")
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

png(file.path(outpath, "corr_heatmap_ICEFREE.png"),
    width = 10, height = 10, units = "in", res = 300)
draw(ht_icefree)
dev.off()

png(file.path(outpath, "corr_heatmap_COASTLINE.png"),
    width = 10, height = 10, units = "in", res = 300)
draw(ht_coastline)
dev.off()

# also draw interactively if running in an R session with a graphics device
# draw(ht_icefree)
# draw(ht_coastline)