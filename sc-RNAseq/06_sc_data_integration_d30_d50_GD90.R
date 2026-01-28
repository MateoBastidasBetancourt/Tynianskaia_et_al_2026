library(readr)
library(dplyr)
library(tibble)
library(pheatmap)
library(tidyverse)
library(pheatmap)

marmoset_pseudocounts_day30 <- read_csv("marmoset_pseudocounts_day30.csv")
marmoset_pseudocounts_day50 <- read_csv("marmoset_pseudocounts_day50.csv")
marmoset_pseudocounts_brains <- read_csv("marmoset_pseudocounts_brains.csv")
################ comparing pseudocells matrixes against each other
# marmoset_pseudocounts_day50
# marmoset_pseudocounts_day30
# marmoset_pseudocounts_brains

# ⃣ Set the rownames to pseudocell IDs and remove first column
day50_mat <- marmoset_pseudocounts_day50 %>%
  column_to_rownames(var = "...1")

day30_mat <- marmoset_pseudocounts_day30 %>%
  column_to_rownames(var = "...1")

brain_mat <- marmoset_pseudocounts_brains %>%
  column_to_rownames(var = "...1")

# find common genes
common_genes <- Reduce(intersect, list(
  colnames(day50_mat),
  colnames(day30_mat),
  colnames(brain_mat)
))

#subset to common genes
day50_mat <- day50_mat[, common_genes]
day30_mat <- day30_mat[, common_genes]
brain_mat <- brain_mat[, common_genes]

# log-transform (if these are counts)
day50_mat <- log1p(day50_mat)
day30_mat <- log1p(day30_mat)
brain_mat <- log1p(brain_mat)

####################From here can operate with the matrix 

library(dplyr)
library(tibble)
library(dplyr)
library(tibble)

library(dplyr)
library(tibble)

aggregate_by_celltype <- function(pseudocount_mat, metadata, cluster_col = "cluster") {
  # ensure both objects have matching rownames
  common_cells <- intersect(rownames(pseudocount_mat), rownames(metadata))
  
  if (length(common_cells) == 0) {
    stop("no matching cell IDs found between pseudocount_mat and metadata.")
  }
  
  pseudocount_mat <- pseudocount_mat[common_cells, , drop = FALSE]
  metadata <- metadata[common_cells, , drop = FALSE]
  
  # add cluster info
  merged <- cbind(metadata[, cluster_col, drop = FALSE], pseudocount_mat)
  colnames(merged)[1] <- "cluster"
  # average expression per cluster
  celltype_avg <- merged %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    tibble::column_to_rownames(var = "cluster")
  
  return(celltype_avg)
}



marmoset_pseudocounts_d30orgs_metadata <- read.csv("/mnt/vast-standard/home/cbastid/u15756/marmoset_pseudocounts_d30orgs_metadata.csv", row.names=1)
marmoset_pseudocounts_d50orgs_metadata <- read.csv("/mnt/vast-standard/home/cbastid/u15756/marmoset_pseudocounts_d50orgs_metadata.csv", row.names=1)
marmoset_pseudocounts_brains_metadata <- read.csv("/mnt/vast-standard/home/cbastid/u15756/marmoset_pseudocounts_brains_metadata.csv", row.names=1)

day30_celltypes <- aggregate_by_celltype(day30_mat, marmoset_pseudocounts_d30orgs_metadata)
day50_celltypes <- aggregate_by_celltype(day50_mat, marmoset_pseudocounts_d50orgs_metadata)
brain_celltypes <- aggregate_by_celltype(brain_mat, marmoset_pseudocounts_brains_metadata)
gc()

add_stage_label <- function(df, stage_label) {
  rownames(df) <- paste0(rownames(df), "_", stage_label)
  return(df)
}

day30_celltypes_labeled <- add_stage_label(day30_celltypes, "d30")
day50_celltypes_labeled <- add_stage_label(day50_celltypes, "d50")
brain_celltypes_labeled <- add_stage_label(brain_celltypes, "GD90")

combined_celltypes <- rbind(day30_celltypes_labeled,
                            day50_celltypes_labeled,
                            brain_celltypes_labeled)

# compute correlation matrix
cor_matrix_celltypes <- cor(t(combined_celltypes), use = "pairwise.complete.obs")
# convert correlation matrix to long format
cor_long <- as.data.frame(as.table(cor_matrix_celltypes)) %>%
  rename(Celltype1 = Var1, Celltype2 = Var2, Correlation = Freq) %>%
  mutate(
    Stage1 = sub(".*_(d\\d+|p\\d+)$", "\\1", Celltype1),
    Stage2 = sub(".*_(d\\d+|p\\d+)$", "\\1", Celltype2)
  ) %>%
  # Keep only correlations between different stages (no within-stage)
  filter(Stage1 != Stage2) %>%
  # Remove self-correlations and duplicates (matrix symmetry)
  filter(as.numeric(factor(Celltype1)) < as.numeric(factor(Celltype2))) %>%
  arrange(desc(Correlation))


# srt by correlation strength
top_cross_stage <- cor_long %>%
  arrange(desc(Correlation)) ###top is between d50 and GD90

pheatmap(cor_matrix_celltypes, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
