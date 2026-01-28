library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(FNN)
library(monocle3)
library(SeuratWrappers)
library(edgeR)

setwd("/mnt/vast-standard/home/cbastid/u15756")
seurat.merged.filtered <- readRDS("seurat.merged.filtered_d50orgs.rds")

seurat.merged.filtered[["percent.mt"]] <- PercentageFeatureSet(
  seurat.merged.filtered,
  pattern = "^MT-"
)

seurat.merged.filtered <- subset(
  seurat.merged.unfiltered,
  subset = nFeature_RNA > 200 & 
    nFeature_RNA < 8000 &
    nCount_RNA > 500 & 
    percent.mt < 10
)

seurat.list <- SplitObject(seurat.merged.filtered, split.by = "Sample")

seurat.list <- lapply(seurat.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = seurat.list)

anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)

seurat.integrated <- IntegrateData(anchorset = anchors)

############# clustering
DefaultAssay(seurat.integrated) <- "integrated"
seurat.integrated <- ScaleData(seurat.integrated)
seurat.integrated <- RunPCA(seurat.integrated)
elbow_plot <- ElbowPlot(seurat.integrated, ndims = 50)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:26)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.45)
seurat.integrated$res_0.45 <- Idents(seurat.integrated)

seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.45)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:26)

cluster.markers_names <- FindAllMarkers(seurat.integrated, only.pos = T, group.by = "seurat_clusters_names")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat.integrated <- CellCycleScoring(seurat.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

######cell annotation. 
##Use CTretrieve function and literature datasets from "seurat_integrated_pipeline_d30.R" script for cell annotation
######################
# Generate marker combination tables
markers_combi_10 <- data.frame(
  cluster = seq(0, n_clusters - 1),
  kanton = CTretrieve(cluster.markers, CT_markers_kanton, 10),
  artegiani = CTretrieve(cluster.markers, CT_markers_artegiani, 10),
  linnarsson = CTretrieve(cluster.markers, CT_markers_linnarsson, 10),
  geschwind = CTretrieve(cluster.markers, CT_markers_geschwind, 10)
)

markers_combi_30 <- data.frame(
  cluster = seq(0, n_clusters - 1),
  kanton = CTretrieve(cluster.markers, CT_markers_kanton, 30),
  artegiani = CTretrieve(cluster.markers, CT_markers_artegiani, 30),
  linnarsson = CTretrieve(cluster.markers, CT_markers_linnarsson, 30),
  geschwind = CTretrieve(cluster.markers, CT_markers_geschwind, 30)
)

markers_combi_50 <- data.frame(
  cluster = seq(0, n_clusters - 1),
  kanton = CTretrieve(cluster.markers, CT_markers_kanton, 50),
  artegiani = CTretrieve(cluster.markers, CT_markers_artegiani, 50),
  linnarsson = CTretrieve(cluster.markers, CT_markers_linnarsson, 50),
  geschwind = CTretrieve(cluster.markers, CT_markers_geschwind, 50)
)

################## celltype assignment

seurat.integrated$seurat_clusters <- as.character(seurat.integrated$seurat_clusters)

seurat.integrated$seurat_clusters_names <- ifelse(seurat.integrated$seurat_clusters %in% c("0"), "EN",
                                                  ifelse(seurat.integrated$seurat_clusters %in% c("3"), "EN-DL",
                                                         ifelse(seurat.integrated$seurat_clusters %in% c("1"), "bRG1",
                                                                ifelse(seurat.integrated$seurat_clusters %in% c("2", "9"), "IN",
                                                                       ifelse(seurat.integrated$seurat_clusters %in% c("8"), "NB",
                                                                              ifelse(seurat.integrated$seurat_clusters %in% c("4"), "tRG",
                                                                                     ifelse(seurat.integrated$seurat_clusters %in% c("5"), "IPC",
                                                                                            ifelse(seurat.integrated$seurat_clusters %in% c("6"), "ChP/EN",
                                                                                                   ifelse(seurat.integrated$seurat_clusters %in% c("7"), "aRG",
                                                                                                          ifelse(seurat.integrated$seurat_clusters %in% c("10"), "bRG3",
                                                                                                                 ifelse(seurat.integrated$seurat_clusters %in% c("12"), "bRG2",
                                                                                                                        ifelse(seurat.integrated$seurat_clusters %in% c("13"), "tRG/ChP",
                                                                                                                               ifelse(seurat.integrated$seurat_clusters %in% c("11"), "Astrocyte", NA)))))))))))))

new_order <- c("aRG", "tRG", "tRG/ChP",  "bRG1", "bRG2", "bRG3", "ChP/EN", "IPC", "NB", "EN", "EN-DL", "IN", "Astrocyte")

seurat.integrated$seurat_clusters_names <- factor(
  seurat.integrated$seurat_clusters_names,
  levels = new_order
)

Idents(seurat.integrated) <- "seurat_clusters_names"

############# PALMD pseudobulk analysis of bRGs
meta <- seurat.integrated@meta.data
meta$cluster <- Idents(seurat.integrated)
meta$group <- paste(meta$Sample, meta$cluster, sep = "_")

# get all layers containing raw counts
count_layers <- grep("^counts", Layers(seurat.integrated[["RNA"]]), value = TRUE)
counts_list <- lapply(count_layers, function(x) seurat.integrated[["RNA"]]@layers[[x]])

#extract gene names
gene_names <- rownames(seurat.integrated[["RNA"]]@features)

# assign them to all count layers
for (x in grep("^counts", Layers(seurat.integrated[["RNA"]]), value = TRUE)) {
  mat <- seurat.integrated[["RNA"]]@layers[[x]]
  rownames(mat) <- gene_names
  seurat.integrated[["RNA"]]@layers[[x]] <- mat
}

# combine
counts <- do.call(cbind, counts_list)
rownames(counts) <- gene_names
# transpose matrix to cells x genes, then aggregate by group
counts_t <- t(as.matrix(counts))               # cells x genes
pb_mat <- rowsum(counts_t, group = meta$group) # rows will be sample_cluster
pb_counts <- t(pb_mat)                         # genes x pseudo-bulk samples

# metadata
pb_meta <- data.frame(
  sample_cluster = rownames(pb_mat)
)
pb_meta$sample  <- sub("_.*", "", pb_meta$sample_cluster)
pb_meta$cluster <- sub(".*_", "", pb_meta$sample_cluster)

#DEG analysis
y <- DGEList(counts = pb_counts)
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + pb_meta$cluster + pb_meta$sample)
colnames(design) <- make.names(colnames(design))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

contrast_aRG <- makeContrasts((pb_meta.clusterbRG1+pb_meta.clusterbRG2+pb_meta.clusterbRG3)/3 - pb_meta.clusteraRG, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast_aRG)
deg_table_aRG <- topTags(qlf, n = Inf)$table

contrast_EN <- makeContrasts((pb_meta.clusterbRG1+pb_meta.clusterbRG2+pb_meta.clusterbRG3)/3 - (pb_meta.clusterEN.DL+pb_meta.clusterEN)/2, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast_EN)
deg_table_Neu <- topTags(qlf, n = Inf)$table

######################pseudocells for data integration
cluster_col <- "seurat_clusters_names"
#extract counts and embeddings from seurat object
count_layers <- grep("^counts", Layers(seurat.integrated[["RNA"]]), value = TRUE)
counts_list <- lapply(count_layers, function(l) {
  LayerData(seurat.integrated[["RNA"]], layer = l)
})
# combine into one sparse matrix (genes × cells)
counts_mat <- do.call(cbind, counts_list)
#ensure column names match cells
colnames(counts_mat) <- colnames(seurat.integrated)
embeddings <- Embeddings(seurat.integrated, reduction = "pca")[, 1:26]         # cells x PCs
meta_df <- seurat.integrated@meta.data

##select target pseudocell size
target_size <- 25
n_cells <- ncol(counts_mat)
n_pseudo <- max(1, floor(n_cells / target_size))  # ~1 pseudocell per 25 cells

set.seed(123)
# pick seed cells
seed_idx <- sample(seq_len(n_cells), n_pseudo)
# assign every cell to its nearest seed in PCA space
nn_to_seed <- get.knnx(data = embeddings[seed_idx, , drop = FALSE],
                       query = embeddings,
                       k = 1)
assignments <- as.vector(nn_to_seed$nn.index)  # length = n_cells, values in 1..n_pseudo
# build sparse membership matrix S (cells x pseudocells)
S <- sparseMatrix(i = seq_len(n_cells), j = assignments, x = 1,
                  dims = c(n_cells, n_pseudo))
#aggregate counts per pseudocell to genes x pseudocells
pseudocell_counts <- counts_mat %*% S
group_sizes <- colSums(S)
pseudocell_counts <- pseudocell_counts %*% Diagonal(x = 1 / as.numeric(group_sizes))
#majority cluster per pseudocell to preserve cluster identity
clusters_char <- as.character(meta_df[[cluster_col]])
cluster_codes <- as.integer(factor(clusters_char))
n_clusters <- max(cluster_codes)

#cells x clusters
Cmat <- sparseMatrix(i = seq_len(n_cells), j = cluster_codes, x = 1,
                     dims = c(n_cells, n_clusters))

# n_pseudo x n_clusters
counts_pc_cluster <- t(S) %*% Cmat
majority_cluster_idx <- max.col(as.matrix(counts_pc_cluster), ties.method = "first")
cluster_levels <- levels(factor(clusters_char))
majority_cluster <- cluster_levels[majority_cluster_idx]

#preserve individual sample identity (3 cell lines d50 cerebral organoids)
sample_char <- as.character(meta_df[["Sample"]])
sample_codes <- as.integer(factor(sample_char))
n_sample <- max(sample_codes)

#cells x samples
Smatrix <- sparseMatrix(i = seq_len(n_cells), j = sample_codes, x = 1,
                        dims = c(n_cells, n_sample))

# n_pseudo x n_samples
counts_pc_samples <- t(S) %*% Smatrix
majority_samples_idx <- max.col(as.matrix(counts_pc_samples), ties.method = "first")
samples_levels <- levels(factor(sample_char))
majority_samples <- samples_levels[majority_samples_idx]

#seurat pseudocell object
pseudocell_obj <- CreateSeuratObject(
  counts = pseudocell_counts,
  meta.data = data.frame(
    cluster = majority_cluster,
    sample = majority_samples,
    size = as.numeric(group_sizes),
    row.names = paste0("pseudo_", seq_len(n_pseudo))
  )
)
colnames(pseudocell_obj) <- rownames(pseudocell_obj@meta.data)

cds_pseudo <- as.cell_data_set(pseudocell_obj)
cell_metadata <- colData(cds_pseudo)
raw_counts <- t(counts(cds_pseudo))
pseudocounts <- as.data.frame(as.matrix(raw_counts)) ##pseudocount matrix

cell_metadata$time <- 50
write.csv(pseudocounts, file = "marmoset_pseudocounts_day50.csv", quote = FALSE)
