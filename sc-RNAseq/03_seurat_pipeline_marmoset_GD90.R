library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(FNN)
library(monocle3)
library(SeuratWrappers)
library(edgeR)
setwd("/mnt/vast-standard/home/cbastid/u15756")
seurat.merged.filtered <- readRDS("seurat.merged.filtered.rds")

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
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:25)

seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.45)
seurat.integrated$res_0.45 <- Idents(seurat.integrated)


seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.45)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:25)

cluster_markers_brains <- FindAllMarkers(seurat.integrated, only.pos = T, logfc.threshold = 0.25, group.by = "seurat_clusters_names") 


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


###### assign celltype names

seurat.integrated$seurat_clusters <- as.character(seurat.integrated$seurat_clusters)


seurat.integrated$seurat_clusters_names <- ifelse(seurat.integrated$seurat_clusters %in% c("0"), "EN-UL",
                                                  ifelse(seurat.integrated$seurat_clusters %in% c("3"), "tRG/bRG",
                                                         ifelse(seurat.integrated$seurat_clusters %in% c("1"), "IPC2",
                                                                ifelse(seurat.integrated$seurat_clusters %in% c("2"), "IN1",
                                                                       ifelse(seurat.integrated$seurat_clusters %in% c("8"), "IN2",
                                                                              ifelse(seurat.integrated$seurat_clusters %in% c("4"), "aRG",
                                                                                     ifelse(seurat.integrated$seurat_clusters %in% c("5"), "NB",
                                                                                            ifelse(seurat.integrated$seurat_clusters %in% c("6"), "MSN",
                                                                                                   ifelse(seurat.integrated$seurat_clusters %in% c("7"), "bRG",
                                                                                                          ifelse(seurat.integrated$seurat_clusters %in% c("10"), "IPC1",
                                                                                                                 ifelse(seurat.integrated$seurat_clusters %in% c("11"), "IPC3",
                                                                                                                        ifelse(seurat.integrated$seurat_clusters %in% c("12"), "OPC/Astrocyte",
                                                                                                                               ifelse(seurat.integrated$seurat_clusters %in% c("13"), "tRG",
                                                                                                                                      ifelse(seurat.integrated$seurat_clusters %in% c("9"), "RBC",
                                                                                                                                             ifelse(seurat.integrated$seurat_clusters %in% c("14"), "ChP/EN",
                                                                                                                                                    ifelse(seurat.integrated$seurat_clusters %in% c("15"), "EN",
                                                                                                                                                           ifelse(seurat.integrated$seurat_clusters %in% c("17"), "EN-DL",
                                                                                                                                                                  ifelse(seurat.integrated$seurat_clusters %in% c("16","18"), "Mic",
                                                                                                                                                                         ifelse(seurat.integrated$seurat_clusters %in% c("19"), "End", NA)))))))))))))))))))


new_order <- c("aRG", "tRG", "tRG/bRG", "ChP/EN", "bRG", "IPC1", "IPC2", "IPC3", "NB", "IN1","IN2", "MSN", "EN", "EN-UL", "EN-DL", "OPC/Astrocyte", "Mic", "RBC", "End")

seurat.integrated$seurat_clusters_names <- factor(
  seurat.integrated$seurat_clusters_names,
  levels = new_order
)

Idents(seurat.integrated) <- "seurat_clusters_names"

####pseudobulk analysis of PALMD expression
# Metadata and grouping
meta <- seurat.integrated@meta.data
meta$cluster <- Idents(seurat.integrated)
meta$group <- paste(meta$Sample, meta$cluster, sep = "_")  # sample-cluster combo

# Extract raw counts (Seurat v4 style)
counts <- seurat.integrated[["RNA"]]@counts

# Transpose (cells x genes) and aggregate counts per group
counts_t <- t(as.matrix(counts))               # cells x genes
pb_mat <- rowsum(counts_t, group = meta$group) # sum counts per sample-cluster group
pb_counts <- t(pb_mat)                         # genes x pseudobulk samples

# Metadata for pseudobulk samples
pb_meta <- data.frame(
  sample_cluster = rownames(pb_mat)
)
pb_meta$sample  <- sub("_.*", "", pb_meta$sample_cluster)
pb_meta$cluster <- sub(".*_", "", pb_meta$sample_cluster)

# Check result
dim(pb_counts)
head(pb_meta)

y <- DGEList(counts = pb_counts)
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + pb_meta$cluster + pb_meta$sample)
colnames(design) <- make.names(colnames(design))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

contrast_aRG <- makeContrasts((pb_meta.clusterbRG+pb_meta.clustertRG.bRG)/2 - pb_meta.clusteraRG, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast_aRG)
topTags(qlf, n = 10)
deg_table_aRG <- topTags(qlf, n = Inf)$table
View(deg_table_aRG)


contrast_EN <- makeContrasts((pb_meta.clusterbRG+pb_meta.clustertRG.bRG)/2 - (pb_meta.clusterEN.DL+pb_meta.clusterEN+pb_meta.clusterChP.EN)/3, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast_EN)
topTags(qlf, n = 10)
deg_table_Neu <- topTags(qlf, n = Inf)$table
View(deg_table_Neu)

######################pseudocells 
cluster_col <- "seurat_clusters_names"
#now extract counts and embeddings 
seurat.merged.filtered$seurat_clusters <- seurat.integrated$seurat_clusters_names
###remove non-ectodermal lineages (not present in cerebral organoids)
Idents(seurat.merged.filtered) <- seurat.merged.filtered$seurat_clusters
seurat.merged.filtered <- subset(
  seurat.merged.filtered,
  idents = setdiff(
    levels(seurat.merged.filtered$seurat_clusters),
    c("RBC", "Mic", "End")
  )
)

seurat.integrated <- subset(
  seurat.integrated,
  idents = setdiff(
    levels(seurat.integrated$seurat_clusters_names),
    c("RBC", "Mic", "End")
  )
)

count_layers <- grep("^counts", Layers(seurat.merged.filtered[["RNA"]]), value = TRUE)

counts_list <- lapply(count_layers, function(l) {
  LayerData(seurat.merged.filtered[["RNA"]], layer = l)
})

#combine into one sparse matrix (genes × cells)
counts_mat <- do.call(cbind, counts_list)
#ensure column names match cells
colnames(counts_mat) <- colnames(seurat.integrated)
embeddings <- Embeddings(seurat.integrated, reduction = "pca")[, 1:25]         # cells x PCs
meta_df <- seurat.integrated@meta.data

#select target pseudocell size
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

#n_pseudo x n_clusters
counts_pc_cluster <- t(S) %*% Cmat
majority_cluster_idx <- max.col(as.matrix(counts_pc_cluster), ties.method = "first")
cluster_levels <- levels(factor(clusters_char))
majority_cluster <- cluster_levels[majority_cluster_idx]

#preserve individual sample identity (3 cj GD90 brain samples)
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
pseudocounts <- as.data.frame(as.matrix(raw_counts))
cell_metadata$time <- 90
write.csv(pseudocounts, file = "marmoset_pseudocounts_brains.csv", quote = FALSE)

