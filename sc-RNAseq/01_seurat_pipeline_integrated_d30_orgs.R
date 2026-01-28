library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FNN)
library(Matrix)

setwd("~/RNAseq/STAR/cellranger/marmoset_2/outs/count")
marmoset_h5 = Read10X_h5("filtered_feature_bc_matrix.h5")
df_metadata <- data.frame(row.names = 1, "cell_names"=colnames(seurat_obj), "cell_line"=c(rep(c("cj4"),15268),rep(c("cj5"),(32481-15268)),rep(c("cj6"),(47734-32481))))
seurat_obj <- CreateSeuratObject(marmoset_h5,min.cells = 3, min.features = 200 , assay = "RNA",names.delim = "-", meta.data = df_metadata, names.field = 2)
remove(marmoset_h5)
remove(df_metadata)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA >500 & percent.mt < 10)
seurat_obj <- SCTransform(seurat_obj)
SaveSeuratRds(seurat_obj,"seur_obj.Rds")

###### human datasets
setwd("~/RNAseq/STAR/cellranger/human_1/outs/count")
#ids <- read.csv(file.path("../", "aggregation.csv"))
human_h5 = Read10X_h5("filtered_feature_bc_matrix.h5")
seurat_obj_human <- CreateSeuratObject(human_h5,min.cells = 3, min.features = 200, assay = "RNA",names.delim = "-", names.field = 2)

df_metadata_human <- data.frame(row.names = 1, "cell_names"=colnames(seurat_obj_human), "cell_line"=c(rep(c("IDP"),23356),rep(c("sc"),(45449-23356)), rep(c("iLonza"), (54981-45449))))

seurat_obj_human <- CreateSeuratObject(human_h5,min.cells = 3, min.features = 200, assay = "RNA",names.delim = "-", meta.data = df_metadata_human, names.field = 2)
remove(human_h5)


seurat_obj_human[["percent.mt"]] <- PercentageFeatureSet(seurat_obj_human, pattern = "^MT")
seurat_obj_human <- subset(seurat_obj_human, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA >500 & percent.mt < 10)
seurat_obj_human <- SCTransform(seurat_obj_human)
SaveSeuratRds(seurat_obj_human,"seurat_obj_human.Rds")
SaveSeuratRds(seurat_obj,"seurat_obj.Rds")
seurat_obj_human <- readRDS("seurat_obj_human.Rds")


##### integrate for species

seurat_obj_combined <- merge(seurat_obj,seurat_obj_human)
df_metadata_combined <- data.frame(row.names = 1, "cell_names"=colnames(seurat_obj_combined), "species"=c(rep(c("Marmoset"),46667),rep(c("Human"),(97742-46667))))
remove(seurat_obj)
remove(seurat_obj_human)
seurat_obj_combined <- AddMetaData(seurat_obj_combined,df_metadata_combined, col.name = "species")
#seurat_obj_combined <- JoinLayers(seurat_obj_combined[["SCT"]])
seurat_obj_combined[["SCT"]] <- split(seurat_obj_combined[["SCT"]], f = seurat_obj_combined$species)
seurat_obj_combined <- SCTransform(seurat_obj_combined, vst.flavor = "v2")
#seurat_obj_combined <- ScaleData(seurat_obj_combined, verbose = FALSE)
seurat_obj_combined <- RunPCA(seurat_obj_combined, npcs = 30, verbose = FALSE)
ElbowPlot(seurat_obj_combined, ndims = 30)
# one-liner to run Integration
seurat_obj_combined <- IntegrateLayers(object = seurat_obj_combined, method = HarmonyIntegration,
                                       orig.reduction = "pca", new.reduction = 'harmony',
                                       assay = "SCT", verbose = FALSE)
seurat_obj_combined <- FindNeighbors(seurat_obj_combined, reduction = "harmony", dims = 1:19) ###ACHTUNG, 19 DIMENSIONS WERE USED!
seurat_obj_combined <- FindClusters(seurat_obj_combined, resolution = 0.27, cluster.name = "harmony_clusters")
seurat_obj_combined$harmony_clusters_res_0.27 <- Idents(seurat_obj_combined)



seurat_obj_combined <- RunUMAP(seurat_obj_combined, dims = 1:19)
##value based on inte4rmediate between elbow and summation'

##object was saved until here, this next step wasn't really saved. Careful!
seurat_obj_combined <- PrepSCTFindMarkers(seurat_obj_combined)
cluster.markers <- FindAllMarkers(seurat_obj_combined, only.pos = T, logfc.threshold = 0.5, group.by = "harmony_clusters") 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_obj_combined <- CellCycleScoring(seurat_obj_combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

################## Annotation:

CT_markers_geschwind <- read.csv("/home/snst1/ZGPZ/cbastid/RNAseq/STAR/gene_markers_geschwind.csv")[2:3]
CT_markers_artegiani <- read.csv("/home/snst1/ZGPZ/cbastid/RNAseq/STAR/gene_markers_artegiani.csv")[1:2]
CT_markers_kanton <- read.csv("/home/snst1/ZGPZ/cbastid/RNAseq/STAR/marker_genes_kanton.csv")[7:8]

CT_markers_linnarsson <- read.csv("/home/snst1/ZGPZ/cbastid/RNAseq/STAR/marker_genes_linnarsson.csv")

CT_markers_linnarsson <- mutate(CT_markers_linnarsson, class_subclass = paste(Class, Subclass, sep = "-"))

# Step 2: Split marker_genes column into individual genes
CT_markers_linnarsson <- separate_rows(CT_markers_linnarsson, MarkerGenes, sep = " ")

# Step 3: Remove non-unique class-subclass + marker assignments
CT_markers_linnarsson <- CT_markers_linnarsson %>%
  group_by(class_subclass, MarkerGenes) %>%
  filter(n() == 1) %>%
  ungroup()
CT_markers_linnarsson <- CT_markers_linnarsson[,c(10,42)]
CT_markers_linnarsson [1] <- toupper(CT_markers_linnarsson$MarkerGenes)

colnames(CT_markers_linnarsson) <- c("gene","cluster")

colnames(CT_markers_geschwind) <- c("gene","cluster")
colnames(CT_markers_artegiani) <- c("gene","cluster")


###This is a custom function that retrieves the most likely cell all clusters belong to, based on literature tables

CTretrieve <- function(x,y,z){
  cluster <- c()
  for(i in 1:length(levels(x$cluster))){
    #top10 <- x %>% filter(cluster %in% levels(x$cluster)[i]) %>% pull(-1)
    top10 <- x %>%
      filter(cluster == levels(x$cluster)[i]) %>%
      arrange(desc(avg_log2FC)) %>%  # Sort by avg_log2FC in descending order
      slice_head(n = as.numeric(z)) %>%         # Select top 10 rows
      pull(-1)  
    top_cluster <- y %>% 
      filter(gene %in% top10) %>% 
      count(cluster) %>% 
      arrange(desc(n)) %>% 
      slice(1) %>% 
      pull(1)
    
    cluster <- c(cluster, if (length(top_cluster) == 0) NA else top_cluster)
    if (length(cluster[i]) == 0){
      cluster = c(cluster,"NA")
    }
  }
  return(cluster)
  
}


#cluster.markers_alt1$cluster <- factor(cluster.markers_alt1$cluster, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
#cluster.markers = cluster.markers_alt1
markers_combi_10 <- data.frame("cluster" = seq(0,16),"kanton"=CTretrieve(cluster.markers, CT_markers_kanton, 10),"artegiani"=CTretrieve(cluster.markers, CT_markers_artegiani, 10), "linnarsson"= CTretrieve(cluster.markers, CT_markers_linnarsson, 10), "geschwind"= CTretrieve(cluster.markers, CT_markers_geschwind, 10))
markers_combi_30 <- data.frame("cluster" = seq(0,16),"kanton"=CTretrieve(cluster.markers, CT_markers_kanton, 30),"artegiani"=CTretrieve(cluster.markers, CT_markers_artegiani, 30), "linnarsson"= CTretrieve(cluster.markers, CT_markers_linnarsson, 30), "geschwind"= CTretrieve(cluster.markers, CT_markers_geschwind, 30))
markers_combi_50 <- data.frame("cluster" = seq(0,16),"kanton"=CTretrieve(cluster.markers, CT_markers_kanton, 50),"artegiani"=CTretrieve(cluster.markers, CT_markers_artegiani, 50), "linnarsson"= CTretrieve(cluster.markers, CT_markers_linnarsson, 50), "geschwind"= CTretrieve(cluster.markers, CT_markers_geschwind, 50))
########Celltype assignment

####NAME ASSIGNMENT OF SEURAT CLUSTERS TO ACTUAL CELL NAMES
seurat_obj_combined <- readRDS("seurat_obj_combined.Rds")

seurat_obj_combined$harmony_clusters <- as.character(seurat_obj_combined$harmony_clusters)

seurat_obj_combined$harmony_clusters_names <- ifelse(seurat_obj_combined$harmony_clusters %in% c("0"), "EN",
                                                     ifelse(seurat_obj_combined$harmony_clusters %in% c("13"), "IN",
                                                            ifelse(seurat_obj_combined$harmony_clusters %in% c("6", "9"), "aRG3",
                                                                   ifelse(seurat_obj_combined$harmony_clusters %in% c("3"), "bRG1", 
                                                                          ifelse(seurat_obj_combined$harmony_clusters %in% c("7"), "bRG2",
                                                                                 ifelse(seurat_obj_combined$harmony_clusters %in% c("4", "14"), "tRG", 
                                                                                        ifelse(seurat_obj_combined$harmony_clusters %in% c("2"), "aRG4",
                                                                                               ifelse(seurat_obj_combined$harmony_clusters %in% c("11"), "NCC",
                                                                                                      ifelse(seurat_obj_combined$harmony_clusters %in% c("5"), "OPC",
                                                                                                             ifelse(seurat_obj_combined$harmony_clusters %in% c("1", "10"), "aRG1",
                                                                                                                    ifelse(seurat_obj_combined$harmony_clusters %in% c("8", "12"), "aRG2",
                                                                                                                           ifelse(seurat_obj_combined$harmony_clusters %in% c("16"), "IPC",
                                                                                                                                  ifelse(seurat_obj_combined$harmony_clusters %in% c("15"), "Others", NA
                                                                                                                                  ))))))))))))))


new_order <- c("aRG1", "aRG2", "aRG3", "aRG4", "tRG", "bRG1", "bRG2", 
               "IPC", "EN", "IN", 
               "OPC", "NCC")

seurat_obj_combined$harmony_clusters_names <- factor(
  seurat_obj_combined$harmony_clusters_names,
  levels = new_order
)

Idents(seurat_obj_combined) <- "harmony_clusters_names"

#############creating df for cell_composition_analysis.R
cell_info <- seurat_obj_combined@meta.data[, c("harmony_clusters_names", 
                                               "cell_line", 
                                               "species", 
                                               "Phase")]

#3D UMAP cluster/species plot for visualization:
set.seed(123) 

seurat_obj_combined <- RunUMAP(
  seurat_obj_combined,
  reduction = "pca",      
  dims = 1:19,                
  n.components = 3,           
  reduction.name = "umap_3d"    
)

#extract 3D UMAP coordinates
umap3d <- Embeddings(seurat_obj_combined, reduction = "umap_3d")[, 1:3]


#get cluster annotation
clusters <- seurat_obj_combined$harmony_clusters_names
species <- seurat_obj_combined$species

#combine into a data frame
df <- data.frame(
  x = umap3d[, 1],
  y = umap3d[, 2],
  z = umap3d[, 3],
  cluster = clusters,
  species = species
)

# Randomly sample 5000 cells (or all if fewer than 5000)
df_sample <- df[sample(nrow(df), min(5000, nrow(df))), ]
write.csv(df_sample, "cells3d_sample.csv", row.names = FALSE) ###plot per cluster and per species
####


######## Pseudocell extraction and species trajectory analysis
#to monocle3
cds_pseudo <- as.cell_data_set(seurat_obj_combined) ####as.cds

colData(cds_pseudo)$orig_cluster <- seurat_obj_combined$harmony_clusters_names
####################################################### Trajectory per species

# subset human and marmoset cells
cds_pseudo <- cds_pseudo[, !(colData(cds_pseudo)$orig_cluster %in% "Others")]
#get cell metadata
cell_meta <- colData(cds_pseudo)
#subset Human
human_cells <- rownames(cell_meta)[cell_meta$species == "Human"]
cds_human <- cds_pseudo[, human_cells]

#subset Marmoset
marmoset_cells <- rownames(cell_meta)[cell_meta$species == "Marmoset"]
cds_marmoset <- cds_pseudo[, marmoset_cells]

#########
cds_human <- preprocess_cds(cds_human, num_dim = 19)

cds_human <- reduce_dimension(cds_human, reduction_method = "UMAP")

cds_human <- cluster_cells(cds_human, resolution = 1e-6)

cds_human <- learn_graph(
  cds_human,
  use_partition = FALSE,
  close_loop = FALSE,
  learn_graph_control = list(
    minimal_branch_len = 10,       
    ncenter = 180,                
    geodesic_distance_ratio = 0.3 
  )
)


cds_human <- order_cells(cds_human, root_cells = colnames(cds_human)[cds_human$ident == "aRG1"])

# Repeat for Marmoset
cds_marmoset <- preprocess_cds(cds_marmoset, num_dim = 19)
cds_marmoset <- reduce_dimension(cds_marmoset, reduction_method = "UMAP")
cds_marmoset <- cluster_cells(cds_marmoset, resolution = 1e-6)
cds_marmoset <- learn_graph(
  cds_marmoset,
  use_partition = FALSE,
  close_loop = FALSE,
  learn_graph_control = list(
    minimal_branch_len = 10,       
    ncenter = 160,               
    geodesic_distance_ratio = 0.3 
  )
)

cds_marmoset <- order_cells(cds_marmoset, root_cells = colnames(cds_marmoset)[cds_marmoset$ident == "aRG1"])

###plot
#Human
plot_cells(
  cds_human,
  color_cells_by = "ident",  
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  cell_size = 0.4,
  group_label_size=5
)

# Marmoset
plot_cells(
  cds_marmoset,
  color_cells_by = "ident",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  cell_size = 0.4,
  group_label_size=4
)

##### Quantitative pseudotime analysis
cds_pseudo_df <- data.frame(
  species = colData(cds_pseudo)$species,
  pseudotime = pseudotime(cds_pseudo), 
  cluster = colData(cds_pseudo)$ident,
  cell_line = colData(cds_pseudo)$cell_line
)

summary_pseudo_df <- cds_pseudo_df %>%
  group_by(cell_line, cluster) %>%
  summarize(
    median_pt = median(pseudotime, na.rm = TRUE),
    mean_pt   = mean(pseudotime, na.rm = TRUE),
    n_cells   = n()
  ) %>%
  arrange(cluster, cell_line)

summary_pseudo_df$species <- ifelse(
  summary_pseudo_df$cell_line %in% c("IDP", "iLonza", "sc"),
  "Human",
  "Marmoset"
)

en_df <- summary_pseudo_df %>%
  filter(cluster %in% c("EN"))##or OPC

library(dplyr)

species_stats <- en_df %>%
  ungroup() %>%       
  group_by(species, cluster) %>%
  summarize(
    mean_pt = mean(median_pt, na.rm = TRUE),
    #sd_pt   = sd(mean_pt, na.rm = TRUE),
    n_lines = n()
  )
#species_stats[1,3] <- sd(c(0.62, 0.86, 1.16))
#species_stats[2,3] <- sd(c(0.97, 0.66, 1.01))

ggplot(species_stats, aes(x = species, y = mean_pt, fill = species)) +
  geom_col(width = 0.6, alpha = 0.8) +
  labs(
    x = "Species",
    y = "Median pseudotime"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

#####PALMD expression
df <- data.frame(
  Expression = FetchData(seurat_obj_combined, vars = "PALMD")[,1],
  Species   = seurat_obj_combined$species # or whatever metadata column you used
)

ggplot(df, aes(x = Species, y = Expression, fill = Species)) +
  geom_violin(trim = TRUE, adjust = 2, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.18) +   # more transparent + smaller
  theme_bw(base_size = 16) +
  labs(y = "Expression Level") +  theme(
    axis.title.y = element_text(size = 25),
    axis.text.x  = element_text(size = 25)
  )+
  theme(
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 20),
    legend.key.size = unit(1.2, "lines")   # optional: makes the color boxes bigger
  )

####
#Pseudocell inference

# Keep only cells that are NOT in "Others"
subset(seurat_obj_combined@meta.data, harmony_clusters_names != "Others")
seurat_obj_combined <- subset(seurat_obj_combined, cells = cells_keep)
remove(cells_keep)
gc()
embeddings <- Embeddings(seurat_obj_combined, reduction = "pca")[, 1:19]
meta_df <- seurat_obj_combined@meta.data

# Extract counts and embeddings from filtered object
counts_mat <- GetAssayData(seurat_obj_combined, assay = "RNA", slot = "counts")

#choose how many pseudocells 
target_size <- 25 
n_cells <- ncol(counts_mat)
n_pseudo <- max(1, floor(n_cells / target_size))  # ~1 pseudocell per 25 neighboring cells

set.seed(123)


# pick seed cells
seed_idx <- sample(seq_len(n_cells), n_pseudo)

# assign every cell to its nearest seed in PCA space
nn_to_seed <- get.knnx(data = embeddings[seed_idx, , drop = FALSE],
                       query = embeddings,
                       k = 1)
assignments <- as.vector(nn_to_seed$nn.index)  # length = n_cells, values in 1..n_pseudo

# build sparse membership matrix S (cells x pseudocells): 1 if cell belongs to that pseudocell
S <- sparseMatrix(i = seq_len(n_cells), j = assignments, x = 1,
                  dims = c(n_cells, n_pseudo))


#aggregate counts per pseudocell (genes x pseudocells)
pseudocell_counts <- counts_mat %*% S
group_sizes <- colSums(S)
pseudocell_counts <- pseudocell_counts %*% Diagonal(x = 1 / as.numeric(group_sizes))
#majority cluster per pseudocell to preserve the cluster identity
cluster_col <- "harmony_clusters_names"
clusters_char <- as.character(meta_df[[cluster_col]])
cluster_codes <- as.integer(factor(clusters_char))
n_clusters <- max(cluster_codes)

# cells x clusters
Cmat <- sparseMatrix(i = seq_len(n_cells), j = cluster_codes, x = 1,
                     dims = c(n_cells, n_clusters))

# counts of clusters within each pseudocell (n_pseudo x n_clusters)
counts_pc_cluster <- t(S) %*% Cmat
majority_cluster_idx <- max.col(as.matrix(counts_pc_cluster), ties.method = "first")
cluster_levels <- levels(factor(clusters_char))
majority_cluster <- cluster_levels[majority_cluster_idx]

#majority species per pseudocell to preserve species identity / extract marmoset pseudocells
species_char <- as.character(meta_df[["species"]])
species_codes <- as.integer(factor(species_char))
n_species <- max(species_codes)

#cells x species
Smatrix <- sparseMatrix(i = seq_len(n_cells), j = species_codes, x = 1,
                        dims = c(n_cells, n_species))

# counts of species within each pseudocell (n_pseudo x n_species)
counts_pc_species <- t(S) %*% Smatrix
majority_species_idx <- max.col(as.matrix(counts_pc_species), ties.method = "first")
species_levels <- levels(factor(species_char))
majority_species <- species_levels[majority_species_idx]

#create pseudocell object with cluster and species information
pseudocell_obj <- CreateSeuratObject(
  counts = pseudocell_counts,
  meta.data = data.frame(
    cluster = majority_cluster,
    species = majority_species,
    size = as.numeric(group_sizes),
    row.names = paste0("pseudo_", seq_len(n_pseudo))
  )
)
colnames(pseudocell_obj) <- rownames(pseudocell_obj@meta.data)
cds_pseudo <- as.cell_data_set(pseudocell_obj)
#get pseudocells coming from marmoset
cell_metadata <- colData(cds_pseudo)
#filter cells from marmoset
marmoset_pseudocells <- rownames(cell_metadata)[cell_metadata$species == "Marmoset"]
# Extract raw pseudo-count matrix
raw_counts <- t(counts(cds_pseudo))

# Subset the counts matrix to only marmoset cells
marmoset_pseudocounts <- raw_counts[marmoset_pseudocells,]
marmoset_pseudocounts <- as.data.frame(as.matrix(marmoset_pseudocounts))
dim(marmoset_pseudocounts)
cell_metadata$time <- 30
cell_metadata_marmoset <- cell_metadata[cell_metadata$species=="Marmoset",]
write.csv(marmoset_pseudocounts, file = "marmoset_pseudocounts_day30_marmoset.csv", quote = FALSE)
