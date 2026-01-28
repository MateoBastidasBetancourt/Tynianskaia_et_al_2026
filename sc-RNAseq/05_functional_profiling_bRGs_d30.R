library(Seurat)
library(MAST)
library(ACTIONetExperiment)
library(scHumanNet)
library(SCINET)
library(SingleCellExperiment)
library(igraph)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
####pseudobulk analysis
setwd("/home/basti/Downloads/pseudo_objects_231224/")
seurat_obj_combined_npcs <- readRDS("IPs.rds")

Idents(seurat_obj_combined_npcs)

pseudo_ifnb <- AggregateExpression(seurat_obj_combined_npcs, return.seurat = T, group.by = c("species", "cell_line", "harmony_clusters"))
print(head(pseudo_ifnb))
print(unique(Idents(pseudo_ifnb)))
# each 'cell' is a donor-condition-celltype pseudobulk profile
Cells(pseudo_ifnb)
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$harmony_clusters, pseudo_ifnb$species, sep = "_")
print(unique(pseudo_ifnb$celltype.stim))
print(unique(pseudo_ifnb$orig.ident))
head(pseudo_ifnb$orig.ident)
head(pseudo_ifnb$cell_line)

Idents(pseudo_ifnb) <- "celltype.stim"

###################################################################

bRG1_pseudo <- FindMarkers(pseudo_ifnb, ident.1 = c("3_Marmoset"), ident.2 = c("3_Human"), logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1)
bRG1_pseudo$gene <- rownames(bRG1_pseudo)
bRG1_pseudo$p_val.adj <- p.adjust(bRG1_pseudo$p_val, method = "BH")

bRG2_pseudo <- FindMarkers(pseudo_ifnb, ident.1 = c("7_Marmoset"), ident.2 = c("7_Human"), logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1)
bRG2_pseudo$gene <- rownames(bRG2_pseudo)
bRG2_pseudo$p_val.adj <- p.adjust(bRG2_pseudo$p_val, method = "BH")


##### differential network analysis (scHumanNet pipeline, Cha et al., 2023)
##### take d30 seurat object / human and marmoset cerebral organoids
seurat_obj_combined <- subset(seurat_obj_combined, subset = harmony_clusters_names %in% c("bRG1", "bRG2"))
labels <- paste(
  seurat_obj_combined$species,
  seurat_obj_combined$harmony_clusters_names,
  sep = "_"
)

meta_celltype_spp <- data.frame(
  species  = as.character(seurat_obj_combined$species),
  celltype = as.character(seurat_obj_combined$harmony_clusters_names),
  labels   = as.character(labels),
  stringsAsFactors = FALSE
)

rownames(meta_celltype_spp) <- colnames(seurat_obj_combined)


# to SingleCellExperiment
sce <- Seurat::as.SingleCellExperiment(seurat_obj_combined)

# ACE reduction
ace <- reduce.ace(sce)
ace = run.ACTIONet(ace = ace)
ace[["Labels"]] <- meta_celltype_spp$labels

ace = compute.cluster.feature.specificity(ace, ace$Labels, "celltype_specificity_scores")
Celltype.specific.networks = run.SCINET.clusters(ace, specificity.slot.name = "celltype_specificity_scores_feature_specificity")

data('HNv3_XC_LLS')
sorted.net.list <- SortAddLLS(Celltype.specific.networks, reference.network = graph.hn3)

strength.list <- GetCentrality(method='degree', net.list = sorted.net.list)


rank.df.final <- CombinePercRank(strength.list)

top.df <- TopHub(rank.df.final, top.n = 50)
head(top.df)
print("top_Hub")

diffPR.df <- DiffPR(rank.df.final, celltypes = 'celltype', condition = 'species', control = 'Human', meta = meta_celltype_spp)
head(diffPR.df)

###OUTPUT
#1 non-parametric
diffPR.df.sig <- FindDiffHub(rank.df.final = rank.df.final, celltypes = 'celltype', condition = 'species', control = 'Human', meta = meta_celltype_spp, net.list=sorted.net.list, q.method='BH', centrality="degree", min.cells= 15)
diffPR.df.sig <- diffPR.df.sig[diffPR.df.sig$pvalue<0.05,]
#2 parametric
diffPR.df.top <- TopDiffHub(diffPR.df, top.percent = 0.05)

######## estimate intersection between both approaches and do GO/KEGG analysis

marker_vector <- list(bRG1_pseudo, bRG2_pseudo)
name_vector <- c("bRG1_pseudo", "bRG2_pseudo")
hub_clusters <- c("bRG1", "bRG2")

library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)

#### with diffPR sig
for (j in 1:length(name_vector)){
  print(name_vector[j])
  print(hub_clusters[j])
  genes_hub  <- subset(diffPR.df.sig,
                       diffPR.df.sig$celltype == hub_clusters[j])
  genes_hub <- subset(genes_hub, genes_hub$diffPR > 0)
  genes_pseudo  <- subset(marker_vector[[j]],
                          marker_vector[[j]]$avg_log2FC > 0.5 & p_val.adj < 0.05)
  #& marker_vector[[j]]$p_val_adj <  0.05)
  #genes_pseudo <- genes_pseudo$gene
  
  common_genes <- intersect(genes_pseudo$gene, genes_hub$gene)
  
  result_up <- clusterProfiler::enrichGO(common_genes,
                                         "org.Hs.eg.db",
                                         keyType = "SYMBOL",
                                         ont = "ALL",
                                         pvalueCutoff = 0.1,
                                         minGSSize = 50)
  result_up <- clusterProfiler::simplify(result_up)
  
  result_up <- result_up@result
  
  kegg_res <- enrichKEGG(gene = bitr(common_genes, fromType="SYMBOL",
                                     toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.05,
                         minGSSize = 50)
  kegg_res <- kegg_res@result
  result_up <- bind_rows(result_up, kegg_res)
  
  var_name <- paste0(paste0(hub_clusters[j], "_integrated_up_"), "modules")
  
  assign(var_name, result_up, envir = .GlobalEnv)
  
  #genes down
  genes_hub  <- subset(diffPR.df.sig,
                       diffPR.df.sig$celltype == hub_clusters[j])
  genes_hub <- subset(genes_hub, genes_hub$diffPR < 0) ###CONTINUE FROM HERE ON.
  #genes_hub <- genes_hub$gene
  
  genes_pseudo  <- subset(marker_vector[[j]],
                          marker_vector[[j]]$avg_log2FC < -0.5 & p_val.adj < 0.05)
  #& marker_vector[[j]]$p_val_adj <  0.05)
  #genes_pseudo <- genes_pseudo$gene
  
  common_genes <- intersect(genes_pseudo$gene, genes_hub$gene)
  
  result_down <- clusterProfiler::enrichGO(common_genes,
                                           "org.Hs.eg.db",
                                           keyType = "SYMBOL",
                                           ont = "ALL",
                                           pvalueCutoff = 0.1,
                                           minGSSize = 50)
  result_down <- clusterProfiler::simplify(result_down)
  
  result_down <- result_down@result
  
  kegg_res <- enrichKEGG(gene = bitr(common_genes, fromType="SYMBOL",
                                     toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.1,
                         minGSSize = 50)
  kegg_res <- kegg_res@result
  result_down <- bind_rows(result_down, kegg_res)
  
  var_name <- paste0(paste0(hub_clusters[j], "_integrated_down_"), "modules")
  assign(var_name, result_down, envir = .GlobalEnv)
  print(common_genes)
  
  ####DEGs
  genes_hub <- subset(diffPR.df.sig,
                      diffPR.df.sig$celltype == hub_clusters[j])
  
  
  #genes <- genes$gene
  genes_pseudo <- subset(
    marker_vector[[j]],
    abs(avg_log2FC) > 0.5 & p_val.adj < 0.05
  )
  
  
  common_genes <- intersect(genes_pseudo$gene, genes_hub$gene)
  
  result_de <- clusterProfiler::enrichGO(common_genes,
                                         "org.Hs.eg.db",
                                         keyType = "SYMBOL",
                                         ont = "ALL",
                                         pvalueCutoff = 0.1,
                                         minGSSize = 50)
  
  result_de <- clusterProfiler::simplify(result_de)
  
  result_de <- result_de@result
  
  kegg_res <- enrichKEGG(gene = bitr(common_genes, fromType="SYMBOL",
                                     toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.1,
                         minGSSize = 50)
  kegg_res <- kegg_res@result
  result_de <- bind_rows(result_de, kegg_res)
  
  var_name <- paste0(paste0(hub_clusters[j], "_integrated_de_"), "modules")
  
  assign(var_name, result_de, envir = .GlobalEnv)
}
