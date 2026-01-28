library(edgeR)
library(dplyr)
library(pheatmap)
####starting from read count table (EGA)
dataset_merged <- read.table("raw_counts_eb.csv", sep = "\t", header = TRUE) #open count matrix
pheno_group <- read.table("pheno_data_EBs.csv", sep = "\t", header = TRUE) #open sample group file in metadata folder
pheno_group$species_stage <- factor(pheno_group$species_stage)

design <- model.matrix(~ 0 + species_stage, data = pheno_group)
dge_obj <- DGEList(counts=dataset_merged, group=pheno_group[,2], genes=rownames(dataset_merged))
identical(pheno_group$sample.id, colnames(dataset_merged)) #See if order in count table is identical than in the metadata
dge_TMM <- calcNormFactors(dge_obj, method="TMM")
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, keep.lib.sizes=FALSE] 
#####
plotMDS(dge_TMM, method="bcv", col = c("blue", "red")[factor(pheno_group$species)], pch = 16, cex = 1.5,xlab = "MDS1", ylab = "MDS2")

#rownames(design) <- pheno_group$sample.id

disp <- estimateDisp(dge_TMM, design, robust=TRUE)
plotBCV(disp)

CONTRASTS <- makeContrasts(HumanvsMarmoset = (species_stagehuman_d3+species_stagehuman_d6)/2-(species_stagemarmoset_d3+species_stagemarmoset_d6)/2, 
                           d6vsd3 = (species_stagehuman_d6+species_stagemarmoset_d6)/2-(species_stagehuman_d3+species_stagemarmoset_d3)/2, 
                           HumanvsMarmoset_NIM = (species_stagehuman_d6-species_stagehuman_d3)-(species_stagemarmoset_d6-species_stagemarmoset_d3) ,
                           d6vsd3_human = (species_stagehuman_d6-species_stagehuman_d3) , 
                           d6vsd3_marmoset = (species_stagemarmoset_d6-species_stagemarmoset_d3),
                           d6humanvsd6marmoset = (species_stagehuman_d6-species_stagemarmoset_d6),
                           d3humanvsd3marmoset  = (species_stagehuman_d3-species_stagemarmoset_d3), levels=design)

for (i in 1:ncol(CONTRASTS)){
  
  contrast_name <- colnames(CONTRASTS)[i]
  current.glmQLFTest <- glmQLFTest(fit, contrast = CONTRASTS[, i]) # perform the test
  assign(contrast_name, current.glmQLFTest) #ssign the result to a variable with the contrast name
}


#####PVALUE CORRECTION AND SORTING GENES ACCORDING TO THEIR PVALUE
result_HumanvsMarmoset <- topTags(HumanvsMarmoset,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d6vsd3 <- topTags(d6vsd3,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_HumanvsMarmoset_NIM <- topTags(HumanvsMarmoset_NIM,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d6vsd3_human <- topTags(d6vsd3_human,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d6vsd3_marmoset <- topTags(d6vsd3_marmoset,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d6humanvsd6marmoset = topTags(d6humanvsd6marmoset,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d3humanvsd3marmoset  = topTags(d3humanvsd3marmoset,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
######################
filter_results <- function(result) {
  filtered_result <- result$table[
    (result$table$logFC > 0.5 | result$table$logFC < -0.5) & result$table$FDR < 0.01, ]
  return(filtered_result)
}
########################
result_HumanvsMarmoset <- filter_results(result_HumanvsMarmoset)
result_HumanvsMarmoset_NIM <- filter_results(result_HumanvsMarmoset_NIM)
result_d6vsd3 <- filter_results(result_d6vsd3)
result_d6vsd3_human <- filter_results(result_d6vsd3_human)
result_d6vsd3_marmoset <- filter_results(result_d6vsd3_marmoset)
result_d6humanvsd6marmoset <- filter_results(result_d6humanvsd6marmoset)
result_d3humanvsd3marmoset <- filter_results(result_d3humanvsd3marmoset)

######### heatmap
library(pheatmap)
library(dplyr)

logCPM <- cpm(dge_TMM, log = TRUE, prior.count = 1)
# Define gene groups

group_m <- c("TUBB3", "NOG", "ZIC1", "SOX3", "DCX", "CENPE", "FABP7", "ASCL1", "NEUROG1", "NEUROG2")
group_non <- c("NES", "GLI3", "VIM", "SFRP1", "ZIC2", "OTX2", "LHX2", "CER1",
               "CHRD", "DKK1", "SIX3", "PAX6", "HES1", "HES5", "INSM1", "NOTCH1", "SOX1", "SOX2")
group_h <- c("FGF2", "EGF", "NLN", "FOXD3", "TP53I3", "CHK1", "LDHA", "MAP3K8", "RRAS2", "FEZF1", "RPL10", "EIF4E1B")

# Combine all into a named group vector
gene_groups <- c(
  setNames(rep("Marmoset up", length(group_m)), group_m),
  setNames(rep("Non-DE", length(group_non)), group_non),
  setNames(rep("Human up", length(group_h)), group_h)
)
##########################################################################3
selected_genes <- c(group_m, group_non, group_h)
genes_in_matrix <- intersect(selected_genes, rownames(logCPM)) # keep those genes present in expression matrix
expr_subset <- logCPM[genes_in_matrix, ] # subset expression matrix

annotation_row <- data.frame(Group = gene_groups[rownames(expr_subset)])
rownames(annotation_row) <- rownames(expr_subset)

ann_colors <- list(
  Group = c(
    "Marmoset up" = "#377EB8",    # Blue
    "Non-DE" = "#4DAF4A",          # Green
    "Human up" = "#E41A1C"
  )
)
library(pheatmap)
pheatmap(expr_subset,
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 8,
         fontsize_row = 8,
         annotation_row = annotation_row,
         annotation_colors = ann_colors)

