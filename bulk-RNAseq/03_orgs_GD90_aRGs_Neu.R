library(dplyr)
library(edgeR)
library(biomaRt)
library(ggplot2)

dataset_aRG <- read.table("raw_counts_AP_orgs.csv", sep = "\t", header = TRUE)
dataset_brains <- read.table("raw_counts_AP_Neu_brains.csv", sep = "\t", header = TRUE)
dataset_combined <- cbind(dataset_brains, dataset_aRG[, 3:ncol(dataset_aRG)])
pheno_group <- read.table("pheno_data_marmoset_orgs&brains.csv", sep = "\t", header = TRUE) #open sample group file in metadata folder

######## Starts the statistical analysis
design <- model.matrix(~ 0 +  origin, data = pheno_group)
dge_obj <- DGEList(counts=dataset_combined[,3:17], group=factor(pheno_group[,2]), genes=dataset_combined$gene_id_final)
identical(pheno_group$stage_line, colnames(dataset_combined[,3:17])) #See if order in count table is identical than in the metadata

dge_TMM <- calcNormFactors(dge_obj, method="upperquartile")
keep <- filterByExpr(dge_TMM, design = design, min.count = 15, min.prop = 0.7, min.total.count = 15)
dge_TMM <- dge_TMM[keep, keep.lib.sizes=FALSE]
rownames(dge_TMM) <- make.unique(dge_TMM$genes$genes)

##################calculate tpm of marker genes
library(biomaRt)

mart <- useEnsembl(biomart = "genes")
mart <- useEnsembl(biomart = "genes", dataset = "cjacchus_gene_ensembl")
gene_list <- dataset_combined$ensembl_gene_id

gene_info <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id", "transcript_length"),
  filters = "ensembl_gene_id",
  values = gene_list,
  mart = mart
)

library(dplyr)

gene_info2 <- gene_info %>%
  group_by(ensembl_gene_id) %>%
  summarise(gene_length = max(transcript_length, na.rm = TRUE))

#match count matrix
gene_length <- gene_info2$gene_length[ match(gene_list, gene_info2$ensembl_gene_id) ]
gene_length_kb <- gene_length / 1000
rpk <- dataset_combined[3:17] / gene_length_kb #RPK

scaling_factors <- colSums(rpk)
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6 #TPM
tpm$gene_id_final <- dataset_combined$gene_id_final
######## plotting markers
#################################
ap_markers <- c("PARD3", "SOX2", "HES1", "HES5", "PROM1", "CENPF", "MKI67",
                "TJP1", "TJP2", "CDH2", "SLC1A3", "CRYAB", "PALS", "NUSAP", "CDH2", "ASPM", "KNL1")

bp_markers <- c("EOMES", "NEUROG2", "INSM1", "BCL11A", "PPP1R17",
                "PTPRZ1", "LIFR", "FEZF2", "PDGFD", "BMP7", "HOPX","TMEM14B")


neu_markers <- c("DCX", "TUBB3", "MAP2", "RBFOX3", "NEUN", "SYN1", "STMN2", "BCL11B", "SATB2", "SLC17A7", "SLC17A6")


ap_genes_found <- ap_markers[ap_markers %in% tpm$gene_id_final]
bp_genes_found <- bp_markers[bp_markers %in% tpm$gene_id_final]
neu_genes_found <- neu_markers[neu_markers %in% tpm$gene_id_final]

ap_score <- colMeans(tpm[tpm$gene_id_final %in% ap_genes_found, c(1:15)])
bp_score <- colMeans(tpm[tpm$gene_id_final %in% bp_genes_found, c(1:15)])
neu_score <- colMeans(tpm[tpm$gene_id_final %in% neu_genes_found, c(1:15)])
### AP org
boxplot(
  list(aRG = ap_score[c(1,3,5)], "bRG/BP" = bp_score[c(1,3,5)], Neu = neu_score[c(1,3,5)]),
  col = c("skyblue", "tomato", "lightgreen"),
  ylab = "Marker expression score (TPM)"
)
##GD90 aRG
boxplot(
  list(aRG = ap_score[7:15], "bRG/BP" = bp_score[7:15], Neu = neu_score[7:15]),
  col = c("skyblue", "tomato", "lightgreen"),
  ylab = "Marker expression score (TPM)"
)
###GD90 neu
boxplot(
  list(aRG = ap_score[c(2,4,6)], "bRG/BP" = bp_score[c(2,4,6)], Neu = neu_score[c(2,4,6)]),
  col = c("skyblue", "tomato", "lightgreen"),
  ylab = "Marker expression score (TPM)"
)

#################################### 
#euclidean distance of FACS samples

# Read the table
go_table <- read.delim("/home/basti/Downloads/GO_term_summary_20250731_115531.txt", header = F ,row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
# Extract, capitalize, and remove duplicates from the Symbol column
go_genes <- unique(toupper(go_table$V2))
# Append to your existing vector
neurodev_genes <- unique(c(neurodev_genes, go_genes))

genes_present <- neurodev_genes[neurodev_genes %in% rownames(logCPM)]
logCPM <- cpm(dge_TMM, log = TRUE, prior.count = 1)
logCPM <- logCPM[genes_present, ]
log_expr_t <- t(logCPM)
#####DO I HAVE TO TRANSPOSE THE MATRIX? SEE IT
dist_matrix <- as.matrix(dist(log_expr_t, method = "euclidean"))
#sample groups
aRG_samples <- grep("aRG", rownames(dist_matrix), value = TRUE)
Neu_samples <- grep("Neu", rownames(dist_matrix), value = TRUE)
d30_samples   <- grep("d30", rownames(dist_matrix), value = TRUE)
d40_samples   <- grep("d40", rownames(dist_matrix), value = TRUE)
d50_samples   <- grep("d50", rownames(dist_matrix), value = TRUE)

# store average distances for each individual brain
avg_distances_aRG <- data.frame(
  aRG = character(),
  stage = character(),
  distance = numeric()
)

avg_distances_Neu <- data.frame(
  Neu = character(),
  stage = character(),
  distance = numeric()
)

for (aRG in aRG_samples) {
  avg_distances_aRG <- rbind(avg_distances_aRG, data.frame(
    aRG = aRG, stage = "d30", distance = mean(dist_matrix[aRG, d30_samples])
  ))
  avg_distances_aRG <- rbind(avg_distances_aRG, data.frame(
    aRG = aRG, stage = "d40", distance = mean(dist_matrix[aRG, d40_samples])
  ))
  avg_distances_aRG <- rbind(avg_distances_aRG, data.frame(
    aRG = aRG, stage = "d50", distance = mean(dist_matrix[aRG, d50_samples])
  ))
}

for (Neu in Neu_samples) {
  avg_distances_Neu <- rbind(avg_distances_Neu, data.frame(
    Neu = Neu, stage = "d30", distance = mean(dist_matrix[Neu, d30_samples])
  ))
  avg_distances_Neu <- rbind(avg_distances_Neu, data.frame(
    Neu = Neu, stage = "d40", distance = mean(dist_matrix[Neu, d40_samples])
  ))
  avg_distances_Neu <- rbind(avg_distances_Neu, data.frame(
    Neu = Neu, stage = "d50", distance = mean(dist_matrix[Neu, d50_samples])
  ))
}

# average across all brain samples for each stage
summary_df_aRG <- avg_distances_aRG %>%
  group_by(stage) %>%
  summarise(
    mean_distance = mean(distance),
    sem = sd(distance) / sqrt(n())
  )

summary_df_Neu <- avg_distances_Neu %>%
  group_by(stage) %>%
  summarise(
    mean_distance = mean(distance),
    sem = sd(distance) / sqrt(n())
  )

####To aRG samples FACS
ggplot(summary_df_aRG, aes(x = stage, y = mean_distance, fill = stage)) +
  geom_col(width = 0.85) +
  geom_errorbar(aes(ymin = mean_distance - sem, ymax = mean_distance + sem),
                width = 0.2, size = 0.7) + coord_cartesian(ylim = c(40, 55)) +
  geom_text(
    aes(label = round(mean_distance, 2), y = mean_distance + 2.5),
    size = 4
  )+
  theme_minimal() +
  labs(title = "Average Euclidean Distance from FACS aRG Samples to Each Organoid Stage",
       y = "Average Euclidean Distance",
       x = "Organoid Stage") +
  theme(
    legend.position = "none",
    text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10))
  ) +
  scale_fill_manual(values = c("d30" = "skyblue", "d40" = "gold", "d50" = "tomato")) + geom_jitter(data = avg_distances_aRG,
                                                                                                   aes(x = stage, y = distance),
                                                                                                   width = 0.15, shape = 21, size = 2,
                                                                                                   fill = "black", alpha = 0.6, inherit.aes = FALSE)

####To neuronal samples FACS
ggplot(summary_df_Neu, aes(x = stage, y = mean_distance, fill = stage)) +
  geom_col(width = 0.85) +
  geom_errorbar(aes(ymin = mean_distance - sem, ymax = mean_distance + sem),
                width = 0.2, size = 0.7) + coord_cartesian(ylim = c(40, 62.5)) +
  geom_text(
    aes(label = round(mean_distance, 2), y = mean_distance + 2.5),
    size = 4
  )+
  theme_minimal() +
  labs(title = "Average Euclidean Distance from FACS Neuronal Samples to Each Organoid Stage",
       y = "Average Euclidean Distance",
       x = "Organoid Stage") +
  theme(
    legend.position = "none",
    text = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10))
  ) +
  scale_fill_manual(values = c("d30" = "skyblue", "d40" = "gold", "d50" = "tomato")) + geom_jitter(data = avg_distances_Neu,
                                                                                                   aes(x = stage, y = distance),
                                                                                                   width = 0.15, shape = 21, size = 2,
                                                                                                   fill = "black", alpha = 0.6, inherit.aes = FALSE)
