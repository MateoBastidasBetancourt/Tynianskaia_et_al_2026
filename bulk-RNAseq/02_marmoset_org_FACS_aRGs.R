library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
###starting from read count table
dataset_aRG <- read.table("raw_counts_AP_orgs.csv", sep = "\t", header = TRUE)
pheno_group <- read.table("pheno_data_aRG_orgs.csv", sep = "\t", header = TRUE) #open sample group file in metadata folder
pheno_group$stage <- factor(pheno_group$stage, levels = c("d30", "d40", "d50"))
pheno_group$line <- factor(pheno_group$line)

design <- model.matrix(~ 0 + stage + line, data = pheno_group)
dge_obj <- DGEList(counts=dataset_aRG[,3:11], group=factor(pheno_group[,2]), genes=dataset_aRG$gene_id_final)
identical(pheno_group$stage_line, colnames(dataset_aRG[,3:11])) #See if order in count table is identical than in the metadata
dge_TMM <- calcNormFactors(dge_obj, method="upperquartile")
keep <- filterByExpr(dge_TMM, design = design, min.count = 15, min.prop = 0.7, min.total.count = 15)
dge_TMM <- dge_TMM[keep, keep.lib.sizes=FALSE]
dim(dge_TMM)
rownames(dge_TMM) <- make.unique(dge_TMM$genes$genes)
disp <- estimateDisp(dge_TMM, design, robust=TRUE)
plotBCV(disp)
plotMDS(dge_TMM)

fit <- glmQLFit(disp, design, robust=TRUE)

###########STATISTICAL TEST - D50 
CONTRASTS <- makeContrasts(d30vsd50 = staged30 -staged50, d30vsd40 = staged30 -staged40, d40vsd50 = staged40 - staged50 ,d30vsd40_50 = staged30 - (staged40+staged50)/2, d30_40vsd50 = (staged30+staged40)/2 - staged50, levels=design)

for (i in 1:ncol(CONTRASTS)){
  
  contrast_name <- colnames(CONTRASTS)[i]
  
  # Perform the test
  current.glmQLFTest <- glmQLFTest(fit, contrast = CONTRASTS[, i])
  
  # Assign the result to a variable with the contrast name
  assign(contrast_name, current.glmQLFTest)
}


#####PVALUE CORRECTION AND SORTING GENES ACCORDING TO THEIR PVALUE
result_d30vsd40 <- topTags(d30vsd40,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d40vsd50 <- topTags(d40vsd50,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d30vsd50 <- topTags(d30vsd50,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d30_40vsd50 <- topTags(d30_40vsd50,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_d30vsd40_50 <- topTags(d30vsd40_50,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")

######################
filter_results <- function(result) {
  filtered_result <- result$table[
    (result$table$logFC > 0.5 | result$table$logFC < -0.5) & result$table$FDR < 0.2, ]
  return(filtered_result)
}

# Apply the function to each topTags object
result_d30vsd40 <- filter_results(result_d30vsd40)
result_d40vsd50 <- filter_results(result_d40vsd50)
result_d30vsd50 <- filter_results(result_d30vsd50)
result_d30_40vsd50 <- filter_results(result_d30_40vsd50)
result_d30vsd40_50 <- filter_results(result_d30vsd40_50)

dfs <- list(result_d30vsd40 = result_d30vsd40, result_d40vsd50 = result_d40vsd50,
            result_d30vsd50 = result_d30vsd50, result_d30_40vsd50 = result_d30_40vsd50, result_d30vsd40_50=result_d30vsd40_50)

### Extract cpm values for all genes (without applying the filter_results function)
log_cpm <- cpm(dge_TMM, log = TRUE)
common_rows <- intersect(rownames(result_d30vsd50), rownames(log_cpm))
log_cpm_all_genes <- log_cpm[common_rows, , drop = FALSE]
result_d30vsd50_all_genes <- cbind(result_d30vsd50, log_cpm_all_genes)

##### SAVE RESULTS
result_d30vsd50_all_genes$d30_av <- rowMeans(result_d30vsd50_all_genes[c(7,10,13)])
result_d30vsd50_all_genes$d40_av <- rowMeans(result_d30vsd50_all_genes[c(8,11,14)])
result_d30vsd50_all_genes$d50_av <- rowMeans(result_d30vsd50_all_genes[c(9,12,15)])
result_d30vsd50_all_genes <- result_d30vsd50_all_genes[ , -c(7:15)]
write.table(result_d30vsd50_all_genes, "result_d30vsd50_cpm_all_genes.txt",row.names = F, col.names = T, quote = F)
###########PLOTS, DECLARE GENES THAT WERE USED FOR THESE PLOTS. PUT IN THE TEXT MB?

result_d30vsd50_cpm_all_genes <- read.csv("~/result_d30vsd50_cpm_all_genes.txt", sep="")
#1: all genes
result_fixed <- result_d30vsd50_cpm_all_genes %>%
  rename(avg_logCPM = logCPM)

result_long <- result_fixed %>%
  pivot_longer(cols = c(d30_av, d40_av, d50_av),
               names_to = "timepoint",
               values_to = "logCPM") %>%
  mutate(timepoint = recode(timepoint,
                            d30_av = "d30",
                            d40_av = "d40",
                            d50_av = "d50"))

result_long <- result_long %>%
  mutate(time_numeric = case_when(
    timepoint == "d30" ~ 30,
    timepoint == "d40" ~ 40,
    timepoint == "d50" ~ 50
  ))
############################## Grouping and plotting genes
genes_group1 <- c("THBS4", "THBS2", "CEMIP", "ALDH1A3")

plot1 <- result_long %>%
  filter(genes %in% genes_group1) %>%
  ggplot(aes(x = timepoint, y = logCPM, group = genes, color = genes)) +
  geom_line(size = 1.8) +
  geom_point(size = 3) +
  labs(title = "Fate decision and delamination",
       x = "Timepoint", y = "logCPM") +
  theme_minimal()+
  theme(
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 25),
    axis.title.y = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
    plot.title = element_text(size = 30),
    axis.text.y = element_text(size = 25), legend.position = c(0.45, 0.5)
  ) + scale_x_discrete(expand = c(0,0)) + coord_cartesian(xlim = c(0.8, 7))  # assuming your x axis is a factor with 3 levels
plot1
