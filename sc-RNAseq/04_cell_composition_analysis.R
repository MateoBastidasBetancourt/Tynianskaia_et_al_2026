library(speckle)
library(ggplot2)
library(dplyr)
cell_info_d30 <- read.table("/home/basti/Downloads/cell_clusters_composition_analysis.tab",sep = "\t", header = T)
cell_info_d50 <- read.table("/home/basti/Downloads/cell_clusters_composition_analysis_d50orgs.tab",sep = "\t", header = T)
cell_info_brains <- read.table("/home/basti/Downloads/cell_clusters_composition_analysis_brains.tab",sep = "\t", header = T)

target_types <- c("aRG1", "aRG2", "aRG3", "aRG4") #insert desired celltypes

# filter the data frame
cell_info <- cell_info[cell_info$harmony_clusters_names %in% target_types, ]

# plot stacked barplots of average celltype proportions
plot_cell_type <- plotCellTypeProps(clusters=cell_info$harmony_clusters_names, sample=cell_info$species) 

#### barplots with interstage differences
cell_info_d30$stage <- "d30"
cell_info_d50$stage <- "d50"
cell_info_brains$stage <- "GD90"

cell_info_d30_marmoset <- cell_info_d30[cell_info_d30$species == "Marmoset",c(1,2,4,5)]
colnames(cell_info_d30_marmoset) <- colnames(cell_info_d50)

cell_info <- rbind(cell_info_d30_marmoset, cell_info_d50, cell_info_brains)


cell_info <- cell_info %>%
  mutate(
    gross_celltype = case_when(
      seurat_clusters_names %in% c("IPC1","IPC2","IPC3", "IPC") ~ "IPC",
      seurat_clusters_names %in% c("bRG1", "bRG2", "bRG3","bRG", "tRG/bRG") ~ "bRG",
      seurat_clusters_names %in% c("aRG", "aRG1", "aRG2", "aRG3", "aRG4") ~ "aRG",
      seurat_clusters_names %in% c("EN-DL", "EN-UL", "EN", "IN", "IN1", "IN2") ~ "Neu",
      TRUE ~ NA_character_  # Optional: set to NA for other clusters
    )
  )

#


prop_df <- cell_info %>%
  dplyr::count(Sample, gross_celltype) %>%
  group_by(Sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


prop_df <- prop_df %>%
  mutate(stage = case_when(
    grepl("^S[1-3]$", Sample) ~ "GD90",
    grepl("^cj", Sample) ~ "d30",
    TRUE ~ "d50"
  ))

summary_df <- prop_df %>%
  group_by(stage, gross_celltype) %>%
  summarise(
    mean_prop = mean(proportion),
    se = sd(proportion) / sqrt(n())
  )

summary_df <- summary_df %>% 
  dplyr::filter(!is.na(gross_celltype))


ggplot(summary_df, aes(x = gross_celltype, y = mean_prop, fill = stage)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black",
    size = 0.3
  ) +
  geom_errorbar(
    aes(ymin = mean_prop - se, ymax = mean_prop + se),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  scale_fill_brewer(palette = "Set1") +  
  labs(
    x = "Cluster",
    y = "Mean Proportion",
    fill = "Stage"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 15)
  )
