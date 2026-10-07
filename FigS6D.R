library(ggplot2)
library("ggpubr")
library(scales)
library("openxlsx")
library(ggtext)

setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig2_Computational/Ave_raw_data")

master_df = read.csv("master_df_XX20250907_correct.csv")
combo_df = master_df[master_df$gene_combo == "Double_gene",]
combo_df$Pathway_type = "Multi"
pathway_info = read.xlsx("FigS6B_raw.xlsx", sheet = 1)

for (i in 1:nrow(combo_df)) {
  cur_gene_X = combo_df$Gene_X[i]
  cur_gene_Y = combo_df$Gene_Y[i]
  path_X_num = pathway_info$Pathway_number[which(pathway_info$Gene == cur_gene_X)[1]]
  path_Y_num = pathway_info$Pathway_number[which(pathway_info$Gene == cur_gene_Y)[1]]
  if (path_X_num == 1 & path_Y_num == 1) {
    X_pathway = pathway_info$Pathway1[which(pathway_info$Gene == cur_gene_X)]
    Y_pathway = pathway_info$Pathway1[which(pathway_info$Gene == cur_gene_Y)]
    if (X_pathway == Y_pathway) {combo_df$Pathway_type[i] = "Same"}
  }
}
table(combo_df$Pathway_type)
write.csv(combo_df, "combo_df_XX20250911.csv")

setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig2_Computational")
types_total_data = read.xlsx("FigS5_raw_data_ATB20250828.xlsx", sheet = 1)
table(types_total_data$S_co)


multi_combo = types_total_data[types_total_data$Pathway_type == "Multi",]
same_combo = types_total_data[types_total_data$Pathway_type == "Same",]

table(multi_combo$E_o)
table(same_combo$E_co)
a = data.frame(table(multi_combo$E_mono), table(same_combo$E_mono))
a
a = a[,-1]
a = a[,-2]
fisher.test(a)

mono_E = types_total_data$E_mono_X
co_E = types_total_data$E_co_X
t.test(mono_E,co_E)
