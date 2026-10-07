library(ggplot2)
library("ggpubr")
library(scales)
library("openxlsx")
library(ggtext)
setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig1_TECAN_screens")
master_exp_FLP = read.xlsx("Total_data_screen.xlsx", sheet = 7)
combo = master_exp_FLP[master_exp_FLP$Condition_type == "Double",]
single = master_exp_FLP[master_exp_FLP$Condition_type == "Single",]

combo$SE_ratio_linear_exp = 0
for (i in 1:nrow(combo)) {
  cur_rep = combo$Replicate[i]
  cur_geneX = combo$DrugX[i]
  cur_geneY = combo$Drug_Y[i]
  single_gene_X_target = single[single$DrugX == cur_geneX & single$Replicate == cur_rep,]$Final_YFP_CFP_Ratio
  single_gene_Y_target = single[single$DrugX == cur_geneY & single$Replicate == cur_rep,]$Final_YFP_CFP_Ratio
  cur_wt_ratio = single[single$DrugX == "No_Drug" & single$Replicate == cur_rep,]$Final_YFP_CFP_Ratio
  
  combo$SE_ratio_linear_exp[i] = single_gene_X_target / cur_wt_ratio * single_gene_Y_target
}
combo$SE_ratio_rel_to_linear = combo$Final_YFP_CFP_Ratio / combo$SE_ratio_linear_exp
write.csv(combo, "combo_exp_SE_ratio.csv")

combo2 = read.xlsx("Total_data_screen.xlsx", sheet = 7)
combo2 = combo2[combo2$Mutualism_Drug_Condition!="GEN+TMP",]
combo2$Mutualism_Drug_Condition <- reorder(combo2$Mutualism_Drug_Condition, combo2$Final_YFP_CFP_Ratio, mean)
ggplot(combo2,aes(x = Mutualism_Drug_Condition, y = Final_YFP_CFP_Ratio, fill = Less_syn)) +
  labs(x = NULL, y = "S:E ratios based on FLP")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_hline(yintercept= 0.38, linetype="dashed", color = "grey")+
  geom_boxplot()+
  geom_jitter()+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_markdown(), 
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )  

mean_SE_ratio_FLP <- combo2 %>%
  group_by(Mutualism_Drug_Condition) %>%
  summarise(mean_SE_rel = mean(Final_YFP_CFP_Ratio),
            median_SE_rel = median(Final_YFP_CFP_Ratio))
write.csv(mean_SE_ratio_FLP,"mean_SE_ratio_FLP.csv")
means = read.csv("mean_SE_ratio_FLP.csv")
wilcox.test(mean_SE_rel~Less_syn,data=means)
wilcox.test(mean_SE_rel~Altered,data=means)


ggplot(combo2,aes(x = Mutualism_Drug_Condition, y = Final_YFP_CFP_Ratio, fill = Altered)) +
  labs(x = NULL, y = "S:E ratios based on FLP")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_hline(yintercept= 0.38, linetype="dashed", color = "grey")+
  geom_boxplot()+
  geom_jitter()+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_markdown(), 
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )  

