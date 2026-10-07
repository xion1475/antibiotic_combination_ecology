setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig2_Computational")
interaction = read.csv("gene_perturbation_drug_interaction_data.csv")
interaction$Type = "No_interaction"
interaction$Type = ifelse(interaction$X < -0.05, "Antagonistic", interaction$Type)
interaction$Type = ifelse(interaction$X > 0.05, "Synergistic", interaction$Type)
write.csv(interaction, "gene_perturbation_drug_interaction_data_XX20260603.csv")
focus_sub_type = interaction[interaction$strain == "S" & interaction$interaction == "monoculture",]
table(focus_sub_type$Type)

setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig2_Computational")
interaction2 = read.xlsx("FigS6_raw_data_ATB20250828.xlsx", sheet = 1)
total_gr = read.xlsx("FigS6_raw_data_ATB20250828.xlsx", sheet = 2)

interaction2$g_diff_S_co_minus_mono = 0
interaction2$mu_diff_S_co_minus_mono = 0
interaction2$g_diff_E_co_minus_mono = 0
interaction2$mu_diff_E_co_minus_mono = 0

for (i in 1:nrow(interaction2)) {
  cur_gene = interaction2$gene[i]
  cur_combo = total_gr[total_gr$gene_xy == cur_gene,]
  cur_geneX_S_diff = cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "coculture",]$gx_norm - cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "monoculture",]$gx_norm
  cur_geneY_S_diff = cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "coculture",]$gy_norm - cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "monoculture",]$gy_norm
  
  cur_geneX_mu_diff = cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "coculture",]$gx - cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "monoculture",]$gx
  cur_geneY_mu_diff = cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "coculture",]$gy - cur_combo[cur_combo$strain == "S" & cur_combo$interaction == "monoculture",]$gy
  interaction2$g_diff_S_co_minus_mono[i] = (cur_geneX_S_diff+cur_geneY_S_diff)/2
  interaction2$mu_diff_S_co_minus_mono[i] = (cur_geneX_mu_diff+cur_geneY_mu_diff)/2
  
  cur_geneX_S_diff = cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "coculture",]$gx_norm - cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "monoculture",]$gx_norm
  cur_geneY_S_diff = cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "coculture",]$gy_norm - cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "monoculture",]$gy_norm
  
  cur_geneX_mu_diff = cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "coculture",]$gx - cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "monoculture",]$gx
  cur_geneY_mu_diff = cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "coculture",]$gy - cur_combo[cur_combo$strain == "E" & cur_combo$interaction == "monoculture",]$gy
  interaction2$g_diff_E_co_minus_mono[i] = (cur_geneX_S_diff+cur_geneY_S_diff)/2
  interaction2$mu_diff_E_co_minus_mono[i] = (cur_geneX_mu_diff+cur_geneY_mu_diff)/2
}

#E coli growth differences
ggplot(interaction2, aes(x = mu_diff_E_co_minus_mono, y = X_co_E_minus_mono_E)) + 
  labs(x = NULL, y = NULL) +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "loess")+
  #xlim(-1,-0.1)+
  geom_hline(yintercept= 0, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) #+ stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

ggplot(interaction2, aes(x=g_diff_E_co_minus_mono, y=X_co_E_minus_mono_E)) + 
  labs(x = NULL, y = NULL) +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "loess")+
  geom_hline(yintercept= 0, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) #+ stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#S. enterica growth differences
ggplot(interaction2, aes(x = mu_diff_S_co_minus_mono, y = X_co_S_minus_mono_S)) + 
  labs(x = NULL, y = NULL) +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "loess")+
  #xlim(-1,-0.1)+
  geom_hline(yintercept= 0, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) #+ stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

ggplot(interaction2, aes(x=g_diff_S_co_minus_mono, y=X_co_S_minus_mono_S)) + 
  labs(x = NULL, y = NULL) +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "loess")+
  geom_hline(yintercept= 0, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) #+ stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")



