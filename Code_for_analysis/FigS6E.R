library(ggplot2)
library("ggpubr")
library(scales)
library("openxlsx")
library(ggtext)
setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig2_Computational/Ave_raw_data")

single_gene_num = 40
gene_pairs_num = 780

master_df = read.csv("master_df_XX20250907_correct.csv")
write.csv(healthy_combo, "healthy_combo_XX20250907.csv")
healthy_combo = read.csv("healthy_combo_XX202509011.csv")
combo_master = master_df[master_df$gene_combo=="Double_gene",]
single_master = master_df[master_df$gene_combo=="Single_gene",]
healthy_combo = combo_master[combo_master$ES_total_yield>2E-7,]
healthy_combo$EX_met_L_e_flux_in_E_co_max = as.numeric(healthy_combo$EX_met_L_e_flux_in_E_co_max)
healthy_combo$EX_met_L_e_flux_in_S_co_max = as.numeric(healthy_combo$EX_met_L_e_flux_in_S_co_max)
master_df$EX_met_L_e_flux_in_S_co_max = as.numeric(master_df$EX_met_L_e_flux_in_S_co_max)
master_df$EX_met_L_e_flux_in_E_co_max = as.numeric(master_df$EX_met_L_e_flux_in_E_co_max)


wt_SE_ratio = 0.346949031
wt_E_co_gal_flux = 0.7719825
wt_E_co_ac_flux = 6.399907e-02
wt_E_co_met_flux = 2.292389e-02
wt_S_co_met_flux = 0.074528004
wt_E_carbon_flux = 4.75989314
wt_S_carbon_flux = 11.6892602
wt_S_co_gal_flux = 1.8965
wt_S_co_ac_flux = 1.572102e-01

healthy_master = master_df[master_df$ES_total_yield>2E-7,]

my_comparisons <- list(c("Combo: less synergistic", "Single"),
                       c("Combo: not less synergistic", "Single"))

figS6Ea = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y=log(SE_ratio,10), fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = log(1,10), linetype="dashed", color = "grey")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_hline(yintercept = log(wt_SE_ratio, 10), linetype = "dashed", color = "red")+
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.6) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "log<sub>10</sub>(S:E ratio)") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Eb = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y=E_carbon_flux_max_total, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_E_carbon_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Carbon excretion flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Ec = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = 10^EX_met_L_e_flux_in_S_co_max, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = 0.077, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.6) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Methionine excretion flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Ed = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = g_S-g_E, fill = E_less_syn_value_YN)) + 
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype="dashed", color = "red")+
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "S:E growth rate (g) diff.") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Ee = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = 10^EX_met_L_e_flux_in_E_co_max, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_E_co_met_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Methionine uptake flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Ef = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = EX_gal_e_flux_in_E_co_max, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_E_co_gal_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Galactose excretion flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Eg = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = EX_ac_e_flux_in_E_co_max, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_E_co_ac_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Acetate excretion flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Eh = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = EX_gal_e_flux_in_S_co_max, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_S_co_gal_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Galactose uptake flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6Ej = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = EX_ac_e_flux_in_S_co_max, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_S_co_ac_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Acetate uptake flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

figS6El = ggplot(healthy_master, aes(x=E_less_syn_value_YN, y = S_carbon_flux_total, fill = E_less_syn_value_YN)) + 
  geom_hline(yintercept = wt_S_carbon_flux, linetype="dashed", color = "red")+
  geom_jitter(width = 0.3, size = 2, alpha = 0.4, color = "grey") +
  geom_boxplot(outlier.shape=NA,width=0.6, alpha = 0.7) +
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#7CAE00"))+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="purple") +
  labs(x = NULL, y = "Carbon uptake flux") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(), 
    axis.text.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2, label = "p.adjust", p.adjust.method = "tukey") # #Wilcoxon P=1.1E-5 

ggarrange(figS6Eb, figS6Ee, figS6Ef, figS6Eg, figS6Ea,
          figS6El, figS6Ec, figS6Eh, figS6Ej, figS6Ed, 
          ncol = 5, nrow = 2) #FigS6E
