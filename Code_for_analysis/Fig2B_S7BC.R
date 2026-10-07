library(ggplot2)
library("ggpubr")
library(scales)
library("openxlsx")
library(ggtext)
setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig2_Computational/Ave_raw_data")

single_gene_num = 40
gene_pairs_num = 780

#Biomass data analysis
biomass_data = read.csv("gene_perturbation_biomass_data.csv")
coculture_biomass = biomass_data[biomass_data$interaction=="coculture",]
gene_list = unique(biomass_data$gene)
E_yield = c(); S_yield = c(); SE_ratio = c()
for (gn in gene_list) {
  current_biomass = coculture_biomass[coculture_biomass$gene==gn,]
  Ey = current_biomass$max_yield[1]
  Sy = current_biomass$max_yield[2]
  
  E_yield = c(E_yield, Ey)
  S_yield = c(S_yield, Sy)
  SE_ratio = c(SE_ratio, Sy/Ey)
}
master_df = data.frame(gene_list, E_yield, S_yield, SE_ratio)
write.csv(master_df, "master_df_XX20250904.csv")
master_df = read.csv("master_df_XX20250904.csv")

#Add interaction data
X_diff_dist = read.xlsx("../FigS6_raw_data_ATB.xlsx", sheet = 1)
master_df$X_co_E_minus_mono_E = NA
master_df$E_less_syn_value_YN = NA
master_df$E_less_syn_type_YN = NA
master_df$gene_combo = "Single_gene"
master_df$ES_total_yield = master_df$E_yield + master_df$S_yield

for (i in (single_gene_num + 1) : (nrow(master_df) - 1)) {
  master_df$gene_combo[i] = "Double_gene"
  gn = master_df$gene_list[i]
  ln_num = which(X_diff_dist$gene == gn)
  master_df$X_co_E_minus_mono_E[i] = X_diff_dist$X_co_E_minus_mono_E[ln_num]
  master_df$E_less_syn_value_YN[i] = X_diff_dist$E_less_syn_value[ln_num]
  master_df$E_less_syn_type_YN[i] = X_diff_dist$E_less_syn_type[ln_num]
}

#Metabolic exchange analysis
fluxd = read.csv("gene_perturbation_flux_data.csv")
omitted_data = read.csv("omitted_mid_log_fluxes.csv")
fluxd <- rbind(fluxd, omitted_data)
#unique_reactions = unique(fluxd$reaction)
#write.table(unique_reactions,"unique_reactions_list.txt",sep="\t")

#Made metabolite information table for each species and metabolite.
metabolites_we_care = c("EX_met__L_e")
for (mblt in metabolites_we_care) {
  metab_list = fluxd[fluxd$reaction==mblt,]
  gal_co = metab_list[metab_list$interaction=="coculture",]
  E_co_flux = gal_co[gal_co$strain=="E",]
  S_co_flux = gal_co[gal_co$strain=="S",]
  
  write.csv(E_co_flux,paste0(mblt, "_flux_in_E_co.csv"))
  write.csv(S_co_flux,paste0(mblt, "_flux_in_S_co.csv"))
}

#Analyze data
species_list = c("E", "S")
for (species in species_list) {
  for (mblt in metabolites_we_care) {
    if (species == "E") {
      species_suffix = "_flux_in_E_co"
    } else {
      species_suffix = "_flux_in_S_co"
    }
    current_flux_file = read.csv(paste0(mblt, species_suffix,".csv"))
    for (gn_num in 1:(gene_pairs_num + single_gene_num)) {
      current_gene = master_df$master_df.gene_list[gn_num]
      current_gene_flux = current_flux_file[current_flux_file$gene==current_gene,]
      current_gene_flux$flux = abs(current_gene_flux$flux)
      master_df[gn_num,paste0(mblt,species_suffix,"_max")] = max(current_gene_flux$flux)
      master_df[gn_num,paste0(mblt,species_suffix,"_median")] = median(current_gene_flux$flux)
      master_df[gn_num,paste0(mblt,species_suffix,"_mean")] = mean(current_gene_flux$flux)
      
      mid_log_flux = current_gene_flux[current_gene_flux$cycle == current_gene_flux$time_to_mid_log,]
      master_df[gn_num,paste0(mblt,species_suffix,"_mid_log")] = mid_log_flux$flux[1]
      master_df[gn_num,paste0(mblt,species_suffix,"_sum")] = sum(current_gene_flux$flux)
    }
  }
}


#Plot data
healthy_combo = read.csv("healthy_combo_XX20260603.csv")
single_master = read.csv("single_master_XX20260603.csv")

wt_SE_ratio = 0.346949031
wt_E_co_gal_flux = 0.7719825
wt_E_co_ac_flux = 6.399907e-02
wt_E_co_met_flux = 2.292389e-02
wt_S_co_met_flux = 0.074528004
wt_E_carbon_flux = 4.75989314
wt_S_carbon_flux = 11.69
 
#Fig2B
healthy_combo$E_carbon_higher_than_single = "No"
healthy_combo$E_carbon_higher_single_flux = 0
healthy_combo$E_carbon_X = 0
wt_target = wt_E_carbon_flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_SE_target = healthy_combo$E_carbon_flux_max_total[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$E_carbon_flux_max_total
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$E_carbon_flux_max_total
  higher_single_target = max(single_gene_X_target, single_gene_Y_target)
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_SE_target / wt_target
  healthy_combo$E_carbon_X[i] = max(1 - combo_rel_target / gene_X_rel_target / gene_Y_rel_target, -1)
  healthy_combo$E_carbon_higher_single_flux[i] = higher_single_target
  if (cur_combo_SE_target > higher_single_target) {healthy_combo$E_carbon_higher_than_single[i] = "Yes"}
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=E_carbon_flux_max_total / E_carbon_linear_exp)) + 
  labs(x = NULL, y = "Carbon excretion relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey1")+
  guides()+
  geom_smooth(method = "lm")+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1), breaks = c(10^-1, 10^0,10^1,10^2))+
  theme(
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.x = element_markdown(), 
    axis.title.y = element_markdown()
  )# + stat_compare_means(method = "wilcox.test")
cor.test(healthy_combo$E_carbon_flux_max_total / healthy_combo$E_carbon_linear_exp, healthy_combo$X_co_E_minus_mono_E,method = "spearman", exact = FALSE)


#S. enterica methionine excretion
healthy_combo$S_met_higher_than_single = "No"
healthy_combo$S_met_higher_single_flux = 0
healthy_combo$S_met_linear_exp = 0
healthy_combo$S_met_X = 0
wt_target = wt_S_co_met_flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_SE_target = healthy_combo$EX_met_L_e_flux_in_S_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = as.numeric(single_master[single_master$gene_list==cur_geneX,]$EX_met_L_e_flux_in_S_co_max)
  single_gene_Y_target = as.numeric(single_master[single_master$gene_list==cur_geneY,]$EX_met_L_e_flux_in_S_co_max)
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  healthy_combo$S_met_X[i] = max(1 - combo_rel_target / gene_X_rel_target / gene_Y_rel_target, -1)
  healthy_combo$S_met_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target * wt_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y = EX_met_L_e_flux_in_S_co_max / S_met_linear_exp))+ 
  labs(x = NULL, y = "Methionine excretion relative to null expectation") +
  geom_point(size = 2, alpha = 0.4)+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1), breaks = c(10^-1, 10^0,10^1,10^2))+
  theme(
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.x = element_markdown(), 
    axis.title.y = element_markdown()
  )+ stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")


#S:E ratio
healthy_combo$SE_ratio_higher_than_single = "No"
healthy_combo$SE_ratio_higher_singles = 0
wt_target = wt_SE_ratio
for (i in 1:nrow(healthy_combo)) {
  cur_combo_SE_ratio = healthy_combo$SE_ratio[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$SE_ratio
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$SE_ratio
  higher_single_target = max(single_gene_X_target, single_gene_Y_target)
  healthy_combo$SE_ratio_higher_singles[i] = higher_single_target
  if (cur_combo_SE_ratio > higher_single_target) {healthy_combo$SE_ratio_higher_than_single[i] = "Yes"}
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=SE_ratio / SE_ratio_linear_exp)) + 
  labs(x = NULL, y = "S:E ratio relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1), breaks = c(10^-1, 10^0,10^1,10^2))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")


#E. coli acetate
healthy_combo$E_ac_linear_exp=0
wt_target = wt_E_co_ac_flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$EX_ac_e_flux_in_E_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$EX_ac_e_flux_in_E_co_max
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$EX_ac_e_flux_in_E_co_max
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$E_ac_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target * wt_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=EX_ac_e_flux_in_E_co_max / E_ac_linear_exp)) + 
  labs(x = NULL, y = "Acetate excretion relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")


#S:E growth rate ratio
single_master$g_gr_ratio = single_master$g_S / single_master$g_E
wt_target = 1
healthy_combo$g_gr_ratio_linear_exp=0

for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$g_gr_ratio[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$g_gr_ratio
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$g_gr_ratio
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$g_gr_ratio_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=g_gr_ratio / g_gr_ratio_linear_exp)) + 
  labs(x = NULL, y = "Growth rate (S:E) ratio relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#E. coli Galactose excretion
healthy_combo$E_gal_linear_exp=0
wt_target = wt_E_co_gal_flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$EX_gal_e_flux_in_E_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$EX_gal_e_flux_in_E_co_max
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$EX_gal_e_flux_in_E_co_max
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$E_gal_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target * wt_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=EX_gal_e_flux_in_E_co_max / E_gal_linear_exp)) + 
  labs(x = NULL, y = "Galactose excretion relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#S. enterica galactose intake
healthy_combo$S_gal_linear_exp=0
wt_target = 1.8958767 #S no_drug galactose intake flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$EX_gal_e_flux_in_S_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$EX_gal_e_flux_in_S_co_max
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$EX_gal_e_flux_in_S_co_max
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$S_gal_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target * wt_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=EX_gal_e_flux_in_S_co_max / S_gal_linear_exp)) + 
  labs(x = NULL, y = "Galactose intake relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#S. enterica acetate intake
healthy_combo$S_ac_linear_exp=0
wt_target = 0.157 #S no_drug acetate intake flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$EX_ac_e_flux_in_S_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$EX_ac_e_flux_in_S_co_max
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$EX_ac_e_flux_in_S_co_max
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$S_ac_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=EX_ac_e_flux_in_S_co_max / S_ac_linear_exp)) + 
  labs(x = NULL, y = "Acetate intake relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#E. coli methionine intake
healthy_combo$E_met_linear_exp=0
wt_target = wt_E_co_met_flux #E no_drug methionine intake flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$EX_met_L_e_flux_in_E_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$EX_met_L_e_flux_in_E_co_max
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$EX_met_L_e_flux_in_E_co_max
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$E_met_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target * wt_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=EX_met_L_e_flux_in_E_co_max / E_met_linear_exp)) + 
  labs(x = NULL, y = "Methionine intake relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#S. enterica methionine excretion
healthy_combo$S_met_linear_exp=0
wt_target = wt_S_co_met_flux #S no_drug methionine excretion flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$EX_met_L_e_flux_in_S_co_max[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$EX_met_L_e_flux_in_S_co_max
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$EX_met_L_e_flux_in_S_co_max
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$S_met_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target * wt_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=EX_met_L_e_flux_in_S_co_max / S_met_linear_exp)) + 
  labs(x = NULL, y = "Methionine excretion relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")

#S carbon intake
healthy_combo$S_carbon_flux_linear_exp=0
wt_target = wt_S_carbon_flux #S no_drug acetate intake flux
for (i in 1:nrow(healthy_combo)) {
  cur_combo_ace_target = healthy_combo$S_carbon_flux_total[i]
  cur_geneX = healthy_combo$Gene_X[i]
  cur_geneY = healthy_combo$Gene_Y[i]
  single_gene_X_target = single_master[single_master$gene_list==cur_geneX,]$S_carbon_flux_total
  single_gene_Y_target = single_master[single_master$gene_list==cur_geneY,]$S_carbon_flux_total
  
  gene_X_rel_target = single_gene_X_target / wt_target
  gene_Y_rel_target = single_gene_Y_target / wt_target
  combo_rel_target = cur_combo_ace_target / wt_target
  healthy_combo$S_carbon_flux_linear_exp[i] = gene_X_rel_target * gene_Y_rel_target
}
ggplot(healthy_combo, aes(x=X_co_E_minus_mono_E, y=S_carbon_flux_total / S_carbon_flux_linear_exp)) + 
  labs(x = NULL, y = "Carbon intake relative to null expectation") +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  geom_hline(yintercept= 1, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(labels = label_log(digits = 1))+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "spearman", cor.coef.name = "rho") # + stat_compare_means(method = "wilcox.test")
