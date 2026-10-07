library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)
library(scales)
library(brms)
library("ggpubr")
library(ggtext)
setwd("/Users/xiongxia/Desktop/Drug_Interaction_Cross_Feeding/zResults_Data/Fig1_TECAN_screens")
X_values_all = read.xlsx("Total_data_screen.xlsx", sheet = 1)

mean_diff_g_single_drug = c(); X_diff = c()
for (i in 1:108) {
  cur_combo = X_values_all$Drug_Combination[i]
  rep_num = X_values_all$Replicate[i]
  
  co_current = X_values_all[X_values_all$Drug_Combination == cur_combo & X_values_all$Replicate == rep_num & X_values_all$Ecology == "Mutualism",]
  mono_current = X_values_all[X_values_all$Drug_Combination == cur_combo & X_values_all$Replicate == rep_num & X_values_all$Ecology == "Monoculture",]
  
  delta_gA = co_current$g_A[1] - mono_current$g_A[1]
  delta_gB = co_current$g_B[1] - mono_current$g_B[1]
  mean_diff_g_single_drug = c(mean_diff_g_single_drug, mean(delta_gA, delta_gB))
  X_diff = c(X_diff, co_current$X[1] - mono_current$X[1])
}

df = data.frame(mean_diff_g_single_drug, X_diff)

mean_X_vals <- X_values_all %>%
  group_by(Drug_Combination, Ecology) %>%
  summarise(mean_gA = mean(g_A, na.rm = TRUE), 
            mean_gB = mean(g_B, na.rm = TRUE),
            mean_X = mean(X, na.rm = TRUE),
            sd_X = sd(X, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(Ecology)

mean_diff_g_single_drug = c(); X_diff = c()
for (i in 1:36) {
  cur_combo = mean_X_vals$Drug_Combination[i]
  co_current = mean_X_vals[mean_X_vals$Drug_Combination == cur_combo & mean_X_vals$Ecology == "Mutualism",]
  mono_current = mean_X_vals[mean_X_vals$Drug_Combination == cur_combo & mean_X_vals$Ecology == "Monoculture",]
  
  delta_gA = co_current$mean_gA[1] - mono_current$mean_gA[1]
  delta_gB = co_current$mean_gB[1] - mono_current$mean_gB[1]
  mean_diff_g_single_drug = c(mean_diff_g_single_drug, mean(delta_gA, delta_gB))
  X_diff = c(X_diff, co_current$mean_X[1] - mono_current$mean_X[1])
}

df = data.frame(X_diff, mean_diff_g_single_drug)
ggplot(df, aes(x = mean_diff_g_single_drug, y = X_diff)) + 
  labs(x = NULL, y = NULL) +
  geom_point(size = 2, alpha = 0.4, color = "grey2")+
  geom_smooth(method = "lm")+
  #xlim(-1,-0.1)+
  geom_hline(yintercept= 0, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = 'black'),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title.y = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  ) + stat_cor(method = "pearson", cor.coef.name = "rho")
