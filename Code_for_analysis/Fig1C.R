library(dplyr)
library(ggplot2)
library(openxlsx)

# ------------------------------
# 1. Make sure ecology has the right order
# ------------------------------
setwd("/Users/xiongxia/Desktop")
df = read.xlsx(""Xiong_et_al_antibiotic_combination_ecology_raw_data.xlsx, sheet = 2)

# ------------------------------
# 2. Rank drug pairs by co-culture X value
#    (here I use the mean co-culture X; change to median if preferred)
# ------------------------------
pair_order <- df %>%
  filter(Ecology == "Mutualism") %>%
  group_by(Drug_Combination) %>%
  summarise(co_mean = mean(X, na.rm = TRUE), .groups = "drop") %>%
  arrange(co_mean) %>%
  pull(Drug_Combination)

df <- df %>%
  mutate(Drug_Combination = factor(Drug_Combination, levels = pair_order))

# ------------------------------
# 3. Wilcoxon test within each pair
# ------------------------------
pvals <- df %>%
  group_by(Drug_Combination) %>%
  summarise(
    p_value = t.test(X ~ Ecology, exact = FALSE)$p.value,
    y_pos   = max(X_g, na.rm = TRUE) * 1.08,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("P = ", signif(p_value, 3))
  )

# ------------------------------
# 4. Plot: one panel per pair
# ------------------------------

ggplot(df, aes(x = Drug_Combination, y = X)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey")+
  geom_boxplot(
    aes(fill = Ecology),
    position = position_dodge(width = 0),
    width = 1.2,
    outlier.shape = NA,
    alpha = 0.5,
    color = "black"
  ) +
  geom_point(aes(shape = Ecology), position = position_jitterdodge(jitter.width = 1, dodge.width = 0),
    size = 2,
    alpha = 0.4,
    color = "black",
    show.legend = TRUE
  ) +
  scale_fill_manual(values = c("Monoculture" = "cyan", "Mutualism" = "yellow")) +
  labs(x = NULL, y = expression(italic(X))) +
  theme_bw() +
  theme(
    #legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

df_mean <- df %>%
  group_by(Drug_Combination, Ecology) %>%
  summarise(mean_X = mean(X, na.rm = TRUE), .groups = "drop")

ggplot(df, aes(x = Drug_Combination, y = X)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "grey") +
  
  geom_boxplot(
    aes(fill = Ecology),
    position = position_dodge(width = 0),
    width = 1.2,
    outlier.shape = NA,
    alpha = 0.5,
    color = "black",
    median.colour = NA   # hide median line
  ) +
  
  # add mean as the middle line
  geom_crossbar(
    data = df_mean,
    aes(x = Drug_Combination, y = mean_X, ymin = mean_X, ymax = mean_X, group = Ecology),
    position = position_dodge(width = 0),
    width = 1.2,
    fatten = 0,
    color = "black"
  ) +
  
  geom_point(
    aes(shape = Ecology),
    position = position_jitterdodge(jitter.width = 1, dodge.width = 0),
    size = 2,
    alpha = 0.4,
    color = "black",
    show.legend = TRUE
  ) +
  
  scale_fill_manual(values = c("Monoculture" = "cyan", "Mutualism" = "yellow")) +
  labs(x = NULL, y = expression(italic(X))) +
  theme_bw() +
  theme(    
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

