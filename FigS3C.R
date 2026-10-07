library(ggplot2)
library(dplyr)
library(openxlsx)
library(tidyr)
library(scales)
library(brms)
setwd("/Users/xiongxia/Desktop")
shift_by_cutoff = read.xlsx("Interaction_count_by_threshold.xlsx", sheet = 4)

# Make sure Yeh_extremes and Ecology are factors
shift_by_cutoff$Yeh_threshold    <- factor(shift_by_cutoff$Yeh_threshold)

# Aggregate counts (in case of duplicate rows)
plot_df <- shift_by_cutoff %>%
  group_by(Yeh_threshold, Interaction_shift_from_mono_to_co) %>%
  summarise(Count = sum(Count), .groups = "drop")

# Plot: stacked by Interaction_type, dodged by Ecology
ggplot(plot_df,
       aes(x = Yeh_threshold,
           y = Count,
           fill = Interaction_shift_from_mono_to_co)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("Toward synergy" = "purple",
                               "Toward antagonism" = "orange",
                               "Stay unchanged"  = "white")) +
  labs(x = "Cutoffs for extreme-based tests",
       y = "Total count",
       fill = "Interaction type") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 10))

binom.test(3, 3, p = 0.5, alternative = "greater")
