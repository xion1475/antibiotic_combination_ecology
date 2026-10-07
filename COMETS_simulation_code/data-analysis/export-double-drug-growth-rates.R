# ATB 
# export growth rates from double drug data

# load libraries
library("tidyverse")
library("readxl")

# type of double drug data
alpha_type = "unique_mono_no_gal" # unique_mono, E_mono, S_mono

# base directory
base = here::here("model-data", "alpha-values")

# source functions
source(here::here("functions", "baranyi-helper-functions.r"))

# growth rates for double drugs
double_drug_base = here::here(base, "double-drug", "biomasses")

double_growth <- data.frame()
double_ratios <- data.frame()
for (name in dir(double_drug_base)){
  
  geneid <- sub("\\.csv$", "", name)
  temp <- read_csv(paste0(double_drug_base, "/", name)) %>% 
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain))
  
 stationary <- temp %>% group_by(strain, gene, interaction) %>% 
    filter(cycle == max(cycle) | cycle == max(cycle) - 1) %>% select(-c(`...1`)) %>% 
    pivot_wider(names_from = cycle, values_from = biomass) %>% mutate(stationary = ifelse(`1999` == `2000`, TRUE, FALSE))
  
  growth <- temp %>% mutate(gene = geneid) %>% 
    group_by(strain, interaction, gene) %>% summarize(
      growth_rate = estimate_log_linear_growth(cycle, biomass),
      .groups = "drop") %>%
    inner_join(., stationary %>% select(-c(`1999`, `2000`)), by = c("strain", "interaction", "gene"))
  
  ratio <- temp %>% filter(interaction == "coculture") %>% group_by(gene, strain, interaction) %>% filter(biomass == max(biomass)) %>% slice_head(n = 1) %>% 
    ungroup() %>% select(biomass, gene, strain, interaction) %>% 
    pivot_wider(names_from = strain, values_from = biomass) %>% mutate(percent_s = S / (E + S)) %>%
    inner_join(., stationary %>% filter(interaction == "coculture") %>% 
                 select(-c(`1999`, `2000`)) %>% pivot_wider(names_from = strain, values_from = stationary) %>% 
                 rename(E_stationary = E, S_stationary = S), by = c("interaction", "gene"))
  
  double_growth <- rbind(double_growth, growth)
  double_ratios <- rbind(double_ratios, ratio)
}

write.csv(double_growth, file = here::here(base, "double-drug", paste0("double_drug_growth_rates_", alpha_type, "tw-media", ".csv")))
write.csv(double_ratios, file = here::here(base, "double-drug", paste0("double_drug_max_yield_", alpha_type,"tw-media", ".csv")))

