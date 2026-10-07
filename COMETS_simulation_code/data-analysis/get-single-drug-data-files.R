# ATB
# 15 Aug 2025
# get growth rate, biomass and flux data for single drugs based on raw COMETS data

# source functions
source(here::here("functions", "baranyi-helper-functions.r"))

base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/single-drug/biomasses"
single_drug_growth_rates_old_log_linear <- data.frame()
single_drug_growth_rates_new_log_linear <- data.frame()
single_drug_biomasses <- data.frame()
for (data in dir(base)){
  biomass <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>% select(-c(`...1`))
  
  growth_rate <- biomass %>% 
    group_by(strain, interaction, gene) %>% summarize(
      growth_rate = estimate_max_growth(cycle, biomass),
      .groups = "drop")
  
  growth_rate_new <- biomass %>% 
    group_by(strain, interaction, gene) %>% summarize(
      growth_rate = estimate_log_linear_growth(cycle, biomass),
      .groups = "drop")
  
  max_yield <- biomass %>% 
    group_by(strain, interaction, gene) %>% filter(biomass == max(biomass)) %>% slice_head(n = 1)
  
  is_stationary <- biomass %>% group_by(strain, interaction, gene) %>% 
    filter(cycle == max(cycle) | cycle == (max(cycle) - 1)) %>% pivot_wider(names_from = cycle, values_from = biomass) %>%
    mutate(is_stationary = ifelse(`749` == `750`, TRUE, FALSE)) %>% select(gene, strain, interaction, is_stationary)
  
  half_max <- biomass %>% inner_join(., max_yield %>% 
                                       rename(max_yield = biomass, time_to_max_yield = cycle), 
                                     by = c("gene", "strain", "interaction")) %>% mutate(half = biomass / max_yield) %>% 
    group_by(gene, strain, interaction) %>% slice_min(abs(half - 0.5), with_ties = FALSE) %>%
    ungroup() %>% rename(half_max = biomass, time_to_half_max = cycle) %>% 
    select(time_to_half_max, half_max, gene, strain, interaction, time_to_max_yield, max_yield)
  
  biomass_stats <- half_max %>% inner_join(., is_stationary, by = c("gene", "interaction", "strain"))
  single_drug_biomasses <- rbind(single_drug_biomasses, biomass_stats)
  single_drug_growth_rates_old_log_linear <- rbind(single_drug_growth_rates_old_log_linear, growth_rate)
  single_drug_growth_rates_new_log_linear <- rbind(single_drug_growth_rates_new_log_linear, growth_rate_new)
}

not_in_stationary <- single_drug_biomasses %>% filter(is_stationary == FALSE)

base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/single-drug/fluxes"
single_drug_fluxes <- data.frame()
for (data in dir(base)){
  fluxes <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>% select(-c(`...1`))
  
  non_zero_fluxes <- fluxes %>% group_by(reaction, gene, strain, interaction) %>% 
    summarize(mean_flux = mean(flux)) %>% filter(mean_flux != 0) %>% pull(reaction) %>% unique()
  
  single_drug_fluxes <- rbind(single_drug_fluxes, fluxes %>% filter(reaction %in% non_zero_fluxes))
}

mid_log_fluxes <- single_drug_fluxes %>% inner_join(., single_drug_biomasses, by = c("gene", "strain", "interaction")) %>% 
  filter((cycle >= time_to_half_max - 10) & (cycle <= time_to_half_max + 10)) %>%
  group_by(strain, interaction, gene, reaction) %>% mutate(mean = mean(flux)) %>% filter(mean != 0) %>%
  select(strain, interaction, gene, cycle, reaction, flux, time_to_half_max)

write.csv(mid_log_fluxes, here::here("model-data", "alpha-values", "single-drug", "single_gene_mid_log_fluxes.csv"))
write.csv(single_drug_biomasses, here::here("model-data", "alpha-values", "single-drug", "single_gene_biomass_summary.csv"))
write.csv(single_drug_growth_rates_old_log_linear, here::here("model-data", "alpha-values", "single-drug", "single_gene_growth_rates_old_log_linear.csv"))
write.csv(single_drug_growth_rates_new_log_linear, here::here("model-data", "alpha-values", "single-drug", "single_gene_growth_rates_new_log_linear.csv"))
