# ATB
# 15 Aug 2025
# get growth rate, biomass and flux data for double drugs based on raw COMETS data

# source functions
source(here::here("functions", "baranyi-helper-functions.r"))

# set base directory
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/double-drug/biomasses"

# loop through directory and pull out biomass and growth rate data
double_drug_growth_rates_old_log_linear <- data.frame()
double_drug_growth_rates_new_log_linear <- data.frame()
double_drug_biomasses <- data.frame()
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
    filter(cycle == max(cycle) | cycle == (max(cycle) - 1)) %>% 
    mutate(index = ifelse(cycle %in% c(10000, 5000), "final", "near_final")) %>% select(-c(cycle)) %>%
    pivot_wider(names_from = index, values_from = biomass) %>%
    mutate(is_stationary = ifelse(near_final == final, TRUE, FALSE)) %>% select(gene, strain, interaction, is_stationary)
  
  half_max <- biomass %>% inner_join(., max_yield %>% 
                                       rename(max_yield = biomass, time_to_max_yield = cycle), 
                                     by = c("gene", "strain", "interaction")) %>% mutate(half = biomass / max_yield) %>% 
    group_by(gene, strain, interaction) %>% slice_min(abs(half - 0.5), with_ties = FALSE) %>%
    ungroup() %>% rename(half_max = biomass, time_to_half_max = cycle) %>% 
    select(time_to_half_max, half_max, gene, strain, interaction, time_to_max_yield, max_yield)
  
  biomass_stats <- half_max %>% inner_join(., is_stationary, by = c("gene", "interaction", "strain"))
  double_drug_biomasses <- rbind(double_drug_biomasses, biomass_stats)
  double_drug_growth_rates_old_log_linear <- rbind(double_drug_growth_rates_old_log_linear, growth_rate)
  double_drug_growth_rates_new_log_linear <- rbind(double_drug_growth_rates_new_log_linear, growth_rate_new)
}

# check that all simulations are in stationary - redo simulations that are not
not_in_stationary <- double_drug_biomasses %>% filter(is_stationary == FALSE)

# update NA values
no_growth <- double_drug_growth_rates_new_log_linear %>% filter(is.na(growth_rate)) %>% pull(gene) %>% unique()

no_luck <- c()
less_stringent <- data.frame()
for (gene in no_growth){
  data = paste0(gene, ".csv")
  biomass <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>% select(-c(`...1`))
  temp <- biomass %>% 
    group_by(strain, interaction, gene) %>% summarize(
      growth_rate = estimate_log_linear_growth(cycle, biomass, r2_thresh = 0.8),
      .groups = "drop")
  if (any(is.na(temp %>% pull(growth_rate)))) {
    no_luck <- c(no_luck, gene)
  }
  less_stringent <- rbind(less_stringent, temp)
}

less_stringent <- less_stringent %>% mutate(growth_rate = ifelse(is.na(growth_rate) & gene %in% no_luck, 0, growth_rate))
double_drug_growth_rates_new_log_linear <- full_join(double_drug_growth_rates_new_log_linear, less_stringent %>% rename(new = growth_rate), by = c("strain", "interaction", "gene")) %>%
  mutate(growth_rate = ifelse(is.na(growth_rate), new, growth_rate))

# set base directory for flux data
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/double-drug/fluxes"
double_drug_fluxes <- data.frame()
for (data in dir(base)){
  fluxes <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>% select(-c(`...1`))
  
  non_zero_fluxes <- fluxes %>% group_by(reaction, gene, strain, interaction) %>% 
    summarize(mean_flux = mean(flux)) %>% filter(mean_flux != 0) %>% pull(reaction) %>% unique()
  
  trimmed_fluxes <- fluxes %>% filter(reaction %in% non_zero_fluxes) %>%
    inner_join(., double_drug_biomasses, by = c("gene", "strain", "interaction")) %>% 
    filter((cycle >= time_to_half_max - 10) & (cycle <= time_to_half_max + 10)) %>%
    group_by(strain, interaction, gene, reaction) %>% mutate(mean = mean(flux)) %>% filter(mean != 0) %>%
    select(strain, interaction, gene, cycle, reaction, flux, time_to_half_max) %>% ungroup()
  
  double_drug_fluxes <- rbind(double_drug_fluxes, trimmed_fluxes)
  
  rm(fluxes, non_zero_fluxes, trimmed_fluxes)
}

# save all files
write.csv(double_drug_fluxes, here::here("model-data", "alpha-values", "double-drug", "double_gene_mid_log_fluxes.csv"))
write.csv(double_drug_biomasses, here::here("model-data", "alpha-values", "double-drug", "double_gene_biomass_summary.csv"))
write.csv(double_drug_growth_rates_new_log_linear %>% select(-c(new)), here::here("model-data", "alpha-values", "double-drug", "double_gene_growth_rates_new_log_linear.csv"))
write.csv(double_drug_growth_rates_old_log_linear, here::here("model-data", "alpha-values", "double-drug", "double_gene_growth_rates_old_log_linear.csv"))
