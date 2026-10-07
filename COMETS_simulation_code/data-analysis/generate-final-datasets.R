# ATB
# 15 Aug 2025
# get final gene perturbation data

# function for gene splitting
split_genes <- function(gene_string, known_genes) {
  matches <- sapply(known_genes, function(g1) {
    g2 <- sub(paste0("^", g1), "", gene_string)
    if (g1 != gene_string && g2 %in% known_genes) {
      return(c(g1, g2))
    } else {
      return(NULL)
    }
  })
  # Filter out NULLs and keep first valid match
  matches <- Filter(Negate(is.null), matches)
  if (length(matches) > 0) {
    return(matches[[1]])
  } else {
    return(c(NA, NA))
  }
}

# source functions
source(here::here("functions", "baranyi-helper-functions.r"))

# import data for no drug
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/no-drug"
data = "WT-biomasses.csv"

no_drug_growth_rates_new_log_linear <- read_csv(here::here(base, data)) %>%
  mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                            strain == "iML1515" ~ "E",
                            TRUE ~ strain)) %>% 
  group_by(strain, interaction, gene) %>% summarize(
    growth_rate = estimate_log_linear_growth(cycle, biomass),
    .groups = "drop") %>%
  mutate(growth_rate = case_when(interaction == "monoculture" & strain == "E" ~ 0.5167917,
                                 interaction == "monoculture" & strain == "S" ~ 0.4225747,
                                 TRUE ~ growth_rate))

no_drug_biomass <- read_csv(here::here(base, data)) %>%
  mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                            strain == "iML1515" ~ "E",
                            TRUE ~ strain)) %>% select(-c(`...1`))

max_yield <- no_drug_biomass %>% 
  group_by(strain, interaction, gene) %>% filter(biomass == max(biomass)) %>% slice_head(n = 1)

is_stationary <- no_drug_biomass %>% group_by(strain, interaction, gene) %>% 
  filter(cycle == max(cycle) | cycle == (max(cycle) - 1)) %>% pivot_wider(names_from = cycle, values_from = biomass) %>%
  mutate(is_stationary = ifelse(`199` == `200`, TRUE, FALSE)) %>% select(gene, strain, interaction, is_stationary)

half_max <- no_drug_biomass %>% inner_join(., max_yield %>% 
                                     rename(max_yield = biomass, time_to_max_yield = cycle), 
                                   by = c("gene", "strain", "interaction")) %>% mutate(half = biomass / max_yield) %>% 
  group_by(gene, strain, interaction) %>% slice_min(abs(half - 0.5), with_ties = FALSE) %>%
  ungroup() %>% rename(half_max = biomass, time_to_half_max = cycle) %>% 
  select(time_to_half_max, half_max, gene, strain, interaction, time_to_max_yield, max_yield)

biomass_stats <- half_max %>% inner_join(., is_stationary, by = c("gene", "interaction", "strain"))

data = "WT-fluxes.csv"

fluxes <- read_csv(here::here(base, data)) %>%
  mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                            strain == "iML1515" ~ "E",
                            TRUE ~ strain)) %>% select(-c(`...1`)) %>%
  filter(grepl("EX_", reaction) == TRUE | grepl("BIOMASS", reaction) == TRUE)

non_zero_fluxes <- fluxes %>% group_by(reaction, gene, strain, interaction) %>% 
  summarize(mean_flux = mean(flux)) %>% filter(mean_flux != 0) %>% pull(reaction) %>% unique()

trimmed_fluxes <- fluxes %>% filter(reaction %in% non_zero_fluxes) %>%
  inner_join(., biomass_stats, by = c("gene", "strain", "interaction")) %>% 
  filter((cycle >= time_to_half_max - 10) & (cycle <= time_to_half_max + 10)) %>%
  group_by(strain, interaction, gene, reaction) %>% mutate(mean = mean(flux)) %>% filter(mean != 0) %>%
  select(strain, interaction, gene, cycle, reaction, flux, time_to_half_max) %>% ungroup()

# make combined data set for biomasses (half max, time to max, etc) - single and double, plus WT!!!!
# cols: strain, interaction, gene, time_to_mid_log, max_yield, time_to_stationary, is_stationary
single_biomasses <- read_csv(here::here("model-data", "alpha-values", "single-drug", "single_gene_biomass_summary.csv"))
double_biomasses <- read_csv(here::here("model-data", "alpha-values", "double-drug", "double_gene_biomass_summary.csv"))

all_biomasses <- rbind(single_biomasses %>% select(-c(`...1`)), double_biomasses %>% select(-c(`...1`)),
      biomass_stats %>% mutate(gene = "no_drug")) %>% 
  select(strain, interaction, gene, time_to_half_max, max_yield, time_to_max_yield, is_stationary) %>%
  rename(time_to_mid_log = time_to_half_max, time_to_stationary = time_to_max_yield)

write.csv(all_biomasses, here::here("model-data", "alpha-values", "gene_perturbation_biomass_data.csv"))

# make combined data set for fluxes - single and double, plus WT!!!!
# cols: strain, interaction, gene, cycle, reaction, flux, time_to_mid_log
single_fluxes <- read_csv(here::here("model-data", "alpha-values", "single-drug", "single_gene_mid_log_fluxes.csv"))
double_fluxes <- read_csv(here::here("model-data", "alpha-values", "double-drug", "double_gene_mid_log_fluxes.csv"))

all_fluxes <- rbind(single_fluxes %>% select(-c(`...1`)), double_fluxes %>% select(-c(`...1`)),
      trimmed_fluxes %>% mutate(gene = "no_drug")) %>% rename(time_to_mid_log = time_to_half_max)

write.csv(all_fluxes, here::here("model-data", "alpha-values", "gene_perturbation_flux_data.csv"))

# make combined data set for interactions
# cols: strain, interaction, gene_combo, gene_x, gene_y, no_drug_growth_rate, gx_norm, gy_norm, gxy_norm, s_norm, class_norm

known_genes <- single_drug_growth_rates %>% pull(gene) %>% unique()

interactions <- double_drug_growth_rates %>% select(c(strain, interaction, gene, growth_rate)) %>%
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>% ungroup() %>%
  select(-split) %>% rename(gxy = growth_rate) %>%
  inner_join(., single_drug_growth_rates %>% rename(gene1 = gene, gx = growth_rate) %>% select(-c(`...1`)), by = c("strain", "interaction", "gene1")) %>%
  inner_join(., single_drug_growth_rates %>% rename(gene2 = gene, gy = growth_rate) %>% select(-c(`...1`)), by = c("strain", "interaction", "gene2")) %>%
  inner_join(., no_drug_growth_rate %>% select(-c(gene)) %>% rename(no_drug_growth_rate = growth_rate), by = c("strain", "interaction")) %>%
  mutate(gx_norm = gx / no_drug_growth_rate,
         gy_norm = gy / no_drug_growth_rate,
         gxy_norm = gxy / no_drug_growth_rate) %>%
  mutate(s_norm = ((gx_norm * gy_norm) - gxy_norm) / (gx_norm * gy_norm)) %>% 
  mutate(class_norm = case_when(s_norm > 0.05 ~ "synergistic",
                                s_norm < -0.05 ~ "antagonistic",
                                TRUE ~ "additive"))

drug_interactions <- interactions %>% rename(gene_x = gene1, gene_y = gene2, gene_xy = gene, g_no_drug = no_drug_growth_rate) %>% 
  rename(s = s_norm, class = class_norm) %>%
  select(strain, interaction, gene_x, gene_y, gene_xy, gx, gy, gxy, g_no_drug, gx_norm, gy_norm, gxy_norm, s, class)

write.csv(drug_interactions, here::here("model-data", "alpha-values", "gene_perturbation_drug_interaction_data.csv"))

