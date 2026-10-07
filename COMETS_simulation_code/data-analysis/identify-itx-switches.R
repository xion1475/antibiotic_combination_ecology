# ATB 
# look at alpha value data
# single drug data

# load libraries
library("tidyverse")
library("readxl")

# base directory
base = here::here("model-data", "alpha-values")

# source functions
source(here::here("functions", "baranyi-helper-functions.r"))

# functions
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


# load no drug data
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/no-drug"

no_drug_growth_rates <- data.frame()
no_drug_ratios <- data.frame()
for (data in c("WT-biomasses-no-gal.csv", "WT-biomasses.csv")){
  run_type <- gsub("^WT-|\\.csv$", "", data)
  temp <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>% 
    group_by(strain, interaction, gene) %>% summarize(
      growth_rate = estimate_log_linear_growth(cycle, biomass),
      .groups = "drop") %>% mutate(type = run_type)
  no_drug_growth_rates <- rbind(no_drug_growth_rates, temp)
  
  temp <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>%
    filter(interaction == "coculture") %>% 
    group_by(gene, strain, interaction) %>% filter(biomass == max(biomass)) %>% slice_head(n = 1) %>% 
    ungroup() %>% select(biomass, gene, strain, interaction) %>% 
    pivot_wider(names_from = strain, values_from = biomass) %>% 
    mutate(percent_s = S / (E + S)) %>% mutate(type = run_type)
    no_drug_ratios <- rbind(no_drug_ratios, temp)
}


# load single drug data
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/single-drug/model-testing"

single_drug_growth_rates <- list()
single_drug_biomasses <- list()
single_drug_ratios <- list()
for (data in dir(base)){
  i = which(data == dir(base))
  run_type <- gsub("^biomasses_|\\.csv$", "", data)
  single_drug_growth_rates[[i]] <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>% 
    group_by(strain, interaction, gene) %>% summarize(
      growth_rate = estimate_log_linear_growth(cycle, biomass),
      .groups = "drop")
  names(single_drug_growth_rates)[[i]] <- run_type
  single_drug_biomasses[[i]] <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain))
  names(single_drug_biomasses)[[i]] <- run_type
  single_drug_ratios[[i]] <- read_csv(here::here(base, data)) %>%
    mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                              strain == "iML1515" ~ "E",
                              TRUE ~ strain)) %>%
    filter(interaction == "coculture") %>% group_by(gene, strain, interaction) %>% filter(biomass == max(biomass)) %>% slice_head(n = 1) %>% 
    ungroup() %>% select(biomass, gene, strain, interaction) %>% 
    pivot_wider(names_from = strain, values_from = biomass) %>% mutate(percent_s = S / (E + S))
  names(single_drug_ratios)[[i]] <- run_type
}

merged_df_single <- reduce(
  names(single_drug_growth_rates),
  function(df1, df2_name) {
    df2 <- single_drug_growth_rates[[df2_name]]
    
    df1 %>%
      left_join(
        df2 %>%
          select(strain, interaction, gene, growth_rate) %>%
          rename(!!df2_name := growth_rate),
        by = c("strain", "interaction", "gene")
      )
  },
  .init = single_drug_growth_rates[[1]] %>%
    select(strain, interaction, gene, growth_rate) %>%
    rename(!!names(single_drug_growth_rates)[1] := growth_rate)
) %>% select(-c("E_mono_no_gal.x")) %>% rename(E_mono_no_gal = E_mono_no_gal.y)

# double drug
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/double-drug/model-testing"

double_drug_growth_rates <- list()
double_drug_ratios <- list()
for (data in dir(base)){
  run_type <- sub("^double_drug_(growth_rates|max_yield)_(.*)\\.csv$", "\\2", data)
  df <- read_csv(here::here(base, data)) %>% select(-c(`...1`))
  if (grepl("growth_rates", data)) {
    double_drug_growth_rates[[length(double_drug_growth_rates) + 1]] <- df
    names(double_drug_growth_rates)[length(double_drug_growth_rates)] <- run_type
  }
  if (grepl("max_yield", data)) {
    double_drug_ratios[[length(double_drug_ratios) + 1]] <- df
    names(double_drug_ratios)[length(double_drug_ratios)] <- run_type
  }
}

merged_df_double <- reduce(
  names(double_drug_growth_rates),
  function(df1, df2_name) {
    df2 <- double_drug_growth_rates[[df2_name]]
    
    df1 %>%
      left_join(
        df2 %>%
          select(strain, interaction, gene, growth_rate) %>%
          rename(!!df2_name := growth_rate),
        by = c("strain", "interaction", "gene")
      )
  },
  .init = double_drug_growth_rates[[1]] %>%
    select(strain, interaction, gene, growth_rate) %>%
    rename(!!names(double_drug_growth_rates)[1] := growth_rate)
) %>% select(-c("E_mono_no_gal.x")) %>% rename(E_mono_no_gal = E_mono_no_gal.y)

# compare drug interactions
known_genes <- merged_df_single %>% pull(gene) %>% unique()

cols_to_check <- c("E_mono_no_gal", "E_mono", "S_mono_no_gal", "S_mono", "unique_mono_no_gal", "unique_mono")

interactions <- data.frame()
for (col in cols_to_check){
  subset_single <- merged_df_single %>% select(c(strain, interaction, gene, {{col}})) %>%
    rename(model_fit = {{col}})
  
  if (grepl("no_gal", col)){
    baseline <- no_drug_growth_rates %>% filter(type == "biomasses-no-gal") %>%
      select(strain, interaction, growth_rate) %>% rename(norm = growth_rate)
  } else{
    baseline <- no_drug_growth_rates %>% filter(type == "biomasses") %>%
      select(strain, interaction, growth_rate) %>% rename(norm = growth_rate)
  }
  
  classes <- merged_df_double %>% select(c(strain, interaction, gene, {{col}})) %>%
    rowwise() %>%
    mutate(split = list(split_genes(gene, known_genes)),
           gene1 = split[[1]],
           gene2 = split[[2]]) %>% ungroup() %>%
    select(-split) %>% rename(gxy = {{col}}) %>%
    inner_join(., subset_single %>% rename(gene1 = gene, gx = model_fit), by = c("strain", "interaction", "gene1")) %>%
    inner_join(., subset_single %>% rename(gene2 = gene, gy = model_fit), by = c("strain", "interaction", "gene2")) %>%
    inner_join(., baseline, by = c("strain", "interaction")) %>%
    mutate(gx_norm = gx / norm,
           gy_norm = gy / norm,
           gxy_norm = gxy / norm) %>%
    mutate(s = ((gx * gy) - gxy) / (gx * gy)) %>% 
    mutate(s_norm = ((gx_norm * gy_norm) - gxy_norm) / (gx_norm * gy_norm)) %>% 
    mutate(s = ifelse((gx == 0 | gy == 0) & gxy == 0, 0, s)) %>%
    mutate(s_norm = ifelse(gx_norm == 0 | gy_norm == 0, 0, s_norm)) %>%
    mutate(class = case_when(s > 0.05 ~ "synergistic",
                             s < -0.05 ~ "antagonistic",
                             TRUE ~ "additive")) %>% 
    mutate(class_norm = case_when(s_norm > 0.05 ~ "synergistic",
                             s_norm < -0.05 ~ "antagonistic",
                             TRUE ~ "additive")) %>% mutate(alpha_type = col)
  
  interactions <- rbind(interactions, classes)
}

baseline <- no_drug_growth_rates %>%
  mutate(galactose = ifelse(grepl("no-gal", type), "- gal", "+ gal"))

merged_df_double %>% pivot_longer(cols = -c(strain, gene, interaction)) %>%
  rename(alpha_type = name) %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  mutate(alpha_type = ifelse(grepl("tw", alpha_type), "unique_mono_tw", alpha_type)) %>%
  inner_join(., baseline, by = c("strain", "interaction")) %>% filter(value == 0) %>%
  group_by(strain, interaction, alpha_type, galactose) %>%
  summarize(n = n()) %>%
  mutate(interaction = factor(interaction, levels = c("monoculture", "coculture"))) %>%
  ggplot(aes(x = galactose, y = n / 435, fill = interaction)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_grid(strain ~ alpha_type) +
  geom_vline(xintercept = 0.0, color = "red", linetype = "dashed") + 
  theme_bw(base_size = 16) + theme(axis.title.x = element_blank()) + ylab("percent of total pairs (n = 435)")

merged_df_single %>% filter(interaction == "coculture") %>% pivot_longer(cols = -c(strain, gene, interaction)) %>%
  rename(alpha_type = name) %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  mutate(alpha_type = ifelse(grepl("tw", alpha_type), "unique_mono_tw", alpha_type)) %>%
  inner_join(., baseline %>% select(-c(type, gene)) %>% rename(norm = growth_rate), by = c("strain", "interaction", "galactose")) %>% 
  ggplot(aes(x = log(value / norm))) + geom_histogram(aes(fill = galactose)) + facet_grid(strain ~ alpha_type) +
  geom_vline(xintercept = 0.0, color = "red", linetype = "dashed") + 
  theme_bw(base_size = 16) + xlab("log(normalized single drug growth rate)") + ylab("number of drugs (n = 30)")

interactions %>% filter(interaction == "coculture") %>% 
  mutate(class_norm = ifelse((gx_norm == 0 | gy_norm == 0) & gxy_norm != 0, "antagonistic", class_norm)) %>%
  mutate(s_norm = ifelse((gx_norm == 0 | gy_norm == 0) & gxy_norm != 0, -0.06, s_norm)) %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(s_norm = ifelse(s_norm < -5, -5, s_norm)) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  ggplot(aes(x = s_norm)) + geom_histogram(aes(fill = galactose)) + facet_grid(strain ~ alpha_type) +
  theme_bw(base_size = 16) + xlab("interaction class value") + ylab("number of gene pairs (n = 435)")

interactions %>% filter(interaction == "coculture") %>% 
  filter(!is.na(s_norm)) %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  ggplot(aes(x = galactose, y = fct_reorder(gene, s_norm), fill = class_norm)) + geom_tile() + facet_grid(alpha_type ~ strain) +
  theme_bw(base_size = 16) + theme(axis.text.y = element_blank(), axis.title.x = element_blank()) + 
  ylab("gene pair") 

interactions %>% 
  mutate(class_norm = ifelse((gx_norm == 0 | gy_norm == 0) & gxy_norm != 0, "antagonistic", class_norm)) %>%
  mutate(s_norm = ifelse((gx_norm == 0 | gy_norm == 0) & gxy_norm != 0, -0.06, s_norm)) %>%
  select(strain, interaction, gene, class_norm, alpha_type) %>%
  group_by(strain, gene, alpha_type) %>% pivot_wider(names_from = interaction, values_from = class_norm) %>%
  ungroup() %>% mutate(switches = ifelse(coculture == monoculture, "no switch", "switches")) %>%
  group_by(strain, alpha_type, switches) %>% summarize(n = n()) %>%
  filter(switches == "switches") %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  ggplot(aes(x = galactose, y = n / 435)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_grid(alpha_type ~ strain) +
  theme_bw(base_size = 16) + theme(axis.title.x = element_blank()) + ylab("percent of gene pairs (n = 435)")


interactions %>% 
  mutate(class_norm = ifelse((gx_norm == 0 | gy_norm == 0) & gxy_norm != 0, "antagonistic", class_norm)) %>%
  mutate(s_norm = ifelse((gx_norm == 0 | gy_norm == 0) & gxy_norm != 0, -0.06, s_norm)) %>%
  select(strain, interaction, gene, class_norm, alpha_type) %>%
  group_by(strain, alpha_type, interaction, class_norm) %>% summarize(n = n()) %>%
  mutate(total = 435) %>% mutate(percent = n / total) %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  filter(galactose == "- gal") %>%
  mutate(interaction = factor(interaction, levels = c("monoculture", "coculture"))) %>%
  ggplot(aes(x = class_norm, y = percent, fill = interaction)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_grid(alpha_type ~ strain) +
  theme_bw(base_size = 16) + theme(axis.title.x = element_blank()) + ylab("percent of gene pairs (n = 435)")
  

interactions %>% inner_join(., double_drug_ratios$unique_mono %>% select(-c(E, S, E_stationary, S_stationary)), by = c("gene", "interaction")) %>% 
  filter(alpha_type == "unique_mono") %>% filter(gx != 0 & gy != 0) %>% 
  select(strain, gene, gene1, gene2, gx_norm, gy_norm, gxy_norm, norm, s_norm, class_norm, percent_s,interaction) %>%
  mutate(baseline_percent = 0.258) %>%
  mutate(relative = (percent_s / baseline_percent)) %>%
  mutate(greater = case_when(relative > 1.05 ~ "more S", 
                             relative < 0.95 ~ "less S",
                             TRUE ~ "no change")) %>%
  group_by(greater, strain, class_norm) %>% 
  summarize(n = n()) %>%
  ggplot(aes(x = greater, y = n / 378, fill = class_norm)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_wrap(~strain)
  










interactions %>% select(strain, interaction, gene, class_norm, alpha_type) %>%
  group_by(strain, gene, alpha_type) %>% pivot_wider(names_from = interaction, values_from = class_norm) %>%
  ungroup() %>% mutate(switches = ifelse(coculture == monoculture, "no switch", "switches")) %>% 
  filter(switches == "switches") %>% 
  mutate(type = case_when(monoculture == "synergistic" & coculture == "antagonistic" ~ "syn to antag",
                          monoculture == "synergistic" & coculture == "additive" ~ "syn to add", 
                          monoculture == "additive" & coculture == "antagonistic" ~ "add to antag",
                          monoculture == "additive" & coculture == "synergistic" ~ "add to syn",
                          monoculture == "antagonistic" & coculture == "synergistic" ~ "antag to syn",
                          monoculture == "antagonistic" & coculture == "additive" ~ "antag to add",
                          TRUE ~ "unknown")) %>%
  group_by(strain, alpha_type, type) %>% summarize(n = n()) %>%
  mutate(galactose = ifelse(grepl("no_gal", alpha_type), "- gal", "+ gal")) %>%
  mutate(alpha_type = ifelse(grepl("no_gal", alpha_type), gsub("_no_gal", "", alpha_type), alpha_type)) %>%
  ggplot(aes(x = galactose, y = n, fill = strain)) + 
  geom_bar(stat = "identity") + facet_grid(alpha_type ~ type) +
  theme_bw(base_size = 16) + theme(axis.title.x = element_blank()) + ylab("number of gene pairs")



merged_df_double %>%
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>%
  select(-split) %>% rename(gxy = model_fit) %>%
  ungroup() %>% 











# does galactose matter for the outcome in single drug conditions?
df_to_test <- single_drug_growth_rates
alpha_type <- "E_mono"
base_names <- unique(sub("_no_gal$", "", names(df_to_test)))

find_differences <- function(df1, df2, shared_cols = c("strain", "gene", "interaction"), value_col = "growth_rate") {
  
  # Rename value columns to distinguish
  df1 <- df1 |> dplyr::rename(biomass_no_gal = {{ value_col }})
  df2 <- df2 |> dplyr::rename(biomass_with_gal = {{ value_col }})
  
  # Full join to combine by shared columns
  joined <- dplyr::full_join(df1, df2, by = shared_cols)
  
  # Filter rows where biomass differs (allowing for small numeric differences)
  differences <- joined |>
    dplyr::filter(abs(biomass_no_gal - biomass_with_gal) > 1e-8 | is.na(biomass_no_gal) | is.na(biomass_with_gal))
  
  return(differences)
}

df1 <- df_to_test[[paste0("biomasses_", alpha_type, "_no_gal")]]
df2 <- df_to_test[[paste0("biomasses_", alpha_type)]]

differences <- find_differences(df1, df2)

differences %>% pivot_longer(cols = c(biomass_no_gal:biomass_with_gal)) %>% 
  mutate(name = ifelse(name == "biomass_no_gal", "no_gal", "gal")) %>%
  ggplot(aes(x = fct_reorder(gene, value), y = value, color = name)) + geom_point(size = 2) + facet_grid(interaction ~ strain) +
  theme_bw(base_size = 16) + coord_flip() + ylab("log-linear growth rate") + xlab("gene") + labs(color = "") + theme(legend.position = "bottom")


merged_df %>%
  ggplot(aes(x = biomasses_E_mono, y = biomasses_S_mono, color = strain)) + geom_point() + facet_wrap(~ interaction)






# LOAD IN DOUBLE DRUG DATA
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/double-drug"

double_growth <- read_csv(here::here(base, "double_drug_growth_rates.csv"))
double_ratios <- read_csv(here::here(base, "double_drug_max_yields.csv"))

# plots
all_combos <- double_growth %>% mutate(model_fit = ifelse(is.na(model_fit), 0, model_fit)) %>%
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>%
  select(-split) %>% select(-c(fit_variable)) %>% rename(gxy = model_fit) %>%
  ungroup() %>% inner_join(., single_growth %>% rename(gene1 = gene, gx = model_fit) %>% select(-c(fit_variable)), by = c("strain", "interaction", "gene1")) %>%
  inner_join(., single_growth %>% rename(gene2 = gene, gy = model_fit) %>% select(-c(fit_variable)), by = c("strain", "interaction", "gene2")) %>%
  mutate(s = ((gx * gy) - gxy) / (gx * gy)) %>% 
  mutate(class = case_when(s > 0.05 ~ "synergistic",
                           s < -0.05 ~ "antagonistic",
                           TRUE ~ "additive")) %>%
  mutate(interaction = ifelse(interaction == "monoculture", "mono", "co")) %>%
  mutate(interaction = factor(interaction, levels = c("mono", "co"))) %>%
  arrange(strain, interaction, class) %>%
  mutate(class = factor(class, levels = unique(class))) %>%
  ggplot(aes(x = interaction, y = gene, fill = class)) +
  geom_tile() + facet_wrap(~strain) +
  theme_bw(base_size = 16) + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                   legend.position = "bottom") +
  labs(fill = "")

bar_charts <- double_growth %>% mutate(model_fit = ifelse(is.na(model_fit), 0, model_fit)) %>%
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>%
  select(-split) %>% select(-c(fit_variable)) %>% rename(gxy = model_fit) %>%
  ungroup() %>% inner_join(., single_growth %>% rename(gene1 = gene, gx = model_fit) %>% select(-c(fit_variable)), by = c("strain", "interaction", "gene1")) %>%
  inner_join(., single_growth %>% rename(gene2 = gene, gy = model_fit) %>% select(-c(fit_variable)), by = c("strain", "interaction", "gene2")) %>%
  mutate(s = ((gx * gy) - gxy) / (gx * gy)) %>% 
  mutate(class = case_when(s > 0.05 ~ "synergistic",
                           s < -0.05 ~ "antagonistic",
                           TRUE ~ "additive")) %>%
  mutate(interaction = ifelse(interaction == "monoculture", "mono", "co")) %>%
  mutate(interaction = factor(interaction, levels = c("mono", "co"))) %>%
  group_by(interaction, strain, class) %>% summarize(n = n()) %>%
  mutate(class = factor(class, levels = unique(class))) %>%
  ggplot(aes(x = interaction, y = n, fill = class)) +
  geom_bar(stat = "identity") + facet_wrap(~strain) +
  theme_bw(base_size = 16) + theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                   legend.position = "bottom") +
  labs(fill = "")

different_per_strain <- double_growth %>% mutate(model_fit = ifelse(is.na(model_fit), 0, model_fit)) %>%
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>%
  select(-split) %>% select(-c(fit_variable)) %>% rename(gxy = model_fit) %>%
  ungroup() %>% inner_join(., single_growth %>% rename(gene1 = gene, gx = model_fit) %>% select(-c(fit_variable)), by = c("strain", "interaction", "gene1")) %>%
  inner_join(., single_growth %>% rename(gene2 = gene, gy = model_fit) %>% select(-c(fit_variable)), by = c("strain", "interaction", "gene2")) %>%
  mutate(s = ((gx * gy) - gxy) / (gx * gy)) %>% 
  mutate(class = case_when(s > 0.05 ~ "synergistic",
                           s < -0.05 ~ "antagonistic",
                           TRUE ~ "additive")) %>%
  mutate(interaction = ifelse(interaction == "monoculture", "mono", "co")) %>%
  mutate(interaction = factor(interaction, levels = c("mono", "co"))) %>% select(strain, interaction, class, gene) %>% 
  pivot_wider(names_from = interaction, values_from = class) %>%
  filter(co != mono) %>% mutate(direction = case_when(mono == "synergistic" & co == "antagonistic" ~ "synergistic to antagonistic",
                                                      co == "synergistic" & mono == "antagonistic" ~ "antagonistic to synergistic",
                                                      mono == "additive" & co == "antagonistic" ~ "additive to antagonistic",
                                                      co == "synergistic" & mono == "additive" ~ "additive to synergistic",
                                                      co == "additive" & mono == "antagonistic" ~ "antagonistic to additive",
                                                      mono == "synergistic" & co == "additive" ~ "synergistic to additive")) %>%
  group_by(strain, direction) %>% summarize(n = n()) %>%
  ggplot(aes(x = strain, y = n)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) + facet_wrap(~direction) +
  theme_bw(base_size = 16) + theme(axis.title.x = element_blank()) + ylab("number of combinations")

# ratio changes
base = "/Users/abisesi/Desktop/PhD/Projects/Antibiotics-Ecology/model-data/alpha-values/no-drug"
wt <- read_csv(here::here(base, "WT-biomasses.csv")) %>%
  mutate(strain = case_when(strain == "STM_v1_0" ~ "S",
                            strain == "iML1515" ~ "E",
                            TRUE ~ strain))

ratio_changes_single_gene <- biomasses %>% filter(interaction == "coculture") %>% 
  group_by(gene, strain, interaction) %>%
  filter(biomass == max(biomass)) %>% slice_head(n = 1) %>% ungroup() %>%
  select(biomass, gene, strain, interaction) %>% pivot_wider(names_from = strain, values_from = biomass) %>%
  mutate(per_s = S / (E + S)) %>%
  mutate(wt = wt %>% rename(baseline = biomass) %>% group_by(strain, interaction) %>% 
           filter(baseline == max(baseline)) %>% slice_head(n = 1) %>% ungroup() %>% select(strain, interaction, baseline) %>% 
           pivot_wider(names_from = strain, values_from = baseline) %>% 
           mutate(per_s_wt = S / (E + S)) %>% filter(interaction == "coculture") %>% pull(per_s_wt)) %>%
  ggplot(aes(x = fct_reorder(gene, per_s), y = log10(per_s))) + 
  geom_point(size = 4) + geom_hline(aes(yintercept = log10(wt)), linetype = "dashed") +
  theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank()) + ylab("log10(percent S)")

ratio_changes_double_gene <- double_ratios %>% 
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>%
  select(-split) %>% 
  ungroup() %>% mutate(wt = wt %>% rename(baseline = biomass) %>% group_by(strain, interaction) %>% 
                         filter(baseline == max(baseline)) %>% slice_head(n = 1) %>% ungroup() %>% select(strain, interaction, baseline) %>% 
                         pivot_wider(names_from = strain, values_from = baseline) %>% 
                         mutate(per_s_wt = S / (E + S)) %>% filter(interaction == "coculture") %>% pull(per_s_wt)) %>%
  ggplot(aes(x = fct_reorder(gene, percent_s), y = log10(percent_s))) + geom_point(size = 4) + 
  theme_bw(base_size = 16) + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()) +
  ylab("log10(percent S)") + geom_hline(aes(yintercept = log10(wt)), linetype = "dashed") + xlab("drug combo (n = 435)")

percents <- double_ratios %>% 
  rowwise() %>%
  mutate(split = list(split_genes(gene, known_genes)),
         gene1 = split[[1]],
         gene2 = split[[2]]) %>%
  select(-split) %>% 
  ungroup() %>% mutate(wt = wt %>% rename(baseline = biomass) %>% group_by(strain, interaction) %>% 
                         filter(baseline == max(baseline)) %>% slice_head(n = 1) %>% ungroup() %>% select(strain, interaction, baseline) %>% 
                         pivot_wider(names_from = strain, values_from = baseline) %>% 
                         mutate(per_s_wt = S / (E + S)) %>% filter(interaction == "coculture") %>% pull(per_s_wt)) %>%
  mutate(direction = case_when(round(percent_s,2) == round(wt,2) ~ "no change",
                               round(percent_s,2) > round(wt,2) ~ "more S",
                               round(percent_s,2) < round(wt,2) ~ "less S")) %>%
  group_by(direction) %>% summarize(n = n()) %>%
  ggplot(aes(x = direction, y = n)) + geom_point(size = 4) + 
  theme_bw(base_size = 16) + 
  ylab("number of combinations") + theme(axis.title.x = element_blank())


# check that curve data for single drugs makes sense
mono_baseline <- biomasses %>% filter(interaction == "monoculture") %>% 
  inner_join(., wt %>% rename(baseline = biomass) %>% group_by(strain, interaction) %>% 
               filter(baseline == max(baseline)) %>% slice_head(n = 1) %>% ungroup() %>% select(strain, interaction, baseline), by = c("strain", "interaction")) %>% 
  select(cycle, biomass, gene, interaction, baseline, strain) %>% ggplot(aes(x = cycle, y = biomass, color = strain)) + 
  geom_line() + facet_wrap(~gene) + geom_hline(aes(yintercept = baseline, color = strain), linetype = "dashed") +
  theme_bw(base_size = 16)

co_baseline <- biomasses %>% filter(interaction == "coculture") %>% 
  inner_join(., wt %>% rename(baseline = biomass) %>% group_by(strain, interaction) %>% 
               filter(baseline == max(baseline)) %>% slice_head(n = 1) %>% ungroup() %>% select(strain, interaction, baseline), by = c("strain", "interaction")) %>% 
  select(cycle, biomass, gene, interaction, baseline, strain) %>% ggplot(aes(x = cycle, y = biomass, color = strain)) + 
  geom_line() + facet_wrap(~gene) + geom_hline(aes(yintercept = baseline, color = strain), linetype = "dashed") +
  theme_bw(base_size = 16)

