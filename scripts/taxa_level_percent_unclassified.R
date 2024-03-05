library(tidyverse)
library(magrittr)
library(vegan)
library(glue)

source("scripts/aux_functions.R")

main_out <- "results/unclassified_plots"
dir.create(main_out, recursive = TRUE)


# Read in  ----------------------------------------------------------------
timed_meta <- get_timed_metadata()


# 1. Read in kaiju tables --------------------------------------------------

shallow_path <- "data/kaiju_shallow"
shallow_phylum_long <- read_tsv(glue("{shallow_path}/kaiju_phylum_summary.tsv"))
shallow_class_long <- read_tsv(glue("{shallow_path}/kaiju_class_summary.tsv"))
shallow_order_long <- read_tsv(glue("{shallow_path}/kaiju_order_summary.tsv"))
shallow_family_long <- read_tsv(glue("{shallow_path}/kaiju_family_summary.tsv"))
shallow_genus_long <- read_tsv(glue("{shallow_path}/kaiju_genus_summary.tsv"))
shallow_species_long <- read_tsv(glue("{shallow_path}/kaiju_species_summary.tsv"))

shallow_taxa_tables <- list(shallow_phylum_long, shallow_class_long, shallow_order_long,
                            shallow_family_long, shallow_genus_long, shallow_species_long)

deep_path <- "data/kaiju_deep/"
phylum_long <- read_tsv(glue("{deep_path}/kaiju_phylum_summary.tsv"))
class_long <- read_tsv(glue("{deep_path}/kaiju_class_summary.tsv"))
order_long <- read_tsv(glue("{deep_path}/kaiju_order_summary.tsv"))
family_long <- read_tsv(glue("{deep_path}/kaiju_family_summary.tsv"))
genus_long <- read_tsv(glue("{deep_path}/kaiju_genus_summary.tsv"))
species_long <- read_tsv(glue("{deep_path}/kaiju_species_summary.tsv"))

deep_taxa_tables <- list(phylum_long,class_long,order_long,family_long,genus_long,species_long)
all_taxa_ranks <- list("phylum","class","order","family","genus","species")

# 2. Define functions ----
filter_unclassified <- function(df) {
  df %<>%
    filter(str_detect(taxon_name, "unclassified$|^cannot|Viruses|[ ;]bacterium[; ]| sp.[ ;]"))  %>%
    mutate(taxon_name = case_when(str_detect(taxon_name, "cannot") ~ "cannot",
                                  TRUE ~ taxon_name))
  return(df)
}

mutate_cols <- function(df, taxa_rank){
  mutate(df, taxa_rank = taxa_rank,
         sample = basename(file) |> str_extract("S00JY-\\d{4}"))
}

# 3. Perform filtering and merge results ----
## Deep samples ----
deep_unclassified_list <- lapply(X = deep_taxa_tables,
                            FUN = filter_unclassified)
deep_unclassified_list <- mapply(mutate_cols, deep_unclassified_list, all_taxa_ranks,
                            SIMPLIFY = FALSE)

deep_unclassified <- bind_rows(deep_unclassified_list)

taxa_factor_levels <-c("phylum",
                       "class",
                       "order",
                       "family",
                       "genus",
                       "species")
# Order taxa level factor for plotting
deep_unclassified$taxa_rank <- ordered(deep_unclassified$taxa_rank,
                                           levels = taxa_factor_levels)

## Shallow samples ----
shallow_unclassified_list <- lapply(X = shallow_taxa_tables,
                            FUN = filter_unclassified)
shallow_unclassified_list <- mapply(mutate_cols, shallow_unclassified_list, all_taxa_ranks,
                            SIMPLIFY = FALSE)

shallow_unclassified <- bind_rows(shallow_unclassified_list)

# Order taxa level factor for plotting
shallow_unclassified$taxa_rank <- ordered(shallow_unclassified$taxa_rank,
                                           levels = taxa_factor_levels)


# 4. Create plots ---------------------------------------------------------

unclassified_components_box <- function(df, subtitle) {
  df |>
    inner_join(timed_meta) |>
    filter(taxon_name %in% c("cannot","unclassified","Viruses")) |>
    ggplot(aes(y = percent, x = Depth)) +
    geom_boxplot() +
    labs(title = "Distribution of % unclassified reads breakdown by soil depth",
         subtitle = subtitle,
         y = "Percent unclassified",
         x = "Soil Depth") +
    geom_jitter(aes(color = Depth), size=1, alpha=0.6)  +
    facet_grid(taxon_name ~ taxa_rank, scales = "free_y")
}

## Create plots ----

# combined
du2 <- deep_unclassified |> mutate(seq_depth = "deep")
su2 <- shallow_unclassified |> mutate(seq_depth = "shallow")
all_depths_unclassified <- bind_rows(du2, su2)
all_depths_unclassified_plot <- all_depths_unclassified |>
  summarize(percent_unclassified = sum(percent), .by = c(sample, taxa_rank, seq_depth)) |>
  inner_join(timed_meta) |>
  ggplot(aes(x = factor(seq_depth, levels=c("shallow","deep")), y = percent_unclassified)) +
  geom_boxplot() +
  geom_jitter(aes(color = Depth), size=1, alpha=0.4) +
  labs(title = "Percentage of unclassified reads at each taxa rank",
       x = "Sequencing depth",
       y = "Percent unclassified",
       color = "Soil depth") +
  facet_grid(.~taxa_rank)
ggsave(glue("{main_out}/unclassified_percents.png"), all_depths_unclassified_plot, width = 8)

uc_components_deep_plot <- unclassified_components_box(deep_unclassified, "All 42 deep metagenomic samples - Timed dataset")
uc_components_shallow_plot <- unclassified_components_box(shallow_unclassified, "All 124 shallow metagenomic samples - Timed dataset")

uc_components_deep_plot
shallow_box

