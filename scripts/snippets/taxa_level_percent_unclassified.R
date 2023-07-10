library(tidyverse)
library(magrittr)
library(vegan)
library(glue)


# Read in and tidy metadata -----------------------------------------------

timed_meta <- read_csv("../metadata/Metadata-IAG-Timed2-v3_2.csv")
# remove all NA rows, remove blanks
timed_meta <- timed_meta[rowSums(is.na(timed_meta)) != ncol(timed_meta),] |>
  rename(sample = `MBI_ID`) |>
  filter(!str_detect(Site, "B\\d"))

timed_meta <- timed_meta |>
  mutate(Month = case_when(Year == "1998" | Year == "2007" ~ Year,
                           TRUE ~ Month))
timed_meta$Month <- ordered(timed_meta$Month, levels = c("0",
                                                         "0.07",
                                                         "0.5",
                                                         "1",
                                                         "3",
                                                         "6",
                                                         "12",
                                                         "18",
                                                         "1998",
                                                         "2007"))
timed_meta$Year <- as.character(timed_meta$Year)


# Read in kaiju tables ----------------------------------------------------

root_path <- "data/deep_Timed/"
phylum_long <- read_tsv(glue("{root_path}/kaiju_phylum_summary.tsv"))
class_long <- read_tsv(glue("{root_path}/kaiju_class_summary.tsv"))
order_long <- read_tsv(glue("{root_path}/kaiju_order_summary.tsv"))
family_long <- read_tsv(glue("{root_path}/kaiju_family_summary.tsv"))
genus_long <- read_tsv(glue("{root_path}/kaiju_genus_summary.tsv"))
species_long <- read_tsv(glue("{root_path}/kaiju_species_summary.tsv"))

all_taxa_tables <- list(phylum_long,class_long,order_long,family_long,genus_long,species_long)

# Define functions ----
filter_unclassified <- function(df) {
  df %<>%
    filter(str_detect(taxon_name, "unclassified|cannot | bacterium;$|^bacterium;$")) %>%
    select(-c(taxon_id, reads)) %>%
    mutate(taxon_name = case_when(str_detect(taxon_name, "cannot") ~ "cannot",
                                  TRUE ~ taxon_name))
  return(df)
}

sum_unclassified <- function(df, rank_string) {
  df_filtered <- df %>% filter_unclassified()
  df_filtered %>% group_by(file) %>% summarise({{rank_string}} := sum(percent))
}

# Perform filtering and merge results ----
phylum_unclass <- sum_unclassified(phylum_long, "phylum")

merged_unclass <- merge(phylum_unclass, sum_unclassified(class_long, "class")) %>%
  merge(sum_unclassified(order_long, "order")) %>%
  merge(sum_unclassified(family_long, "family")) %>%
  merge(sum_unclassified(genus_long, "genus")) %>%
  merge(sum_unclassified(species_long, "species"))

# Species has a lower total unclassified than genus?

merged_unclass_long <- merged_unclass %>%
  pivot_longer(cols = !file, names_to = "taxa_rank",values_to = "percent_unclassified") |>
  mutate(sample = basename(file) |> str_extract("S00JY-\\d{4}"))

# Order taxa level factor for plotting ----
merged_unclass_long$taxa_rank <- ordered(merged_unclass_long$taxa_rank,
                                           levels = c("phylum",
                                                      "class",
                                                      "order",
                                                      "family",
                                                      "genus",
                                                      "species"))

# Plot violin distributions of % unclassified across taxa levels
v_plot <- merged_unclass_long |>
  left_join(timed_meta) |>
  ggplot(aes(x = taxa_rank, y = percent_unclassified)) +
  geom_violin() +
  geom_jitter(aes( color = Depth), size=1, alpha=0.4) +
  labs(title = "Distribution of % unclassified reads at each taxa rank",
       subtitle = "All 42 deep metagenomic reads - Timed dataset",
       x = "Taxa rank")

v_plot

