library(ggplot2)
library(tidyverse)
library(vegan)
library(glue)

source("scripts/aux_functions.R")

main_out <- "results_feb16/ancombc"
plot_out <- glue("{main_out}/DA_heatmaps")

timed_meta <- get_timed_metadata()



# Heatmaps - Taxa --------------------------------------------

# Retrieve ANCOMBC2 results for the Month factor for a given taxa level. Reads
# in all of the shallow subset results, as well as the deep result table.
get_combined_month_DA_results <- function(taxa_level) {
  files_month_ancom <- list.files(path = glue("{main_out}/shallow_sequencing/{taxa_level}/"),
                                  pattern = ".*month.*.csv", full.names = TRUE)
  combined_results <- map(files_month_ancom, ~read_csv(.x)
                          %>% mutate(filename = basename(.x))) |>
    bind_rows()

  deep_phylum <- read_csv(glue("{main_out}/deep_sequencing/{taxa_level}/da-table_month-continuous.csv")) |>
    mutate(filename = "da-table_month-continuous_deep.csv") # rename in table to ID as deeply sequenced

  cr2 <- combined_results |> bind_rows(deep_phylum) |>
    mutate(subset = str_extract(filename, ".*month-continuous(.*)\\.csv", group = 1) |>
             str_remove("^_"),
           subset = ifelse(subset == "", "shallow", subset))
}

# formats the Process factor pairwise comparisons: creates shorthands and makes
# it a factor for proper plotting
format_process_comparisons <- function(results) {
  comparison_factor_order <- c("Sieved_Fresh" = "S/F",
                               "Air-dried_Fresh" = "A/F",
                               "Ground_Fresh" = "G/F",
                               "Air-dried_Sieved" = "A/S",
                               "Ground_Sieved" = "G/S",
                               "Ground_Air-dried" = "G/A")
  results2 <- results |>
    mutate(comparison_code = comparison_factor_order[comparison])
  results2$comparison_code <- fct(results2$comparison_code,
                                  levels=comparison_factor_order)
  results2
}

# Calculates average taxa count from shallow count table and uses these values
# to categorize taxa in common, mid, or rare, with cutoffs depending on taxa
# level.
format_taxa_results <- function(results, counts, taxa_level,
                                metafactor = c("process","month")) {
  metafactor <- rlang::arg_match(metafactor)
  rare_cutoff <- if (taxa_level == "phylum") 0.001 else 0.001
  common_cutoff <- if(taxa_level == "phylum") 0.01 else 0.005

  avg_percents <- counts$filtered |>
    summarize(avg_percent = mean(classified_percent), .by = c(taxon_id, taxon_name))

  cr2 <- results |>
    filter(diff & passed.ss) |>
    mutate(occurrence_tally = n(), .by = taxon) |>
    left_join(avg_percents, by = join_by(taxon == taxon_name)) |>
    mutate(rarity = case_when(avg_percent < rare_cutoff ~ "rare",
                              avg_percent > common_cutoff ~ "common",
                              TRUE ~ "mid"))

  cr2$rarity <- factor(cr2$rarity,
                       levels=c("common","mid","rare"))

  if (metafactor == "month") {
    cr2 <- cr2 |>
      mutate(facet_group = case_when(str_detect(subset, "0-15") ~ "0-15",
                                     str_detect(subset, "15-30") ~ "15-30",
                                     TRUE ~ subset))
    cr2$facet_group <- factor(cr2$facet_group,
                              levels=c("shallow","0-15","15-30", "deep"))
  }

  cr2
}

# Given a numeric table column, provides a ggplot color scale with 2 gradients,
# sharply distinct at 0
rescale_2gradient <- function(fc) {
  scale_fill_gradientn(
    colors = c("#8F001F", "rosybrown1","aquamarine", "#1F8F00"),
    values = scales::rescale(c(min(fc), 0, 0.00000001, max(fc)))
  )
}

## Phylum ----
phyla_counts <- get_processed_taxonomy("shallow", "phylum", timed_meta)

### Month factor ----
phyla_month <- get_combined_month_DA_results("phylum") |>
  format_taxa_results(phyla_counts, "phylum", "month")

phy_month_plot_tab <- phyla_month |>
  filter(occurrence_tally > 1) # keep only phyla in >1 subset

phy_plot_month_gradient <- phy_month_plot_tab |>
  ggplot(aes(x = subset, y = fct_reorder(taxon, avg_percent))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(rarity~facet_group, scales = "free", space = "free") +
  labs(title="Differentially abundant phyla for factor Month",
       subtitle="Mean abundance: common (>1%), rare (<0.1%); significant results only",
       y="Phylum") +
  geom_tile(aes(fill = log2fc*18), color="black") +
  rescale_2gradient(phy_month_plot_tab$log2fc) +
  labs(fill=str_wrap("Log2 fold-change (18 months)",16))
phy_plot_month_gradient
ggsave(glue("{plot_out}/month_phyla_heatmap.png"), phy_plot_month_gradient, width = 8, height = 6)

### Process factor ----

phyla_process <- read_csv(glue("{main_out}/shallow_sequencing/phylum/da-table_Process.csv")) |>
  format_taxa_results(phyla_counts, "phylum", "process") |>
  format_process_comparisons()

phy_plot_process <- phyla_process |>
  ggplot(aes(x = "", y = fct_reorder(taxon, avg_percent))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(rarity~comparison_code, drop=FALSE, scales = "free", space = "free_y") +
  labs(title="Differentially abundant phyla for factor Process",
       subtitle="Mean abundance: common (>1%), rare (<0.1%); significant results only",
       y="Phylum",
       x="Comparison") +
  geom_tile(aes(fill = log2fc)) +
  rescale_2gradient(phyla_process$log2fc)+
  labs(fill="Log2 fold-change")
phy_plot_process
ggsave(glue("{plot_out}/process_phyla_heatmap.png"), phy_plot_process, width = 10, height = 6)

## Genus ----
genus_counts <- get_processed_taxonomy("shallow", "genus", timed_meta)

### Month factor ----
genus_month <- get_combined_month_DA_results("genus") |>
  format_taxa_results(genus_counts, "genus", "month")

g_month_plot_tab <- genus_month |>
  filter(occurrence_tally > 1) |>  # keep only phyla in >1 subset
  filter(rarity == "common")

g_plot_month <- g_month_plot_tab |>
  ggplot(aes(x = subset, y = fct_reorder(taxon, avg_percent))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8)) +
  facet_grid(rarity~facet_group, scales = "free", space = "free") +
  labs(title="Differentially abundant genera for factor Month",
       subtitle="significant results only; common (mean abundance >0.05%) genera only",
       y="Genus") +
  geom_tile(aes(fill = log2fc*18), color="black") +
  rescale_2gradient(g_month_plot_tab$log2fc) +
  labs(fill=str_wrap("Log2 fold-change (18 months)",16))
g_plot_month
ggsave(glue("{plot_out}/month_genus_heatmap.png"), g_plot_month, width = 8, height = 6)

### Process factor ----

genus_process <- read_csv(glue("{main_out}/shallow_sequencing/genus/da-table_Process.csv")) |>
  format_taxa_results(genus_counts, "genus", "process") |>
  format_process_comparisons()

g_process_plot_tab <- genus_process |>
  filter(occurrence_tally > 1) |> # keep only phyla in >1 subset
  filter(rarity == "common")
g_process_plot_tab$rarity <- g_process_plot_tab$rarity |>
  as.character() |>
  fct(levels=c("common")) # refactor with only common for proper plotting

g_plot_process <- g_process_plot_tab |>
  ggplot(aes(x = "", y = fct_reorder(taxon, avg_percent))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(rarity~comparison_code, drop=FALSE, scales = "free", space = "free_y") +
  labs(title="Differentially abundant genera for factor Process",
       subtitle = "significant results only; common (mean abundance >0.05%) genera only",
       y="genus",
       x="Comparison") +
  geom_tile(aes(fill = log2fc)) +
  rescale_2gradient(g_process_plot_tab$log2fc) +
  labs(fill="Log2 fold-change")
g_plot_process
ggsave(glue("{plot_out}/process_genus_heatmap.png"), g_plot_process, width = 10, height = 6)

# Heatmaps - Woltka ------------------------------------------

## Month factor ----
pw_counts <- get_processed_pathways(timed_meta)

pw_month_results <- read_csv(glue("{main_out}/deep_sequencing/woltka_pathways/da-table_month-continuous.csv"))
pw_avg_pathway_counts <- pw_counts$long |>
  summarise(avg_pw_reads = mean(reads), .by = pathway_name)

pw_month_plot_tab <- pw_month_results |>
  left_join(pw_avg_pathway_counts) |>
  filter(diff & passed.ss)

pw_month_plot <- pw_month_plot_tab |>
  ggplot(aes(x = "", y = fct_reorder(pathway_name, avg_pw_reads))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(fct_reorder(factor(pathway_group), -avg_pw_reads, .fun = mean)~., scales = "free", space = "free", switch = "y") +
  labs(title=str_wrap("Differentially abundant Woltka pathways for factor Month", 40),
       subtitle = "significant results only",
       y="Pathway group & pathway",
       x=element_blank()) +
  geom_tile(aes(fill = log2fc*18)) +
  rescale_2gradient(pw_month_plot_tab$log2fc) +
  labs(fill=str_wrap("Log2 fold-change (18 months)",16)) +
  theme(strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside") +
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 25))
pw_month_plot
ggsave(glue("{plot_out}/month_woltka_pathway_heatmap.png"), pw_month_plot, width = 8, height = 9)

## Process factor ----
pw_process_results <- read_csv(glue("{main_out}/deep_sequencing/woltka_pathways/da-table_Process.csv")) |>
  format_process_comparisons()

pw_process_plot_tab <- pw_process_results |>
  left_join(pw_avg_pathway_counts) |>
  filter(diff & passed.ss)

pw_process_plot <- pw_process_plot_tab |>
  ggplot(aes(x = "", y = fct_reorder(pathway_name, avg_pw_reads))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_grid(fct_reorder(factor(pathway_group),-avg_pw_reads, .fun = mean)~comparison_code,
             scales = "free", space = "free", switch = "y") +
  labs(title=str_wrap("Differentially abundant Woltka pathways for factor Process", 40),
       subtitle = "significant results only",
       y="Pathway group & pathway",
       x=element_blank()) +
  geom_tile(aes(fill = log2fc), color="black") +
  rescale_2gradient(pw_process_plot_tab$log2fc) +
  labs(fill=str_wrap("Log2 fold-change (18 months)",16)) +
  theme(strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside")
pw_process_plot
ggsave(glue("{plot_out}/process_woltka_pathway_heatmap.png"), pw_process_plot, width = 12, height = 20)
