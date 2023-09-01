library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(indicspecies)

source("scripts/aux_functions.R")

# 0.1 Define globals/output paths --------------------------------------------------
run_shallow_seqdepth <- FALSE # FALSE means run with deep seqdepth
taxa_level <- "phylum"

# main_out <- "results/shallow/"
main_out <- ifelse(run_shallow_seqdepth, "results/shallow/", "results/deep/")
taxa_out <- glue("{main_out}/{taxa_level}/")

# Create output directories
lapply(c(main_out, taxa_out), dir.create, showWarnings = FALSE)


# 0.2 Read in metadata ----------------------------------------------------
timed_meta <- get_timed_metadata()


# ---- 1. Read in and filter taxonomy data ----------------------------------------
taxa_long <- if(run_shallow_seqdepth){
  read_tsv(glue("data/shallow_taxonomy/kaiju_{taxa_level}_summary.tsv"))
} else {
  read_tsv(glue("data/deep_Timed/kaiju_{taxa_level}_summary.tsv"))
}

# Reformat some columns
taxa_1_long <- taxa_long |>
  mutate(sample = str_extract(file, "S00JY-\\d{4}")) |>
  rename("taxon_lineage" = taxon_name) |>
  mutate(taxon_lineage = str_replace(taxon_lineage, "cellular organisms;", ""))

# Remove pre-2020 and non-Timed experiment samples, filter unclassified taxa,
# remove suspect sample 0597, and apply 0.01% per-sample threshold
taxa_2_filtered <- taxa_1_long |>
  filter(sample %in% filter(timed_meta, Year == 2020)$sample) |>
  filter(!str_detect(taxon_lineage, "unclassified$|^cannot|Viruses|[ ;]bacterium[; ]")) |>
  filter(sample != "S00JY-0597") |>
  filter(reads > sum(reads)*0.0001, .by = sample)

# Calculate relative abundance
taxa_3_relabund <- taxa_2_filtered |>
  mutate(rel_abund = reads/sum(reads), .by = sample)

# Wide relative abundance for following analyses
taxa_4_wide_relabund <- taxa_3_relabund |>
  select(-c(file, taxon_id, percent, reads)) |>
  pivot_wider(names_from = taxon_lineage, values_from = rel_abund)


# 2. Perform NMDS, ANOSIM and Indicator Species Analyses plots --------------------
anosim_result_file <- glue("{taxa_out}/nmds_ANOSIM.txt")
cat(glue("ANOSIM Results for taxa level {str_to_upper(taxa_level)}\n\n\n"),
    file = anosim_result_file)

for (foi in names(select(timed_meta, !matches("sample|Site|Year")))) {
  # Skip the metadata factor if there is only 1 distinct value in it after filtering
  if(length(unique(filter(timed_meta, sample %in% taxa_4_wide_relabund$sample)[[foi]])) < 2) { next }

  random_seed <- 86752155
  formatted_foi <- prep_foi_data(taxa_4_wide_relabund, timed_meta, foi)

  anosim_results <- run_anosim(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
  save_anosim(anosim_result_file, anosim_results$ano_string, foi)

  nmds_results <- run_nmds(formatted_foi$matrix, formatted_foi$metadata, random_seed)
  save_nmds(nmds_results$nmds_scores, anosim_results$ano_string, foi, taxa_out, taxa_level)

  if (foi != "Month") {
    indic_results <- run_indicspecies(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
    save_indicspecies(indic_results$indic_obj, indic_results$indic_table, foi, taxa_out)
  }

  # Run analyses separately for different soil depths for certain factors
  if (!foi %in% c("sample","Site","Year", "Depth")){
    for (cur_depth in c("0-15", "15-30")) {
      filtered_meta <- filter(timed_meta, Depth == cur_depth)

      formatted_foi <- prep_foi_data(taxa_4_wide_relabund, filtered_meta, foi)

      anosim_results <- run_anosim(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
      save_anosim(anosim_result_file, anosim_results$ano_string, foi, cur_depth = cur_depth)

      nmds_results <- run_nmds(formatted_foi$matrix, formatted_foi$metadata, random_seed)
      save_nmds(nmds_results$nmds_scores, anosim_results$ano_string, foi, taxa_out, taxa_level, cur_depth = cur_depth)

      if (foi != "Month") {
        indic_results <- run_indicspecies(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
        save_indicspecies(indic_results$indic_obj, indic_results$indic_table, foi, taxa_out, cur_depth = cur_depth)
      }
    }
  }
}


# 3. Alpha diversity ------------------------------------------------------

calc_diversity_df <- function(d){
  sample_col <- select(d, "sample")
  sample_data <- select(d, -"sample")
  observed_richness <- specnumber(sample_data)
  invsimpson <- diversity(sample_data, index="invsimpson")
  shannon <- diversity(sample_data, index="shannon")
  div_df <- data.frame(
    ID = sample_col,
    Observed_Richness = observed_richness,
    Inv_Simpson = invsimpson,
    Shannon = shannon
  )
  return(div_df)
}

# runs and saves kruskal wallis on shannon index for both
# replicates and treatments
alpha_KW_test <- function(d, alpha_index, meta_factor, outfile) {
  heading <- glue("{alpha_index} Index by {meta_factor} - Kruskal Wallis Results:\n")
  factor_alpha <- stats::kruskal.test(stats::formula(glue("{alpha_index}~{meta_factor}")), d)
  cat(heading, file = outfile, append = T)
  utils::capture.output(factor_alpha, file = outfile, append = T)
  factor_alpha
}

# plot alpha diversity
plot_alpha <- function(d, alpha, group_var, taxa_rank, KW_result) {
  index <- sym(alpha)
  grp <- sym(group_var)
  alpha_plot <- ggplot(d, aes(x = !!grp, y = !!index)) +
    geom_boxplot() +
    geom_jitter(aes(color = Season), alpha = 0.7) +
    labs(title = glue("Alpha Diversity: {KW_result$data.name}"),
         subtitle = glue("Taxonomic rank: {taxa_rank}   KW-test p-val: {signif(KW_result$p.value, 3)}"),
         x = group_var)
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
make_all_alpha_plots <- function(data, taxa_rank, outdir) {
  alpha_div_data <- calc_diversity_df(data)
  alpha_div_data <- merge(alpha_div_data, timed_meta, by.x = "sample", by.y = "sample")

  utils::write.csv(alpha_div_data,
                   file = glue("{outdir}/alpha_div_data.csv"),
                   row.names = F)


  alpha_stats_filepath <- glue("{outdir}/alpha_div_stat_tests.txt")
  file.create(alpha_stats_filepath)

  meta_vars <- list("Month", "Process", "month_group")
  alpha_vars <- list("Observed_Richness", "Shannon", "Inv_Simpson")
  for (mv in meta_vars){
    for (av in alpha_vars) {
      if (mv == "Month" | mv == "month_group") {
        plot_data <- filter(alpha_div_data, Process == "Ground")
      } else if (mv == "Process") {
        plot_data <- filter(alpha_div_data, Month == "0")
      } else {
        plot_data <- alpha_div_data
      }
      KW_result <- alpha_KW_test(plot_data, av, mv, alpha_stats_filepath)
      plot <- plot_alpha(plot_data, av, mv, taxa_rank, KW_result)
      ggsave(plot = plot, filename = glue("{outdir}/alpha_div_{mv}_{av}.png"), bg = "white",
             width = 5.6)
    }
  }
}

# Use read counts, not relative abundance
taxa_wide_reads <- taxa_3_relabund |>
  select(-c(file, taxon_id, percent, rel_abund)) |>
  pivot_wider(names_from = taxon_lineage, values_from = reads)

taxa_wide_reads[is.na(taxa_wide_reads)] <- 0

make_all_alpha_plots(taxa_wide_reads, taxa_level, taxa_out)

