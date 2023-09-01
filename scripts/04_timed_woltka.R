library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(indicspecies)

source("scripts/aux_functions.R")

outdir <- "results/woltka/rpk"
dir.create(outdir)

# 0.2 Read in metadata ----------------------------------------------------
timed_meta <- get_timed_metadata()


# 1. Read in format KO gene count data ------------------------------------

ko <- read_tsv("data/ko_rpk_filtered-0-001p.tsv") |>
  relocate("Name") |>
  rename("name" = Name,
         "ko_gene_ID" = `#FeatureID`)

# Pivot longer and get relative abundance
ko_long <- ko |>
  pivot_longer(cols = starts_with("S00JY"),
               names_to = "sample",
               values_to = "reads") |>
  mutate(relative_abundance = reads/sum(reads), .by = sample,
         sample = str_extract(sample,"S00JY-\\d\\d\\d\\d"))

# Pivot for samples as rows, genes as columns
ko_pivot <- ko_long |>
  pivot_wider(id_cols = sample,
              names_from = ko_gene_ID,
              values_from = reads)

# Pivot for samples as rows, genes as columns
rel_ko_pivot <- ko_long |>
  pivot_wider(id_cols = sample,
              names_from = ko_gene_ID,
              values_from = relative_abundance)

# FILTER out the weirdo sample
rel_ko_pivot <- rel_ko_pivot |>
  filter(sample != "S00JY-0597" & sample %in% filter(timed_meta, Year == 2020)$sample)

anosim_result_file <- glue("{outdir}/nmds_ANOSIM.txt")
cat("ANOSIM Results for WOLTKA functional profiling\n\n\n",
    file = anosim_result_file)
for (foi in names(select(timed_meta, !matches("sample|Site|Year")))) {
  # Skip the metadata factor if there is only 1 distinct value in it after filtering
  if(length(unique(filter(timed_meta, sample %in% rel_ko_pivot$sample)[[foi]])) < 2) { next }

  # results <- run_ano_nmds_indic(rel_ko_pivot, timed_meta, foi)

  random_seed <- 86752155
  formatted_foi <- prep_foi_data(rel_ko_pivot, timed_meta, foi)

  anosim_results <- run_anosim(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
  save_anosim(anosim_result_file, anosim_results$ano_string, foi)

  nmds_results <- run_nmds(formatted_foi$matrix, formatted_foi$metadata, random_seed)
  save_nmds(nmds_results$nmds_scores, anosim_results$ano_string, foi, outdir)

  if (foi != "Month") {
    indic_results <- run_indicspecies(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
    save_indicspecies(indic_results$indic_obj, indic_results$indic_table, foi, outdir)
  }

  # Run analyses separately for different soil depths for certain factors
  if (!foi %in% c("sample","Site","Year", "Depth")){
    for (cur_depth in c("0-15", "15-30")) {
      filtered_meta <- filter(timed_meta, Depth == cur_depth)
      # results_cur_depth <- run_ano_nmds_indic(rel_ko_pivot, filtered_meta, foi)

      formatted_foi <- prep_foi_data(rel_ko_pivot, filtered_meta, foi)

      anosim_results <- run_anosim(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
      save_anosim(anosim_result_file, anosim_results$ano_string, foi, cur_depth = cur_depth)

      nmds_results <- run_nmds(formatted_foi$matrix, formatted_foi$metadata, random_seed)
      save_nmds(nmds_results$nmds_scores, anosim_results$ano_string, foi, outdir, cur_depth = cur_depth)

      if (foi != "Month") {
        indic_results <- run_indicspecies(formatted_foi$matrix, formatted_foi$metadata, foi, random_seed)
        save_indicspecies(indic_results$indic_obj, indic_results$indic_table, foi, outdir, cur_depth = cur_depth)
      }
    }
  }
}


# Pathway Analysis --------------------------------------------------------



