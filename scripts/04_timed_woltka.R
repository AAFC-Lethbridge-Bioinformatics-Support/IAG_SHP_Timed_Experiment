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

anosim_result_file <- glue("{outdir}/anosim.txt")
cat("ANOSIM Results for WOLTKA functional profiling\n\n\n",
    file = anosim_result_file)
for (foi in names(select(timed_meta, !matches("sample|Site|Year")))) {
  # Skip the metadata factor if there is only 1 distinct value in it after filtering
  if(length(unique(filter(timed_meta, sample %in% rel_ko_pivot$sample)[[foi]])) < 2) { next }
  results <- run_ano_nmds_indic(rel_ko_pivot, timed_meta, foi)

  # Save analysis results
  cat(glue("FACTOR: {foi}\n"), file = anosim_result_file, append = T)
  utils::capture.output(
    results$anosim,
    file = anosim_result_file, append = T)


  nmds_plot <- plot_nmds_by_factor(
    df = results$nmds_scores,
    meta_factor = foi,
    dataset_name = "Timed",
    shape_factor = "Depth",
    polygon = TRUE,
    save_image = TRUE,
    save_path = glue("{outdir}/nmds_{foi}.png"),
    subtitle_append = results$ano_string
  )
  nmds_plot

  if(!is.null(results$indic_table) & !is.null(results$indic)){
    write_csv(results$indic_table, glue("{outdir}/indicspecies_{foi}.csv"))

    sink(glue("{outdir}/indicspecies_{foi}.txt"))
    summary(results$indic)
    sink()
  }

  # Run analyses separately for different soil depths for certain factors
  if (!foi %in% c("sample","Site","Year", "Depth")){
    for (cur_depth in c("0-15", "15-30")) {
      filtered_meta <- filter(timed_meta, Depth == cur_depth)
      results_cur_depth <- run_ano_nmds_indic(rel_ko_pivot, filtered_meta, foi)

      # Save analysis results
      cat(glue("FACTOR: {foi} {cur_depth}\n"), file = anosim_result_file, append = T)
      utils::capture.output(
        results_cur_depth$anosim,
        file = anosim_result_file, append = T)

      nmds_plot <- plot_nmds_by_factor(
        df = results_cur_depth$nmds_scores,
        meta_factor = foi,
        dataset_name = "Timed",
        shape_factor = "Depth",
        polygon = TRUE,
        save_image = TRUE,
        save_path = glue("{outdir}/nmds_{foi}_{cur_depth}.png"),
        subtitle_append = paste(results_cur_depth$ano_string,
                                glue("Soil depth: {cur_depth} (cm)"))
      )
      nmds_plot

      if(!is.null(results_cur_depth$indic_table) & !is.null(results_cur_depth$indic)){
        write_csv(results_cur_depth$indic_table, glue("{outdir}/indicspecies_{foi}_{cur_depth}.csv"))

        sink(glue("{outdir}/indicspecies_{foi}_{cur_depth}.txt"))
        summary(results_cur_depth$indic)
        sink()
      }
    }
  }
}


# Pathway Analysis --------------------------------------------------------



