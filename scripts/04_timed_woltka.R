library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(indicspecies)

source("scripts/aux_functions.R")

outdir <- "results/woltka/rpk/filtered/"
dir.create(outdir)

# 0.2 Read in metadata ----------------------------------------------------
timed_meta <- get_timed_metadata()


sample_read_counts <- read.csv("data/qcd_read_count.txt", header = FALSE)

# Add forward and reverse read counts together
sample_read_counts <- sample_read_counts |>
  rename("sample" = V1) |>
  summarise(reads = sum(V2), .by = sample) |>
  mutate(sample = str_extract(sample,"S00JY-\\d\\d\\d\\d"))

# Read protein count data
ko <- read_tsv("data/ko_rpk.tsv")
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

# ko_filtered <- read_tsv("data/ko_rpk_filtered-0-001p.tsv") |>
#   relocate("Name") |>
#   rename("name" = Name,
#          "ko_gene_ID" = `#FeatureID`)
#
# # Pivot longer and get relative abundance
# ko_long_filtered <- ko_filtered |>
#   pivot_longer(cols = starts_with("S00JY"),
#                names_to = "sample",
#                values_to = "reads") |>
#   mutate(relative_abundance = reads/sum(reads), .by = sample,
#          sample = str_extract(sample,"S00JY-\\d\\d\\d\\d"))

# write.csv(ko_long, "data/woltka_timed_protein_ko_rpk_LONG.csv", row.names = FALSE)

# Applying a threshold
# ko_long <- ko_long |>
#   mutate(avg_gene_reads = mean(reads), .by = ko_gene_ID) |>
#   filter(avg_gene_reads > 50)

sample_read_counts <- ko_long |>
  mutate(total_gene_reads = sum(reads), .by = sample) |>
  select(sample, total_gene_reads) |>
  distinct() |> merge(sample_read_counts) |>
  mutate(percent_mapped_to_gene = total_gene_reads/reads*100)
# write.csv(sample_read_counts, "data/timed_deep_samples_read_count.csv", row.names = FALSE)

sample_read_counts |> ggplot(aes(x = "", y = percent_mapped_to_gene)) +
  geom_violin() +
  geom_jitter()

# Pivot for samples as rows, genes as columns
ko_pivot <- ko_long |>
  pivot_wider(id_cols = sample,
              names_from = ko_gene_ID,
              values_from = reads)

# Now for relative abundance values
rel_ko_pivot <- ko_long |>
  pivot_wider(id_cols = sample,
              names_from = ko_gene_ID,
              values_from = relative_abundance)

# FILTER out the weirdo sample
rel_ko_pivot <- rel_ko_pivot |>
  filter(sample != "S00JY-0597")

anosim_result_file <- glue("{outdir}/anosim.txt")
cat("ANOSIM Results for WOLTKA functional profiling\n\n\n",
    file = anosim_result_file)
for (foi in names(select(timed_meta, !matches("sample|Site|Year")))) {
  filtered_meta <- timed_meta |>
    filter(sample %in% rel_ko_pivot$sample & Year == 2020)

  # specific filters for Month and Process factors
  if (foi == "Month"){
    filtered_meta <- filtered_meta |> filter(Process == "Ground")
  } else if (foi == "Process") {
    filtered_meta <- filtered_meta |> filter(Month == "0")
  }

  # filter for the appropriate subset of samples
  ko_foi <- rel_ko_pivot |>
    filter(sample %in% filtered_meta$sample)

  # a meta table with the same sample order used for the abundance matrix
  foi_metadata <- merge(ko_foi$sample, filtered_meta,
                        by.x = "x", by.y = "sample") |>
    rename(sample = x)
  # Skip the metadata factor if there is only 1 distinct value in it after filtering
  if(length(unique(foi_metadata[[foi]])) < 2) { next }

  rel_abund_matrix_foi <- select(ko_foi, -sample) |>
    as.matrix()
  # Replace NAs with 0
  rel_abund_matrix_foi[is.na(rel_abund_matrix_foi)] <- 0

  set.seed(1234)
  ano <- anosim(rel_abund_matrix_foi, foi_metadata[[foi]], distance = "bray", permutations = 999)

  # Write result to file
  cat(glue("FACTOR: {foi}\n"), file = anosim_result_file, append = T)
  utils::capture.output(
    ano,
    file = anosim_result_file, append = T)

  subtitle <- glue("   ANOSIM p-val: {ano$signif} R-stat: {ano$statistic}")

  ## Perform NMDS ----

  nmds = metaMDS(rel_abund_matrix_foi, distance = "bray")
  plot(nmds)

  data.scores = as.data.frame(scores(nmds)$sites)
  data.scores$sample <- ko_foi$sample

  # Merge NMDS results with metadata ----
  data.scores <- merge(data.scores, foi_metadata)
  p <- plot_nmds_by_factor(
    df = data.scores,
    meta_factor = foi,
    dataset_name = "Timed",
    shape_factor = "Depth",
    polygon = TRUE,
    save_image = TRUE,
    save_path = glue("{outdir}/{foi}_nmds_by_depth.png"),
    subtitle_append = subtitle
  )


  # Perform indicspecies analysis ----
  # indicspecies wrapper function
  run_indicspecies <- function(rel_abund_matrix, metadata_tb, foi){
    indic_by_foi = multipatt(rel_abund_matrix, metadata_tb[[foi]], func = "r.g",
                             control = how(nperm=9999))

    output_tbl <- indic_by_foi$sign |>
      mutate("significant" = p.value <= 0.05) |>
      rownames_to_column("taxon_lineage") |>
      arrange(p.value, -stat)

    write_csv(output_tbl, glue("{outdir}/indicspecies_{foi}.csv"))

    sink(glue("{outdir}/indicspecies_{foi}.txt"))
    summary(indic_by_foi)
    sink()
  }
  if (foi == "Month") {
    foi_metadata <- foi_metadata |>
      mutate(month_group = case_when(Month < 1 ~ "<1 month",
                                     Month <= 6 ~ "1-6 months",
                                     Month <= 18 ~ "12-18 months"))
    meta_factor <- "month_group"
  } else {
    meta_factor <- foi
  }

  run_indicspecies(rel_abund_matrix_foi, foi_metadata, meta_factor)
}

# Run Month separately for separate soil depths
for (foi in c("Month", "Season")) {
  for (soil_depth in c("0-15", "15-30")) {

    filtered_meta <- timed_meta |>
      filter(Depth == soil_depth) |>
      filter(sample %in% rel_ko_pivot$sample & Year == 2020)

    # specific filters for Month
    if (foi == "Month"){
      filtered_meta <- filtered_meta |> filter(Process == "Ground")
    }

    # filter for the appropriate subset of samples
    ko_foi <- rel_ko_pivot |>
      filter(sample %in% filtered_meta$sample)

    # a meta table with the same sample order used for the abundance matrix
    foi_metadata <- merge(ko_foi$sample, filtered_meta,
                          by.x = "x", by.y = "sample") |>
      rename(sample = x)
    # Skip the metadata factor if there is only 1 distinct value in it after filtering
    if(length(unique(foi_metadata[[foi]])) < 2) { next }

    rel_abund_matrix_foi <- select(ko_foi, -sample) |>
      as.matrix()
    # Replace NAs with 0
    rel_abund_matrix_foi[is.na(rel_abund_matrix_foi)] <- 0

    set.seed(1234)
    ano <- anosim(rel_abund_matrix_foi, foi_metadata[[foi]], distance = "bray", permutations = 999)

    # Write result to file
    cat(glue("FACTOR: {foi} Soil-depth: {soil_depth}\n"), file = anosim_result_file, append = T)
    utils::capture.output(
      ano,
      file = anosim_result_file, append = T)

    subtitle <- glue("   Soil-depth: {soil_depth} ANOSIM p-val: {ano$signif} R-stat: {ano$statistic}")

    ## Perform NMDS ----

    nmds = metaMDS(rel_abund_matrix_foi, distance = "bray")
    plot(nmds)

    data.scores = as.data.frame(scores(nmds)$sites)
    data.scores$sample <- ko_foi$sample

    # Merge NMDS results with metadata ----
    data.scores <- merge(data.scores, foi_metadata)
    p <- plot_nmds_by_factor(
      df = data.scores,
      meta_factor = foi,
      dataset_name = "Timed",
      shape_factor = "Season",
      polygon = TRUE,
      save_image = TRUE,
      save_path = glue("{outdir}/{foi}_nmds_by_depth_{soil_depth}.png"),
      subtitle_append = subtitle
    )


    # Perform indicspecies analysis ----
    # indicspecies wrapper function
    run_indicspecies <- function(rel_abund_matrix, metadata_tb, foi){
      indic_by_foi = multipatt(rel_abund_matrix, metadata_tb[[foi]], func = "r.g",
                               control = how(nperm=9999))

      output_tbl <- indic_by_foi$sign |>
        mutate("significant" = p.value <= 0.05) |>
        rownames_to_column("taxon_lineage") |>
        arrange(p.value, -stat)

      write_csv(output_tbl, glue("{outdir}/indicspecies_{foi}_{soil_depth}.csv"))

      sink(glue("{outdir}/indicspecies_{foi}_{soil_depth}.txt"))
      summary(indic_by_foi)
      sink()
    }
    if (foi == "Month") {
      foi_metadata <- foi_metadata |>
        mutate(month_group = case_when(Month < 1 ~ "<1 month",
                                       Month <= 6 ~ "1-6 months",
                                       Month <= 18 ~ "12-18 months"))
      meta_factor <- "month_group"
    } else {
      meta_factor <- foi
    }

    run_indicspecies(rel_abund_matrix_foi, foi_metadata, meta_factor)
  }
}


