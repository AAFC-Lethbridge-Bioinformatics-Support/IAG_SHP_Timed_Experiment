library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(indicspecies)

source("scripts/aux_functions.R")

# 0.1 Define globals/output paths --------------------------------------------------
taxa_level <- "genus"

main_out <- "results/shallow/"
taxa_out <- glue("{main_out}/{taxa_level}")
nmds_out <- glue("{taxa_out}/NMDS")
indicspecies_out <- glue("{taxa_out}/indicspecies")
alphadiv_out <- glue("{taxa_out}/alpha_diversity")

create_dir_if_nonexistant <- function(path) {
  ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
}

lapply(c(main_out, taxa_out, nmds_out, indicspecies_out, alphadiv_out),
       create_dir_if_nonexistant)

# ---- 0.2 Read in and tidy timed metadata -----------------------------------------
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

# how many NA's are there in each column?
sapply(timed_meta, function(x) sum(length(which(is.na(x)))))

# ---- 1. Read in and filter taxonomy data ----------------------------------------
taxa_long_deep <- read_tsv(glue("data/deep_Timed/kaiju_{taxa_level}_summary.tsv"))
taxa_long_shallow <- read_tsv(glue("data/shallow_taxonomy/kaiju_{taxa_level}_summary.tsv"))

# reformat some column aspects
taxa_long <- taxa_long_shallow |>
  mutate(sample = str_extract(file, "S00JY-\\d{4}")) |>
  rename("taxon_lineage" = taxon_name) |>
  mutate(taxon_lineage = str_replace(taxon_lineage, "cellular organisms;", ""))

# remove non-Timed experiment samples, filter unclassified, apply read
# thresholds at 0.01% on per-sample basis
taxa_long_filtered <- taxa_long |>
  filter(sample %in% timed_meta$sample) |>
  filter(!str_detect(taxon_lineage, "unclassified$|^cannot|Viruses|[ ;]bacterium[; ]")) |>
  mutate(sample_threshold = sum(reads)*0.0001, .by = sample) |>
  filter(reads > sample_threshold) |>
  select(-sample_threshold)

# recalculate percentage so that it adds to 100 after removing unclassified and
# applying the threshold
taxa_longf_total <- taxa_long_filtered |>
  mutate(recalculated_percent = reads/sum(reads), .by = sample)

# remove suspect sample and pre-2020 samples
taxa_longf_total <- taxa_longf_total |>
  filter(sample != "S00JY-0597") |>
  filter(sample %in% filter(timed_meta, Year == 2020)$sample)

# 2. Perform NMDS, ANOSIM and Indicator Species Analyses plots --------------------
anosim_result_file <- glue("{nmds_out}/anosim_{taxa_level}.txt")
cat(glue("ANOSIM Results for taxa level {str_to_upper(taxa_level)}\n\n\n"),
    file = anosim_result_file)

for (foi in names(select(timed_meta, !matches("sample|Site|Year")))) {

  filtered_meta <- timed_meta

  # specific filters for Month and Process factors
  if (foi == "Month"){
    filtered_meta <- filtered_meta |> filter(Process == "Ground")
  } else if (foi == "Process") {
    filtered_meta <- filtered_meta |> filter(Month == "0")
  }

  # filter for the appropriate subset of samples
  taxa_long_foi <- taxa_longf_total |>
    filter(sample %in% filtered_meta$sample)

  # Pivot and make proportional for proper ANOSIM format
  taxa_wide_foi <- taxa_long_foi |>
    select(-c(file, taxon_id, percent, reads)) |>
    pivot_wider(names_from = taxon_lineage, values_from = recalculated_percent)

  # a meta table with the same sample order used for the abundance matrix
  foi_metadata <- select(taxa_wide_foi, sample) |>
    left_join(filtered_meta)

  # Skip the metadata factor if there is only 1 distinct value in it after filtering
  if(length(unique(foi_metadata[[foi]])) < 2) { next }

  rel_abund_matrix_foi <- select(taxa_wide_foi, -sample) |>
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

  nmds = metaMDS(rel_abund_matrix_foi, distance = "bray", weakties = FALSE)
  plot(nmds)

  data.scores = as.data.frame(scores(nmds)$sites)
  data.scores$sample <- taxa_wide_foi$sample

  # Merge NMDS results with metadata ----
  data.scores <- left_join(data.scores, foi_metadata)
  nmds_plot <- plot_nmds_by_factor(
    df = data.scores,
    meta_factor = foi,
    dataset_name = "Timed",
    taxa_level = taxa_level,
    shape_factor = "Depth",
    polygon = TRUE,
    save_image = TRUE,
    save_path = glue("{nmds_out}/{foi}_{taxa_level}_nmds_by_depth.png"),
    subtitle_append = subtitle
  )
  nmds_plot

  # Perform indicspecies analysis ----
  # indicspecies wrapper function
  run_indicspecies <- function(rel_abund_matrix, metadata_tb, foi){
    indic_by_foi = multipatt(rel_abund_matrix, metadata_tb[[foi]], func = "r.g",
                             control = how(nperm=9999))

    output_tbl <- indic_by_foi$sign |>
      mutate("significant" = p.value <= 0.05) |>
      rownames_to_column("taxon_lineage") |>
      arrange(p.value, -stat)

    write_csv(output_tbl, glue("{indicspecies_out}/indicspecies_{foi}.csv"))

    sink(glue("{indicspecies_out}/indicspecies_{foi}.txt"))
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



# NMDS split by soil depth ------------------------------------------------
anosim_result_file <- glue("{nmds_out}/anosim_{taxa_level}_by_soil_depth.txt")
cat(glue("ANOSIM Results for taxa level {str_to_upper(taxa_level)}\n\n\n"),
    file = anosim_result_file)

for (foi in names(select(timed_meta, !matches("sample|Site|Year")))) {
  for (cur_depth in c("0-15", "15-30")) {
    filtered_meta <- timed_meta |> filter(Depth == cur_depth)

    # specific filters for Month and Process factors
    if (foi == "Month"){
      filtered_meta <- filtered_meta |> filter(Process == "Ground")
    } else if (foi == "Process") {
      filtered_meta <- filtered_meta |> filter(Month == "0")
    }

    # filter for the appropriate subset of samples
    taxa_long_foi <- taxa_longf_total |>
      filter(sample %in% filtered_meta$sample)

    # Pivot and make proportional for proper ANOSIM format
    taxa_wide_foi <- taxa_long_foi |>
      select(-c(file, taxon_id, percent, reads)) |>
      pivot_wider(names_from = taxon_lineage, values_from = recalculated_percent)

    # a meta table with the same sample order used for the abundance matrix
    foi_metadata <- merge(taxa_wide_foi$sample, filtered_meta,
                          by.x = "x", by.y = "sample") |>
      rename(sample = x)
    if(length(unique(foi_metadata[[foi]])) < 2) { next }

    rel_abund_matrix_foi <- select(taxa_wide_foi, -sample) |>
      as.matrix()
    # Replace NAs with 0
    rel_abund_matrix_foi[is.na(rel_abund_matrix_foi)] <- 0

    set.seed(1234)
    ano <- anosim(rel_abund_matrix_foi, foi_metadata[[foi]], distance = "bray", permutations = 999)

    # Write result to file
    cat(glue("FACTOR: {foi} {cur_depth}\n"), file = anosim_result_file, append = T)
    utils::capture.output(
      ano,
      file = anosim_result_file, append = T)
    subtitle <- glue(" Soil depth: {cur_depth} (cm)   ANOSIM p-val: {ano$signif} R-stat: {ano$statistic}")

    # Perform NMDS ----

    nmds = metaMDS(rel_abund_matrix_foi, distance = "bray")
    plot(nmds)

    data.scores = as.data.frame(scores(nmds)$sites)
    data.scores$sample <- taxa_wide_foi$sample

    # Merge NMDS results with metadata ----
    data.scores <- merge(data.scores, foi_metadata)
    plot_nmds_by_factor(
      df = data.scores,
      meta_factor = foi,
      dataset_name = "Timed",
      taxa_level = taxa_level,
      shape_factor = "Season",
      polygon = TRUE,
      save_image = TRUE,
      save_path = glue("{nmds_out}/{foi}_{taxa_level}_nmds_{cur_depth}.png"),
      subtitle_append = subtitle
    )

    run_indicspecies <- function(rel_abund_matrix, metadata_tb, foi){
      indic_by_foi = multipatt(rel_abund_matrix, metadata_tb[[foi]], func = "r.g",
                               control = how(nperm=9999))

      output_tbl <- indic_by_foi$sign |>
        mutate("significant" = p.value <= 0.05) |>
        rownames_to_column("taxon_lineage") |>
        arrange(p.value, -stat)

      write_csv(output_tbl, glue("{indicspecies_out}/indicspecies_{foi}_{cur_depth}.csv"))

      sink(glue("{indicspecies_out}/indicspecies_{foi}_{cur_depth}.txt"))
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
}

# plot alpha diversity
plot_alpha <- function(d, alpha, group_var) {
  index <- sym(alpha)
  grp <- sym(group_var)
  alpha_plot <- ggplot(d, aes(x = !!grp, y = !!index)) +
    geom_boxplot() +
    geom_jitter(aes(color = Season), alpha = 0.7) +
    labs(title = "Alpha Diversity",
         subtitle = glue("Metric: {alpha}"),
         x = group_var)
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
make_all_alpha_plots <- function(data, taxa_rank, outdir) {
  alpha_div_data <- calc_diversity_df(data)
  alpha_div_data <- merge(alpha_div_data, timed_meta, by.x = "sample", by.y = "sample")
  alpha_div_data$Month <- factor(alpha_div_data$Month,
                                 levels = c("0", "0.5", "0.07", "1", "3", "6", "12", "18", "1998", "2007"))
  # filter out Blanks
  alpha_div_data <- filter(alpha_div_data, Month != "Blank")

  utils::write.csv(alpha_div_data,
                   file = glue("{outdir}/alpha_div_data_{taxa_rank}.csv"),
                   row.names = F)


  alpha_stats_filepath <- glue("{outdir}/alpha_stats.txt")
  file.create(alpha_stats_filepath)

  meta_vars <- list("Month", "Process")
  alpha_vars <- list("Observed_Richness", "Shannon", "Inv_Simpson")
  for (mv in meta_vars){
    for (av in alpha_vars) {
      if (mv == "Month") {
        plot_data <- filter(alpha_div_data, Process == "Ground")
      } else if (mv == "Process") {
        plot_data <- filter(alpha_div_data, Month == "0")
      } else {
        plot_data <- alpha_div_data
      }
      alpha_KW_test(plot_data, av, mv, alpha_stats_filepath)
      plot <- plot_alpha(plot_data, av, mv)
      ggsave(plot = plot, filename = glue("{outdir}/alpha_div_{mv}_{av}_{taxa_rank}.png"), bg = "white",
             width = 7.6)
    }
  }
}

taxa_wide_reads <- taxa_longf_total |>
  select(-c(file, taxon_id, percent, recalculated_percent)) |>
  pivot_wider(names_from = taxon_lineage, values_from = reads)

taxa_wide_reads[is.na(taxa_wide_reads)] <- 0

make_all_alpha_plots(taxa_wide_reads, taxa_level, alphadiv_out)

