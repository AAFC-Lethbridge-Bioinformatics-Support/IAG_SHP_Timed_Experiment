pkgs = c("tidyverse", "vegan", "glue", "indicspecies")
lapply(pkgs, library, character.only = TRUE)

source("scripts/aux_functions.R")


# 0.1 Define globals/output paths --------------------------------------------------
main_out <- "results/deep/timed/phyloflash"
nmds_out <- glue("{main_out}/NMDS")
indicspecies_out <- glue("{main_out}/indicspecies")
alphadiv_out <- glue("{main_out}/alpha_diversity")

create_dir_if_nonexistant <- function(path) {
  ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
}

lapply(c(main_out, nmds_out, indicspecies_out, alphadiv_out),
       create_dir_if_nonexistant)


# 0.2 Read in metadata ----------------------------------------------------
timed_meta <- get_timed_metadata()


# ---- 1. Read in and filter taxonomy data ----------------------------------------
pf_long = read.csv("pipeline_outputs/deep_Timed/timed_phyloflash_NTUabundances_concatenated.csv", header = FALSE)

pf_long <- pf_long |>
  rename(
    taxa_lineage = V1,
    count = V2,
    sample = V3
  )

# apply 0.01% threshold, per sample
pf_long_filtered <- pf_long |>
  mutate(sample_threshold = sum(count)*0.0001, .by = sample) |>
  filter(count > sample_threshold)

# calculate relative abundance after filter
pf_longf_total <- pf_long_filtered |>
  mutate(total_sample_count = sum(count),
         rel_abund = count/total_sample_count*100,
         .by = sample)

pf_wide <- pf_longf_total |>
  pivot_wider(id_cols = sample,
              names_from = taxa_lineage,
              values_from = rel_abund)

pf_wide[is.na(pf_wide)] <- 0


# 2. Perform NMDS, ANOSIM and Indicator Species Analyses plots --------------------
anosim_result_file <- glue("{nmds_out}/anosim.txt")
cat(glue("ANOSIM Results for Phyloflash \n\n\n"),
    file = anosim_result_file)

for (foi in names(select(timed_meta, !matches("sample|Site|Season|Year")))) {

  filtered_meta <- timed_meta |>
    filter(sample %in% pf_wide$sample)

  # specific filters for Month and Process factors
  if (foi == "Month"){
    filtered_meta <- filtered_meta |> filter(Process == "Ground")
  } else if (foi == "Process") {
    filtered_meta <- filtered_meta |> filter(Month == "0")
  }

  # filter for the appropriate subset of samples
  taxa_wide_foi <- pf_wide |>
    filter(sample %in% filtered_meta$sample)

  # a meta table with the same sample order used for the abundance matrix
  foi_metadata <- merge(taxa_wide_foi$sample, filtered_meta,
                        by.x = "x", by.y = "sample") |>
    rename(sample = x)
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

  nmds = metaMDS(rel_abund_matrix_foi, distance = "bray")
  plot(nmds)

  data.scores = as.data.frame(scores(nmds)$sites)
  data.scores$sample <- taxa_wide_foi$sample

  # Merge NMDS results with metadata
  data.scores <- merge(data.scores, foi_metadata)
  plot_nmds_by_factor(
    df = data.scores,
    meta_factor = foi,
    dataset_name = "Timed",
    shape_factor = "Season",
    polygon = TRUE,
    save_image = TRUE,
    save_path = glue("{nmds_out}/{foi}_nmds.png"),
    subtitle_append = subtitle
  )

  run_indicspecies <- function(rel_abund_matrix, metadata_tb, foi){
    indic_by_foi = multipatt(rel_abund_matrix, metadata_tb[[foi]], func = "r.g",
                             control = how(nperm=9999),
                             duleg = TRUE)

    output_tbl <- indic_by_foi$sign |>
      mutate("significant" = p.value <= 0.05) |>
      rownames_to_column("taxon_lineage") |>
      arrange(p.value, -stat)

    write_csv(output_tbl, glue("{indicspecies_out}/indicspecies_{foi}.csv"))

    sink(glue("{indicspecies_out}/indicspecies_{foi}.txt"))
    summary(indic_by_foi)
    sink()
  }

  run_indicspecies(rel_abund_matrix_foi, foi_metadata, foi)
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
make_all_alpha_plots <- function(data, outdir) {

  alpha_div_data <- calc_diversity_df(data)
  alpha_div_data <- merge(alpha_div_data, timed_meta, by.x = "sample", by.y = "sample")
  alpha_div_data$Month <- factor(alpha_div_data$Month,
                                 levels = c("0", "0.5", "0.07", "1", "3", "6", "12", "18", "1998", "2007"))
  # filter out Blanks
  alpha_div_data <- filter(alpha_div_data, Month != "Blank")

  utils::write.csv(alpha_div_data,
                   file = glue("{outdir}/alpha_div_data.csv"),
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
      ggsave(plot = plot, filename = glue("{outdir}/alpha_div_{mv}_{av}.png"), bg = "white")
    }
  }
}

make_all_alpha_plots(pf_wide, alphadiv_out)
