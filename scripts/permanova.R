# Load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)

source("scripts/aux_functions.R")

# 0.1 Define globals/output paths --------------------------------------------------
run_shallow_seqdepth <- TRUE # FALSE means run with deep sequencing depth
seq_depth <- ifelse(run_shallow_seqdepth, "shallow", "deep")
taxa_level <- "phylum"

n_permutations <- 9999

main_out <- glue("results/{seq_depth}")
main_out <- glue("{main_out}/permanova/")
taxa_out <- glue("{main_out}/{taxa_level}/")

# Create output directories
lapply(c(main_out, taxa_out), dir.create, recursive = TRUE, showWarnings = FALSE)


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

# Wide relative abundance for following analyses
taxa_3_wide <- taxa_2_filtered |>
  select(-c(file, taxon_id, percent)) |>
  pivot_wider(names_from = taxon_lineage, values_from = reads)

taxa_matrix <- taxa_3_wide |>
  column_to_rownames("sample") |>
  as.matrix()
taxa_matrix[is.na(taxa_matrix)] <- 0


# 2. Create components of phyloseq object ---------------------------------
taxa_otu_table <- otu_table(taxa_matrix, taxa_are_rows=FALSE)
timed_meta$month_group <- as.character(timed_meta$month_group)
proj_sample_data <- timed_meta |>
  filter(sample %in% taxa_3_wide$sample) |>
  column_to_rownames("sample") |>
  sample_data()

# Create phyloseq object
base_pseq <- phyloseq(taxa_otu_table, proj_sample_data)
# Pick relative abundances (compositional)
pseq.rel <- transform_sample_counts(base_pseq, function(x) x / sum(x))



# Run analyses ------------------------------------------------------------

fcs <- if(run_shallow_seqdepth) {
  c("month_group","POS","Season","Depth", "Process")
} else {
  c("month_group","Season","Depth", "Process")
}

summary_filename <- glue("{main_out}/permanova_results_summary.csv")
results_summary <- if (file.exists(summary_filename)){
  read.csv(summary_filename)
} else {
  data.frame(meta_factor=character(),
             taxa_level=character(),
             sequence_depth=character(),
             soil_depth=character(),
             season=character(),
             homogeneity_pval=double(),
             permanova_pval=double(),
             permanova_R2=double())
}

results_summary <- analyses_wrap(pseq = pseq.rel,
                                 fcs = fcs,
                                 base_outfile = "allfactors") %>%
  keep(!names(.) %in% c("Process", "month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth)

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Depth=="0-15"),
                                 fcs = fcs[fcs!="Depth"],
                                 base_outfile = "allfactors_0-15") %>%
  keep(!names(.) %in% c("Process", "month_group", "Depth")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Depth=="15-30"),
                                 fcs = fcs[fcs!="Depth"],
                                 base_outfile = "allfactors_15-30") %>%
  keep(!names(.) %in% c("Process", "month_group", "Depth")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30")

# Process factor
# Month==0 is required for proper Process comparison. This leaves too few
# samples to model factor interactions (e.g. POS*Process)
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Month==0),
                                 fcs = fcs[fcs!="month_group"],
                                 base_outfile = "t0",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth)

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Month==0&Depth=="0-15"),
                                 fcs = fcs[fcs!="month_group"&fcs!="Depth"],
                                 base_outfile = "t0_0-15",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Month==0&Depth=="15-30"),
                                 fcs = fcs[fcs!="month_group"&fcs!="Depth"],
                                 base_outfile = "t0_15-30",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30")


# Month factor
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"),
                                 fcs = fcs[fcs!="Process"],
                                 base_outfile = "ground") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth)

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="0-15"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                 base_outfile = "ground0-15") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground0-15spring") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     season = "Spring")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground0-15fall") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     season = "Fall")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                 base_outfile = "ground15-30") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&Season=="Spring"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground15-30spring") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     season = "Spring")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&Season=="Fall"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground15-30Fall") %>%
  keep(names(.) %in% c("month_group")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     season = "Fall")

write_csv(results_summary, file = glue("{main_out}/permanova_results_summary.csv"))


analyses_wrap <- function(pseq,
                          fcs,
                          base_outfile,
                          factor_interaction = TRUE) {
  factor_results_list <- list()
  PERMANOVA_result_file <- glue("{taxa_out}/PERMANOVA_{base_outfile}_{taxa_level}.txt")
  cat(glue("PERMANOVA via vegan::adonis2 for taxa level {str_to_upper(taxa_level)}\n\n\n"),
      file = PERMANOVA_result_file)

  otu <- otu_table(pseq)
  meta <- as(sample_data(pseq), "data.frame")
  rownames(meta) <- sample_names(pseq)
  dist <- vegdist(otu)

  homog_test <- homogeneity_test(dist, meta, fcs, PERMANOVA_result_file)
  adonis_result <- run_adonis(otu, meta, fcs, PERMANOVA_result_file, factor_interaction)

  for (fc in fcs) {
    fc_homog_pval <- homog_test[[fc]]$`Pr(>F)`[1]
    subtitle <- glue("\n\nHomogeneity of variance p.val: {signif(fc_homog_pval, 4)}")
    fc_adonis <- adonis_result$aov.tab[fc,]
    fc_adonis_R2 <- fc_adonis$R2[1]
    fc_adonis_pval <- fc_adonis$`Pr(>F)`[1]
    subtitle <- glue("{subtitle}\nAdonis PERMANOVA test p.val: {signif(fc_adonis_pval, 4)}  R2: {signif(fc_adonis_R2, 4)}")
    plot_nmds(dist, meta, fc,
              subtitle=subtitle,
              save_image=TRUE,
              save_path=glue("{taxa_out}/nmds_{base_outfile}_{fc}.png"))
    factor_results_list[[fc]] <- list(homog_pval = fc_homog_pval,
                                      adonis_pval = fc_adonis_pval,
                                      adonis_R2 = fc_adonis_R2)
  }
  # plot_top_coefficients(adonis_result, meta, paste0("_",base_outfile,"_"))
  factor_results_list
}

add_row_to_summary <- function(result_list,
                               summary_file,
                               taxa_level,
                               seq_depth,
                               soil_depth="All",
                               season="All"){
  for (fc in names(result_list)){
    summary_file <- summary_file |>
      add_row(meta_factor = fc,
              taxa_level = taxa_level,
              sequence_depth = seq_depth,
              soil_depth = soil_depth,
              season = season,
              homogeneity_pval = result_list[[fc]]$homog_pval,
              permanova_pval = result_list[[fc]]$adonis_pval,
              permanova_R2 = result_list[[fc]]$adonis_R2)
  }
  summary_file
}



# Analysis functions ------------------------------------------------------

plot_nmds <- function(dist, meta, fc, save_image=FALSE, save_path, ...){
  nmds_obj <- metaMDS(dist)
  nmds_scores <- scores(nmds_obj) |>
    as.data.frame() |>
    rownames_to_column("sample")
  # Merge NMDS results with metadata
  meta_temp <- meta |> rownames_to_column("sample")
  nmds_scores <- inner_join(meta_temp, nmds_scores)
  centroid <- nmds_scores |> summarize(NMDS1 = mean(NMDS1),
                                       NMDS2 = mean(NMDS2),
                                       .by = fc)
  # Ellipse won't work if a group has < 4 items. Check for this
  small_groups <- any(table(nmds_scores[[fc]]) < 4)

  fc_eval <- sym(fc)
  nmds_plot <- plot_nmds_by_factor(
    df = nmds_scores,
    meta_factor = fc,
    dataset_name = "Timed",
    taxa_level = taxa_level,
    shape_factor = "Depth",
    polygon = small_groups,
    ...
  )

  if (!small_groups) {
    nmds_plot <- nmds_plot + stat_ellipse(aes(color = !!fc_eval), show.legend = FALSE)
  }
  nmds_plot <- nmds_plot +
    labs(subtitle = paste0(nmds_plot$labels$subtitle, "\nNMDS stress: ", signif(nmds_obj$stress, 4))) +
    # centroid
    geom_point(data=centroid, aes(fill = !!fc_eval), size = 5, shape=21, stroke = 1.5, alpha=0.8)
  browser()
  if (save_image){
    ggsave(plot = nmds_plot,
           filename = save_path,
           bg = "white",
           width = 10, height = 7.8)
  }
  nmds_plot
}

run_adonis <- function(otu, meta, fcs, outfile=NULL, factor_interaction=TRUE){
  formula_char <- ifelse(factor_interaction, "*", "+")
  formula <- as.formula(glue("otu ~ {str_flatten(fcs, collapse = formula_char)}"))
  permanova <- adonis(formula = formula,
                      data = meta,
                      permutations = n_permutations,
                      method = "bray")
  if (!is.null(outfile)){
    sink(outfile, append = TRUE)
    print(permanova$aov.tab)
    sink()
  }
  permanova
}

homogeneity_test <- function(dist, meta, fcs, outfile=NULL){
  test_results <- list()
  for (fc in fcs){
    if(n_distinct(meta[[fc]]) > 1) {
      bd <- anova(betadisper(dist, meta[[fc]]))
      test_results[[fc]] <- bd
      if(!is.null(outfile)){
        cat("\nANOVA check on homogeneity of", fc, "\n", file = outfile, append = T)
        sink(outfile, append = TRUE)
        print(bd)
        sink()
      }
    }
  }
  test_results
}

plot_top_coefficients <- function(adonis_result, meta, filestring){
  for (fc in rownames(adonis_result$coefficients)){
    # Trim digits off the end to see if factor matches those of metadata
    if (str_extract(fc, "(.*)(?=\\d+)") %in% colnames(meta)){
      factor_coef <- coefficients(adonis_result)[fc,]
      df <- data.frame(taxa = names(factor_coef),
                       coef = factor_coef,
                       abs_coef = abs(factor_coef)) |>
        mutate(taxa = str_extract(taxa,"(?<=;)[^;]*(?=;$)"))
      bplt <- slice_max(df, order_by = abs_coef, n = 20) |>
        ggplot(aes(x = reorder(taxa, coef), y = coef)) +
        geom_bar(stat = "identity") +
        coord_flip() +  # Horizontal bars
        labs(
          title = paste("Top 20 coef in factor:", fc),
          x = "Taxa",
          y = "Coefficient"
        ) +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8))
      ggsave(glue("{taxa_out}/top_coefs_{filestring}_{fc}.png"), plot = bplt, width = 6, height = 4)
    }
  }
}


