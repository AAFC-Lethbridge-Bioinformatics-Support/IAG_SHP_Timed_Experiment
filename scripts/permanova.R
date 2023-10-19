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
  mutate(taxon_lineage = str_extract(taxon_lineage, ";([^;]+);$", group = 1)) |>
  filter(sample != "S00JY-0597") |>
  filter(reads > sum(reads)*0.0001, .by = sample)

# Wide table for following analyses
taxa_3_wide <- taxa_2_filtered |>
  select(-c(file, taxon_id, percent)) |>
  pivot_wider(names_from = taxon_lineage, values_from = reads)

taxa_matrix <- taxa_3_wide |>
  column_to_rownames("sample") |>
  as.matrix()
taxa_matrix[is.na(taxa_matrix)] <- 0


# 2. Create components of phyloseq object ---------------------------------
taxa_otu_table <- otu_table(taxa_matrix, taxa_are_rows=FALSE)

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
  c("month_continuous","POS","Season","Depth", "Process")
} else {
  c("month_continuous","Season","Depth", "Process")
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
             pos=character(),
             homogeneity_pval=double(),
             permanova_pval=double(),
             permanova_R2=double())
}

results_summary <- analyses_wrap(pseq = pseq.rel,
                                 fcs = fcs,
                                 base_outfile = "allfactors") %>%
  # process and month are treated separately in later analysis calls
  keep(!names(.) %in% c("Process", "month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth)

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Depth=="0-15"),
                                 fcs = fcs[fcs!="Depth"],
                                 base_outfile = "allfactors_0-15") %>%
  keep(!names(.) %in% c("Process", "month_continuous", "Depth")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Depth=="15-30"),
                                 fcs = fcs[fcs!="Depth"],
                                 base_outfile = "allfactors_15-30") %>%
  keep(!names(.) %in% c("Process", "month_continuous", "Depth")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30")

# Process factor
# Month==0 is required for proper Process comparison. This leaves too few
# samples to model factor interactions (e.g. POS*Process)
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Month==0),
                                 fcs = fcs[fcs!="month_continuous"],
                                 base_outfile = "t0",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     note = "timepoint 0 samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Month==0&Depth=="0-15"),
                                 fcs = fcs[fcs!="month_continuous"&fcs!="Depth"],
                                 base_outfile = "t0_0-15",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     note = "timepoint 0 samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Month==0&Depth=="15-30"),
                                 fcs = fcs[fcs!="month_continuous"&fcs!="Depth"],
                                 base_outfile = "t0_15-30",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     note = "timepoint 0 samples only")


# Month factor
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"),
                                 fcs = fcs[fcs!="Process"],
                                 base_outfile = "ground") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="0-15"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                 base_outfile = "ground0-15") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground0-15spring") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     season = "Spring",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground0-15fall") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     season = "Fall",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                 base_outfile = "ground15-30") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     note = "ground samples only")

# split deep soils by POS, instead of season
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&POS=="M"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                 base_outfile = "ground15-30M") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     pos = "M",
                     note = "ground samples only")
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&POS=="L"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                 base_outfile = "ground15-30L") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     pos = "L",
                     note = "ground samples only")
results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&POS=="U"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                 base_outfile = "ground15-30U") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     pos = "U",
                     note = "ground samples only")


results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&Season=="Spring"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground15-30spring") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     season = "Spring",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq.rel, Process=="Ground"&Depth=="15-30"&Season=="Fall"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground15-30Fall") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     season = "Fall",
                     note = "ground samples only")

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
    subtitle <- ""
    if (fc %in% names(homog_test)){
      fc_homog_pval <- homog_test[[fc]]$`Pr(>F)`[1]
      subtitle <- glue("\n\nHomogeneity of variance p.val: {signif(fc_homog_pval, 4)}")
    } else {
      browser() # does this condition ever happen?
    }
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
  factor_results_list
}

add_row_to_summary <- function(result_list,
                               summary_file,
                               taxa_level,
                               seq_depth,
                               soil_depth="All",
                               season="All",
                               pos="All",
                               note=""){
  for (fc in names(result_list)){
    summary_file <- summary_file |>
      add_row(meta_factor = fc,
              taxa_level = taxa_level,
              sequence_depth = seq_depth,
              soil_depth = soil_depth,
              season = season,
              pos = pos,
              note = note,
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

  # Ellipse won't work if a group has < 4 items or if variable is continuous
  small_groups <- any(table(nmds_scores[[fc]]) < 4)
  polygon <- small_groups && !is.numeric(nmds_scores[[fc]])

  fc_eval <- sym(fc)
  nmds_plot <- plot_nmds_by_factor(
    df = nmds_scores,
    meta_factor = fc,
    dataset_name = "Timed",
    taxa_level = taxa_level,
    shape_factor = "Depth",
    polygon = polygon,
    ...
  )

  if (!is.numeric(nmds_scores[[fc]])) {
    if (!small_groups) { # ellipse
      nmds_plot <- nmds_plot + stat_ellipse(aes(color = !!fc_eval), show.legend = FALSE)
    }
    centroid <- nmds_scores |> summarize(NMDS1 = mean(NMDS1),
                                         NMDS2 = mean(NMDS2),
                                         .by = fc)
    nmds_plot <- nmds_plot +
      geom_point(data=centroid, aes(fill = !!fc_eval), size = 5, shape=21, stroke = 1.5, alpha=0.8)
  } else {
    my.envfit <- envfit(nmds_obj, meta, permutations = 9999)
    env.scores <- as.data.frame(scores(my.envfit, display = "vectors"))
    env.scores <- cbind(env.scores, env.variables = rownames(env.scores))
    env.scores <- cbind(env.scores, pval = my.envfit$vectors$pvals)

    scaled_arrow <- scale_arrow(end_x = env.scores$NMDS1,
                                end_y = env.scores$NMDS2,
                                minx = min(nmds_scores$NMDS1),
                                maxx = max(nmds_scores$NMDS1),
                                miny = min(nmds_scores$NMDS2),
                                maxy = max(nmds_scores$NMDS2))
    nmds_plot <- nmds_plot +
      geom_segment(aes(x = 0, xend = scaled_arrow["x"], y = 0, yend = scaled_arrow["y"]),
                   arrow = arrow(length = unit(0.25, "cm")), color = "grey10", lwd=0.3) +
      labs(subtitle = paste0(nmds_plot$labels$subtitle, "\nArrow envfit p.val: ", signif(env.scores$pval, 4)))
  }
  nmds_plot <- nmds_plot +
    labs(subtitle = paste0(nmds_plot$labels$subtitle, "\nNMDS stress: ", signif(nmds_obj$stress, 4)))

  if (save_image){
    ggsave(plot = nmds_plot,
           filename = save_path,
           bg = "white",
           width = 10, height = 7.8)
  }
  nmds_plot
}

scale_arrow <- function(end_x, end_y, minx, maxx, miny, maxy) {

  if (end_x < maxx && end_x > minx && end_y < maxy && end_y > miny){
    return(c(x = end_x, y = end_y))
  }

  scale_factor <- c(maxx/end_x,
                    minx/end_x,
                    maxy/end_y,
                    miny/end_y) |>
    abs() |>
    Filter(\(x) x <= 1, x=_) |>
    min()
  a_scaled <- c(x = end_x * scale_factor, y = end_y * scale_factor)
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

