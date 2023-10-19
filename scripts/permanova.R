# Load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)
library(glue)

source("scripts/aux_functions.R")
source("scripts/permanova_aux.R")

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

