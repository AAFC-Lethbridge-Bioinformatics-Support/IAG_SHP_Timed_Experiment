# Load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)
library(glue)

source("scripts/aux_functions.R")
source("scripts/permanova_aux.R")


# 0.1 Define globals/output paths -----------------------------------------

seq_depth <- "shallow" # shallow or deep
taxa_level <- "phylum" # phylum or genus

n_permutations <- 9999

main_out <- glue("results/{seq_depth}")
main_out <- glue("{main_out}/permanova/")
taxa_out <- glue("{main_out}/{taxa_level}/")

# Create output directories
lapply(c(main_out, taxa_out), dir.create, recursive = TRUE, showWarnings = FALSE)


# 0.2 Read in data and metadata -------------------------------------------

taxa_counts <- get_processed_taxonomy(seq_depth, taxa_level)
timed_meta <- get_timed_metadata()
# Create phyloseq object
pseq <- create_phyloseq(taxa_counts, timed_meta)


# 0.3 Analysis prep -------------------------------------------------------

fcs <- if(seq_depth == "shallow") {
  c("month_continuous","Process","POS","Season","Depth")
} else if(seq_depth == "deep") {
  c("month_continuous","Process","Season","Depth")
} else { stop("Sequence depth should be specified as 'shallow' or 'deep'")}

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


# 1. Run analyses ---------------------------------------------------------

results_summary <- analyses_wrap(pseq = pseq$rel,
                                 fcs = fcs,
                                 base_outfile = "allfactors") %>%
  # process and month are treated separately in later analysis calls
  keep(!names(.) %in% c("Process", "month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth)

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Depth=="0-15"),
                                 fcs = fcs[fcs!="Depth"],
                                 base_outfile = "allfactors_0-15") %>%
  keep(!names(.) %in% c("Process", "month_continuous", "Depth")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Depth=="15-30"),
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
results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Month==0),
                                 fcs = fcs[fcs!="month_continuous"],
                                 base_outfile = "t0",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     note = "timepoint 0 samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Month==0&Depth=="0-15"),
                                 fcs = fcs[fcs!="month_continuous"&fcs!="Depth"],
                                 base_outfile = "t0_0-15",
                                 factor_interaction = FALSE) %>%
  keep(names(.) %in% c("Process")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     note = "timepoint 0 samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Month==0&Depth=="15-30"),
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
results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"),
                                 fcs = fcs[fcs!="Process"],
                                 base_outfile = "ground") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="0-15"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                 base_outfile = "ground0-15") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground0-15spring") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     season = "Spring",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground0-15fall") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "0-15",
                     season = "Fall",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                 base_outfile = "ground15-30") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     note = "ground samples only")

# split deep soils by POS, instead of season
results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&POS=="M"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                 base_outfile = "ground15-30M") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     pos = "M",
                     note = "ground samples only")
results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&POS=="L"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                 base_outfile = "ground15-30L") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     pos = "L",
                     note = "ground samples only")
results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&POS=="U"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                 base_outfile = "ground15-30U") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     pos = "U",
                     note = "ground samples only")


results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&Season=="Spring"),
                                 fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                 base_outfile = "ground15-30spring") %>%
  keep(names(.) %in% c("month_continuous")) |>
  add_row_to_summary(summary_file = results_summary,
                     taxa_level = taxa_level,
                     seq_depth = seq_depth,
                     soil_depth = "15-30",
                     season = "Spring",
                     note = "ground samples only")

results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&Season=="Fall"),
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

