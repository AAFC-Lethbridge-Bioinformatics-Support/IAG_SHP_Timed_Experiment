# Differential abundance via ANCOM-BC2
#
# Identifies differentially abundant units (taxa or genes) across groups for
# categorical factors (Process) or directionally across a continuous variable.
# Applies analysis/plotting to specific subsets in accordance with PERMANOVA
# subsetting.
#
# Shallow and deep sequencing is analysed at the phylum and genus levels. Woltka
# functional profiling (KEGG KO counts) are also analysed. Deep sequencing (and
# thus woltka counts, which were produced from deep sequencing), has fewer
# samples and therefore cannot tolerate subsetting well due to low sample count.


library(ggplot2)
library(tidyverse)
library(dplyr)
library(vegan)
library(phyloseq)
library(ANCOMBC)
library(glue)

source("scripts/aux_functions.R")
source("scripts/ancombc2_aux.R")


# Run ANCOMBC2 analyses for month ------------------------------------------

main_out <- glue("results/ancombc")
# Rename since '_' is a functional character in ancombc result column names
timed_meta <- get_timed_metadata() |>
  rename(month.continuous = "month_continuous")

# Returns an appropriate plot subtitle based on the taxa level/woltka level
get_subtitle <- function(level = c("phylum", "genus", "woltka_KO", "woltka_pathways")) {
  level <- rlang::arg_match(level)
  if (level == "woltka_KO"){
    "Functional profile (KO)"
  } else if (level == "woltka_pathways") {
    "Functional profile (pathways)"
  } else {
    glue("Taxa level: {level}")
  }
}

# Provides an initial vector of metadata variables to be used in ANCOMBC2
# formula, depending on seq_depth and the main factor of interest
get_initial_factors <- function(seq_depth = c("shallow", "deep"),
                                main_factor = c("month", "process")) {
  main_factor <- rlang::arg_match(main_factor)
  if (main_factor == "month") {

    fcs <- c("month.continuous", "Depth", "Season")
  } else {
    fcs <- c("Process.character", "Depth", "Season")
  }
  if (seq_depth == "shallow") {
    # Shallow has multiple POS groups, deep doesn't
    fcs[4] <- "POS"
  }
  fcs
}

# Wrapper for all the month ancombc2 analysis calls
call_month_analyses <- function(pseq, factors, subtitle, taxa_out, seq_depth) {
  ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"),
                     factors = factors,
                     subtitle = subtitle,
                     outdir = taxa_out)

  # Shallow has more samples and can withstand subsetting
  if (seq_depth == "shallow") {
    factors <- c("month.continuous", "Season", "POS")
    # split by depth
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_0-15")
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_15-30")

    factors <- c("month.continuous", "POS")
    # Split shallow soil by season
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_0-15_Spring")
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_0-15_Fall")

    factors <- c("month.continuous", "Season")
    # Split deep soil by POS
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="U"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_15-30_POS-U")
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="M"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_15-30_POS-M")
    ancombc_wrap_month(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="L"),
                       factors = factors,
                       subtitle = subtitle,
                       outdir = taxa_out,
                       filename_append = "_15-30_POS-L")
  }
}

for (seq_depth in c("shallow", "deep")) {
  depth_out <- glue("{main_out}/{seq_depth}_sequencing/")
  dir.create(depth_out, recursive = TRUE)

  for (unit_level in c("phylum", "genus", "woltka_KO", "woltka_pathways")) {
    # Woltka is only for deep sequences
    if (seq_depth == "shallow" && str_detect(unit_level, "woltka")) next

    taxa_out <- glue("{depth_out}/{unit_level}/")
    dir.create(taxa_out)

    subtitle <- get_subtitle(unit_level)

    counts <- if (unit_level == "woltka_KO"){
      get_processed_KO(timed_meta)
    } else if (unit_level == "woltka_pathways") {
      get_processed_pathways(timed_meta)
    } else {
      get_processed_taxonomy(seq_depth, unit_level, timed_meta)
    }

    # Metavariable: Month ----
    # pseq <- create_phyloseq(counts, timed_meta)
    # month_factors <- get_initial_factors(seq_depth, "month")
    # call_month_analyses(pseq, month_factors, subtitle, taxa_out, seq_depth)

    # Metavariable: Process ----
    # Change Process to character because ANCOMBC2 doesn't work best with ordered factor
    timed_meta$Process.character <- as.character(timed_meta$Process)
    pseq <- create_phyloseq(counts, timed_meta)
    process_factors <- get_initial_factors(seq_depth, "process")

    ancombc_wrap_process(pseq = subset_samples(pseq$base, Month == 0),
                         factors = process_factors,
                         subtitle = subtitle,
                         outdir = taxa_out)

    # Add KO/pathway names to written csv result files
    if (unit_level == "woltka_KO") {
      csv_paths <- Sys.glob(glue("{taxa_out}/*.csv"))
      lapply(csv_paths, \(x) {
        csv <- read_csv(x)
        joined <- left_join(csv, counts$gene_key, by = join_by(taxon == ko_gene_ID))
        write_csv(joined, x)
      })
    }
    if (unit_level == "woltka_pathways") {
      csv_paths <- Sys.glob(glue("{taxa_out}/*.csv"))
      lapply(csv_paths, \(x) {
        csv <- read_csv(x)
        joined <- left_join(csv, counts$pathway_key, by = join_by(taxon == pathwayID))
        write_csv(joined, x)
      })
    }
  }
}

