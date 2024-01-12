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

for (seq_depth in c("shallow", "deep")) {
  depth_out <- glue("{main_out}/{seq_depth}_sequencing/")
  dir.create(depth_out, recursive = TRUE)

  for (taxa_level in c("phylum", "genus", "woltka")) {
    if (seq_depth == "shallow" && taxa_level == "woltka") next
    taxa_out <- glue("{depth_out}/{taxa_level}/")
    dir.create(taxa_out)

    subtitle <- if (taxa_level == "woltka"){
      "Functional profile (KO)"
    } else {
      glue("Taxa level: {taxa_level}")
    }

    counts <- if (taxa_level == "woltka") {
      get_processed_KO(timed_meta)
    } else {
      get_processed_taxonomy(seq_depth, taxa_level, timed_meta)
    }

    # Metavariable: Month ----
    pseq <- create_phyloseq(counts, timed_meta)

    factors <- if (seq_depth == "shallow") {
      c("month.continuous", "Depth", "Season", "POS")
    } else {
      c("month.continuous", "Depth", "Season")
    }

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

    # Metavariable: Process ----
    # Change Process to character because ANCOMBC2 doesn't work best with ordered factor
    timed_meta$Process.character <- as.character(timed_meta$Process)
    pseq <- create_phyloseq(counts, timed_meta)

    factors <- if (seq_depth == "shallow") {
      c("Process.character", "Depth", "Season", "POS")
    } else {
      c("Process.character", "Depth", "Season")
    }
    ancombc_wrap_process(pseq = subset_samples(pseq$base, Month == 0),
                         factors = factors,
                         subtitle = subtitle,
                         outdir = taxa_out)

    if (taxa_level == "woltka") {
      # Add KO names to written csv result files
      csv_paths <- Sys.glob(glue("{taxa_out}/*.csv"))
      lapply(csv_paths, \(x) {
        csv <- read_csv(x)
        joined <- left_join(csv, counts$gene_key, by = join_by(taxon == ko_gene_ID))
        write_csv(joined, x)
      })
    }
  }
}
