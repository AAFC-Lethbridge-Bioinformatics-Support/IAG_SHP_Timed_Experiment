library(ggplot2)
library(tidyverse)
library(dplyr)
library(vegan)
library(phyloseq)
library(ANCOMBC)
library(glue)

source("scripts/aux_functions.R")
source("scripts/ancombc2_aux.R")


seq_depth <- "shallow" # shallow or deep
taxa_level <- "genus" # phylum or genus

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
    pseq <- create_phyloseq(counts, timed_meta)

    factors <- if (seq_depth == "shallow") {
      c("month.continuous", "Depth", "Season", "POS")
    } else {
      c("month.continuous", "Depth", "Season")
    }

    ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"),
                 formula = str_flatten(factors, collapse = "+"),
                 meta = timed_meta,
                 subtitle = subtitle,
                 outdir = taxa_out)

    # Shallow has more samples and can withstand subsetting
    if (seq_depth == "shallow") {
      factors <- c("month.continuous", "Season", "POS")
      # split by depth
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_0-15")
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_15-30")

      factors <- c("month.continuous", "POS")
      # Split shallow soil by season
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_0-15_Spring")
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_0-15_Fall")

      factors <- c("month.continuous", "Season")
      # Split deep soil by POS
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="U"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_15-30_POS-U")
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="M"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_15-30_POS-M")
      ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="L"),
                   formula = str_flatten(factors, collapse = "+"),
                   meta = timed_meta,
                   subtitle = subtitle,
                   outdir = taxa_out,
                   filename_append = "_15-30_POS-L")
    }
  }
}


