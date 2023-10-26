library(ggplot2)
library(tidyverse)
library(dplyr)
library(vegan)
library(phyloseq)
library(ANCOMBC)

source("scripts/aux_functions.R")
source("scripts/ancombc2_aux.R")

# 0.1 Define globals/output paths -----------------------------------------

seq_depth <- "shallow" # shallow or deep
taxa_level <- "genus" # phylum or genus

main_out <- glue("results/{seq_depth}")
main_out <- glue("{main_out}/ancombc/")
taxa_out <- glue("{main_out}/{taxa_level}/")

# Create output directories
lapply(c(main_out, taxa_out), dir.create, recursive = TRUE, showWarnings = FALSE)


# 0.2 Read in data and metadata -------------------------------------------

# Rename since '_' is a functional character in ancombc result column names
timed_meta <- get_timed_metadata() |>
  rename(month.continuous = "month_continuous")
taxa_counts <- get_processed_taxonomy(seq_depth, taxa_level, timed_meta)
# Create phyloseq object
pseq <- create_phyloseq(taxa_counts, timed_meta)


# Run ANCOMBC2 analyses for month------------------------------------------

ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"),
             formula = str_flatten(c("month.continuous", "Depth", "Season", "POS"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out)

# split by depth
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"),
             formula = str_flatten(c("month.continuous", "Season", "POS"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_0-15")
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"),
             formula = str_flatten(c("month.continuous", "Season", "POS"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_15-30")

# Split shallow depth by season
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
             formula = str_flatten(c("month.continuous", "POS"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_0-15_Spring")
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
             formula = str_flatten(c("month.continuous", "POS"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_0-15_Fall")

# Split deeper depth by season
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&Season=="Spring"),
             formula = str_flatten(c("month.continuous", "POS"), collapse = "+"),
             taxa_level = taxa_level,
             filename_append = "_15-30_Spring")
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&Season=="Fall"),
             formula = str_flatten(c("month.continuous", "POS"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_15-30_Fall")


# Split deeper depth by POS
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="U"),
             formula = str_flatten(c("month.continuous", "Season"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_15-30_POS-U")
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="M"),
             formula = str_flatten(c("month.continuous", "Season"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_15-30_POS-M")
ancombc_wrap(pseq = subset_samples(pseq$base, Process=="Ground"&Depth=="15-30"&POS=="L"),
             formula = str_flatten(c("month.continuous", "Season"), collapse = "+"),
             meta = timed_meta,
             taxa_level = taxa_level,
             taxa_out = taxa_out,
             filename_append = "_15-30_POS-L")

