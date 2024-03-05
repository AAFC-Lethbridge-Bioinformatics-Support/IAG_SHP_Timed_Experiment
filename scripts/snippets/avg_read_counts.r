# Get the basic read count stats of shallow and deep reads of kaiju output.
# Include only samples that are used in analysis (e.g. remove Blank samples)
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)
library(glue)

source("scripts/aux_functions.R")
tm <- get_timed_metadata()
shallow <- get_processed_taxonomy(sequence_depth = "shallow", taxa_level = "phylum",meta = tm)

raw_summary_shallow <- shallow$raw |> summarise(reads = sum(reads), .by = file) |>
  mutate(sample = str_extract(file, "S00JY-\\d{4}"),
         seq_depth = "shallow")

deep <- get_processed_taxonomy(sequence_depth = "deep", taxa_level = "phylum",meta = tm)

raw_summary_deep <- deep$raw |> summarise(reads = sum(reads), .by = file) |>
  mutate(sample = str_extract(file, "S00JY-\\d{4}"),
         seq_depth = "deep")


raw_summary <- full_join(raw_summary_shallow, raw_summary_deep) |>
  left_join(tm)

raw_summary_shallow |>
  left_join(tm) |>
  ggplot(aes(x = '', y = reads)) +
  geom_boxplot() +
  geom_jitter(color = "blue")

raw_summary_deep |>
  left_join(tm) |>
  ggplot(aes(x = '', y = reads)) +
  geom_boxplot() +
  geom_jitter(color = "blue")

raw_summary_deep$reads |> summary()
raw_summary_shallow$reads |> summary()
