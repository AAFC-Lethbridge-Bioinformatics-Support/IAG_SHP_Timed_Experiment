library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(indicspecies)

source("scripts/aux_functions.R")

# 0.1 Define globals/output paths ----------------------------------------------
taxa_level <- "phylum"

main_out <- "results/combined_seq_depth_NMDS"


create_dir_if_nonexistant <- function(path) {
  ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
}

lapply(c(main_out), create_dir_if_nonexistant)


# ---- 0.2 Read in and tidy timed metadata -------------------------------------
timed_meta <- read_csv("../metadata/Metadata-IAG-Timed2-v3_2.csv")
# remove all NA rows, remove blanks
timed_meta <- timed_meta[rowSums(is.na(timed_meta)) != ncol(timed_meta),] |>
  rename(sample = `MBI_ID`) |>
  filter(!str_detect(Site, "B\\d"))

timed_meta %<>% mutate(Month = case_when(Year == "1998" | Year == "2007" ~ Year,
                                         TRUE ~ Month))
timed_meta$Month <- ordered(timed_meta$Month, levels = c("0",
                                                         "0.07",
                                                         "0.5",
                                                         "1",
                                                         "3",
                                                         "6",
                                                         "12",
                                                         "18",
                                                         "1998",
                                                         "2007"))
timed_meta$Year <- as.character(timed_meta$Year)


# ---- 1. Read in and filter taxonomy data -------------------------------------
files <- c(glue("data/deep_Timed/kaiju_{taxa_level}_summary.tsv"),
           glue("../pipeline_outputs/shallow_719_all_datasets/kaiju_{taxa_level}_summary.tsv"))
taxa_long <- read_tsv(files, id = "sequencing_depth")


# reformat some column aspects
taxa_long <- taxa_long |>
  mutate(sample = str_extract(file, "S00JY-\\d{4}"),
         sequencing_depth = ifelse(str_detect(sequencing_depth, "deep_Timed"),
                                   "deep",
                                   "shallow")
         ) |>
  rename("taxon_lineage" = taxon_name) |>
  mutate(taxon_lineage = str_replace(taxon_lineage, "cellular organisms;", ""))

# filter from this set if you want only the 42 samples with deep sequencing
deep_sample_list <- taxa_long |> filter(sequencing_depth == "deep") |>
  distinct(sample)

taxa_long_filtered <- taxa_long |>
  filter(sample %in% timed_meta$sample) |>
  filter(!str_detect(taxon_lineage, "unclassified$|^cannot|Viruses|[ ;]bacterium[; ]")) |>
  mutate(sample_threshold = sum(reads)*0.0001, .by = c(sequencing_depth, sample)) |>
  filter(reads > sample_threshold)

sample_thresholds <- taxa_long_filtered |>
  distinct(sequencing_depth, sample, sample_threshold)

taxa_long_filtered <- select(taxa_long_filtered, -sample_threshold)

# recalculate percentage so that it adds to 100 after removing unclassified and
# applying the threshold
taxa_longf_total <- taxa_long_filtered |>
  summarize(total_reads = sum(reads), .by = c(sequencing_depth, sample)) |>
  merge(taxa_long_filtered, y=_) |>
  mutate(recalculated_percent = reads/total_reads*100)


# 2. Perform NMDS, ANOSIM ------------------------------------------------------
anosim_result_file <- glue("{main_out}/anosim_combined_seqdepth_all_samples.txt")
cat(glue("ANOSIM Results for taxa level {str_to_upper(taxa_level)}\n\n\n"),
    file = anosim_result_file,
    append = T)


filtered_meta <- timed_meta

# filter for the appropriate subset of samples
taxa_long_foi <- taxa_longf_total |>
  filter(sample %in% filtered_meta$sample)

# Pivot and make proportional for proper ANOSIM format
taxa_wide_foi <- taxa_longf_total |>
  select(-c(file, taxon_id, total_reads, percent, reads)) |>
  pivot_wider(names_from = taxon_lineage, values_from = recalculated_percent)

# a meta table with the same sample order used for the abundance matrix
foi_metadata <- merge(taxa_wide_foi$sample, filtered_meta,
                      by.x = "x", by.y = "sample") |>
  rename(sample = x)

rel_abund_matrix_foi <- select(taxa_wide_foi, -c(sequencing_depth,sample)) |>
  as.matrix()
# Replace NAs with 0
rel_abund_matrix_foi[is.na(rel_abund_matrix_foi)] <- 0

set.seed(1234)
ano <- anosim(rel_abund_matrix_foi, taxa_wide_foi$sequencing_depth, distance = "bray", permutations = 999)

# Write result to file
cat(glue("FACTOR: sequencing depth\n"), file = anosim_result_file, append = T)
utils::capture.output(
  ano,
  file = anosim_result_file, append = T)
subtitle <- glue("   ANOSIM p-val: {ano$signif} R-stat: {ano$statistic}")

## Perform NMDS ----

nmds = metaMDS(rel_abund_matrix_foi, distance = "bray")
plot(nmds)

data.scores = as.data.frame(scores(nmds)$sites)
data.scores$sample <- taxa_wide_foi$sample
data.scores$sequencing_depth <- taxa_wide_foi$sequencing_depth

## Merge NMDS results with metadata ----
data.scores_meta <- merge(data.scores, foi_metadata) |> distinct()
plot_nmds_by_factor(
  df = data.scores_meta,
  meta_factor = "sequencing_depth",
  dataset_name = "Timed",
  taxa_level = taxa_level,
  shape_factor = "Depth",
  polygon = TRUE,
  save_image = TRUE,
  save_path = glue("{main_out}/sequencing_depth_{taxa_level}_all_samples_nmds.png"),
  subtitle_append = subtitle
)


# Procrustes Analysis -----------------------------------------------------

# Only run this if you've filtered the shallow sequencing to contain only the 42
# samples that have deep sequencing (i.e. 1 to 1 sample correspondance)


## Shallow vs Deep ---------------------------------------------------------

# SHALLOW
taxa_wide_shallow <- taxa_wide_foi |>
  filter(sample %in% deep_sample_list$sample) |>
  filter(sequencing_depth == "shallow")
rel_abund_matrix_shallow <- select(taxa_wide_shallow, -c(sequencing_depth,sample)) |>
  as.matrix()
# Replace NAs with 0
rel_abund_matrix_shallow[is.na(rel_abund_matrix_shallow)] <- 0

nmds_shallow = metaMDS(rel_abund_matrix_shallow, distance = "bray")
plot(nmds_shallow)

data.scores_shallow = as.data.frame(scores(nmds_shallow)$sites)
data.scores_shallow$sample <- taxa_wide_shallow$sample
data.scores_shallow$sequencing_depth <- taxa_wide_shallow$sequencing_depth

# DEEP
taxa_wide_deep <- taxa_wide_foi |> filter(sequencing_depth == "deep")
rel_abund_matrix_deep <- select(taxa_wide_deep, -c(sequencing_depth,sample)) |>
  as.matrix()
# Replace NAs with 0
rel_abund_matrix_deep[is.na(rel_abund_matrix_deep)] <- 0

nmds_deep = metaMDS(rel_abund_matrix_deep, distance = "bray")
plot(nmds_deep)

data.scores_deep = as.data.frame(scores(nmds_deep)$sites)
data.scores_deep$sample <- taxa_wide_deep$sample
data.scores_deep$sequencing_depth <- taxa_wide_deep$sequencing_depth


pcs <- protest(nmds_shallow,nmds_deep)

procrustes_test_file <- glue("{main_out}/procrustes_shallow_vs_deep_{taxa_level}.txt")
cat(glue("Procrustes test (vegan::protest) Results for taxa level {str_to_upper(taxa_level)}\n\n\n"),
    file = procrustes_test_file, append = T)
sink(procrustes_test_file, append = T)
pcs
sink()

png(glue("{main_out}/procrustes_shallow_vs_deep_{taxa_level}.png"))
plot(pcs,
     main = "Procrustes errors - shallow vs deep sequencing")
dev.off()

# Phyloflash procrustes ---------------------------------------------------


# Phyloflash
pf_long = read.csv("data/deep_Timed/timed_phyloflash_NTUabundances_concatenated.csv", header = FALSE)

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

#NMDS plot
pf_wide <- pf_longf_total |>
  pivot_wider(id_cols = sample,
              names_from = taxa_lineage,
              values_from = rel_abund)

pf_wide[is.na(pf_wide)] <- 0

nmds_pf <-  pf_wide |> column_to_rownames(var = "sample") |> as.matrix() |>
  metaMDS(distance = "bray")
data.scores_pf <- as.data.frame(scores(nmds)$sites) |>
  rownames_to_column("sample")


pcs <- protest(nmds_deep,nmds_pf)
pcs
png(glue("{main_out}/procrustes_kaiju_vs_phyloflash_{taxa_level}.png"))
plot(pcs,
     main = "Procrustes errors - kaiju vs phyloflash")
dev.off()
plot(pcs)
