# For PERMANOVA and NMDS to compare shallow and deep sequencing depths
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)
library(glue)

source("scripts/aux_functions.R")
source("scripts/permanova_aux.R")

# 0.1 Define globals/output paths ----------------------------------------------

main_out <- "results/combined_seq_depth_NMDS2"
dir.create(main_out)


# 0.2 Read in metadata ----------------------------------------------------
timed_meta <- get_timed_metadata()


# ---- 1. Read in and filter taxonomy data -------------------------------------
for (taxa_level in c("phylum", "genus")) {
  shallow_taxa <- get_processed_taxonomy("shallow", taxa_level, timed_meta)
  shallow_to_merge <- shallow_taxa$filtered |>
    mutate("sequencing_depth" = "shallow",
           sample = str_c(sample, "_shallow"))
  deep_taxa <- get_processed_taxonomy("deep", taxa_level, timed_meta)
  deep_to_merge <- deep_taxa$filtered |>
    mutate("sequencing_depth" = "deep",
           sample = str_c(sample, "_deep"))
  taxa <- bind_rows(shallow_to_merge, deep_to_merge)

  # Build a metadata table with separate entries for shallow and deep samples,
  # even though they have all the same metadata attributes. This allows the
  # workflow to distinguish shallow/deep samples with the same ID
  shallow_meta <- timed_meta |>
    mutate("sequencing_depth" = "shallow",
           sample = str_c(sample, "_shallow")) |>
    filter(sample %in% shallow_to_merge$sample)
  deep_meta <- timed_meta |>
    mutate("sequencing_depth" = "deep",
           sample = str_c(sample, "_deep")) |>
    filter(sample %in% deep_to_merge$sample)
  meta <- bind_rows(shallow_meta, deep_meta)


  taxa_wide <- taxa |>
    select(-c(file, taxon_id, percent)) |>
    pivot_wider(names_from = taxon_name, values_from = reads)

  taxa_matrix <- select(taxa_wide, -c(sequencing_depth)) |>
    column_to_rownames("sample") |>
    as.matrix()
  taxa_matrix[is.na(taxa_matrix)] <- 0

  # Build phyloseq obj
  # format needed for create_phyloseq function
  count_holder <- list()
  count_holder$wide <- taxa_wide
  count_holder$matrix <- taxa_matrix
  pseq <- create_phyloseq(count_holder, meta)

  analysis_res <- analyses_wrap(pseq$rel,
                                c("sequencing_depth"),
                                main_out,outfile_extension = taxa_level,
                                taxa_level)
}
