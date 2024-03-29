library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(phyloseq)

color_palette <- c("#89C5DA", "#DA5724", "#689030", "#CE50CA", "#3F4921",
                   "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", "#D7C1B1",
                   "#D14285", "#AD6F3B", "#CD9BCD", "#74D944", "#6DDE88", "#652926", "#7FDCC0",
                   "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#38333E", "#599861")

# Read in Timed metadata and return a reformatted version
get_timed_metadata <- function(path = "../metadata/Metadata-IAG-Timed2-v3_2.csv") {
  message("Loading metadata at path ", path)

  timed_meta <- read_csv(path)
  # remove rows with all NA, remove blanks
  timed_meta <- timed_meta[rowSums(is.na(timed_meta)) != ncol(timed_meta),] |>
    rename(sample = "MBI_ID") |>
    filter(!str_detect(Site, "B\\d")) # Blank samples denoted by B#
  timed_meta$Month <- as.numeric(timed_meta$Month)
  timed_meta <- timed_meta |>
    mutate(Month = case_when(Year == "1998" | Year == "2007" ~ as.numeric(Year),
                             TRUE ~ Month))
  timed_meta$Month <- ordered(timed_meta$Month,
                              levels = c("0",
                                         "0.07",
                                         "0.5",
                                         "1",
                                         "3",
                                         "6",
                                         "12",
                                         "18",
                                         "1998",
                                         "2007"))

  timed_meta$month_continuous <- as.numeric(levels(timed_meta$Month))[timed_meta$Month]
  timed_meta$Year <- as.character(timed_meta$Year)
  timed_meta$Process <- ordered(timed_meta$Process,
                                levels = c("Fresh",
                                           "Sieved",
                                           "Air-dried",
                                           "Ground"))
  timed_meta
}

# Given the sequencing depth and taxa level, read taxonomic classification data,
# process and format the data, then return a named list with several formats:
# raw counts, filtered counts, wide filtered counts, and wide filtered counts as
# a matrix (NAs set to 0)
get_processed_taxonomy <- function(sequence_depth, taxa_level, meta, min_filter=TRUE){
  taxa_long <- if(sequence_depth == "shallow"){
    read_tsv(glue("data/kaiju_shallow/kaiju_{taxa_level}_summary.tsv"))
  } else if(sequence_depth == "deep"){
    read_tsv(glue("data/kaiju_deep/kaiju_{taxa_level}_summary.tsv"))
  } else {
    stop("Error: sequence depth should be specified as 'shallow' or 'deep'")
  }

  # Reformat and filter some columns.
  # Remove pre-2020 samples, filter unclassified taxa, remove suspect sample
  # 0597
  taxa_filtered_samples <- taxa_long |>
    mutate(sample = str_extract(file, "S00JY-\\d{4}")) |>
    filter(sample != "S00JY-0597") |>
    filter(sample %in% filter(meta, Year == 2020)$sample)

  taxa_1_filtered <- taxa_filtered_samples |>
    filter(!str_detect(taxon_name, "unclassified$|^cannot|Viruses|[ ;]bacterium[; ]"))

    # Shorten lineage to the name at the relevant taxa level
  if (taxa_level == "phylum") {
    taxa_1_filtered <- taxa_1_filtered |>
      mutate(taxon_name = str_extract(taxon_name, ";([^;]+);$", group = 1))
  } else if (taxa_level == "genus") {
    # Handle clashes where the same name is present >once despite distinct
    # lineage (much more likely for genera)
    clashing_taxa <- taxa_1_filtered |>
      distinct(taxon_name) |>
      mutate(taxon_name = str_extract(taxon_name, ";([^;]+);$", group = 1)) |>
      count(taxon_name) |> filter(n > 1) |> pull(taxon_name)
    taxa_1_filtered <- taxa_1_filtered |>
      mutate(lineage_tokens = taxon_name |> str_replace(";$", "") |> str_split(";"),
             final_2_tokens = lapply(lineage_tokens, tail, 2),
             taxon_name = ifelse(lapply(final_2_tokens, last) %in% clashing_taxa,
                                 str_c(lapply(final_2_tokens, last), " (", lapply(final_2_tokens, first), ")"),
                                 lapply(final_2_tokens, last)) |> unlist()) |>
      select(-c(lineage_tokens, final_2_tokens))
  }

  taxa_1_filtered <- taxa_1_filtered |>
    mutate(classified_percent = reads/sum(reads), .by = sample)

  # apply 0.01% per-sample threshold
  if (min_filter) {
    taxa_1_filtered <- taxa_1_filtered |>
      filter(reads > sum(reads)*0.0001, .by = sample)
  }

  taxa_2_wide <- taxa_1_filtered |>
    select(-c(file, taxon_id, percent, classified_percent)) |>
    pivot_wider(names_from = taxon_name, values_from = reads)

  taxa_matrix <- taxa_2_wide |>
    column_to_rownames("sample") |>
    as.matrix()
  taxa_matrix[is.na(taxa_matrix)] <- 0

  list(raw = taxa_filtered_samples,
       filtered = taxa_1_filtered,
       wide = taxa_2_wide,
       matrix = taxa_matrix)
}

# Read the woltka KO data, process and format it, then return a named list with
# several formats: long (tidy) counts, wide (genes as columns), and wide matrix
get_processed_KO <- function(meta) {
  ko <- read_tsv("data/woltka_results/ko_filtered-0-01p.tsv") |>
    relocate("Name") |>
    rename("name" = Name,
           "ko_gene_ID" = `#FeatureID`)

  ko_long <- ko |>
    pivot_longer(cols = starts_with("S00JY"),
                 names_to = "sample",
                 values_to = "reads")

  # FILTER out the odd sample and keep only 2020 samples
  ko_long_filtered <- ko_long |>
    filter(sample != "S00JY-0597" & sample %in% filter(meta, Year == 2020)$sample)

  gene_names <- ko_long_filtered |> distinct(name, ko_gene_ID)

  # Pivot for samples as rows, genes as columns
  ko_wide <- ko_long_filtered |>
    pivot_wider(id_cols = sample,
                names_from = ko_gene_ID,
                values_from = reads)

  ko_matrix <- ko_wide |>
    column_to_rownames("sample") |>
    as.matrix()
  ko_matrix[is.na(ko_matrix)] <- 0

  list(long = ko_long_filtered,
       wide = ko_wide,
       matrix = ko_matrix,
       gene_key = gene_names)

}

# Read the woltka pathway data, process and format it, then return a named list
# with several formats: long (tidy) counts, wide (pathways as columns), and wide
# matrix
get_processed_pathways <- function(meta) {

  # Names from woltka workflow kegg_query.py
  pw_names <- read_tsv("data/woltka_results/pathway_name.txt",
                       col_names = FALSE)
  colnames(pw_names) <- c("pathway_ID","pathway_name")
  # Pathway ID to pathway class/groups from woltka workflow kegg_query.py
  pw_groups <- read_tsv("data/woltka_results/pathway-to-class.txt",
                        col_names = FALSE)
  colnames(pw_groups) <- c("pathway_ID","pathway_class", "pathway_group")
  # pathway-to-class.txt misses global/overiew group, so add this manually
  pw_global_overview_group <- read_tsv("data/woltka_results/global_and_overview_pathways.txt")
  pw_groups <- full_join(pw_groups, pw_global_overview_group)

  pw_names_full <- left_join(pw_names, pw_groups)

  pw_counts <- read_tsv("data/woltka_results/pathways_filtered-0-01p.tsv") |>
    rename("pathway_ID" = `#FeatureID`)

  # add pathway names and groups
  pw_counts <- pw_counts |>
    left_join(pw_names_full)

  # Pivot longer and get relative abundance
  pw_long <- pw_counts |>
    pivot_longer(cols = starts_with("S00JY"),
                 names_to = "sample",
                 values_to = "reads")

  # FILTER out the odd sample and keep only 2020 samples
  pw_long_filtered <- pw_long |>
    filter(sample != "S00JY-0597" & sample %in% filter(meta, Year == 2020)$sample)

  pathway_key <- pw_long_filtered |> distinct(pathway_ID, pathway_name, pathway_group, pathway_class)

  # Pivot for samples as rows
  pw_wide <- pw_long_filtered |>
    pivot_wider(id_cols = sample,
                names_from = pathway_ID,
                values_from = reads)

  pw_matrix <- pw_wide |>
    column_to_rownames("sample") |>
    as.matrix()
  pw_matrix[is.na(pw_matrix)] <- 0

  list(long = pw_long_filtered,
       wide = pw_wide,
       matrix = pw_matrix,
       pathway_key = pathway_key)
}

# Given the taxa data outputted from get_processed_taxonomy() and the metadata
# from get_timed_metadata(), return a base and relative abundance phyloseq
# object
create_phyloseq <- function(taxa_counts, metadata){
  taxa_otu_table <- otu_table(taxa_counts$matrix, taxa_are_rows=FALSE)

  proj_sample_data <- metadata |>
    filter(sample %in% taxa_counts$wide$sample) |>
    column_to_rownames("sample") |>
    sample_data()

  # Create phyloseq object
  pseq_base <- phyloseq(taxa_otu_table, proj_sample_data)
  # Transform for relative abundances (compositional)
  pseq_rel <- transform_sample_counts(pseq_base, function(x) x / sum(x))

  list(base = pseq_base,
       rel = pseq_rel)
}
