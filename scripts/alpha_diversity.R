# Alpha diversity testing and visualization for Month and Process factors
#
# Calculates alpha diversity using Inverse Simpson, Shannon, and Observed
# Richness measurements. The Process factor is categorical, is visualized with
# box plots, and tested with KW and Dunn tests. The Month factor is continuous,
# is visualized as scatter plots with linear fits, and tested via linear
# modelling.
#
# Shallow and deep sequencing is analysed at the phylum and genus levels. Woltka
# functional profiling (KEGG KO counts) are also analysed.

library(tidyverse)
library(vegan)
library(ggplot2)
library(glue)
library(dunn.test)

source("scripts/aux_functions.R")
source("scripts/alpha_diversity_aux.R")

main_out <- glue("results/alpha_diversity_unfiltered")
timed_meta <- get_timed_metadata()

# Analysis loops for taxa counts
for (seq_depth in c("shallow", "deep")) {
    depth_out <- glue("{main_out}/{seq_depth}_sequencing/")
  for (taxa_level in c("phylum", "genus")) {
    taxa_out <- glue("{depth_out}/{taxa_level}/")

    # Create output directories
    lapply(c(main_out, depth_out, taxa_out), dir.create, recursive = TRUE, showWarnings = FALSE)

    taxa_counts <- get_processed_taxonomy(seq_depth, taxa_level, timed_meta, min_filter = FALSE)
    pseq <- create_phyloseq(taxa_counts, timed_meta)

    make_all_alpha_plots(pseq = pseq$base,
                         outdir = taxa_out,
                         subtitle=glue("Taxonomic rank: {taxa_level}"))
  }
}

# Alpha div on functional profiles
#
# Phyloseq alpha-div estimate function doesn't work on fractional values
# (non-integer), which Woltka counts are due to normalization for gene length
# (RPK). Additionally, alpha diversity measurements should be taken with
# unfiltered count data, but a minimum filter was applied at the Woltka analysis

# woltka_out <- glue("{main_out}/woltka/")
# dir.create(woltka_out)
#
# ko_counts <- get_processed_KO(timed_meta)
# pseq <- create_phyloseq(ko_counts, timed_meta)
#
# make_all_alpha_plots(pseq$base,
#                      outdir = woltka_out,
#                      subtitle = "Functional profile (KO)")
