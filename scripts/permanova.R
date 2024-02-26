# Beta diversity testing and visualization via PERMANOVA and NMDS
#
# Calculates Bray-Curtis distances on count data and uses those values to test
# for statistically separate groups within metadata factors. PERMANOVA models
# how much each included factor contributes to the data's variance. These values
# were used to determine which factors to subset by to investigate the
# contribution of Month, which has a weak R2 but significant pvalue.
#
# Shallow and deep sequencing is analysed at the phylum and genus levels. Woltka
# functional profiling (KEGG KO counts) are also analysed. Deep sequencing (and
# thus woltka counts, which were produced from deep sequencing), has fewer
# samples and therefore cannot tolerate subsetting well due to low sample count.

library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)
library(glue)

source("scripts/aux_functions.R")
source("scripts/permanova_aux.R")

# Setup -------------------------------------------------------------------
main_out <- glue("results/permanova")
timed_meta <- get_timed_metadata()

summary_filename <- glue("{main_out}/permanova_results_summary.csv")
results_summary <- data.frame(meta_factor=character(),
                              unit_level=character(),
                              sequence_depth=character(),
                              subset=character(),
                              homogeneity_pval=double(),
                              permanova_pval=double(),
                              permanova_R2=double())


# 1. Run analyses ---------------------------------------------------------
for (seq_depth in c("shallow", "deep")) {
  depth_out <- glue("{main_out}/{seq_depth}_sequencing/")
  dir.create(depth_out, recursive = TRUE)

  fcs <- if(seq_depth == "shallow") {
    c("month_continuous","Process","POS","Season","Depth")
  } else if(seq_depth == "deep") {
    # Deep sequencing has only POS=="M"
    c("month_continuous","Process","Season","Depth")
  } else { stop("Sequence depth should be specified as 'shallow' or 'deep'")}

  for (unit_level in c("phylum", "genus", "woltka_KO", "woltka_pathways")) {
    # Woltka is only for deep sequences
    if (seq_depth == "shallow" &&  str_detect(unit_level, "woltka")) next

    taxa_out <- glue("{depth_out}/{unit_level}/")
    dir.create(taxa_out)

    counts <- if (unit_level == "woltka_KO") {
      get_processed_KO(timed_meta)
    } else if (unit_level == "woltka_pathways") {
      get_processed_pathways(timed_meta)
    } else {
      get_processed_taxonomy(seq_depth, unit_level, timed_meta)
    }
    pseq <- create_phyloseq(counts, timed_meta)

    # All factors (no filters)
    results_summary <- analyses_wrap(pseq = pseq$rel,
                                     fcs = fcs,
                                     outdir = taxa_out,
                                     unit_level = unit_level,
                                     outfile_extension = "all") %>%
      # exclude results of key factors because they require specific filters
      keep(!names(.) %in% c("Process", "month_continuous")) |>
      add_row_to_summary(summary_file = results_summary,
                         unit_level = unit_level,
                         seq_depth = seq_depth,
                         subset = "all")
    # Process factor
    # Month==0 is required for proper Process comparison across processes. This
    # leaves too few samples to model factor interactions (e.g. POS*Process)
    results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Month==0),
                                     fcs = fcs[fcs!="month_continuous"],
                                     outdir = taxa_out,
                                     unit_level = unit_level,
                                     outfile_extension = "t0",
                                     factor_interaction = FALSE) |>
      add_row_to_summary(summary_file = results_summary,
                         unit_level = unit_level,
                         seq_depth = seq_depth,
                         subset = "t0")

    # Month factor
    # Process=="Ground" is required for proper comparison across time
    results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"),
                                     fcs = fcs[fcs!="Process"],
                                     outdir = taxa_out,
                                     unit_level = unit_level,
                                     outfile_extension = "ground") |>
      add_row_to_summary(summary_file = results_summary,
                         unit_level = unit_level,
                         seq_depth = seq_depth,
                         subset = "ground")

    # Shallow has more samples and can withstand subsetting
    if (seq_depth == "shallow") {
      # Month factor - subset by soil depth
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="0-15"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground0-15") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground0-15")
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground15-30") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground15-30")

      # Subset shallow soils by Season
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="0-15"&Season=="Spring"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground0-15spring") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground0-15spring")
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="0-15"&Season=="Fall"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="Season"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground0-15fall") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground0-15fall")

      # Subset deep soils by POS
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&POS=="M"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground15-30M") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground15-30M")
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&POS=="L"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground15-30L") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground15-30L")
      results_summary <- analyses_wrap(pseq = subset_samples(pseq$rel, Process=="Ground"&Depth=="15-30"&POS=="U"),
                                       fcs = fcs[fcs!="Process"&fcs!="Depth"&fcs!="POS"],
                                       outdir = taxa_out,
                                       unit_level = unit_level,
                                       outfile_extension = "ground15-30U") |>
        add_row_to_summary(summary_file = results_summary,
                           unit_level = unit_level,
                           seq_depth = seq_depth,
                           subset = "ground15-30U")
    }
  }
}

write_csv(results_summary, file = glue("{main_out}/permanova_results_summary.csv"))

