# Mantel test and dissimilarity over time with linear model
#
# To investigate whether communities differ significantly over time, and to
# quantify the change, if linear model fits significantly
#
# Shallow and deep sequencing is analysed at the phylum and genus levels. Woltka
# functional profiling (KEGG KO counts) are also analysed.

library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)
library(glue)

source("scripts/aux_functions.R")
source("scripts/mantel_aux.R")

main_out <- glue("results/mantel")
timed_meta <- get_timed_metadata()

for (seq_depth in c("shallow", "deep")) {
  depth_out <- glue("{main_out}/{seq_depth}_sequencing/")
  for (taxa_level in c("phylum", "genus", "woltka")) {
    # woltka only applies to deep seqs
    if (seq_depth == "shallow" && taxa_level == "woltka") next

    taxa_out <- glue("{depth_out}/{taxa_level}/")
    lapply(c(main_out, depth_out, taxa_out), dir.create, recursive = TRUE, showWarnings = FALSE)

    subtitle <- if (taxa_level == "woltka"){
      "Functional profile (KO)"
    } else {
      glue("Taxa level: {taxa_level}")
    }
    counts_obj <- if(taxa_level == "woltka") {
      get_processed_KO(timed_meta)
    } else {
      get_processed_taxonomy(seq_depth, taxa_level, timed_meta)
    }
    pseq <- create_phyloseq(counts_obj, timed_meta)

    count_meta_table <- pseq$rel |> otu_table() |> as.data.frame() |>
      rownames_to_column("sample") |>
      inner_join(timed_meta, y = _) |>
      filter(Process == "Ground") |> column_to_rownames("sample")

    counts <- count_meta_table[,10:ncol(count_meta_table)]
    counts[is.na(counts)] <- 0


    for (dist_method in c("bray", "jaccard")){
      # Get distance matrices for taxa counts and Month
      # Rarefy counts via avgdist
      dist.abund <- vegdist(counts, method = dist_method)
      dist.month <- vegdist(count_meta_table["month_continuous"], method = "euclidean")

      mtl_all <- mantel(dist.abund, dist.month, method = "spearman", permutations = 9999, na.rm = T)
      capture.output(mtl_all,
                     file = glue("{taxa_out}/mantel_{dist_method}_v_time.txt"))

      # Corresponding linear model test
      combined_dists <- data.frame(abundance_dist = as.vector(dist.abund),
                                   month_dist = as.vector(dist.month))
      l_fit <- lm(abundance_dist ~ month_dist, data = combined_dists) |> summary()
      capture.output(l_fit,
                     file = glue("{taxa_out}/linear_fit_{dist_method}_v_time.txt"))

      subtitle2 <- glue("{subtitle}\nR2: {signif(l_fit$r.squared,4)} P.value: {signif(l_fit$coefficients[2,4],4)}")

      mm_lm <- make_distance_lm_plot(dist.abund, dist.month, dist_method,
                                     subtitle = subtitle2)
      ggsave(glue("{taxa_out}/linear_fit_{dist_method}_v_time.png"), plot = mm_lm)
    }
  }
}
