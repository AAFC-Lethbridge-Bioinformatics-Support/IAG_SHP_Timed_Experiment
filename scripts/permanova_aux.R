# This script contains functions used solely by the permanova.R script

# Load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)

analyses_wrap <- function(pseq,
                          fcs,
                          base_outfile,
                          factor_interaction = TRUE) {
  factor_results_list <- list()
  PERMANOVA_result_file <- glue("{taxa_out}/PERMANOVA_{base_outfile}_{taxa_level}.txt")
  cat(glue("PERMANOVA via vegan::adonis2 for taxa level {str_to_upper(taxa_level)}\n\n\n"),
      file = PERMANOVA_result_file)

  otu <- otu_table(pseq)
  meta <- as(sample_data(pseq), "data.frame")
  rownames(meta) <- sample_names(pseq)
  dist <- vegdist(otu)

  homog_test <- homogeneity_test(dist, meta, fcs, PERMANOVA_result_file)
  adonis_result <- run_adonis(otu, meta, fcs, PERMANOVA_result_file, factor_interaction)

  for (fc in fcs) {
    subtitle <- ""
    if (fc %in% names(homog_test)){
      fc_homog_pval <- homog_test[[fc]]$`Pr(>F)`[1]
      subtitle <- glue("\n\nHomogeneity of variance p.val: {signif(fc_homog_pval, 4)}")
    } else {
      browser() # does this condition ever happen?
    }
    fc_adonis <- adonis_result$aov.tab[fc,]
    fc_adonis_R2 <- fc_adonis$R2[1]
    fc_adonis_pval <- fc_adonis$`Pr(>F)`[1]
    subtitle <- glue("{subtitle}\nAdonis PERMANOVA test p.val: {signif(fc_adonis_pval, 4)}  R2: {signif(fc_adonis_R2, 4)}")
    plot_nmds(dist, meta, fc,
              subtitle=subtitle,
              save_image=TRUE,
              save_path=glue("{taxa_out}/nmds_{base_outfile}_{fc}.png"))
    factor_results_list[[fc]] <- list(homog_pval = fc_homog_pval,
                                      adonis_pval = fc_adonis_pval,
                                      adonis_R2 = fc_adonis_R2)
  }
  factor_results_list
}

add_row_to_summary <- function(result_list,
                               summary_file,
                               taxa_level,
                               seq_depth,
                               soil_depth="All",
                               season="All",
                               pos="All",
                               note=""){
  for (fc in names(result_list)){
    summary_file <- summary_file |>
      add_row(meta_factor = fc,
              taxa_level = taxa_level,
              sequence_depth = seq_depth,
              soil_depth = soil_depth,
              season = season,
              pos = pos,
              note = note,
              homogeneity_pval = result_list[[fc]]$homog_pval,
              permanova_pval = result_list[[fc]]$adonis_pval,
              permanova_R2 = result_list[[fc]]$adonis_R2)
  }
  summary_file
}



# Analysis functions ------------------------------------------------------

plot_nmds <- function(dist, meta, fc, save_image=FALSE, save_path, ...){
  nmds_obj <- metaMDS(dist)
  nmds_scores <- scores(nmds_obj) |>
    as.data.frame() |>
    rownames_to_column("sample")
  # Merge NMDS results with metadata
  meta_temp <- meta |> rownames_to_column("sample")
  nmds_scores <- inner_join(meta_temp, nmds_scores)

  # Ellipse won't work if a group has < 4 items or if variable is continuous
  small_groups <- any(table(nmds_scores[[fc]]) < 4)
  polygon <- small_groups && !is.numeric(nmds_scores[[fc]])

  fc_eval <- sym(fc)
  nmds_plot <- plot_nmds_by_factor(
    df = nmds_scores,
    meta_factor = fc,
    dataset_name = "Timed",
    taxa_level = taxa_level,
    shape_factor = "Depth",
    polygon = polygon,
    ...
  )

  if (!is.numeric(nmds_scores[[fc]])) {
    if (!small_groups) { # ellipse
      nmds_plot <- nmds_plot + stat_ellipse(aes(color = !!fc_eval), show.legend = FALSE)
    }
    centroid <- nmds_scores |> summarize(NMDS1 = mean(NMDS1),
                                         NMDS2 = mean(NMDS2),
                                         .by = fc)
    nmds_plot <- nmds_plot +
      geom_point(data=centroid, aes(fill = !!fc_eval), size = 5, shape=21, stroke = 1.5, alpha=0.8)
  } else {
    my.envfit <- envfit(nmds_obj, meta, permutations = 9999)
    env.scores <- as.data.frame(scores(my.envfit, display = "vectors"))
    env.scores <- cbind(env.scores, env.variables = rownames(env.scores))
    env.scores <- cbind(env.scores, pval = my.envfit$vectors$pvals)

    scaled_arrow <- scale_arrow(end_x = env.scores$NMDS1,
                                end_y = env.scores$NMDS2,
                                minx = min(nmds_scores$NMDS1),
                                maxx = max(nmds_scores$NMDS1),
                                miny = min(nmds_scores$NMDS2),
                                maxy = max(nmds_scores$NMDS2))
    nmds_plot <- nmds_plot +
      geom_segment(aes(x = 0, xend = scaled_arrow["x"], y = 0, yend = scaled_arrow["y"]),
                   arrow = arrow(length = unit(0.25, "cm")), color = "grey10", lwd=0.3) +
      labs(subtitle = paste0(nmds_plot$labels$subtitle, "\nArrow envfit p.val: ", signif(env.scores$pval, 4)))
  }
  nmds_plot <- nmds_plot +
    labs(subtitle = paste0(nmds_plot$labels$subtitle, "\nNMDS stress: ", signif(nmds_obj$stress, 4)))

  if (save_image){
    ggsave(plot = nmds_plot,
           filename = save_path,
           bg = "white",
           width = 10, height = 7.8)
  }
  nmds_plot
}

scale_arrow <- function(end_x, end_y, minx, maxx, miny, maxy) {

  if (end_x < maxx && end_x > minx && end_y < maxy && end_y > miny){
    return(c(x = end_x, y = end_y))
  }

  scale_factor <- c(maxx/end_x,
                    minx/end_x,
                    maxy/end_y,
                    miny/end_y) |>
    abs() |>
    Filter(\(x) x <= 1, x=_) |>
    min()
  a_scaled <- c(x = end_x * scale_factor, y = end_y * scale_factor)
}

run_adonis <- function(otu, meta, fcs, outfile=NULL, factor_interaction=TRUE){
  formula_char <- ifelse(factor_interaction, "*", "+")
  formula <- as.formula(glue("otu ~ {str_flatten(fcs, collapse = formula_char)}"))

  permanova <- adonis(formula = formula,
                      data = meta,
                      permutations = n_permutations,
                      method = "bray")

  if (!is.null(outfile)){
    sink(outfile, append = TRUE)
    print(permanova$aov.tab)
    sink()
  }
  permanova
}

homogeneity_test <- function(dist, meta, fcs, outfile=NULL){
  test_results <- list()
  for (fc in fcs){
    if(n_distinct(meta[[fc]]) > 1) {
      bd <- anova(betadisper(dist, meta[[fc]]))
      test_results[[fc]] <- bd
      if(!is.null(outfile)){
        cat("\nANOVA check on homogeneity of", fc, "\n", file = outfile, append = T)
        sink(outfile, append = TRUE)
        print(bd)
        sink()
      }
    }
  }
  test_results
}

