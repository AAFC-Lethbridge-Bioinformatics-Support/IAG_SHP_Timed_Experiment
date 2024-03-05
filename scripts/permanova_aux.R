# This script contains functions used solely by the permanova.R script

# Load libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)

analyses_wrap <- function(pseq,
                          fcs,
                          outdir,
                          outfile_extension,
                          unit_level,
                          factor_interaction = TRUE,
                          rarefaction = FALSE) {
  factor_results_list <- list()
  PERMANOVA_result_file <- glue("{outdir}/PERMANOVA_{outfile_extension}.txt")
  adonis_text_header_extension <- if (unit_level == "woltka_KO") {
    glue("for Woltka functional profiling (KEGG KO counts)")
  } else if (unit_level == "woltka_pathways") {
    glue("for Woltka functional profiling (KEGG pathways)")
  } else {
    glue("for taxa level {str_to_upper(unit_level)}")
  }
  adonis_text_header <- glue("PERMANOVA via vegan::adonis2 {adonis_text_header_extension}\n\n", .trim = FALSE)
  cat(adonis_text_header, file = PERMANOVA_result_file)

  otu <- otu_table(pseq)
  meta <- as(sample_data(pseq), "data.frame")
  rownames(meta) <- sample_names(pseq)

  if (rarefaction) {
    smallest_sample_count <- otu |> rowSums() |> min()
    if (smallest_sample_count == 1) {
      stop("Relative abundance values detected. Rarefaction option must be used with raw counts. Use pseq$base instead of pseq$rel.")
    }
    # Rarefy counts via avgdist
    dist <- otu |> as.data.frame() |>
      avgdist(dmethod = "bray", sample = smallest_sample_count)
  } else {
    dist <- vegdist(otu)
  }

  adonis_result <- run_adonis(otu, meta, fcs, factor_interaction)

  sink(PERMANOVA_result_file, append = TRUE)
  print(adonis_result)
  sink()

  for (fc in fcs) {
    subtitle <- if (unit_level == "woltka_KO") {
      "Woltka functional profiling (KEGG KO counts)"
    } else if (unit_level == "woltka_pathways") {
      "Woltka functional profiling (KEGG pathways)"
    } else {
      glue("Taxonomic rank: {unit_level}")
    }

    homog_test <- anova(betadisper(dist, meta[[fc]]))
    cat("\nANOVA check on homogeneity of", fc, "\n", file = PERMANOVA_result_file, append = T)
    sink(PERMANOVA_result_file, append = TRUE)
    print(homog_test)
    sink()
    homog_pval <- homog_test$`Pr(>F)`[1]
    subtitle <- glue("{subtitle}\nHomogeneity of variance p.val: {signif(homog_pval, 4)}")

    fc_adonis <- adonis_result[fc,]
    fc_adonis_R2 <- fc_adonis$R2[1]
    fc_adonis_pval <- fc_adonis$`Pr(>F)`[1]
    subtitle <- glue("{subtitle}\nAdonis PERMANOVA test p.val: {signif(fc_adonis_pval, 4)}  R2: {signif(fc_adonis_R2, 4)}")
    plot_nmds(dist, meta, fc,
              subtitle=subtitle,
              save_image=TRUE,
              save_path=glue("{outdir}/nmds_{outfile_extension}_{fc}.png"))
    factor_results_list[[fc]] <- list(homog_pval = homog_pval,
                                      adonis_pval = fc_adonis_pval,
                                      adonis_R2 = fc_adonis_R2)
  }
  factor_results_list
}

add_row_to_summary <- function(result_list,
                               summary_file,
                               unit_level,
                               seq_depth,
                               subset) {
  for (fc in names(result_list)){
    summary_file <- summary_file |>
      add_row(meta_factor = fc,
              unit_level = unit_level,
              sequence_depth = seq_depth,
              subset = subset,
              homogeneity_pval = result_list[[fc]]$homog_pval,
              permanova_pval = result_list[[fc]]$adonis_pval,
              permanova_R2 = result_list[[fc]]$adonis_R2)
  }
  summary_file
}



# Analysis functions ------------------------------------------------------

plot_nmds <- function(dist, meta, fc, save_image=FALSE, save_path, subtitle, ...){
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
    shape_factor = "Depth",
    polygon = polygon,
    subtitle = subtitle,
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

run_adonis <- function(otu, meta, fcs, factor_interaction=TRUE, n_permutations = 9999, ...){
  formula_char <- ifelse(factor_interaction, "*", "+")
  formula <- as.formula(glue("otu ~ {str_flatten(fcs, collapse = formula_char)}"))

  permanova <- adonis2(formula = formula,
                       data = meta,
                       permutations = n_permutations,
                       method = "bray",
                       ...)
}

# Save and return an NMDS plot
plot_nmds_by_factor <- function(df,
                                meta_factor,
                                dataset_name = "",
                                facet_factor = "",
                                shape_factor = "",
                                polygon = FALSE,
                                remove_NA = TRUE,
                                save_image = FALSE,
                                save_path = "nmds.png",
                                subtitle = "",
                                appended_naming = FALSE) {

  if (remove_NA) { df %<>% drop_na(meta_factor) }

  plot <- create_nmds_plot_by_factor(
    original_df = df,
    meta_factor = meta_factor,
    dataset_name = dataset_name,
    facet_factor = facet_factor,
    shape_factor = shape_factor,
    polygon = polygon,
    plot_subtitle = subtitle
  )

  if (save_image) {
    if (appended_naming) {
      if (polygon) save_path %<>% str_replace(".png", "_polygons.png")
      if (!remove_NA) save_path %<>% str_replace(".png", "_with_NA.png")
    }
    if (facet_factor != "") {
      save_path %<>% str_replace(".png", glue("_by_{facet_factor}.png"))
    }
    tryCatch({
      ggsave(plot = plot,
             filename = save_path,
             bg = "white",
             width = 10, height = 7.8)
    },
    error = function(cond) { message(cond) },
    finally = {})
  }

  return(plot)
}


create_nmds_plot_by_factor <- function(original_df,
                                       meta_factor,
                                       dataset_name,
                                       facet_factor,
                                       shape_factor,
                                       polygon,
                                       plot_subtitle) {

  df <- original_df
  meta_var <- sym(meta_factor)

  is_faceted <- facet_factor != ""
  if (is_faceted) { facet_var <- sym(facet_factor) }
  is_shapes <- shape_factor != ""
  if (is_shapes) {
    shapes_var <- sym(shape_factor)
  }

  plot_title <- ifelse(dataset_name == "",
                       glue("NMDS ordination by factor {meta_factor}"),
                       glue("NMDS ordination of {dataset_name} samples by factor {meta_factor}"))
  if (is_faceted) {
    plot_title %<>% paste(glue("across {facet_factor}"))
  }

  plot <- df %>%
    ggplot(aes(x = NMDS1, y = NMDS2), na.rm = TRUE) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(),
          legend.key=element_blank()) +
    theme_bw() +
    labs(x = "NMDS1", y = "NMDS2", colour = meta_factor,
         title = plot_title,
         subtitle = plot_subtitle)

  if (is_faceted) {
    plot <- plot + facet_wrap(vars(!!facet_var))
  }

  if (is.numeric(df[[meta_factor]])) {
    plot <- plot +
      scale_fill_viridis_c() +
      scale_color_viridis_c()
  } else {
    plot <- plot +
      scale_color_manual(values=color_palette) +
      scale_fill_manual(values=color_palette)
  }
  if (polygon) {
    if (is_faceted) {
      hull <- group_by(df, !!facet_var, !!meta_var)
    } else {
      hull <- group_by(df, !!meta_var)
    }
    hull <- hull %>%
      slice(grDevices::chull(NMDS1, NMDS2))
    plot <- plot +
      geom_polygon(data = hull, alpha = 0.1, aes(color = !!meta_var), fill = NA, linewidth = 1.5)
    # plot <- plot + stat_ellipse(mapping = aes(color = !!meta_var), size = 2, alpha = 0.5, type = "t")
  }

  if (is_shapes) {
    plot <- plot + geom_point(mapping = aes(fill = !!meta_var, shape = !!shapes_var),
                              size = 3, alpha = 0.7, na.rm = TRUE) +
      scale_shape_manual(values = c(21,24,23,22))
  } else {
    plot <- plot + geom_point(mapping = aes(fill = !!meta_var),
                              size = 3, alpha = 0.7, na.rm = TRUE, shape = 21)
  }

  plot
}
