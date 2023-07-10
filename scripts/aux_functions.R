library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(roxygen2)

color_palette <- c("#89C5DA", "#DA5724", "#689030", "#CE50CA", "#3F4921",
                   "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", "#D7C1B1",
                   "#D14285", "#AD6F3B", "#CD9BCD", "#74D944", "#6DDE88", "#652926", "#7FDCC0",
                   "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#38333E", "#599861")

#' Save and/or return an NMDS plot
#'
#' @param df a dataframe with NMDS1, NMDS2, and metadata factor columns
#' @param meta_factor string, a metadata factor to plot on the NMDS graph. Must
#'   correspond to a column in `df`
#' @param dataset_name string, optional, a dataset name to be inserted into plot
#'   title
#' @param taxa_level string, optional, a taxa_level to be inserted into plot
#'   subtitle
#' @param facet_factor string, optional, a metadata factor with which to facet
#'   the NMDS plot. Alters the file name by default.
#' @param polygon boolean, whether to include polygons to enclose the groups of
#'   the metadata factor. Defaults to FALSE
#' @param remove_NA boolean, whether to remove rows in which the `meta_factor`
#'   is NA. Defaults to TRUE
#' @param save_image boolean, whether to save the plot as a `.png` file.
#'   Defaults to FALSE
#' @param save_path string, the file path at which to save the `.png` file, if
#'   `save_image` is set to TRUE
#' @param appended_naming boolean, whether to append indicators to the output
#'   image file name base on other arguments. E.g. `_polygons` appended if
#'   `polygon` argument is set to TRUE. Meant for distinguishing files when
#'   multiple versions are to be saved.
plot_nmds_by_factor <- function(df,
                           meta_factor,
                           dataset_name = "",
                           taxa_level = "",
                           facet_factor = "",
                           shape_factor = "",
                           polygon = FALSE,
                           remove_NA = TRUE,
                           save_image = FALSE,
                           save_path = "nmds.png",
                           subtitle_append = "",
                           appended_naming = FALSE) {

  if (remove_NA) { df %<>% drop_na(meta_factor) }

  plot <- create_nmds_plot_by_factor(
    original_df = df,
    meta_factor = meta_factor,
    dataset_name = dataset_name,
    taxa_level = taxa_level,
    facet_factor = facet_factor,
    shape_factor = shape_factor,
    polygon = polygon,
    subtitle_append = subtitle_append
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


create_nmds_plot_by_factor <-
  function(original_df,
           meta_factor,
           dataset_name,
           taxa_level,
           facet_factor,
           shape_factor,
           polygon,
           subtitle_append) {

  df <- original_df
  meta_var <- sym(meta_factor)

  is_faceted <- facet_factor != ""
  if (is_faceted) { facet_var <- sym(facet_factor) }
  is_shapes <- shape_factor != ""
  if (is_shapes) {
    shapes_var <- sym(shape_factor)
  }

  # if numeric, make into factor by splitting range of values into n groups
  if (is.numeric(df[[meta_factor]])) {
    df <- bin_continuous_factor(df, meta_factor, 5)
  }

  plot_title <- ifelse(dataset_name == "",
                       glue("Taxonomic ordination by factor {meta_factor}"),
                       glue("Taxonomic ordination of {dataset_name} samples by factor {meta_factor}"))
  if (is_faceted) {
    plot_title %<>% paste(glue("across {facet_factor}"))
  }
  plot_subtitle <- ifelse(taxa_level == "", "",
                          glue("Taxonomic rank: {taxa_level}"))
  plot_subtitle <- paste(plot_subtitle, subtitle_append)

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

  if (is.factor(df[[meta_factor]]) & nlevels(df[[meta_factor]]) <= 5) {
    plot <- plot +
      scale_fill_viridis_d() +
      scale_color_viridis_d()
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


  return(plot)
  }

bin_continuous_factor <- function(meta_table, factor_of_interest, n_bins=5) {
  still_trying <- TRUE
  while(still_trying) {
    tryCatch({
      meta_table[[factor_of_interest]] <- as.factor(cut_number(meta_table[[factor_of_interest]], n=n_bins))
      still_trying <- FALSE
    },
    error = function(cond) {
      message(paste(cond, "Retrying with 1 less bin."))
    },
    finally = { n_bins <- n_bins - 1})
  }

  return(meta_table)
}
