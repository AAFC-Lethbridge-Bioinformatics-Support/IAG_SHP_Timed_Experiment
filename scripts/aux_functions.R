library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(roxygen2)

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
    rename(sample = `MBI_ID`) |>
    filter(!str_detect(Site, "B\\d")) # Blank samples denoted by B#
  timed_meta$Month <- as.numeric(timed_meta$Month)
  timed_meta <- timed_meta |>
    mutate(Month = case_when(Year == "1998" | Year == "2007" ~ as.numeric(Year),
                             TRUE ~ Month),
           month_group = case_when(Month < 1 ~ "<1 month",
                                   Month <= 6 ~ "1-6 months",
                                   Month <= 18 ~ "12-18 months"))
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
  timed_meta$month_group <- ordered(timed_meta$month_group,
                                    levels = c("<1 month",
                                               "1-6 months",
                                               "12-18 months"))
  timed_meta$month_continuous <- as.numeric(levels(timed_meta$Month))[timed_meta$Month]
  timed_meta$Year <- as.character(timed_meta$Year)
  timed_meta
}

prep_foi_data <- function(count_table, metadata, foi){
  # specific filters for Month and Process factors
  if (foi == "Month" | foi == "month_group"){
    foi_metadata <- filter(metadata , Process == "Ground")
  } else if (foi == "Process") {
    foi_metadata <- filter(metadata, Month == "0")
  } else {
    foi_metadata <- metadata
  }

  # Ensure metadata doesn't have samples that aren't in count table and vice versa
  foi_metadata <- filter(foi_metadata, sample %in% count_table$sample)
  foi_filtered <- filter(count_table, sample %in% foi_metadata$sample)

  # Ensure metadata sample order is identical to that of count table
  foi_metadata <- select(foi_filtered, sample) |>
    left_join(foi_metadata)

  foi_matrix <- foi_filtered |>
    column_to_rownames(var = "sample") |>
    as.matrix()
  foi_matrix[is.na(foi_matrix)] <- 0

  # ASSERTION
  all(rownames(foi_matrix) == foi_metadata$sample) |>
    assertthat::assert_that(msg = "Count table and metadata table do not have identical sample columns after processing. They must match in value and order.")

  list(matrix = foi_matrix,
       metadata = foi_metadata
  )
}

run_anosim <- function(count_matrix, metadata, foi, random_seed){
  set.seed(random_seed)
  ano_obj <- anosim(count_matrix, metadata[[foi]], distance = "bray", permutations = 999)

  ano_string <- glue("ANOSIM p-val: {ano_obj$signif} R-stat: {ano_obj$statistic}")

  list(ano_obj = ano_obj,
       ano_string = ano_string)
}

save_anosim <- function(anosim_result_file, ano_string, foi, cur_depth=NULL){
  header <- if (is.null(cur_depth)){
    cat("\n", file = anosim_result_file, append = T)
    glue("FACTOR: {foi}")
  } else {
    glue("FACTOR: {foi} at soil depth {cur_depth}")
  }
  cat(header, "\n", file = anosim_result_file, append = T)
  cat(ano_string, "\n", file = anosim_result_file, append = T)
}

run_nmds <- function(count_matrix, metadata, random_seed){
  set.seed(random_seed)
  nmds_obj <- metaMDS(count_matrix, distance = "bray")

  nmds_scores <- as.data.frame(scores(nmds_obj)$sites) |>
    rownames_to_column("sample")
  # Merge NMDS results with metadata
  nmds_scores <- left_join(nmds_scores, metadata)

  list(nmds_obj = nmds_obj,
       nmds_scores = nmds_scores)
}

save_nmds <- function(nmds_scores, ano_string, foi, outdir, taxa_level=NULL, cur_depth=NULL){
  if (is.null(cur_depth)) {
    save_path <- glue("{outdir}/nmds_{foi}.png")
    subtitle <- ano_string
  } else {
    save_path <- glue("{outdir}/nmds_{foi}_{cur_depth}.png")
    subtitle <- paste(ano_string,
                      glue("Soil depth: {cur_depth} (cm)"))
  }
  nmds_plot <- plot_nmds_by_factor(
    df = nmds_scores,
    meta_factor = foi,
    dataset_name = "Timed",
    taxa_level = taxa_level,
    shape_factor = "Depth",
    polygon = TRUE,
    save_image = TRUE,
    save_path = save_path,
    subtitle_append = subtitle
  )
  nmds_plot
}

run_indicspecies <- function(count_matrix, metadata, foi, random_seed){
  set.seed(random_seed)
  indic_obj = multipatt(count_matrix, metadata[[foi]], func = "r.g",
                           control = how(nperm=9999))
  # Format significant indicspecies results in table
  indic_table <- indic_obj$sign |>
    mutate("significant" = p.value <= 0.05) |>
    rownames_to_column("feature") |>
    arrange(p.value, -stat)

  list(indic_obj = indic_obj,
       indic_table = indic_table)
}

save_indicspecies <- function(indic_obj, indic_table, foi, outdir, cur_depth=NULL){
  base_name <- glue("{outdir}/indicspecies_{foi}")
  if (!is.null(cur_depth)) {
    base_name <- glue("{base_name}_{cur_depth}")
  }

  write_csv(indic_table, glue("{base_name}.csv"))

  sink(glue("{base_name}.txt"))
  summary(indic_obj)
  sink()
}

# Run ANOSIM, NMDS, and indicspecies analyses. These are grouped because of
# common input format
run_ano_nmds_indic <- function(count_table, metadata, foi){
  random_seed <- 86752155

  # specific filters for Month and Process factors
  if (foi == "Month" | foi == "month_group"){
    foi_metadata <- filter(metadata , Process == "Ground")
  } else if (foi == "Process") {
    foi_metadata <- filter(metadata, Month == "0")
  } else {
    foi_metadata <- metadata
  }

  # Skip the metadata factor if there is only 1 distinct value in it after filtering
  # if(length(unique(foi_metadata[[foi]])) < 2) { next }

  # Ensure metadata doesn't have samples that aren't in count table and vice versa
  foi_metadata <- filter(foi_metadata, sample %in% count_table$sample)
  foi_filtered <- filter(count_table, sample %in% foi_metadata$sample)

  # Ensure metadata sample order is identical to that of count table
  foi_metadata <- select(foi_filtered, sample) |>
    left_join(foi_metadata)

  rel_abund_matrix_foi <- foi_filtered |>
    column_to_rownames(var = "sample") |>
    as.matrix()
  rel_abund_matrix_foi[is.na(rel_abund_matrix_foi)] <- 0

  # ASSERTION
  all(rownames(rel_abund_matrix_foi) == foi_metadata$sample) |>
    assertthat::assert_that(msg = "Count table and metadata table do not have identical sample columns. They must match in value and order.")

  set.seed(random_seed)
  ano_foi <- anosim(rel_abund_matrix_foi, foi_metadata[[foi]], distance = "bray", permutations = 999)

  ano_string_foi <- glue("ANOSIM p-val: {ano_foi$signif} R-stat: {ano_foi$statistic}")

  # Perform NMDS ----
  set.seed(random_seed)
  nmds_foi <-  metaMDS(rel_abund_matrix_foi, distance = "bray")
  plot(nmds_foi)

  nmds_scores_foi <- as.data.frame(scores(nmds_foi)$sites) |>
    rownames_to_column("sample")
  # Merge NMDS results with metadata
  nmds_scores_foi <- left_join(nmds_scores_foi, foi_metadata)

  if (foi == "Month") {
    # Don't run indicspecies for month. It will be handled with month_group
    return(list(anosim = ano_foi,
                ano_string = ano_string_foi,
                nmds = nmds_foi,
                nmds_scores = nmds_scores_foi,
                indic = NULL,
                indic_table = NULL))
  }

  # Perform indicspecies analysis  ----
  set.seed(random_seed)
  indic_by_foi = multipatt(rel_abund_matrix_foi, foi_metadata[[foi]], func = "r.g",
                           control = how(nperm=9999))
  # Format significant indicspecies results in table
  indic_foi_tbl <- indic_by_foi$sign |>
    mutate("significant" = p.value <= 0.05) |>
    rownames_to_column("feature") |>
    arrange(p.value, -stat)

  # Return relevant analysis outputs
  list(anosim = ano_foi,
       ano_string = ano_string_foi,
       nmds = nmds_foi,
       nmds_scores = nmds_scores_foi,
       indic = indic_by_foi,
       indic_table = indic_foi_tbl)
}

# Save and return an NMDS plot
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

  plot_title <- ifelse(dataset_name == "",
                       glue("NMDS ordination by factor {meta_factor}"),
                       glue("NMDS ordination of {dataset_name} samples by factor {meta_factor}"))
  if (is_faceted) {
    plot_title %<>% paste(glue("across {facet_factor}"))
  }
  plot_subtitle <- ifelse(is.null(taxa_level), "",
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
  } else if (is.numeric(df[[meta_factor]])) {
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
