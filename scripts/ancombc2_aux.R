# This script contains functions used solely by the ancombc2.R script

# Load libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(vegan)
library(phyloseq)
library(ANCOMBC)

# Used as argument to ANCOMBC2 for parallel computing
N_CPU <- 4

volc_plot <- function(df, plot_title="", subtitle="") {
  # determine x range
  x_default <- 1.35
  x_data <- if (any(df$diff)) {
    max(abs(filter(df, diff)$log2fc), na.rm = TRUE)
  } else {
    max(abs(df$log2fc), na.rm = TRUE)
  }
  # Use 20% above the greatest absolute, signficant (if any) point to limit
  # extreme outliers below significance
  x_use <- x_data*1.2

  # determine y range
  y_default <- 1.40
  y_data <- -log10(min(df$q))
  y_use <- max(y_default, y_data) + 0.10

  # plot
  da_plot <- df |>
    ggplot(aes(x = log2fc,
               y = -log10(q),
               color = direction,
               shape = passed.ss)
    ) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c(
      "Increase" = "blue",
      "No change" = "black",
      "Decrease" = "red"
    )) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4)) +
    xlim(0 - x_use, 0 + x_use) +
    ylim(0, y_use) +
    geom_hline(
      yintercept = -log10(0.05),
      lty = 4,
      col = "black",
      lwd = 0.8
    ) +
    labs(x = "log2(fold change)",
         y = "-log10(adj. p-value)",
         color = "Direction",
         shape = "Pseudo-count sensitive",
         title = plot_title,
         subtitle = subtitle) +
    theme(legend.position = "right")

  return(da_plot)
}


plot_sigdif_taxa <- function(ancom_result, plot_var, sig_taxa) {
  if (!plot_var %in% colnames(ancom_result$meta)) stop(glue("Specified metadata {plot_var} variable not in metadata table"))
  # Why does the bias correct log table have NAs?
  sig_log_table <- ancom_result$bias_correct_log_table |>
    rownames_to_column("taxon_name") |>
    pivot_longer(cols = !taxon_name, names_to = "sample", values_to = "bias.corrected.log.reads") |>
    filter(taxon_name %in% sig_taxa) |>
    mutate(taxon_name = str_replace(taxon_name, "Candidatus", "Ca")) |>
    inner_join(ancom_result$meta)

  plot <- sig_log_table |>
    ggplot(aes(x = .data[[plot_var]], y = bias.corrected.log.reads)) +
    theme_bw()

  # Numeric factors will make scatter plots, rather than box plots
  plot <- if (is.numeric(sig_log_table[[plot_var]])) {
    plot +
      geom_point() +
      geom_smooth(method = "lm")
  } else {
    plot +
      geom_boxplot() +
      geom_jitter(color = "darkgrey", alpha = 0.7) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  plot +
    facet_wrap(taxon_name ~ ., scales = "free_y") +
    labs(title = glue("Significantly different taxa across {plot_var}"))
}

# split sigdif DA taxa into groups of 25 and save separate plots (for legibility)
save_multi_da_plots <- function(ancom_result, plot_var, sig_taxa, outdir, filename_append){
  chunks <- split(sig_taxa, ceiling(seq_along(sig_taxa)/25))
  for (nm in names(chunks)) {
    da_plot <- plot_sigdif_taxa(ancom_result = ancom_result,
                                plot_var = plot_var,
                                sig_taxa = chunks[[nm]])
    padded_name <- str_pad(nm, width = max(nchar(names(chunks))),
                           side = "left", pad = "0")
    ggsave(da_plot,
           filename = glue("{outdir}/da-plot_{plot_var}{filename_append}_{padded_name}.png"),
           width = 10.5,
           height = 8)
  }
}


# Returns ANCOMBC2 non-pairwise results in tidy format, pivoting and creating a
# column that specifies metavariable.
tidy_ancombc_res <- function(ancom_result) {
    ancom_tidy <- ancom_result$res |>
      rename_with(.cols = starts_with("passed_ss"), .fn = str_replace, "_", ".") |>
      pivot_longer(cols = !taxon, names_sep = "_", names_to = c(".value", "metavariable")) |>
      mutate(log2fc = log2(exp(lfc)),
             direction = ifelse(q < 0.05,
                                ifelse(log2fc > 0,
                                       "Increase",
                                       "Decrease"),
                                "No Change"))

}

# Returns ANCOMBC2 pairwise results in tidy format, pivoting and creating a
# column that specifies the pairwise comparison. Providing the metadata table
# (accessible by ancom_result$meta) will add the missing "default" group to the
# comparison column values. This remedies how ancombc2 drops the name of one
# group from the column names.
tidy_ancombc_res_pair <- function(ancom_result, meta_var) {
  ancom_tidy <- ancom_result$res_pair |>
    rename_with(.cols = starts_with("passed_ss"), .fn = str_replace, "_", ".") |>
    rename_with(.fn = str_replace, .cols = !taxon, "_", ">") |>
    pivot_longer(cols = !taxon, names_sep = ">", names_to = c(".value", "comparison")) |>
    mutate(comparison = str_replace_all(comparison, meta_var, ""),
           log2fc = log2(exp(lfc)),
           direction = ifelse(q < 0.05,
                              ifelse(log2fc > 0,
                                     "Increase",
                                     "Decrease"),
                              "No Change"))

  if (!is.null(ancom_result$meta)) {
    all_vals <- ancom_result$meta[[meta_var]] |> unique()
    result_vals <- ancom_tidy$comparison |> unique()
    missing_val <- setdiff(all_vals, result_vals)
    if (length(missing_val) > 1) stop("More than 1 missing value detected")

    ancom_tidy |>
      mutate(comparison = ifelse(!str_detect(comparison, "_"),
                                 paste0(comparison, "_", missing_val),
                                 comparison))
  } else {
    message("Metadata not attached to ancombc results (ancom_result$meta). Missing group name has not been appended to comparison column.")
    ancom_tidy
  }
}

ancombc2_wrap <- function(pseq,
                          meta_var,
                          factors,
                          outdir,
                          filename_append="",
                          use_saved_object = TRUE,
                          method = c("default", "pairwise")) {
  method <- rlang::arg_match(method)
  ancom_object_path <- glue("{outdir}/ancom_{meta_var}_result{filename_append}.rds")

  if (file.exists(ancom_object_path) && use_saved_object) {
    return(readRDS(ancom_object_path))
  }

  start_time <- Sys.time()
  if (method == "default") {
    ancom_result <- ancombc2(data = pseq,
                             fix_formula = str_flatten(factors, collapse = "+"),
                             p_adj_method = "BH",
                             n_cl = N_CPU,
                             pseudo_sens = TRUE)
  } else if (method == "pairwise") {
    ancom_result <- ancombc2(data = pseq,
                             fix_formula = str_flatten(factors, collapse = "+"),
                             p_adj_method = "BH",
                             n_cl = N_CPU,
                             pseudo_sens = TRUE,
                             group = meta_var,
                             pairwise = TRUE)
  }
  end_time <- Sys.time()
  message(glue("ANCOMBC2 ran for {end_time - start_time} mins. {meta_var}{filename_append}"))
  # Package metadata with ancombc result
  ancom_result$meta <- sample_data(pseq) |> as_tibble(rownames="sample")
  saveRDS(ancom_result, ancom_object_path)

  ancom_result
}


ancombc_wrap_process <- function(pseq, factors, outdir, subtitle="", filename_append="") {
  ancom_result <- ancombc2_wrap(
    pseq = pseq,
    meta_var = "Process.character",
    factors = factors,
    outdir = outdir,
    method = "pairwise"
  )
  ancom_tidy <- tidy_ancombc_res_pair(ancom_result,meta_var = "Process.character")

  # Format ANCOMBC2 results to be consistent with experiment design
  ancom_tidy |> arrange(-diff, -passed.ss, comparison, q) |>
    mutate(comparison = case_when(comparison == "Sieved_Ground" ~ "Ground_Sieved",
                                  comparison == "Sieved_Air-dried" ~ "Air-dried_Sieved",
                                  comparison == "Fresh_Air-dried" ~ "Air-dried_Fresh",
                                  TRUE ~ comparison),
           log2fc = case_when(comparison %in% c("Ground_Sieved", "Air-dried_Sieved", "Air-dried_Fresh") ~ -log2fc,
                              TRUE ~ log2fc),
           lfc = case_when(comparison %in% c("Ground_Sieved", "Air-dried_Sieved", "Air-dried_Fresh") ~ -lfc,
                           TRUE ~ lfc),
           W = case_when(comparison %in% c("Ground_Sieved", "Air-dried_Sieved", "Air-dried_Fresh") ~ -W,
                         TRUE ~ W),
           direction = case_when(comparison %in% c("Ground_Sieved", "Air-dried_Sieved", "Air-dried_Fresh") & direction == "Decrease" ~ "Increase",
                                 comparison %in% c("Ground_Sieved", "Air-dried_Sieved", "Air-dried_Fresh") & direction == "Increase" ~ "Decrease",
                                 TRUE ~ direction)) |>
    write_csv(file = glue("{outdir}/da-table_Process{filename_append}.csv"))

  if(any(ancom_tidy$diff)) {
    sig_taxa <- ancom_tidy |> filter(diff) |> arrange(taxon) |> pull(taxon) |> unique()
    # Use "Process" instead of "Process.character" to order groups in plots
    if (length(sig_taxa) > 25) {
      save_multi_da_plots(ancom_result = ancom_result,
                          plot_var = "Process",
                          sig_taxa = sig_taxa,
                          outdir = outdir,
                          filename_append = filename_append)
    } else {
      da_plot <- plot_sigdif_taxa(ancom_result = ancom_result,
                                  plot_var = "Process",
                                  sig_taxa = sig_taxa)
      ggsave(da_plot,
             filename = glue("{outdir}/da-plot_Process{filename_append}.png"),
             width = 9,
             height = 7)
    }

  } else {
    message("ANCOMBC2: No significantly different taxa were found.")
  }
}

# outputs: volcano plot, sigdif LM plots, sigdif table
ancombc_wrap_month <- function(pseq,
                               factors,
                               outdir,
                               subtitle="",
                               filename_append="") {
  ancom_result <- ancombc2_wrap(
    pseq = pseq,
    meta_var = "month-continuous",
    factors = factors,
    outdir = outdir,
    filename_append = filename_append,
    method = "default"
  )

  ancom_tidy <- tidy_ancombc_res(ancom_result)

  # pull out month results
  month_only <- ancom_tidy |> filter(metavariable == "month.continuous")
  # save csv
  month_only |> select(-metavariable) |>
    arrange(-diff, -passed.ss, q) |>
    write_csv(file = glue("{outdir}/da-table_month-continuous{filename_append}.csv"))

  # make and save volcano plot
  volc_plot(month_only, "Differential abundance across Month",
            subtitle = subtitle) |>
    ggsave(filename = glue("{outdir}/volcplot_month-continuous{filename_append}.png"))

  if(any(month_only$diff)) {
    # make and save differentially abundant taxa linear model plots
    sig_taxa <- month_only |> filter(diff) |> arrange(taxon) |> pull(taxon)
    if (length(sig_taxa) > 25) {
      save_multi_da_plots(ancom_result, "month.continuous", sig_taxa, outdir, filename_append)
    } else {
      da_plot <- plot_sigdif_taxa(ancom_result = ancom_result,
                                  plot_var = "month.continuous",
                                  sig_taxa = sig_taxa)
      ggsave(da_plot,
             filename = glue("{outdir}/da-plot_month-continuous{filename_append}.png"),
             width = 9,
             height = 7)
    }

  } else {
    message("ANCOMBC2: No significantly different taxa were found.")
  }
}
