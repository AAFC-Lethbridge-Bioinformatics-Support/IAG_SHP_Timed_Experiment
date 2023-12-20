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

sigdif_taxa_plot <- function(sigdif_taxa, ancom_result, meta) {
  # Why does the bias correct log table have NAs?
  sig_log_table <- ancom_result$bias_correct_log_table |>
    rownames_to_column("taxon_name") |>
    pivot_longer(cols = !taxon_name, names_to = "sample", values_to = "bias.corrected.log.reads") |>
    filter(taxon_name %in% sigdif_taxa) |>
    mutate(taxon_name = str_replace(taxon_name, "Candidatus", "Ca")) |>
    inner_join(meta)

  sig_log_table |>
    ggplot(aes(x = month.continuous, y = bias.corrected.log.reads)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(taxon_name ~ ., scales = "free_y") +
    theme_bw() +
    labs(title = "Significantly different taxa across Month")
}


# outputs: volcano plot, sigdif LM plots, sigdif table
ancombc_wrap <- function(pseq, formula, meta, outdir, subtitle="", filename_append="", use_saved_object=TRUE) {
  ancom_object_path <- glue("{outdir}/ancom_month_result{filename_append}.rds")
  ancom_result <- if (file.exists(ancom_object_path) && use_saved_object) {
    readRDS(ancom_object_path)
  } else {
    start_time <- Sys.time()
    ancom_result <- ancombc2(data = pseq,
             fix_formula = formula,
             p_adj_method = "BH",
             n_cl = N_CPU,
             pseudo_sens = TRUE)
    end_time <- Sys.time()
    message(glue("ANCOMBC2 ran for {end_time - start_time} mins. {filename_append}"))
    saveRDS(ancom_result, ancom_object_path)
    ancom_result
  }



  ancom_tidy <- format_ancombc_results(ancom_result)

  # pull out month results
  month_only <- ancom_tidy |> filter(metavariable == "month.continuous")

  # save csv
  month_only |> select(taxon, log2fc, q, diff, direction,
                       passed.ss) |>
    arrange(-diff, -passed.ss, q) |>
    write_csv(file = glue("{outdir}/da-table_month-continuous{filename_append}.csv"))

  # make and save volcano plot
  v_plot <- volc_plot(month_only, "Differential abundance across Month",
                      subtitle = subtitle)
  ggsave(v_plot,
         filename = glue("{outdir}/volcplot_month-continuous{filename_append}.png"))

  if(any(month_only$diff)) {
    # make and save differentially abundant taxa linear model plots
    sig_taxa <- month_only |> filter(diff == TRUE) |> arrange(taxon) |> pull(taxon)
    if (length(sig_taxa > 25)) {
      save_multi_da_plots(sig_taxa, ancom_result, meta, outdir, filename_append)
    } else {
      da_plot <- sigdif_taxa_plot(sigdif_taxa = sig_taxa,
                                  ancom_result = ancom_result,
                                  meta = timed_meta)
      ggsave(da_plot,
             filename = glue("{outdir}/da-plot_month-continuous{filename_append}.png"),
             width = 9,
             height = 7)
    }

  } else {
    message("ANCOMBC2: No significantly different taxa were found.")
  }
}

# split sigdif DA taxa into groups of 25 and save separate plots (for legibility)
save_multi_da_plots <- function(sig_taxa, ancom_result, meta, outdir, filename_append){
  chunks <- split(sig_taxa, ceiling(seq_along(sig_taxa)/25))
  for (nm in names(chunks)) {
    da_plot <- sigdif_taxa_plot(sigdif_taxa = chunks[[nm]],
                                ancom_result = ancom_result,
                                meta = meta)
    padded_name <- str_pad(nm, width = max(nchar(names(chunks))),
                           side = "left", pad = "0")
    ggsave(da_plot,
           filename = glue("{outdir}/da-plot_month-continuous{filename_append}_{padded_name}.png"),
           width = 9,
           height = 7)
  }
}

format_ancombc_results <- function(ancom_result) {
  # make tidy format
  ancom_tidy <- ancom_result$res |>
    rename_with(.cols = starts_with("passed_ss"), .fn = str_replace, "_", ".") |>
    pivot_longer(cols = !taxon, names_sep = "_", names_to = c(".value", "metavariable"))

  # add log2, direction, and pseudo-sensitivity (boolean) columns
  ancom_tidy |>
    mutate(log2fc = log2(exp(lfc)),
           direction = ifelse(q < 0.05,
                              ifelse(log2fc > 0,
                                     "Increase",
                                     "Decrease"),
                              "No Change"))
}
