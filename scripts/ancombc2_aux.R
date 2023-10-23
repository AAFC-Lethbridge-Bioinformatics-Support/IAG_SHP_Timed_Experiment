# This script contains functions used solely by the ancombc2.R script

# Load libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(vegan)
library(phyloseq)
library(ANCOMBC)


volc_plot <- function(df, plot_title="", subtitle="") {
  # determine x range
  x_default <- 1.35
  x_data <- max(abs(df$log2fc))
  # x_use <- max(x_default, x_data) + 0.15
  x_use <- x_data + 0.15

  # determine y range
  y_default <- 1.40
  y_data <- -log10(min(df$q))
  y_use <- max(y_default, y_data) + 0.10

  # plot
  da_plot <- df |>
    ggplot(aes(x = log2fc,
               y = -log10(q),
               color = direction,
               shape = pseudo.sensitive)
    ) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c(
      "Increase" = "blue",
      "No change" = "black",
      "Decrease" = "red"
    )) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4)) +
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

sigdif_taxa_plot <- function(sigdif_taxa, counts, meta) {
  # use matrix because NAs have been set to 0 for accurate lm plotting
  sig_counts <- counts$matrix |>
    as.data.frame() |>
    rownames_to_column("sample") |>
    pivot_longer(cols = -sample, names_to = "taxon_name", values_to = "reads") |>
    filter(taxon_name %in% sigdif_taxa) |>
    mutate(taxon_name = str_replace(taxon_name, "Candidatus", "Ca")) |>
    inner_join(meta)

  sig_counts |>
    ggplot(aes(x = month.continuous, y = reads)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(taxon_name ~ ., scales = "free_y") +
    theme_bw() +
    labs(title = "Significantly different taxa across Month")
}


# outputs: volcano plot, sigdif LM plots, sigdif table
ancombc_wrap <- function(pseq, formula, filename_append="") {
  # month as continuous variable
  ancom_result <- ancombc2(data = pseq,
                           fix_formula = formula,
                           p_adj_method = "BH",
                           n_cl = 4,
                           pseudo_sens = TRUE)
  ancom_tidy <- format_ancombc_results(ancom_result)

  # pull out month results
  month_only <- ancom_tidy |> filter(metavariable == "month.continuous")

  # save csv
  month_only |> select(taxon, log2fc, q, diff, direction,
                       pseudo.sensitivity.p, pseudo.sensitive) |>
    arrange(-diff, -pseudo.sensitive, q) |>
    write_csv(file = glue("{taxa_out}/da-table_month-continuous{filename_append}.csv"))

  # make and save volcano polt
  v_plot <- volc_plot(month_only, "Differential abundance across Month",
                      subtitle = glue("Taxa level: {taxa_level}"))
  ggsave(v_plot,
         filename = glue("{taxa_out}/volcplot_month-continuous{filename_append}.png"))

  if(any(month_only$diff)) {
    # make and save differentially abundant taxa linear model plot
    da_plot <- month_only |> filter(diff == TRUE) |> pull(taxon) |>
      sigdif_taxa_plot(counts = taxa_counts, meta = timed_meta)
    ggsave(da_plot,
           filename = glue("{taxa_out}/da-plot_month-continuous{filename_append}.png"),
           width = 9,
           height = 7)
  } else {
    messsage("ANCOMBC2: No significantly different taxa were found.")
  }
}

format_ancombc_results <- function(ancom_result) {
  # make tidy format
  ancom2_tidy <- ancom_result$res |>
    pivot_longer(cols = !taxon, names_sep = "_", names_to = c(".value", "metavariable"))
  # merge pseudo sensitivity test results
  ancom3_tidy <- ancom_result$pseudo_sens_tab |>
    pivot_longer(cols = !taxon,  names_to = "metavariable", values_to = "pseudo.sensitivity.p") |>
    inner_join(ancom2_tidy)

  # add log2, direction, and pseudo-sensitivity (boolean) columns
  ancom4_tidy <- ancom3_tidy |>
    mutate(log2fc = log2(exp(lfc)),
           pseudo.sensitive = pseudo.sensitivity.p > 0.05,
           direction = ifelse(q < 0.05,
                              ifelse(log2fc > 0,
                                     "Increase",
                                     "Decrease"),
                              "No Change"))
}
