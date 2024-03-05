# This script contains functions used solely by the alpha_diversity.R script

library(tidyverse)
library(vegan)
library(ggplot2)
library(glue)
library(dunn.test)

# dunn.test is run for pairwise testing of categorical factors.
# Test results are saved to text file and all sig. results are appended as a
# string and returned for use as plot subtitles
run_dunn_test <- function(d, alpha_index, meta_factor, outfile) {
  subtitle <- ""
  dt <- dunn.test(d[[alpha_index]], g=d[[meta_factor]], method = "bh")
  cat("\tKW test was significant. Running pairwise Dunn test\n", file = outfile, append = T)
  if (all(dt$P.adjusted > 0.05)) {
    cat("\tNo pairwise tests were <0.05 after p adjustment.\n\n", file = outfile, append = T)
    return(subtitle)
  }
  for (i in 1:length(dt$comparisons)) {
    comparison <- dt$comparisons[i]
    adj.p <- dt$P.adjusted[i]
    if (adj.p < 0.05) {
      comparison_string <- glue("Comparison ({comparison}) Adj. P-value: {signif(adj.p, 4)}\n", .trim=FALSE)
      cat(paste0("\t",comparison_string), file = outfile, append = T)
      subtitle <- paste0(subtitle, comparison_string)
    }
  }
  subtitle |> str_trim()
}

# runs and saves kruskal wallis test, returns KW test result object
alpha_KW_test <- function(d, alpha_index, meta_factor, outfile) {
  heading <- glue("{alpha_index} Index by {meta_factor} - Kruskal Wallis Results:\n")
  factor_alpha <- stats::kruskal.test(stats::formula(glue("{alpha_index}~{meta_factor}")), d)
  cat(heading, file = outfile, append = T)
  utils::capture.output(factor_alpha, file = outfile, append = T)

  factor_alpha
}

# Creates plot of alpha div metric for a given metadata factor. Bar plot if
# categorical, scatter plot if numeric. Returns ggplot object.
plot_alpha <- function(d, alpha, x_var, subtitle) {
  isContinuous <- is.numeric(d[[x_var]])

  alpha_plot <- ggplot(d, aes(x = .data[[x_var]], y = .data[[alpha]])) +
    labs(title = glue("Alpha Diversity: {alpha}"),
         subtitle = subtitle,
         x = x_var) +
    theme_bw()
  if (isContinuous) {
    alpha_plot <- alpha_plot +
      geom_point(alpha = 0.7) +
      geom_smooth(color = "black", method = "lm", se=FALSE)
  } else {
    alpha_plot <- alpha_plot+
      geom_boxplot() +
      geom_jitter(aes(color = Season), alpha = 0.7)
  }
  alpha_plot
}

# Wrapper for running KW test, running dunn.test if KW test was significant,
# and creates a plot. Returns KW test result object and ggplot in a named list
categorical_alpha_diversity <- function(plot_data,
                                        alpha_var,
                                        meta_var,
                                        filepath,
                                        subtitle,
                                        outdir,
                                        save_plot=TRUE) {
  kw_result <- alpha_KW_test(plot_data, alpha_var, meta_var, filepath)

  kw_result$subtitle <- ""
  if (kw_result$p.value < 0.05 && n_distinct(plot_data[[meta_var]]) > 2) {
    print(glue("Signif KW test on factor {meta_var} by index {alpha_var}. P.val {kw_result$p.value}"))
    dt_subtitle <- run_dunn_test(plot_data, alpha_var, meta_var, filepath)
    kw_result$subtitle <- dt_subtitle
  }

  subtitle <- glue("{subtitle}\nKW-test p-val: {signif(kw_result$p.value, 3)}\n{kw_result$subtitle}")
  plot <- plot_alpha(plot_data, alpha_var, meta_var, subtitle)
  if (save_plot) {
    ggsave(plot = plot,
           filename = glue("{outdir}/alpha_div_{meta_var}_{alpha_var}.png"),
           bg = "white",
           width = 5.6)
  }
  list(kw_result = kw_result,
       boxplot = plot)
}

# Runs lm on a given alpha diversity measurement by time, appends it to text
# output file, and returns the lm summary object
get_lm <- function(plot_data, alpha_var, filepath, subtitle){
  lm_fit <- lm(formula = glue("{alpha_var} ~ month_continuous"),
               data = plot_data) |> summary()
  sink(file = filepath, append = TRUE)
  print(glue("Linear fit of {alpha_var} over time\n{subtitle}"))
  print(lm_fit)
  sink()
  lm_fit
}

# Linear model and scatter plot with fitted line for factor Month. Returns a
# named list with the linear model object and ggplot
linear_fit_alpha_diversity <- function(plot_data,
                                       alpha_var,
                                       filepath,
                                       subtitle,
                                       outdir,
                                       filename_append,
                                       save_plot=TRUE) {
  lm_fit <- get_lm(plot_data, alpha_var, filepath, subtitle)

  subtitle <- glue("{subtitle}\nR2: {signif(lm_fit$r.squared,4)} P.value: {signif(lm_fit$coefficients[2,4],4)}")

  plot <- plot_alpha(plot_data, alpha_var, "month_continuous", subtitle)
  if (save_plot) {
    ggsave(plot = plot,
           filename = glue("{outdir}/linear_fit_month_{filename_append}.png"),
           bg = "white",
           width = 5.6)
  }
  list(lm_fit = lm_fit,
       scatter_plot = plot)
}

# Calculate alpha diversity of a phyloseq object, merge results with sample info
# (metadata) table, save to file and return the merged table
calc_alpha_diversity <- function(pseq, outdir) {
  alpha_div_data <- estimate_richness(pseq) |>
    rownames_to_column("sample") |>
    mutate(sample = str_replace(sample, "\\.", "-"))

  metadata <- sample_data(pseq) |> as_tibble(rownames="sample")
  alpha_div_data <- merge(alpha_div_data, metadata, by.x = "sample", by.y = "sample")
  utils::write.csv(alpha_div_data,
                   file = glue("{outdir}/alpha_diversity_data.csv"),
                   row.names = F)
  alpha_div_data
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
make_all_alpha_plots <- function(pseq, outdir, subtitle) {
  alpha_div_data <- calc_alpha_diversity(pseq, outdir)
  alpha_vars <- list("Observed", "Chao1", "Shannon", "InvSimpson")

  make_month_plots(alpha_div_data, alpha_vars, subtitle, outdir)
  make_process_plots(alpha_div_data, alpha_vars, subtitle, outdir)
}

# Wrapper to run linear modelling and plotting for the (continuous) Month
# metadata factor for each given alpha div measurement. Applies
# analysis/plotting to specific subsets in accordance with PERMANOVA subsetting.
# Creates saves a plot for each alpha div measurement, faceted to include all
# subsets.
make_month_plots <- function(alpha_div_data, alpha_vars, subtitle, outdir) {
  lm_month_filepath <- glue("{outdir}/linear_fit_month.txt")
  file.create(lm_month_filepath)

  plot_data <- filter(alpha_div_data, Process == "Ground")

  for (av in alpha_vars) {
    continuous_month_plots <- list()
    continuous_month <- linear_fit_alpha_diversity(plot_data, av,
                                                   lm_month_filepath,
                                                   subtitle,
                                                   outdir,
                                                   filename_append = av,
                                                   save_plot = FALSE)
    continuous_month_plots[[av]] <- continuous_month$scatter_plot

    depth_tables <- split(plot_data, plot_data$Depth)
    for (d_name in names(depth_tables)) {
      cur_depth_table <- depth_tables[[d_name]]
      subtitle_depth <- glue("{subtitle}\nSoil depth: {d_name}")
      continuous_month <- linear_fit_alpha_diversity(cur_depth_table, av,
                                                     lm_month_filepath,
                                                     subtitle_depth,
                                                     outdir,
                                                     filename_append = glue("{av}_{d_name}"),
                                                     save_plot = FALSE)
      continuous_month_plots[[glue("{av}_{d_name}")]] <- continuous_month$scatter_plot

      # split shallow soils by season
      if (d_name == "0-15") {
        season_tables <- split(cur_depth_table, cur_depth_table$Season)
        for (s_name in names(season_tables)) {
          subtitle_season <- glue("{subtitle_depth} Season: {s_name}")
          continuous_month <- linear_fit_alpha_diversity(season_tables[[s_name]], av,
                                                         lm_month_filepath,
                                                         subtitle_season,
                                                         outdir,
                                                         filename_append = glue("{av}_{d_name}_{s_name}"),
                                                         save_plot = FALSE)
          continuous_month_plots[[glue("{av}_{d_name}_{s_name}")]] <- continuous_month$scatter_plot
        }
      } else if (d_name == "15-30") {
        # split deep soils by POS
        POS_tables <- split(cur_depth_table, cur_depth_table$POS)
        for (p_name in names(POS_tables)) {
          subtitle_POS <- glue("{subtitle_depth} POS: {p_name}")
          continuous_month <- linear_fit_alpha_diversity(POS_tables[[p_name]], av,
                                                         lm_month_filepath,
                                                         subtitle_POS,
                                                         outdir,
                                                         filename_append = glue("{av}_{d_name}_{p_name}"),
                                                         save_plot = FALSE)
          continuous_month_plots[[glue("{av}_{d_name}_{p_name}")]] <- continuous_month$scatter_plot
        }
      }
    }
    facet_scatterplot <- gridExtra::grid.arrange(grobs=continuous_month_plots)
    ggsave(plot = facet_scatterplot, filename = glue("{outdir}/alpha_div_Month_{av}.png"), bg = "white",
           width = 15, height = 15)
  }
}

# Wrapper to run stat testing and plotting for the (categorical) Process
# metadata factor for each given alpha div measurement. Creates saves a plot for
# each alpha div measurement, faceted to include all subsets.
make_process_plots <- function(alpha_div_data, alpha_vars, subtitle, outdir) {
  alpha_stats_filepath <- glue("{outdir}/alpha_div_stat_tests.txt")
  file.create(alpha_stats_filepath)
  plot_data <- filter(alpha_div_data, Month == "0")

  box_process_plots <- list()
  for (av in alpha_vars) {
    # Box plot and KW test
    categorical <- categorical_alpha_diversity(plot_data, av, "Process",
                                               alpha_stats_filepath, subtitle,
                                               outdir,
                                               save_plot = FALSE)
    box_process_plots[[av]] <- categorical$boxplot
  }

  facet_scatterplot <- gridExtra::grid.arrange(grobs=box_process_plots)
  ggsave(plot = facet_scatterplot, filename = glue("{outdir}/alpha_div_Process.png"), bg = "white",
         width = 9, height = 9)
}
