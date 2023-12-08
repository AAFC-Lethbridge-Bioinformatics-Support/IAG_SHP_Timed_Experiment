library(tidyverse)
library(vegan)
library(ggplot2)
library(glue)
library(dunn.test)

calc_diversity_df <- function(d){
  sample_col <- select(d, "sample")
  sample_data <- select(d, -"sample")
  observed_richness <- specnumber(sample_data)
  invsimpson <- diversity(sample_data, index="invsimpson")
  shannon <- diversity(sample_data, index="shannon")
  div_df <- data.frame(
    ID = sample_col,
    Observed_Richness = observed_richness,
    Inv_Simpson = invsimpson,
    Shannon = shannon
  )
  return(div_df)
}

# If the global KW test is sig, run dunn.test and report the sig pairwise results
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

# runs and saves kruskal wallis on shannon index for both
# replicates and treatments
alpha_KW_test <- function(d, alpha_index, meta_factor, outfile) {

  heading <- glue("{alpha_index} Index by {meta_factor} - Kruskal Wallis Results:\n")
  factor_alpha <- stats::kruskal.test(stats::formula(glue("{alpha_index}~{meta_factor}")), d)
  cat(heading, file = outfile, append = T)
  utils::capture.output(factor_alpha, file = outfile, append = T)

  factor_alpha$subtitle <- ""
  if (factor_alpha$p.value < 0.05 && n_distinct(d[[meta_factor]]) > 2) {
    print(glue("Signif KW test on factor {meta_factor} by index {alpha_index}. P.val {factor_alpha$p.value}"))
    dt_subtitle <- run_dunn_test(d, alpha_index, meta_factor, outfile)
    factor_alpha$subtitle <- dt_subtitle
  }

  factor_alpha
}

plot_alpha <- function(d, alpha, x_var, subtitle) {
  isContinuous <- is.numeric(d[[x_var]])

  index <- sym(alpha)
  cont <- sym(x_var)
  alpha_plot <- ggplot(d, aes(x = !!cont, y = !!index)) +
    labs(title = glue("Alpha Diversity: {alpha}"),
         subtitle = subtitle,
         x = x_var) +
    theme_bw()
  if (isContinuous) {
    message(glue("{x_var} is being treated as a continuous variable for plotting"))
    alpha_plot <- alpha_plot +
      geom_point(aes(color = Season), alpha = 0.7) +
      geom_smooth(aes(group=Season, color = Season),method="lm",se=FALSE) +
      geom_smooth(color = "black", method = "lm", se=FALSE)
  } else {
    message(glue("{x_var} is being treated as a categorical variable for plotting"))
    alpha_plot <- alpha_plot+
      geom_boxplot() +
      geom_jitter(aes(color = Season), alpha = 0.7)
  }
  alpha_plot
}

# KW and Dunn tests and grouped boxplot
categorical_alpha_diversity <- function(plot_data,
                                        alpha_var,
                                        meta_var,
                                        filepath,
                                        subtitle,
                                        outdir) {
  kw_result <- alpha_KW_test(plot_data, alpha_var, meta_var, filepath)

  subtitle <- glue("{subtitle}\nKW-test p-val: {signif(kw_result$p.value, 3)}\n{kw_result$subtitle}")
  plot <- plot_alpha(plot_data, alpha_var, meta_var, subtitle)
  ggsave(plot = plot,
         filename = glue("{outdir}/alpha_div_{meta_var}_{alpha_var}.png"),
         bg = "white",
         width = 5.6)

  list(kw_result = kw_result,
       boxplot = plot)
}

# Run lm on Month 3 times: Spring, Fall, and Both Seasons
run_month_lm_seasons <- function(plot_data, alpha_var, filepath){
  lm_fit <- lm(formula = glue("{alpha_var} ~ month_continuous"),
               data = plot_data) |> summary()
  sink(file = filepath, append = TRUE)
  print(glue("Linear fit of {alpha_var} over time\n"))
  print(lm_fit)
  sink()

  lm_fit_spring <- lm(formula = glue("{alpha_var} ~ month_continuous"),
                      data = filter(plot_data), Season == "Spring") |> summary()
  sink(file = filepath, append = TRUE)
  print(glue("Season==Spring Linear fit of {alpha_var} over time\n"))
  print(lm_fit_spring)
  sink()
  lm_fit_fall <- lm(formula = glue("{alpha_var} ~ month_continuous"),
                    data = filter(plot_data), Season == "Fall") |> summary()
  sink(file = filepath, append = TRUE)
  print(glue("Season==Fall Linear fit of {alpha_var} over time\n"))
  print(lm_fit_fall)
  sink()

  list(all = lm_fit,
       spring = lm_fit_spring,
       fall = lm_fit_fall)
}

# Linear model and scatter plot with fitted line for factor Month
linear_fit_alpha_diversity <- function(plot_data,
                                       alpha_var,
                                       filepath,
                                       subtitle,
                                       outdir) {
  lm_fits <- run_month_lm_seasons(plot_data, alpha_var, filepath)

  for (fit_name in names(lm_fits)) {
    fit <- lm_fits[[fit_name]]
    subtitle <- glue("{subtitle}\nSeason: {fit_name} R2: {signif(fit$r.squared,4)} P.value: {signif(fit$coefficients[2,4],4)}")
  }

  plot <- plot_alpha(plot_data, alpha_var, "month_continuous", subtitle)
  ggsave(plot = plot,
         filename = glue("{outdir}/linear_fit_month_{alpha_var}.png"),
         bg = "white",
         width = 5.6)
  list(lm_fits = lm_fits,
       scatter_plot = plot)
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
make_all_alpha_plots <- function(data, outdir, subtitle) {
  alpha_stats_filepath <- glue("{outdir}/alpha_div_stat_tests.txt")
  file.create(alpha_stats_filepath)
  lm_month_filepath <- glue("{outdir}/linear_fit_month.txt")
  file.create(lm_month_filepath)

  alpha_div_data <- calc_diversity_df(data)
  alpha_div_data <- merge(alpha_div_data, timed_meta, by.x = "sample", by.y = "sample")
  utils::write.csv(alpha_div_data,
                   file = glue("{outdir}/alpha_diversity_data.csv"),
                   row.names = F)

  continuous_month_plots <- list()
  box_process_plots <- list()

  meta_vars <- list("Month", "Process")
  alpha_vars <- list("Observed_Richness", "Shannon", "Inv_Simpson")
  for (mv in meta_vars){
    if (mv == "Month") {
      plot_data <- filter(alpha_div_data, Process == "Ground")
    } else if (mv == "Process") {
      plot_data <- filter(alpha_div_data, Month == "0")
    } else {
      plot_data <- alpha_div_data
    }
    for (av in alpha_vars) {
      if (mv == "Month"){
        continuous_month <- linear_fit_alpha_diversity(plot_data, av,
                                                       lm_month_filepath,
                                                       subtitle,
                                                       outdir)
        continuous_month_plots[[av]] <- continuous_month$scatter_plot
        next
      } else {
        # Box plot and KW test
        categorical <- categorical_alpha_diversity(plot_data, av, mv,
                                                   alpha_stats_filepath, subtitle,
                                                   outdir)
        box_process_plots[[av]] <- categorical$boxplot
      }
    }

    if (mv == "Month") {
      facet_scatterplot <- gridExtra::grid.arrange(grobs=continuous_month_plots, ncol=3)
      ggsave(plot = facet_scatterplot, filename = glue("{outdir}/alpha_div_Month_linear.png"), bg = "white",
             width = 15, height = 5.0)
    } else if (mv == "Process") {
      facet_scatterplot <- gridExtra::grid.arrange(grobs=box_process_plots, ncol=3)
      ggsave(plot = facet_scatterplot, filename = glue("{outdir}/alpha_div_Process.png"), bg = "white",
             width = 15, height = 5.0)
    }
  }
}
