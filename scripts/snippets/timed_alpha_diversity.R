library(tidyverse)
library(magrittr)
library(vegan)
library(glue)

color_palette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921",
                   "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", "#D7C1B1",
                   "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0",
                   "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#38333E", "#599861")

# ---- Read in and tidy timed metadata ----
timed_meta <- read_csv("metadata/Metadata-IAG-Timed2-v3_1.csv")

# Metadata corrections
timed_meta %<>% mutate(
  Month = case_when(Year == "1998" | Year == "2007" ~ Year,
                    TRUE ~ Month)
)

# how many NA's are there in each column?
sapply(timed_meta, function(x) sum(length(which(is.na(x)))))

# ---- Read in taxonomy data ----
taxa_level <- "phylum"
threshold <- ifelse(taxa_level == "phylum",
                    50,
                    ifelse(taxa_level == "genus",
                           10,
                           NULL))
taxa_long <- read_tsv(glue("pipeline_outputs/deep_Timed/kaiju_{taxa_level}_summary.tsv"))
taxa_long %<>% mutate(sample = str_extract(file, "S00JY-\\d{4}")) %>%
  filter(sample %in% timed_meta$`MBI_ID`)
# How many unique taxa?
length(unique(taxa_long$taxon_name))

taxa_long_reads <- taxa_long %>% select(-c(file, percent, taxon_id))

# replace NA reads with 0
taxa_long_reads$reads[is.na(taxa_long_reads$reads)] <- 0

if (!is.na(threshold)) {
  # Filter out counts below threshold
  taxa_long_reads %<>% filter(reads >= threshold)
}
# How many unique taxa above cutoff?
length(unique(taxa_long_reads$taxon_name))

taxa_wide_reads <- taxa_long_reads %>%
  pivot_wider(names_from = taxon_name, values_from = reads)

taxa_wide_reads[is.na(taxa_wide_reads)] <- 0


# Alpha Diversity ---------------------------------------------------------
# calc alpha metrics
# assumes data in tidy format and no extra columns
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

# runs and saves kruskal wallis on shannon index for both
# replicates and treatments
alpha_KW_test <- function(d, alpha_index, meta_factor, outfile) {
  heading <- glue("{alpha_index} Index by {meta_factor} - Kruskal Wallis Results:\n")
  factor_alpha <- stats::kruskal.test(stats::formula(glue("{alpha_index}~{meta_factor}")), d)
  cat(heading, file = outfile, append = T)
  utils::capture.output(factor_alpha, file = outfile, append = T)
}

# plot alpha diversity
plot_alpha <- function(d, alpha, group_var) {
  index <- sym(alpha)
  grp <- sym(group_var)
  alpha_plot <- ggplot(d, aes(x = !!grp, y = !!index)) +
    geom_boxplot() +
    geom_jitter(aes(color = Season), alpha = 0.7) +
    labs(title = "Alpha Diversity",
         subtitle = glue("Metric: {alpha}"),
         x = group_var)
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
make_all_alpha_plots <- function(data, taxa_rank) {
  outdir <- glue("plots/timed/alpha_diversity/{taxa_rank}/")

  alpha_div_data <- calc_diversity_df(data)
  alpha_div_data <- merge(alpha_div_data, timed_meta, by.x = "sample", by.y = "MBI_ID")
  alpha_div_data$Month <- factor(alpha_div_data$Month,
                                 levels = c("0", "0.5", "0.07", "1", "3", "6", "12", "18", "1998", "2007"))
  # filter out Blanks
  alpha_div_data <- filter(alpha_div_data, Month != "Blank")

  utils::write.csv(alpha_div_data,
                   file = glue("{outdir}/alpha_div_data_{taxa_rank}.csv"),
                   row.names = F)


  alpha_stats_filepath <- glue("{outdir}/alpha_stats.txt")
  file.create(alpha_stats_filepath)

  meta_vars <- list("Month", "Process")
  alpha_vars <- list("Observed_Richness", "Shannon", "Inv_Simpson")
  for (mv in meta_vars){
    for (av in alpha_vars) {
      if (mv == "Month") {
        plot_data <- filter(alpha_div_data, Process == "Ground")
      } else if (mv == "Process") {
        plot_data <- filter(alpha_div_data, Month == "0")
      } else {
        plot_data <- alpha_div_data
      }
      alpha_KW_test(plot_data, av, mv, alpha_stats_filepath)
      plot <- plot_alpha(plot_data, av, mv)
      ggsave(plot = plot, filename = glue("{outdir}/alpha_div_{mv}_{av}_{taxa_rank}.png"), bg = "white")
    }
  }
}

make_all_alpha_plots(taxa_wide_reads, taxa_level)
