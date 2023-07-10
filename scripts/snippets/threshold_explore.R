library(tidyverse)
library(magrittr)
library(vegan)
library(glue)
library(indicspecies)


# Read in files -----------------------------------------------------------

ko_long <- read_csv("data/woltka_timed_protein_ko_rpk_LONG.csv")
sample_read_counts <- read_csv("data/timed_deep_samples_read_count.csv")
timed_meta <- read_csv("../metadata/Metadata-IAG-Timed2-v3_2.csv")
# remove all NA rows, remove blanks
timed_meta <- timed_meta[rowSums(is.na(timed_meta)) != ncol(timed_meta),] |>
  rename(sample = `MBI_ID`) |>
  filter(!str_detect(Site, "B\\d"))

timed_meta <- timed_meta |>
  mutate(Month = case_when(Year == "1998" | Year == "2007" ~ Year,
                           TRUE ~ Month))
timed_meta$Month <- ordered(timed_meta$Month, levels = c("0",
                                                         "0.07",
                                                         "0.5",
                                                         "1",
                                                         "3",
                                                         "6",
                                                         "12",
                                                         "18",
                                                         "1998",
                                                         "2007"))
timed_meta$Year <- as.character(timed_meta$Year)

meta_w_counts <- merge(timed_meta, sample_read_counts)


# Exploration -------------------------------------------------------------

meta_w_counts |> ggplot(aes(x = total_gene_reads, y = percent_mapped_to_gene)) +
  geom_point(aes(color = Depth)) +
  geom_smooth(method = "lm", aes(color = Depth), se=FALSE) +
  scale_color_manual(values=color_palette) +
  scale_fill_manual(values=color_palette)

tm <- ko_long |>
  summarise(avg_gene_count = mean(reads),
            total_gene_count = sum(reads),
            .by = ko_gene_ID)

threshold_gene_nums <- tibble(
  # thresh = c(5, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
  thresh = seq(1,100, by=5)
) |>
  mutate(n_genes_less_than = sapply(thresh, function(x) {
    tmp <- filter(tm, total_gene_count < x)
    nrow(tmp)
  }))

ggplot(threshold_gene_nums, aes(x=as.factor(thresh), y = n_genes_less_than)) +
  geom_col()


ko_long |> summarise(avg_gene_count = mean(reads),
                     total_gene_count = sum(reads),
                     .by = ko_gene_ID)

threshold_gene_nums <- tibble(
    # thresh = c(5, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
    thresh = seq(1,10000, by=100)
  ) |>
  mutate(n_genes_less_than = sapply(thresh, function(x) {
    tmp <- filter(tm, avg_gene_count < x)
    nrow(tmp)
  }),
  percent_genes_below = n_genes_less_than/7422)
threshold_gene_nums$n_genes_less_than

ggplot(threshold_gene_nums, aes(x=as.factor(thresh), y = n_genes_less_than/7422)) +
  geom_col()

tm |> ggplot(aes(x = "", y = log(avg_gene_count))) +
  geom_boxplot() +
  geom_jitter()

tm |> count(avg_gene_count<50)
