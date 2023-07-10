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

# how many NA's are there in each column?
sapply(timed_meta, function(x) sum(length(which(is.na(x)))))

# ---- Read in taxonomy data ----
taxa_level <- "phylum"
taxa_long <- read_tsv(glue("pipeline_outputs/deep_Timed/kaiju_summaries/kaiju_{taxa_level}_summary.tsv"))
taxa_long %<>% mutate(sample = str_extract(file, "S00JY-\\d{4}")) %>%
  filter(sample %in% timed_meta$`MBI_ID`)

taxa_read_avg <- taxa_long %>% group_by(taxon_name) %>%
  summarize(avg_reads = mean(reads)) %>%
  arrange(-avg_reads)
taxa_pc_avg <- taxa_long %>% group_by(taxon_name) %>%
  summarize(avg_percent = mean(percent)) %>%
  arrange(-avg_percent)

hist(taxa_read_avg$avg_reads,
     breaks = 50,
     main = glue("Histogram of taxa at the {taxa_level} level by average read count"),
     xlab = "Read count buckets")
hist(taxa_read_avg$avg_reads[taxa_read_avg$avg_reads < 4000],
     breaks = 50,
     main = glue("Histogram of taxa at the {taxa_level} level by average read count"),
     xlab = "Read count buckets")
hist(taxa_read_avg$avg_reads[taxa_read_avg$avg_reads < 100],
     breaks = 50,
     main = glue("Histogram of taxa at the {taxa_level} level by average read count"),
     xlab = "Read count buckets")
hist(taxa_read_avg$avg_reads[taxa_read_avg$avg_reads < 50],
     breaks = 25,
     main = glue("Histogram of taxa at the {taxa_level} level by average read count"),
     xlab = "Read count buckets")
hist(taxa_read_avg$avg_reads[taxa_read_avg$avg_reads < 10],
     breaks = 10,
     main = glue("Histogram of taxa at the {taxa_level} level by average read count"),
     xlab = "Read count buckets")


# Threshold for phylum at 50 reads?
# Threshold for genus at 10 reads?
