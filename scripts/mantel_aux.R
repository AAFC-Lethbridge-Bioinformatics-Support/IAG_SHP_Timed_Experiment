library(tidyverse)
library(vegan)
library(glue)

make_distance_lm_plot <- function(dist.abund, dist.month, dist_method,
                                  subtitle="", log_transform=FALSE){
  y_label <- if (dist_method == "bray") {
    "Bray-Curtis Dissimilarity"
  } else {
    "Jaccard Dissimilarity"
  }
  point_alpha <- if (seq_depth == "deep") 0.03 else 0.05

  mm <- data.frame(abundance_dist = as.vector(dist.abund),
                   month_dist = as.vector(dist.month)) |>
    mutate(avg_dist_by_month = ifelse(log_transform,
                                      mean(log(abundance_dist)),
                                      mean(abundance_dist)),
           .by = month_dist)
  g <- mm |> ggplot(aes(x = month_dist))
  if (log_transform) {
    g <- mm |> ggplot(aes(y = log(abundance_dist), x = month_dist)) +
      geom_point(size = 3, alpha = point_alpha) +
      geom_smooth(method = "lm", colour = "blue")
  } else {
    g <- mm |> ggplot(aes(y = abundance_dist, x = month_dist)) +
      geom_point(size = 3, alpha = point_alpha) +
      geom_smooth(method = "lm", colour = "blue")
  }

  g <- g +
    geom_point(aes(y = avg_dist_by_month), color = "red") +
    labs(x = "Difference in sample time (month)",
         y = y_label,
         title = "Sample dissimilarity over time",
         subtitle = subtitle) +
    theme(axis.text.x = element_text(face = "bold",colour = "black", size = 12),
          axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
          axis.title= element_text(face = "bold", size = 14, colour = "black"),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"))
}
