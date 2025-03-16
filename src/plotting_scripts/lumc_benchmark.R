library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(forcats)
library(cowplot)

compute_and_replace_with_row_averages <- function(data) {
  rows12 <- data %>%
    filter(id %in% c(1, 2)) %>%
    group_by(wuhan_ab, ec_tool, hrt) %>%
    summarise(
      id = 12,
      f1_score = mean(f1_score),
      precision = mean(precision),
      recall = mean(recall),
      split = "combined"
    )

  rows89 <- data %>%
    filter(id %in% c(8, 9)) %>%
    group_by(wuhan_ab, ec_tool, hrt) %>%
    summarise(
      id = 89,
      f1_score = mean(f1_score),
      precision = mean(precision),
      recall = mean(recall),
      split = "combined"
    )

  data <- data %>% filter(!id %in% c(1, 2, 8, 9))
  data <- rbind(data, rows12)
  data <- rbind(data, rows89)

  return(data)
}

select_best_ec_tool_per_hrt_tool <- function(data) {

  selection <- data %>%
    group_by(hrt, ec_tool, split) %>%
    mutate(result_num = length(unique(id))) %>% 
    reframe(f1_score = mean(f1_score), result_num = max(result_num)) %>%
    group_by(hrt) %>%
    arrange(desc(result_num), desc(f1_score), .by_group = TRUE) %>%
    slice(1)

  return(selection)
}


data_file_name <- "data/results.tsv"

data <- read.csv(data_file_name, sep = "\t")

data <- data %>%
  tidyr::separate(sample_name, into = c("wuhan_percentage", "ec_tool", "hrt", "split"), sep = "-") %>%
  tidyr::separate(wuhan_percentage, into = c("id", "wuhan_ab"), sep = "_")

data$id <- as.numeric(data$id)

selection_best <- select_best_ec_tool_per_hrt_tool(data)

df_best <- data %>%
  dplyr::group_by(id, wuhan_ab, hrt) %>%
  semi_join(selection_best, by = c("ec_tool", "hrt", "split"))

df_best <- compute_and_replace_with_row_averages(df_best)

facet_labels <- c(
  "f1_score" = "F1-score",
  "recall" = "Recall",
  "precision" = "Precision"
)

df_long <- df_best %>%
  pivot_longer(
    cols = c(`f1_score`, recall, precision),
    names_to = "Metric",
    values_to = "Value"
  )

df_long$wuhan_ab <- as.numeric(df_long$wuhan_ab)

df_long <- df_long[order(df_long$wuhan_ab), ]

df_long$Metric <- factor(df_long$Metric, levels = c("f1_score", "recall", "precision"))

plot_f1 <- ggplot(df_long %>% filter(Metric == "f1_score"), aes(x = wuhan_ab, y = Value, color = hrt, group = hrt)) +
  geom_line(aes(alpha = 1), show.legend = FALSE) +
  geom_point(aes(shape = hrt), size = 2, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) + # Custom y-axis
  facet_wrap(~Metric, scales = "free_y", labeller = labeller(Metric = facet_labels)) +
  theme_minimal()

plot_recall <- ggplot(df_long %>% filter(Metric == "recall"), aes(x = wuhan_ab, y = Value, color = hrt, group = hrt)) +
  geom_line(aes(alpha = 1), show.legend = FALSE) +
  geom_point(aes(shape = hrt), size = 2, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) + # Custom y-axis
  facet_wrap(~Metric, scales = "free_y", labeller = labeller(Metric = facet_labels)) +
  theme_minimal()

plot_precision <- ggplot(df_long %>% filter(Metric == "precision"), aes(x = wuhan_ab, y = Value, color = hrt, group = hrt)) +
  geom_line(aes(alpha = 1), show.legend = FALSE) +
  geom_point(aes(shape = hrt), size = 2, alpha = 0.8) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) + # Custom y-axis
  facet_wrap(~Metric, scales = "free_y", labeller = labeller(Metric = facet_labels)) +
  theme_minimal()


combined_plot <- (plot_f1 + plot_recall + plot_precision) +
  plot_layout(
    guides = "collect"
  ) &
  plot_annotation(tag_levels = "A") &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
  )

combined_plot <- wrap_elements(panel = combined_plot) +
  labs(tag = "Wuhan percentage (%)") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )

combined_plot &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
  )
