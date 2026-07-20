# produces figure for comparing first intron pair splicing order distribution to all other introns

library(dplyr)
library(ggplot2)
library(readr)

# pooled splicing order data with one row per intron pair, including coordinates and fraction_downstream
df <- read_tsv("merged/0407/splicing_order_pooled.tsv")

df <- df %>%
  group_by(gene_id) %>%
  mutate(is_first_pair = (intron1_start == min(intron1_start))) %>%
  ungroup()

# split distributions into first intron pairs vs all other pairs and calculate summary statistics
summary_by_position <- df %>%
  group_by(is_first_pair) %>%
  summarise(
    n             = n(),
    mean_fd       = mean(fraction_downstream, na.rm = TRUE),
    sem           = sd(fraction_downstream, na.rm = TRUE) / sqrt(n()),
    pct_above_0.5 = mean(fraction_downstream > 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(is_first_pair, "First intron pair", "All other pairs"))

print(summary_by_position)

# create figure

df <- df %>%
  mutate(
    pair_type = ifelse(is_first_pair, "First intron pair", "All other pairs"),
    pair_type = factor(pair_type, levels = c("First intron pair", "All other pairs"))
)
mean_lines <- data.frame(
  mean_fd    = c(0.544, 0.497),
  pair_type  = c("First intron pair", "All other pairs")
)

p <- ggplot(df, aes(x = fraction_downstream, fill = pair_type)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins     = 80,
    alpha    = 0.5,
    position = "identity"
  ) +
  geom_vline(
    data      = mean_lines,
    aes(xintercept = mean_fd, color = pair_type),
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 0.5,
             linetype   = "dashed",
             color      = "black",
             linewidth  = 0.8) +
  scale_fill_manual(
    values = c("First intron pair" = "#E05C5C",
               "All other pairs"   = "#4B8BBE"),
    name   = NULL
  ) +
  scale_color_manual(
    values = c("First intron pair" = "#E05C5C",
               "All other pairs"   = "#4B8BBE"),
    name   = NULL
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "Fraction Downstream Spliced First",
    y = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position  = c(0.02, 0.98),
        legend.justification = c(0, 1))

ggsave("figures/first_intron_check.pdf", p, width = 7, height = 7)
ggsave("figures/first_intron_check.png", p, width = 7, height = 7,
       dpi = 200, bg = "white")

# one-sample t-test comparing mean fraction_downstream for first intron pairs to 0.5
no_first <- df %>% filter(!is_first_pair) %>% pull(fraction_downstream)
t_result <- t.test(no_first, mu = 0.5, alternative = "greater")
cat(sprintf(
  "\nWithout first intron pairs:\n  n = %s\n  mean = %.4f\n  p (vs 0.5) = %.2e\n",
  format(length(no_first), big.mark = ","),
  mean(no_first),
  t_result$p.value
))