library(ggplot2)

# ----------
# READ IN DATA FOR THE SEQUENCE COUNTS, 'plot_seq_count_variables.R'
# ----------

size_compare <-
  ggplot(
    counts,
    aes(
      x = model,
      y = norm,
      fill = model
    )
  ) +
  geom_bar(
    stat = 'identity',
    width = 0.8,
    position = 'dodge'
  ) +
  geom_text(
    aes(
      label = count
    ),
    vjust = -0.4,
    color = 'black',
    position = position_dodge(0.9),
    size = 5
  ) +
  theme(
    plot.title = element_text(
      size = 20
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 15
    ),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(
      size = 15
    ),
    strip.text.x = element_text(
      size = 15
    ),
    strip.text.y = element_text(
      size = 15
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  ylim(0, 1) +
  labs(
    title = plot_title,
    caption = plot_caption
  ) +
  scale_fill_brewer(
    palette = 'GnBu'
  ) +
  facet_grid(
    rows = vars(id_f),
    cols = vars(type)
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
size_compare
dev.off()
