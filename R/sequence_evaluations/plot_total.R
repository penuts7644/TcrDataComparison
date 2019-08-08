library(ggplot2)

# ----------
# READ IN DATA FOR PRODUCTIVE, UNPRODUCTIVE AND ALL MODELS, 'plot_total_variables.R'
# ----------

eval_compare <-
  ggplot(
    data = models,
    aes(
      x = X,
      y = Y,
      color = type
    )
  ) +
  geom_point(
    size = 2,
    alpha = 0.4,
  ) +
  geom_smooth(
    size = 1,
    method = lm,
    se = FALSE,
    fullrange = TRUE,
    show.legend = FALSE
  ) +
  geom_label(
    aes(
      x = 0,
      y = 1,
      label = round(corr.NT, digits = 4),
      col = 'NT'
    ),
    hjust = 0,
    vjust = 1,
    size = 6,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_label(
    aes(
      x = 0,
      y = 0.8,
      label = round(corr.AA, digits = 4),
      col = 'AA'
    ),
    hjust = 0,
    vjust = 1,
    size = 6,
    fontface = "bold",
    show.legend = FALSE
  ) +
  theme(
    plot.title = element_text(
      size = 20
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 15
    ),
    legend.text = element_text(
      size = 15
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1
    ),
    axis.text.y=element_text(
      size = 11
    ),
    strip.text.x = element_text(
      size = 15
    ),
    strip.text.y = element_text(
      size = 15
    ),
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal'
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(
    title = plot_title,
    caption = plot_caption
  ) +
  scale_color_manual(
    values = c('#1f78b4', '#d95f02')
  ) +
  facet_grid(
    rows = vars(name),
    cols = vars(combination)
  )

jpeg(output_filename, width = 3000, height = 4000, res = 300)
eval_compare
dev.off()
