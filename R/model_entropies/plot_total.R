library(ggplot2)

# ----------
# READ IN DATA FOR PRODUCTIVE, UNPRODUCTIVE AND ALL MODELS, 'plot_total_variables.R'
# ----------

entr_compare <-
  ggplot(
    models,
    aes(
      x = L1,
      y = value
    )
  ) +
  geom_boxplot(
    alpha = 0.5,
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    height = 0,
    width = .2,
    size = 4,
    aes(color = L2),
    alpha = 0.8
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
    axis.text.x = element_text(
      size = 15,
    ),
    axis.title.y = element_text(
      size = 15
    ),
    axis.text.y=element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  ylim(0, 1) +
  labs(
    title = plot_title,
    caption = plot_caption,
    y = plot_y
  ) +
  scale_color_manual(
    values=c('#af8dc3', '#7fbf7b') # PRGn
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
entr_compare
dev.off()
