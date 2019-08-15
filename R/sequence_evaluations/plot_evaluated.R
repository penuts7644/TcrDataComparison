library(GGally)

# ----------
# READ IN DATA FOR GIVEN PRODUCTIVE, UNPRODUCTIVE OR ALL MODELS, 'plot_evaluated_variables.R'
# ----------

upper_plot_fn <- function(data, mapping, ...){
  ggally_cor(
    data = data,
    mapping = mapping,
    method = 'spearman',
    use = 'pairwise',
    alignPercent = 1,
    size = 6
  ) +
  theme(
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(
      linetype = "dashed",
      colour = "gray",
      fill = NA
    )
  ) +
  scale_color_manual(
    values = c('#1f78b4', '#d95f02')
  )
}

lower_plot_fn <- function(data, mapping, ...){
  ggplot(
    data = data,
    mapping = mapping
  ) +
  geom_point(
    size = 2,
    alpha = 0.4
  ) +
  geom_smooth(
    size = 1,
    method = lm,
    se = FALSE,
    fullrange = TRUE
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  scale_color_manual(
    values = c('#1f78b4', '#d95f02')
  )
}

legend_plot <- ggplot(data = models) +
  geom_point(
    aes(x = `0`, y = `1`, color = type)
  ) +
  theme(
    legend.text = element_text(
      size = 15
    ),
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal'
  ) +
  scale_color_manual(
    values = c('#1f78b4', '#d95f02')
  )

eval_compare <-
  ggpairs(
    models,
    columns = 2:ncol(models),
    mapping = aes(color = type),
    axisLabels = 'show',
    upper = list(continuous = upper_plot_fn),
    diag = list(continuous = 'blankDiag'),
    lower = list(continuous = lower_plot_fn),
    progress = TRUE,
    legend = grab_legend(legend_plot)
  ) +
  theme(
    plot.title = element_text(
      size = 20
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 15
    ),
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 11
    ),
    strip.text.x = element_text(
      size = 15
    ),
    strip.text.y = element_text(
      size = 15
    ),
    legend.position = 'bottom',
    strip.placement = 'outside'
  ) +
  labs(
    title = plot_title,
    caption = plot_caption
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
eval_compare
dev.off()
