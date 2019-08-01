library(ggplot2)
library(GGally)
library(reshape2)
library(gtable)
library(gridExtra)

# ----------
# READ IN DATA FOR GIVEN PRODUCTIVE, UNPRODUCTIVE OR ALL MODELS, 'plot_evaluated_variables.R'
# ----------

# Combine the model dataframes together on the sequence index and type
models <- na.omit(do.call('cbind', list(
  model_0, model_1[2], model_2[2], model_3[2],
  model_4[2], model_combined[2], model_default[2]
)))
rm(model_0, model_1, model_2, model_3, model_4, model_combined, model_default)

# Rename the sequence type column and values
names(models) <- c('type', '0', '1', '2', '3', '4', 'combined', 'default')
models$type <- as.character(models$type)
models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
models$type <- as.factor(models$type)

# Apply simple min-max normalize to scale NT and AA values accordingly
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
models[models$type == 'NT', -1] <- apply(models[models$type == 'NT', -1], 2, normalize)
models[models$type == 'AA', -1] <- apply(models[models$type == 'AA', -1], 2, normalize)

# Function for plotting the upper half of the graph
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
    values = c('#af8dc3', '#7fbf7b'), # PRGn
  )
}

# Function for plotting the lower half of the graph
lower_plot_fn <- function(data, mapping, ...){
  ggplot(
    data = data,
    mapping = mapping
  ) +
  geom_point(
    size = 1,
    alpha = 0.4,
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
    values = c('#af8dc3', '#7fbf7b') # PRGn
  )
}

# Plot one graph to collect the legend from
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
    values = c('#af8dc3', '#7fbf7b'), # PRGn
  )

# Create the ggpairs combined graph and print
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
    strip.placement = 'outside'
  ) +
  labs(
    title = plot_title,
    caption = plot_caption
  )

# Write out our plot png
jpeg(output_filename, width = 4000, height = 4000, res = 300)
eval_compare
dev.off()
