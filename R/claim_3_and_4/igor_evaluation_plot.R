# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function(y, x, combination) {
  newY <- melt(y[order(y$row_id), ][, 2:3], variable.name = 'type', value.name = 'Y')
  newX <- melt(x[order(x$row_id), ][, 2:3], variable.name = 'type', value.name = 'X')
  model <- na.omit(cbind(newY, newX['X']))
  model$type <- as.character(model$type)
  model$type[model$type == 'nt_pgen_estimate'] <- 'NT'
  model$type[model$type == 'aa_pgen_estimate'] <- 'AA'
  model$type <- as.factor(model$type)
  model$corr.NT[model$type == 'NT'] <- c(cor(model[model$type == 'NT', 2], model[model$type == 'NT', 3], method = 'spearman'), rep(NA, nrow(model[model$type == 'NT', ]) - 1))
  model$corr.AA[model$type == 'AA'] <- c(cor(model[model$type == 'AA', 2], model[model$type == 'AA', 3], method = 'spearman'), rep(NA, nrow(model[model$type == 'AA', ]) - 1))
  model$combination <- combination
  names(model)[1] <- "Sequence type"
  model <- model[sample(nrow(model)), ]
  return (model)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'Pgen score'
plot_x <- 'Pgen score'
output_filename <- '~/Downloads/claim_3_and_4/igor_evaluation_plot.png'

# ----------
# MODEL DATA
# ----------
cc <- c(NA, "NULL", NA, rep("NULL", 4), rep(NA, 2))

model_dejong <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/dejong/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_emerson <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/emerson/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_peakman <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/peakman/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_igor <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/igor/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_b_i <- process_model(model_dejong, model_igor, 'project 1 (left) - IGoR (bottom)')
model_e_i <- process_model(model_emerson, model_igor, 'project 2 (left) - IGoR (bottom)')
model_p_i <- process_model(model_peakman, model_igor, 'project 3 (left) - IGoR (bottom)')
models <- as.data.frame(do.call("rbind", list(model_b_i, model_e_i, model_p_i)))
rm(cc, model_dejong, model_emerson, model_peakman, model_igor, model_b_i, model_e_i, model_p_i)

# ----------
# MAKING THE PLOTS
# ----------
eval_compare <-
  ggplot(
    data = models,
    aes(
      x = X,
      y = Y,
      color = `Sequence type`,
      linetype = `Sequence type`
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = 'gray',
    size = 0.2
  ) +
  geom_point(
    size = 1,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  geom_smooth(
    size = 1.2,
    method = lm,
    se = FALSE,
    fullrange = TRUE
  ) +
  geom_label(
    aes(
      x = 0,
      y = ((max(models$Y) / 5) * 4.9),
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
      y = ((max(models$Y) / 5) * 3.3),
      label = round(corr.AA, digits = 4),
      col = 'AA'
    ),
    hjust = 0,
    vjust = 1,
    size = 6,
    fontface = "bold",
    show.legend = FALSE
  ) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(
      size = 18
    ),
    legend.text = element_text(
      size = 15
    ),
    axis.title.x = element_text(
      size = 18
    ),
    axis.text.x = element_text(
      size = 13,
      angle = 60,
      hjust = 1
    ),
    axis.title.y = element_text(
      size = 18
    ),
    axis.text.y = element_text(
      size = 13
    ),
    strip.text.x = element_text(
      size = 15
    ),
    strip.text.y = element_text(
      size = 15
    ),
    legend.position = 'top',
    legend.direction = 'horizontal'
  ) +
  scale_x_sqrt(
    name = plot_x,
    limits = c(0, max(pmax(models$Y, models$X)))
  ) +
  scale_y_sqrt(
    name = plot_y,
    limits = c(0, max(pmax(models$Y, models$X)))
  ) +
  scale_color_manual(
    values = c('#ca0020', '#0571b0')
  ) +
  facet_wrap(
    vars(combination)
  )

jpeg(output_filename, width = 4000, height = 2000, res = 300)
eval_compare
dev.off()
