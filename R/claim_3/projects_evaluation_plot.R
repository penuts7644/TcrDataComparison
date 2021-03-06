# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)
CORRELATION_METHOD <- 'pearson' # or 'spearman'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function(y, x, combination) {
  newY <- melt(y, id.vars = 'row_id', variable.name = 'type', value.name = 'Y')
  newX <- melt(x, id.vars = 'row_id', variable.name = 'type', value.name = 'X')
  model <- merge(newX, newY, by=c('row_id', 'type'))
  model$type <- as.character(model$type)
  model$type[model$type == 'nt_pgen_estimate'] <- 'NT'
  model$type[model$type == 'aa_pgen_estimate'] <- 'AA'
  model$type <- as.factor(model$type)
  model <- model[!(model$row_id %in% model[is.na(model[, 3:4]), 'row_id']), ][, 2:4]
  model$corr.NT[model$type == 'NT'] <- c(
    round(cor(model[model$type == 'NT', 2], model[model$type == 'NT', 3], method = CORRELATION_METHOD), digits = 2),
    rep(NA, nrow(model[model$type == 'NT', ]) - 1)
  )
  model$corr.AA[model$type == 'AA'] <- c(
    round(cor(model[model$type == 'AA', 2], model[model$type == 'AA', 3], method = CORRELATION_METHOD), digits = 2),
    rep(NA, nrow(model[model$type == 'AA', ]) - 1)
  )
  model$combination <- combination
  names(model)[1] <- "Sequence type"
  model <- model[sample(nrow(model)), ]
  return (model)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'Generation probability score'
plot_x <- 'Generation probability score'
output_filename <- paste('~/Downloads/claim_3/study_evaluation_plot_', CORRELATION_METHOD, '.png', sep = '')

# ----------
# MODEL DATA
# ----------
cc <- c(NA, "NULL", NA, rep("NULL", 4), rep(NA, 2))

model_emerson <- data.frame(read.table('~/Downloads/claim_3/evaluations/emerson/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_dejong <- data.frame(read.table('~/Downloads/claim_3/evaluations/dejong/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_gomez <- data.frame(read.table('~/Downloads/claim_3/evaluations/gomez/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_e_j <- process_model(model_emerson, model_dejong, 'study 1 (left) - study 2 (bottom)')
model_e_g <- process_model(model_emerson, model_gomez, 'study 1 (left) - study 3 (bottom)')
model_j_g <- process_model(model_dejong, model_gomez, 'study 2 (left) - study 3 (bottom)')
models <- as.data.frame(do.call("rbind", list(model_e_j, model_e_g, model_j_g)))
rm(cc, model_dejong, model_emerson, model_gomez, model_e_j, model_e_g, model_j_g)

# ----------
# MAKING THE PLOTS
# ----------
eval_compare <-
  ggplot(
    data = models,
    aes(
      x = X,
      y = Y,
      color = `Sequence type`
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
    alpha = 0.4
  ) +
  geom_label(
    aes(
      x = 0,
      y = ((max(models$Y) / 5) * 4.9),
      label = corr.NT,
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
      y = ((max(models$Y) / 5) * 4),
      label = corr.AA,
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
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    )
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
