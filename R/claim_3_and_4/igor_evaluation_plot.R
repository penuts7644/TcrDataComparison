# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(x, y, combination) {
  newX <- melt(x, variable.name = 'type', value.name = 'X')
  newY <- melt(y, variable.name = 'type', value.name = 'Y')
  model <- na.omit(cbind(newX, newY['Y']))
  model$type <- as.character(model$type)
  model$type[model$type == 'nt_pgen_estimate'] <- 'NT'
  model$type[model$type == 'aa_pgen_estimate'] <- 'AA'
  model$type <- as.factor(model$type)
  model[model$type == 'NT', -1] <- apply(model[model$type == 'NT', -1], 2, normalize)
  model[model$type == 'AA', -1] <- apply(model[model$type == 'AA', -1], 2, normalize)
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
plot_y <- 'Pgen score (normalized)'
plot_x <- 'Pgen score (normalized)'
output_filename <- '~/Downloads/claim_3_and_4/igor_evaluation_plot.png'

# ----------
# MODEL DATA
# ----------
cc <- c(NA, rep("NULL", 6), rep(NA, 2))

model_brusko <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/brusko/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_emerson <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/emerson/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_peakman <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/peakman/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_igor <- data.frame(read.table('~/Downloads/claim_3_and_4/evaluations/igor/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
model_b_i <- process_model(model_brusko, model_igor, 'project 1 (X) - IGoR (Y)')
model_e_i <- process_model(model_emerson, model_igor, 'project 2 (X) - IGoR (Y)')
model_p_i <- process_model(model_peakman, model_igor, 'project 3 (X) - IGoR (Y)')
models <- as.data.frame(do.call("rbind", list(model_b_i, model_e_i, model_p_i)))
rm(cc, model_brusko, model_emerson, model_peakman, model_igor, model_b_i, model_e_i, model_p_i)

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
      shape = `Sequence type`,
      linetype = `Sequence type`
    )
  ) +
  geom_point(
    size = 1.6,
    alpha = 0.2,
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
  xlim(0, 1) +
  ylim(0, 1) +
  labs(
    y = plot_y,
    x = plot_x
  ) +
  scale_color_brewer(
    palette = 'Dark2'
  ) +
  facet_wrap(
    vars(combination)
  )

jpeg(output_filename, width = 3000, height = 4000, res = 300)
eval_compare
dev.off()
