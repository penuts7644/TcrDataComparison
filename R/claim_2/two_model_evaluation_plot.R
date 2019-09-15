# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function() {
  cc1 <- c(NA, "NULL", rep(NA, 7))
  cc2 <- c(NA, rep("NULL", 2), rep(NA, 4))
  model_2_pgen <- data.frame(read.table('~/Downloads/claim_2/evaluations/emerson/subject_2/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
  model_2_cdr3 <- data.frame(read.table('~/Downloads/claim_2/models/emerson/subject_2/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
  model_2_cdr3$origin.y <- 'subject 3'
  model_2_pgen <- merge(model_2_pgen, model_2_cdr3, by = c('rearrangement', 'amino_acid', 'v_gene_choice', 'j_gene_choice'), all.x=TRUE)
  model_2_pgen <- model_2_pgen[!duplicated(model_2_pgen$row_id), ]
  model_2_pgen[is.na(model_2_pgen$origin.y), 'origin.y'] <- 'other'
  model_3_pgen <- data.frame(read.table('~/Downloads/claim_2/evaluations/emerson/subject_3/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
  model_3_cdr3 <- data.frame(read.table('~/Downloads/claim_2/models/emerson/subject_3/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
  model_3_cdr3$origin.x <- 'control 1'
  model_3_pgen <- merge(model_3_pgen, model_3_cdr3, by = c('rearrangement', 'amino_acid', 'v_gene_choice', 'j_gene_choice'), all.x = TRUE)
  model_3_pgen <- model_3_pgen[!duplicated(model_3_pgen$row_id), ]
  model_3_pgen[is.na(model_3_pgen$origin.x), 'origin.x'] <- 'other'
  model_2 <- melt(model_2_pgen[order(model_2_pgen$row_id), ][, 6:8])
  model_3 <- melt(model_3_pgen[order(model_3_pgen$row_id), ][, 6:8])
  models <- na.omit(cbind(model_2[, c(1, 3)], model_3[, c(1, 3)]))
  models[models$origin.y == 'other' & models$origin.x == 'control 1', 'origin.y'] <- 'control 1'
  # models[models$origin.y == 'subject 3' & models$origin.x == 'control 1', 'origin.y'] <- 'both'
  models[models$origin.y == 'subject 3' & models$origin.x == 'control 1', 'origin.y'] <- 'other'
  models <- models[, c(1, 2, 4)]
  names(models) <- c('Sequence origin', 'Y', 'X')
  models <- models[models$`Sequence origin` != 'other', ]
  models <- models[sample(nrow(models)), ]
  return (models)
}

# ----------
# PLOT VARIABLES
# ----------
models <- process_model()
plot_y <- 'Pgen score (subject 3)'
plot_x <- 'Pgen score (control 1)'
output_filename <- '~/Downloads/claim_2/two_model_evaluation_plot.png'

# ----------
# MAKING THE PLOTS
# ----------
eval_compare <-
  ggplot(
    data = models,
    aes(
      x = X,
      y = Y,
      color = `Sequence origin`
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
    alpha = 0.8
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
    # legend.direction = 'vertical'
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
    values = c('#1b9e77', '#d95f02')
  )

jpeg(output_filename, width = 2000, height = 2000, res = 300)
eval_compare
dev.off()
