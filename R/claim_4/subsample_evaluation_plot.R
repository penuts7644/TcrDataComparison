# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
MODELS <- c('subject_0', 'subject_1', 'subject_2', 'subject_3', 'subject_4')
CORRELATION_METHOD <- 'pearson' # or 'spearman'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function(model) {
  cc <- c(NA, "NULL", NA, rep("NULL", 4), rep(NA, 2))
  model_all <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/all/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_50000 <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/50000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_10000 <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/10000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_5000 <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/5000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_1000 <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/1000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_500 <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/500/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_100 <- data.frame(read.table(paste('~/Downloads/claim_4/evaluations/', model, '/100/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_all <- melt(model_all, id.vars = 'row_id')
  model_50000 <- melt(model_50000, id.vars = 'row_id')
  model_10000 <- melt(model_10000, id.vars = 'row_id')
  model_5000 <- melt(model_5000, id.vars = 'row_id')
  model_1000 <- melt(model_1000, id.vars = 'row_id')
  model_500 <- melt(model_500, id.vars = 'row_id')
  model_100 <- melt(model_100, id.vars = 'row_id')
  names(model_all) <- c('row_id', 'type', 'all')
  names(model_50000) <- c('row_id', 'type', '50000')
  names(model_10000) <- c('row_id', 'type', '10000')
  names(model_5000) <- c('row_id', 'type', '5000')
  names(model_1000) <- c('row_id', 'type', '1000')
  names(model_500) <- c('row_id', 'type', '500')
  names(model_100) <- c('row_id', 'type', '100')
  models <- Reduce(function(x, y) merge(x, y, by=c('row_id', 'type')), list(model_all, model_50000, model_10000, model_5000, model_1000, model_500, model_100))
  rm(cc, model_all, model_50000, model_10000, model_5000, model_1000, model_500, model_100)
  models$type <- as.character(models$type)
  models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
  models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
  models$type <- as.factor(models$type)
  models <- models[!(models$row_id %in% models[is.na(models[, 3:9]), 'row_id']), ][, 2:9]
  models <- models[sample(nrow(models)), ]
  names(models)[1] <- "Sequence type"
  return (models)
}

upper_plot_fn <- function(data, mapping, ...){
  ggally_cor(
    data = data,
    mapping = mapping,
    method = CORRELATION_METHOD,
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
    values = c('#ca0020', '#0571b0')
  )
}

lower_plot_fn <- function(data, mapping, ...){
  ggplot(
    data = data,
    mapping = mapping
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
  scale_x_sqrt(
     limits = c(0, max(pmax(data[, 2], data[, 3], data[, 4], data[, 5], data[, 6])))
   ) +
   scale_y_sqrt(
     limits = c(0, max(pmax(data[, 2], data[, 3], data[, 4], data[, 5], data[, 6])))
   ) +
  scale_color_manual(
    values = c('#ca0020', '#0571b0')
  )
}

for (model in MODELS) {

  # ----------
  # EXTRACTING THE DATA
  # ----------
  models <- process_model(model)
  plot_y <- 'Generation probability score'
  plot_x <- 'Generation probability score'
  output_filename <- paste('~/Downloads/claim_4/subsample_evaluation_plot_', model, '_', CORRELATION_METHOD, '.png', sep = '')

  # ----------
  # MAKING THE PLOTS
  # ----------
  legend_plot <- ggplot(data = models) +
    geom_point(
      aes(
        x = `100`,
        y = `500`,
        color = `Sequence type`
      ),
      size = 3,
      alpha = 1
    ) +
    theme_bw() +
    theme(
      legend.title = element_text(
        size = 18
      ),
      legend.text = element_text(
        size = 15
      ),
      legend.position = 'top',
      legend.direction = 'horizontal'
    ) +
    scale_color_manual(
      values = c('#ca0020', '#0571b0')
    )

  eval_compare <-
    ggpairs(
      models,
      columns = 2:ncol(models),
      mapping = aes(
        color = `Sequence type`
      ),
      axisLabels = 'show',
      upper = list(continuous = upper_plot_fn),
      diag = list(continuous = 'blankDiag'),
      lower = list(continuous = lower_plot_fn),
      progress = TRUE,
      legend = grab_legend(legend_plot)
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
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
      legend.position = 'top'
    ) +
    labs(
      y = plot_y,
      x = plot_x
    )

  jpeg(output_filename, width = 4000, height = 4000, res = 300)
  print(eval_compare)
  dev.off()
}
