# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
PROJECTS <- c('dejong', 'emerson', 'gomez')
CORRELATION_METHOD <- 'pearson' # or 'spearman'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function(project) {
  cc <- c(NA, "NULL", NA, rep("NULL", 4), rep(NA, 2))
  model_0 <- data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_0/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_1 <- data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_1/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_2 <- data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_2/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_3 <- data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_3/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_4 <- data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_4/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_0 <- melt(model_0, id.vars = 'row_id')
  model_1 <- melt(model_1, id.vars = 'row_id')
  model_2 <- melt(model_2, id.vars = 'row_id')
  model_3 <- melt(model_3, id.vars = 'row_id')
  model_4 <- melt(model_4, id.vars = 'row_id')
  names(model_0) <- c('row_id', 'type', 'subject 1')
  names(model_1) <- c('row_id', 'type', 'subject 2')
  names(model_2) <- c('row_id', 'type', 'subject 3')
  names(model_3) <- c('row_id', 'type', 'control 1')
  names(model_4) <- c('row_id', 'type', 'control 2')
  models <- Reduce(function(x, y) merge(x, y, by=c('row_id', 'type')), list(model_0, model_1, model_2, model_3, model_4))
  rm(cc, model_0, model_1, model_2, model_3, model_4)
  models$type <- as.character(models$type)
  models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
  models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
  models$type <- as.factor(models$type)
  models <- models[!(models$row_id %in% models[is.na(models[, 3:7]), 'row_id']), ][, 2:7]
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

for (project in PROJECTS) {

  # ----------
  # EXTRACTING THE DATA
  # ----------
  models <- process_model(project)
  plot_y <- 'Generation probability score'
  plot_x <- 'Generation probability score'
  output_filename <- paste('~/Downloads/claim_2/subject_evaluation_plot_', project, '_', CORRELATION_METHOD, '.png', sep = '')

  # ----------
  # MAKING THE PLOTS
  # ----------
  legend_plot <- ggplot(data = models) +
    geom_point(
      aes(
        x = `subject 1`,
        y = `subject 2`,
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
      x = plot_x,
      y = plot_y
    )

  jpeg(output_filename, width = 4000, height = 4000, res = 300)
  print(eval_compare)
  dev.off()
}
