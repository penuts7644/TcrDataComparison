# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
MODELS <- c('subject_0', 'subject_1', 'subject_2', 'subject_3', 'subject_4')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function(model) {
  cc <- c(NA, "NULL", NA, rep("NULL", 4), rep(NA, 2))
  model_all <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/all/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_50000 <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/50000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_10000 <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/10000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_5000 <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/5000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_1000 <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/1000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_500 <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/500/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_100 <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/100/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_all <- melt(model_all[order(model_all$row_id), ][, 2:3])
  model_50000 <- melt(model_50000[order(model_50000$row_id), ][, 2:3])
  model_10000 <- melt(model_10000[order(model_10000$row_id), ][, 2:3])
  model_5000 <- melt(model_5000[order(model_5000$row_id), ][, 2:3])
  model_1000 <- melt(model_1000[order(model_1000$row_id), ][, 2:3])
  model_500 <- melt(model_500[order(model_500$row_id), ][, 2:3])
  model_100 <- melt(model_100[order(model_100$row_id), ][, 2:3])
  models <- na.omit(do.call('cbind', list(
    model_all, model_50000[2], model_10000[2], model_5000[2], model_1000[2], model_500[2], model_100[2]
  )))
  rm(cc, model_all, model_50000, model_10000, model_5000, model_1000, model_500, model_100)
  names(models) <- c('type', 'all', '50000', '10000', '5000', '1000', '500', '100')
  models$type <- as.character(models$type)
  models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
  models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
  models$type <- as.factor(models$type)
  names(models)[1] <- "Sequence type"
  models <- models[sample(nrow(models)), ]
  return (models)
}

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
  geom_smooth(
    size = 1.2,
    method = lm,
    se = FALSE,
    fullrange = TRUE
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
  plot_y <- 'Pgen score'
  plot_x <- 'Pgen score'
  output_filename <- paste('~/Downloads/claim_5/subsample_evaluation_plot_', model, '.png', sep = '')

  # ----------
  # MAKING THE PLOTS
  # ----------
  legend_plot <- ggplot(data = models) +
    geom_smooth(
      aes(
        x = `100`,
        y = `500`,
        color = `Sequence type`,
        linetype = `Sequence type`
      ),
      size = 1.2,
      method = lm,
      se = FALSE,
      fullrange = TRUE
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
        color = `Sequence type`,
        linetype = `Sequence type`
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
