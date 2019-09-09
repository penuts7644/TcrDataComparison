# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
PROJECTS <- c('dejong', 'emerson', 'peakman')

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
  model_0 <- melt(model_0[order(model_0$row_id), ][, 2:3])
  model_1 <- melt(model_1[order(model_1$row_id), ][, 2:3])
  model_2 <- melt(model_2[order(model_2$row_id), ][, 2:3])
  model_3 <- melt(model_3[order(model_3$row_id), ][, 2:3])
  model_4 <- melt(model_4[order(model_4$row_id), ][, 2:3])
  models <- na.omit(do.call('cbind', list(
    model_0, model_1[2], model_2[2], model_3[2], model_4[2]
  )))
  rm(cc, model_0, model_1, model_2, model_3, model_4)
  names(models) <- c('type', 'subject 1', 'subject 2', 'subject 3', 'control 1', 'control 2')
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
    values = c('#ca0020', '#000000')
  )
}

lower_plot_fn <- function(data, mapping, ...){
  ggplot(
    data = data,
    mapping = mapping
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
    scale_color_manual(
      values = c('#ca0020', '#000000')
    )
}

for (project in PROJECTS) {

  # ----------
  # EXTRACTING THE DATA
  # ----------
  models <- process_model(project)
  plot_y <- 'Pgen score'
  plot_x <- 'Pgen score'
  output_filename <- paste('~/Downloads/claim_2/subject_evaluation_plot_', project, '.png', sep = '')

  # ----------
  # MAKING THE PLOTS
  # ----------
  legend_plot <- ggplot(data = models) +
    geom_smooth(
      aes(
        x = `subject 1`,
        y = `subject 2`,
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
      values = c('#ca0020', '#000000')
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
