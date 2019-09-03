# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
PROJECTS <- c('brusko', 'emerson', 'peakman')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(project) {
  cc <- c(NA, rep("NULL", 6), rep(NA, 2))
  model_0 <- melt(data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_0/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_1 <- melt(data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_1/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_2 <- melt(data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_2/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_3 <- melt(data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_3/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_4 <- melt(data.frame(read.table(paste('~/Downloads/claim_2/evaluations/', project, '/subject_4/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  models <- na.omit(do.call('cbind', list(
    model_0, model_1[2], model_2[2], model_3[2], model_4[2]
  )))
  rm(cc, model_0, model_1, model_2, model_3, model_4)
  names(models) <- c('type', '0', '1', '2', '3', '4')
  models$type <- as.character(models$type)
  models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
  models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
  models$type <- as.factor(models$type)
  models[models$type == 'NT', -1] <- apply(models[models$type == 'NT', -1], 2, normalize)
  models[models$type == 'AA', -1] <- apply(models[models$type == 'AA', -1], 2, normalize)
  names(models)[1] <- "Sequence type"
  names(models)[2] <- "disease 1"
  names(models)[3] <- "disease 2"
  names(models)[4] <- "disease 3"
  names(models)[5] <- "control 1"
  names(models)[6] <- "control 2"
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
  scale_color_brewer(
    palette = 'Dark2'
  )
}

lower_plot_fn <- function(data, mapping, ...){
  ggplot(
    data = data,
    mapping = mapping
  ) +
    geom_point(
      size = 1.6,
      alpha = 0.2
    ) +
    geom_smooth(
      size = 1.2,
      method = lm,
      se = FALSE,
      fullrange = TRUE
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    scale_color_brewer(
      palette = 'Dark2'
    )
}

for (project in PROJECTS) {

  # ----------
  # EXTRACTING THE DATA
  # ----------
  models <- process_model(project)
  plot_y <- 'Pgen score (normalized)'
  plot_x <- 'Pgen score (normalized)'
  output_filename <- paste('~/Downloads/claim_2/subject_evaluation_plot_', project, '.png', sep = '')

  # ----------
  # MAKING THE PLOTS
  # ----------
  legend_plot <- ggplot(data = models) +
    geom_smooth(
      aes(
        x = `disease 1`,
        y = `disease 2`,
        color = `Sequence type`,
        linetype = `Sequence type`
      ),
      size = 1.2,
      method = lm,
      se = FALSE,
      fullrange = TRUE
    ) +
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
    scale_color_brewer(
      palette = 'Dark2'
    )

  eval_compare <-
    ggpairs(
      models,
      columns = 2:ncol(models),
      mapping = aes(
        color = `Sequence type`,
        shape = `Sequence type`,
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
