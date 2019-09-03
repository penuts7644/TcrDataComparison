# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
MODELS <- c('subject_0', 'subject_1', 'subject_2', 'subject_3', 'subject_4')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(model) {
  cc <- c(NA, rep("NULL", 6), rep(NA, 2))
  model_all <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/all/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_100000 <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/100000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_50000 <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/50000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_10000 <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/10000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_5000 <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/5000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_1000 <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/1000/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  model_500 <- melt(data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', model, '/500/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  models <- na.omit(do.call('cbind', list(
    model_all, model_100000[2], model_50000[2], model_10000[2], model_5000[2], model_1000[2], model_500[2]
  )))
  rm(cc, model_all, model_100000, model_50000, model_10000, model_5000, model_1000, model_500)
  names(models) <- c('type', 'all', '100000', '50000', '10000', '5000', '1000', '500')
  models$type <- as.character(models$type)
  models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
  models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
  models$type <- as.factor(models$type)
  models[models$type == 'NT', -1] <- apply(models[models$type == 'NT', -1], 2, normalize)
  models[models$type == 'AA', -1] <- apply(models[models$type == 'AA', -1], 2, normalize)
  names(model)[1] <- "Sequence type"
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

for (model in MODELS) {
  
  # ----------
  # EXTRACTING THE DATA
  # ----------
  models <- process_model(model)
  plot_y <- 'Pgen score (normalized)'
  plot_x <- 'Pgen score (normalized)'
  output_filename <- paste('~/Downloads/claim_5/subsample_evaluation_plot_', model, '.png', sep = '')
  
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
