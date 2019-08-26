# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
MODELS <- c('evaluate_0', 'evaluate_1', 'evaluate_2', 'evaluate_3', 'evaluate_4', 'evaluate_combined')
TYPES <- c('productive', 'unproductive', 'all')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(type, model) {
  cc <- c(NA, rep("NULL", 6), rep(NA, 2))
  if (type == 'productive') {
    model_all <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/all/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_100000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/100000/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_50000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/50000/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_10000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/10000/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_5000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/5000/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/1000/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_500 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/500/evaluated_CDR3/', model, '/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  } else if (type == 'unproductive') {
    model_all <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/all/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_100000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/100000/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_50000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/50000/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_10000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/10000/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_5000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/5000/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/1000/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_500 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/500/evaluated_CDR3/', model, '/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  } else if (type == 'all') {
    model_all <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/all/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_100000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/100000/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_50000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/50000/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_10000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/10000/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_5000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/5000/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1000 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/1000/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_500 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/500/evaluated_CDR3/', model, '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  }
  models <- na.omit(do.call('cbind', list(
    model_all, model_100000[2], model_50000[2], model_10000[2],
    model_5000[2], model_1000[2], model_500[2]
  )))
  rm(cc, model_all, model_100000, model_50000, model_10000, model_5000, model_1000, model_500)
  names(models) <- c('type', 'all', '100000', '50000', '10000', '5000', '1000', '500')
  models$type <- as.character(models$type)
  models$type[models$type == 'nt_pgen_estimate'] <- 'NT'
  models$type[models$type == 'aa_pgen_estimate'] <- 'AA'
  models$type <- as.factor(models$type)
  models[models$type == 'NT', -1] <- apply(models[models$type == 'NT', -1], 2, normalize)
  models[models$type == 'AA', -1] <- apply(models[models$type == 'AA', -1], 2, normalize)
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
      values = c('#1f78b4', '#d95f02')
    )
}

lower_plot_fn <- function(data, mapping, ...){
  ggplot(
    data = data,
    mapping = mapping
  ) +
    geom_point(
      size = 2,
      alpha = 0.4
    ) +
    geom_smooth(
      size = 1,
      method = lm,
      se = FALSE,
      fullrange = TRUE
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    scale_color_manual(
      values = c('#1f78b4', '#d95f02')
    )
}

for (model in MODELS) {
  for (type in TYPES) {
  
    # ----------
    # EXTRACTING THE DATA
    # ----------
    models <- process_model(type, model)
    plot_title <- 'Normalized AA and NT sequence pgen compared between single model from various subsamples'
    output_filename <- paste('~/Downloads/process_5_files/evaluated_seqs_', type, '_rplot_', model, '.png', sep = '')
    
    # ----------
    # MAKING THE PLOTS
    # ----------
    legend_plot <- ggplot(data = models) +
      geom_point(
        aes(x = `0`, y = `1`, color = type)
      ) +
      theme(
        legend.text = element_text(
          size = 15
        ),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.direction = 'horizontal'
      ) +
      scale_color_manual(
        values = c('#1f78b4', '#d95f02')
      )
    
    eval_compare <-
      ggpairs(
        models,
        columns = 2:ncol(models),
        mapping = aes(color = type),
        axisLabels = 'show',
        upper = list(continuous = upper_plot_fn),
        diag = list(continuous = 'blankDiag'),
        lower = list(continuous = lower_plot_fn),
        progress = TRUE,
        legend = grab_legend(legend_plot)
      ) +
      theme_bw(
        plot.title = element_text(
          size = 20
        ),
        axis.text.x = element_text(
          size = 11,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(
          size = 11
        ),
        strip.text.x = element_text(
          size = 15
        ),
        strip.text.y = element_text(
          size = 15
        ),
        legend.position = 'bottom',
        strip.placement = 'outside'
      ) +
      labs(
        title = plot_title
      )
    
    jpeg(output_filename, width = 4000, height = 4000, res = 300)
    eval_compare
    dev.off()
  }
}
