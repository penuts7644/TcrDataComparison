# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
SUBSET_IDS <- c('all', '100000', '50000', '10000', '5000', '1000', '500')
TYPES <- c('productive', 'unproductive', 'all')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(type, subset) {
  cc <- c(NA, rep("NULL", 6), rep(NA, 2))
  if (type == 'productive') {
    model_0 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_0/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_1/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_2 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_2/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_3 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_3/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_4 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_4/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_combined <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_combined/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_default <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_igor/pgen_estimate_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  } else if (type == 'unproductive') {
    model_0 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_0/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_1/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_2 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_2/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_3 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_3/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_4 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_4/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_combined <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_combined/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_default <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_igor/pgen_estimate_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  } else if (type == 'all') {
    model_0 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_0/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_1/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_2 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_2/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_3 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_3/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_4 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_4/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_combined <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_combined/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_default <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_igor/pgen_estimate_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  }
  models <- na.omit(do.call('cbind', list(
    model_0, model_1[2], model_2[2], model_3[2],
    model_4[2], model_combined[2], model_default[2]
  )))
  rm(cc, model_0, model_1, model_2, model_3, model_4, model_combined, model_default)
  names(models) <- c('type', '0', '1', '2', '3', '4', 'combined', 'default')
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

for (subset in SUBSET_IDS) {
  for (type in TYPES) {

    # ----------
    # EXTRACTING THE DATA
    # ----------
    models <- process_model(type, subset)
    plot_title <- 'Normalized AA and NT sequence pgen compared between models'
    output_filename <- paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluated_seqs_', type, '_rplot_', subset, '.png', sep = '')

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
