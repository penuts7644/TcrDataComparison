# ----------
# GENERAL
# ----------
library(ggplot2)
library(reshape2)
SUBSET_IDS <- c('all', '100000', '50000', '10000', '5000', '1000', '500')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

for (subset in SUBSET_IDS) {
  
  # ----------
  # PLOT VARIABLES
  # ----------
  plot_title <- 'Normalized subject-specific complete model entropies'
  plot_y <- 'combined KL divergence (complete models) (bits)'
  output_filename <- paste('~/Downloads/process_5_files/', subset, '/model_entropies/entropy_total_rplot_', subset, '.png', sep = '')
  
  # ----------
  # MODEL DATA
  # ----------
  productive <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_total.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  unproductive <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_total.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  all <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_total.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  
  models <- melt(list(
    'productive' = productive,
    'unproductive' = unproductive,
    'all' = all
  ))
  models <- models[models$Var1 != models$Var2, ]
  models <- models[models$Var1 != 'combined', ]
  models$L2 <- 'subject - subject'
  models$L2[models$Var2 == 'combined'] <- 'subject - combined'
  models <- models[order(models$L2, decreasing = TRUE), ]
  models[, 3] <- apply(models[3], 2, normalize)
  rm(productive, unproductive, all)
  
  # ----------
  # MAKING THE PLOTS
  # ----------
  entr_compare <-
    ggplot(
      models,
      aes(
        x = L1,
        y = value
      )
    ) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_jitter(
      height = 0,
      width = .2,
      size = 4,
      aes(color = L2),
      alpha = 0.8
    ) +
    theme_bw(
      plot.title = element_text(
        size = 20
      ),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        size = 15,
      ),
      axis.title.y = element_text(
        size = 15
      ),
      axis.text.y = element_text(
        size = 11
      ),
      legend.text = element_text(
        size = 15
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.direction = "horizontal"
    ) +
    ylim(0, 1) +
    labs(
      title = plot_title,
      y = plot_y
    ) +
    scale_color_manual(
      values = c('#1f78b4', '#d95f02')
    )
  
  jpeg(output_filename, width = 4000, height = 4000, res = 300)
  entr_compare
  dev.off()
}
