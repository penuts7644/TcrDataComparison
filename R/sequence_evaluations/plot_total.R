# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)
SUBSET_IDS <- c('all', '100000', '50000', '10000', '5000', '1000', '500')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(x, y, combination, name) {
  newX <- melt(x, variable.name = 'type', value.name = 'X')
  newY <- melt(y, variable.name = 'type', value.name = 'Y')
  model <- na.omit(cbind(newX, newY['Y']))
  model$type <- as.character(model$type)
  model$type[model$type == 'nt_pgen_estimate'] <- 'NT'
  model$type[model$type == 'aa_pgen_estimate'] <- 'AA'
  model$type <- as.factor(model$type)
  model[model$type == 'NT', -1] <- apply(model[model$type == 'NT', -1], 2, normalize)
  model[model$type == 'AA', -1] <- apply(model[model$type == 'AA', -1], 2, normalize)
  model$corr.NT[model$type == 'NT'] <- c(cor(model[model$type == 'NT', 2], model[model$type == 'NT', 3], method = 'spearman'), rep(NA, nrow(model[model$type == 'NT', ]) - 1))
  model$corr.AA[model$type == 'AA'] <- c(cor(model[model$type == 'AA', 2], model[model$type == 'AA', 3], method = 'spearman'), rep(NA, nrow(model[model$type == 'AA', ]) - 1))
  model$combination <- combination
  model$name <- name
  model <- model[sample(nrow(model)), ]
  return (model)
}

for (subset in SUBSET_IDS) {
  
  # ----------
  # PLOT VARIABLES
  # ----------
  plot_title <- 'Normalized AA and NT sequence pgen compared for each model version'
  output_filename <- paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluated_seqs_total_rplot_', subset, '.png', sep = '')
  
  # ----------
  # MODEL DATA
  # ----------
  cc <- c(NA, rep("NULL", 6), rep(NA, 2))
  
  model_0_productive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_0/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_0_unproductive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_0/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_0_all <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_0/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_0_p_a <- process_model(model_0_productive, model_0_all, 'productive - all', '0')
  model_0_p_u <- process_model(model_0_productive, model_0_unproductive, 'productive - unproductive', '0')
  model_0_a_u <- process_model(model_0_all, model_0_unproductive, 'all - unproductive', '0')
  model_0 <- as.data.frame(do.call("rbind", list(model_0_p_a, model_0_p_u, model_0_a_u)))
  rm(model_0_productive, model_0_unproductive, model_0_all, model_0_p_a, model_0_p_u, model_0_a_u)
  
  model_1_productive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_1/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_1_unproductive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_1/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_1_all <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_1/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_1_p_a <- process_model(model_1_productive, model_1_all, 'productive - all', '1')
  model_1_p_u <- process_model(model_1_productive, model_1_unproductive, 'productive - unproductive', '1')
  model_1_a_u <- process_model(model_1_all, model_1_unproductive, 'all - unproductive', '1')
  model_1 <- as.data.frame(do.call("rbind", list(model_1_p_a, model_1_p_u, model_1_a_u)))
  rm(model_1_productive, model_1_unproductive, model_1_all, model_1_p_a, model_1_p_u, model_1_a_u)
  
  model_2_productive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_2/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_2_unproductive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_2/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_2_all <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_2/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_2_p_a <- process_model(model_2_productive, model_2_all, 'productive - all', '2')
  model_2_p_u <- process_model(model_2_productive, model_2_unproductive, 'productive - unproductive', '2')
  model_2_a_u <- process_model(model_2_all, model_2_unproductive, 'all - unproductive', '2')
  model_2 <- as.data.frame(do.call("rbind", list(model_2_p_a, model_2_p_u, model_2_a_u)))
  rm(model_2_productive, model_2_unproductive, model_2_all, model_2_p_a, model_2_p_u, model_2_a_u)
  
  model_3_productive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_3/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_3_unproductive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_3/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_3_all <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_3/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_3_p_a <- process_model(model_3_productive, model_3_all, 'productive - all', '3')
  model_3_p_u <- process_model(model_3_productive, model_3_unproductive, 'productive - unproductive', '3')
  model_3_a_u <- process_model(model_3_all, model_3_unproductive, 'all - unproductive', '3')
  model_3 <- as.data.frame(do.call("rbind", list(model_3_p_a, model_3_p_u, model_3_a_u)))
  rm(model_3_productive, model_3_unproductive, model_3_all, model_3_p_a, model_3_p_u, model_3_a_u)
  
  model_4_productive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_4/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_4_unproductive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_4/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_4_all <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_4/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_4_p_a <- process_model(model_4_productive, model_4_all, 'productive - all', '4')
  model_4_p_u <- process_model(model_4_productive, model_4_unproductive, 'productive - unproductive', '4')
  model_4_a_u <- process_model(model_4_all, model_4_unproductive, 'all - unproductive', '4')
  model_4 <- as.data.frame(do.call("rbind", list(model_4_p_a, model_4_p_u, model_4_a_u)))
  rm(model_4_productive, model_4_unproductive, model_4_all, model_4_p_a, model_4_p_u, model_4_a_u)
  
  model_combined_productive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_combined/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_combined_unproductive <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_combined/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_combined_all <- data.frame(read.table(paste('~/Downloads/process_5_files/', subset, '/evaluated_CDR3/evaluate_combined/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
  model_combined_p_a <- process_model(model_combined_productive, model_combined_all, 'productive - all', 'combined')
  model_combined_p_u <- process_model(model_combined_productive, model_combined_unproductive, 'productive - unproductive', 'combined')
  model_combined_a_u <- process_model(model_combined_all, model_combined_unproductive, 'all - unproductive', 'combined')
  model_combined <- as.data.frame(do.call("rbind", list(model_combined_p_a, model_combined_p_u, model_combined_a_u)))
  rm(model_combined_productive, model_combined_unproductive, model_combined_all, model_combined_p_a, model_combined_p_u, model_combined_a_u)
  
  models <- do.call('rbind', list(model_0, model_1, model_2, model_3, model_4, model_combined))
  rm(cc, model_0, model_1, model_2, model_3, model_4, model_combined)
  
  # ----------
  # MAKING THE PLOTS
  # ----------
  eval_compare <-
    ggplot(
      data = models,
      aes(
        x = X,
        y = Y,
        color = type
      )
    ) +
    geom_point(
      size = 2,
      alpha = 0.4,
    ) +
    geom_smooth(
      size = 1,
      method = lm,
      se = FALSE,
      fullrange = TRUE,
      show.legend = FALSE
    ) +
    geom_label(
      aes(
        x = 0,
        y = 1,
        label = round(corr.NT, digits = 4),
        col = 'NT'
      ),
      hjust = 0,
      vjust = 1,
      size = 6,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_label(
      aes(
        x = 0,
        y = 0.8,
        label = round(corr.AA, digits = 4),
        col = 'AA'
      ),
      hjust = 0,
      vjust = 1,
      size = 6,
      fontface = "bold",
      show.legend = FALSE
    ) +
    theme_bw(
      plot.title = element_text(
        size = 20
      ),
      legend.text = element_text(
        size = 15
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
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
      legend.title = element_blank(),
      legend.direction = 'horizontal'
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = plot_title
    ) +
    scale_color_manual(
      values = c('#1f78b4', '#d95f02')
    ) +
    facet_grid(
      rows = vars(name),
      cols = vars(combination)
    )
  
  jpeg(output_filename, width = 3000, height = 4000, res = 300)
  eval_compare
  dev.off()
}
