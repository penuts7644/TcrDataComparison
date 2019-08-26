# ----------
# GENERAL
# ----------
library(ggplot2)
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
  if (type == 'productive') {
    V <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_V.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    V_trim_3 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_V_trim_3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_D.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_3 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_D_trim_3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_5 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_D_trim_5.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_J.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J_trim_5 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_J_trim_5.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_VD <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_insert_length_VD.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_DJ <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_insert_length_DJ.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_VD <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_dinuc_markov_VD.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_DJ <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/productive/calc_model_entropy_dinuc_markov_DJ.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  } else if (type == 'unproductive') {
    V <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_V.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    V_trim_3 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_V_trim_3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_D.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_3 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_D_trim_3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_5 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_D_trim_5.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_J.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J_trim_5 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_J_trim_5.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_VD <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_insert_length_VD.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_DJ <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_insert_length_DJ.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_VD <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_dinuc_markov_VD.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_DJ <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/unproductive/calc_model_entropy_dinuc_markov_DJ.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  } else if (type == 'all') {
    V <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_V.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    V_trim_3 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_V_trim_3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_D.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_3 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_D_trim_3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_5 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_D_trim_5.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_J.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J_trim_5 <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_J_trim_5.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_VD <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_insert_length_VD.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_DJ <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_insert_length_DJ.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_VD <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_dinuc_markov_VD.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_DJ <- data.matrix(read.table(paste('~/Downloads/process_5_files/', subset, '/model_entropies/all/calc_model_entropy_dinuc_markov_DJ.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  }
  models <- melt(list(
    'V-gene' = V,
    'D-gene' = D,
    'J-gene' = J,
    'V-gene trim (3)' = V_trim_3,
    'D-gene trim (3)' = D_trim_3,
    'D-gene trim (5)' = D_trim_5,
    'J-gene trim (5)' = J_trim_5,
    'VD insert length' = insert_length_VD,
    'DJ insert length' = insert_length_DJ,
    'DinucMarkov VD' = dinuc_markov_VD,
    'DinucMarkov DJ' = dinuc_markov_DJ
  ))
  rm(V, D, J, V_trim_3, D_trim_3, D_trim_5, J_trim_5, insert_length_VD, insert_length_DJ, dinuc_markov_VD, dinuc_markov_DJ)
  models <- models[models$Var1 != models$Var2, ]
  models <- models[models$Var1 != 'combined', ]
  models$L2 <- 'subject - subject'
  models$L2[models$Var2 == 'combined'] <- 'subject - combined'
  models <- models[order(models$L2, decreasing = TRUE), ]
  models[, 3] <- apply(models[3], 2, normalize)
  return (models)
}

for (subset in SUBSET_IDS) {
  for (type in TYPES) {

    # ----------
    # EXTRACTING THE DATA
    # ----------
    models <- process_model(type, subset)
    plot_title <- 'Normalized subject-specific event-level entropies'
    plot_y <- 'combined KL divergence (bits)'
    output_filename <- paste('~/Downloads/process_5_files/', subset, '/model_entropies/entropy_event_', type, '_rplot_', subset, '.png', sep = '')

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
          angle = 45,
          hjust = 1
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
}
