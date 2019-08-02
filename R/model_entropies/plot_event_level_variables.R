# ----------
# GENERAL
# ----------
library(reshape2)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(type) {
  if (type == 'productive') {
    V <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_V.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    V_trim_3 <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_V_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_D.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_3 <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_D_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_5 <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_D_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_J.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J_trim_5 <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_J_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_VD <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_insert_length_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_DJ <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_insert_length_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_VD <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_dinuc_markov_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_DJ <- data.matrix(read.table('~/Downloads/model_entropies/productive/calc_model_entropy_dinuc_markov_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  } else if (type == 'unproductive') {
    V <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_V.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    V_trim_3 <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_V_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_D.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_3 <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_D_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_5 <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_D_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_J.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J_trim_5 <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_J_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_VD <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_insert_length_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_DJ <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_insert_length_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_VD <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_dinuc_markov_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_DJ <- data.matrix(read.table('~/Downloads/model_entropies/unproductive/calc_model_entropy_dinuc_markov_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
  } else if (type == 'all') {
    V <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_V.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    V_trim_3 <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_V_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_D.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_3 <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_D_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    D_trim_5 <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_D_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_J.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    J_trim_5 <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_J_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_VD <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_insert_length_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    insert_length_DJ <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_insert_length_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_VD <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_dinuc_markov_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
    dinuc_markov_DJ <- data.matrix(read.table('~/Downloads/model_entropies/all/calc_model_entropy_dinuc_markov_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
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
  models$L2 <- 'subject'
  models$L2[models$Var2 == 'combined'] <- 'combined'
  models <- models[order(models$L2, decreasing = TRUE), ]
  models[, 3] <- apply(models[3], 2, normalize)
  return (models)
}

# ----------
# PRODUCTIVE
# ----------
models <- process_model('productive')
plot_title <- 'Normalized subject-specific event-level entropies'
plot_caption <- 'Subject-specific models compared against other subject-specific and combined models that were trained using productive sequences\nfrom each dataset.'
plot_y <- 'combined KL divergence'
output_filename <- '~/Downloads/entropy_event_productive_rplot.png'

# ----------
# UNPRODUCTIVE
# ----------
models <- process_model('unproductive')
plot_title <- 'Normalized subject-specific event-level entropies'
plot_caption <- 'Subject-specific models compared against other subject-specific and combined models that were trained using unproductive sequences\nfrom each dataset.'
plot_y <- 'combined KL divergence'
output_filename <- '~/Downloads/entropy_event_unproductive_rplot.png'

# ----------
# ALL
# ----------
models <- process_model('all')
plot_title <- 'Normalized subject-specific event-level entropies'
plot_caption <- 'Subject-specific models compared against other subject-specific and combined models that were trained using all sequences (productive\nand unproductive) from each dataset.'
plot_y <- 'combined KL divergence'
output_filename <- '~/Downloads/entropy_event_all_rplot.png'
