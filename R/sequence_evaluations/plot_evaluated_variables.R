# ----------
# GENERAL
# ----------
library(reshape2)
SUBSET_ID <- 'all' # Change this value to '100000', '50000', '10000', '5000', '1000' or '500'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(type) {
  cc <- c(NA, rep("NULL", 6), rep(NA, 2))
  if (type == 'productive') {
    model_0 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_0/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_1/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_2 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_2/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_3 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_3/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_4 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_4/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_combined <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_combined/pgen_estimate_productive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_default <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_igor/pgen_estimate_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  } else if (type == 'unproductive') {
    model_0 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_0/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_1/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_2 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_2/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_3 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_3/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_4 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_4/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_combined <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_combined/pgen_estimate_unproductive_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_default <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_igor/pgen_estimate_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
  } else if (type == 'all') {
    model_0 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_0/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_1 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_1/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_2 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_2/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_3 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_3/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_4 <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_4/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_combined <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_combined/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
    model_default <- melt(data.frame(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluate_igor/pgen_estimate_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
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

# ----------
# PRODUCTIVE
# ----------
models <- process_model('productive')
plot_title <- 'Normalized AA and NT sequence pgen compared between models'
plot_caption <- 'Subject-specific, combined and default models compared against each other. The models (except default) were trained using the\nproductive sequences. All models evaluated the same combined data file that is also used for training the combined model.\nSpearman correlation method is used for the upper-right plots.'
output_filename <- paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluated_seqs_productive_rplot_', SUBSET_ID, '.png', sep = '')

# ----------
# UNPRODUCTIVE
# ----------
models <- process_model('unproductive')
plot_title <- 'Normalized AA and NT sequence pgen compared between models'
plot_caption <- 'Subject-specific, combined and default models compared against each other. The models (except default) were trained using the\nunproductive sequences. All models evaluated the same combined data file that is also used for training the combined model.\nSpearman correlation method is used for the upper-right plots.'
output_filename <- paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluated_seqs_unproductive_rplot_', SUBSET_ID, '.png', sep = '')

# ----------
# ALL
# ----------
models <- process_model('all')
plot_title <- 'Normalized AA and NT sequence pgen compared between models'
plot_caption <- 'Subject-specific, combined and default models compared against each other. The models (except default) were trained using all\nsequences (productive and unproductive). All models evaluated the same combined data file that is also used for training the\ncombined model. Spearman correlation method is used for the upper-right plots.'
output_filename <- paste('~/Downloads/process_5_files/', SUBSET_ID, '/evaluated_CDR3/evaluated_seqs_all_rplot_', SUBSET_ID, '.png', sep = '')
