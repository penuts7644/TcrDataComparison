# ----------
# GENERAL
# ----------
library(reshape2)
SUBSET_ID <- 'all'


# ----------
# PLOT VARIABLES
# ----------
plot_title <- 'Normalized subject-specific complete model entropies'
plot_caption <- 'Subject-specific models compared against other subject-specific and combined models that were trained using productive,\nunproductive or all sequences from each dataset.'
plot_y <- 'combined KL divergence\n(complete models)'
output_filename <- paste('~/Downloads/process_5_files/', SUBSET_ID, '/model_entropies/entropy_total_rplot_', SUBSET_ID, '.png', sep = '')

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# ----------
# MODEL DATA
# ----------
productive <- data.matrix(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/model_entropies/productive/calc_model_entropy_total.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
unproductive <- data.matrix(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/model_entropies/unproductive/calc_model_entropy_total.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))
all <- data.matrix(read.table(paste('~/Downloads/process_5_files/', SUBSET_ID, '/model_entropies/all/calc_model_entropy_total.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE))

models <- melt(list(
  'productive' = productive,
  'unproductive' = unproductive,
  'all' = all
))
models <- models[models$Var1 != models$Var2, ]
models <- models[models$Var1 != 'combined', ]
models$L2 <- 'subject'
models$L2[models$Var2 == 'combined'] <- 'combined'
models <- models[order(models$L2, decreasing = TRUE), ]
models[, 3] <- apply(models[3], 2, normalize)
rm(productive, unproductive, all)
