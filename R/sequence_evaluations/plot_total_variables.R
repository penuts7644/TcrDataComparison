# ----------
# PLOT VARIABLES
# ----------
plot_title <- 'Normalized AA and NT sequence pgen compared for each model version'
plot_caption <- 'Subject-specific and combined models compared against their respective trained model versions.\nEach model has been trained with either productive, unproductive or all (productive and\nunproductive) sequences. All models evaluated the same combined data file that is also used for\ntraining the combined model. Spearman correlation method is used for computing the correlations.'
output_filename <- '~/Downloads/evaluated_seqs_total_rplot.png'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

process_model <- function(x, y, combination, name) {
  model <- na.omit(cbind(x, y))
  colnames(model) <- c('NT.x', 'AA.x', 'NT.y', 'AA.y')
  model$NT.corr <- cor(model$NT.x, model$NT.y, method = 'spearman')
  model$AA.corr <- cor(model$AA.x, model$AA.y, method = 'spearman')
  model[, 1:4] <- apply(model[1:4], 2, normalize)
  model$combination <- combination
  model$name <- name
  return (model)
}

# ----------
# MODEL DATA
# ----------
cc <- c(NA, rep("NULL", 6), rep(NA, 2))  # For AA and NT comparison

model_0_productive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_0/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_0_unproductive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_0/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_0_all <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_0/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_0_p_a <- process_model(model_0_productive, model_0_all, 'productive - all', '0')
model_0_p_u <- process_model(model_0_productive, model_0_unproductive, 'productive - unproductive', '0')
model_0_a_u <- process_model(model_0_all, model_0_unproductive, 'all - unproductive', '0')
model_0 <- as.data.frame(do.call("rbind", list(model_0_p_a, model_0_p_u, model_0_a_u)))
rm(model_0_productive, model_0_unproductive, model_0_all, model_0_p_a, model_0_p_u, model_0_a_u)

model_1_productive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_1/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_1_unproductive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_1/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_1_all <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_1/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_1_p_a <- process_model(model_1_productive, model_1_all, 'productive - all', '1')
model_1_p_u <- process_model(model_1_productive, model_1_unproductive, 'productive - unproductive', '1')
model_1_a_u <- process_model(model_1_all, model_1_unproductive, 'all - unproductive', '1')
model_1 <- as.data.frame(do.call("rbind", list(model_1_p_a, model_1_p_u, model_1_a_u)))
rm(model_1_productive, model_1_unproductive, model_1_all, model_1_p_a, model_1_p_u, model_1_a_u)

model_2_productive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_2/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_2_unproductive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_2/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_2_all <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_2/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_2_p_a <- process_model(model_2_productive, model_2_all, 'productive - all', '2')
model_2_p_u <- process_model(model_2_productive, model_2_unproductive, 'productive - unproductive', '2')
model_2_a_u <- process_model(model_2_all, model_2_unproductive, 'all - unproductive', '2')
model_2 <- as.data.frame(do.call("rbind", list(model_2_p_a, model_2_p_u, model_2_a_u)))
rm(model_2_productive, model_2_unproductive, model_2_all, model_2_p_a, model_2_p_u, model_2_a_u)

model_3_productive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_3/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_3_unproductive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_3/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_3_all <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_3/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_3_p_a <- process_model(model_3_productive, model_3_all, 'productive - all', '3')
model_3_p_u <- process_model(model_3_productive, model_3_unproductive, 'productive - unproductive', '3')
model_3_a_u <- process_model(model_3_all, model_3_unproductive, 'all - unproductive', '3')
model_3 <- as.data.frame(do.call("rbind", list(model_3_p_a, model_3_p_u, model_3_a_u)))
rm(model_3_productive, model_3_unproductive, model_3_all, model_3_p_a, model_3_p_u, model_3_a_u)

model_4_productive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_4/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_4_unproductive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_4/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_4_all <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_4/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_4_p_a <- process_model(model_4_productive, model_4_all, 'productive - all', '4')
model_4_p_u <- process_model(model_4_productive, model_4_unproductive, 'productive - unproductive', '4')
model_4_a_u <- process_model(model_4_all, model_4_unproductive, 'all - unproductive', '4')
model_4 <- as.data.frame(do.call("rbind", list(model_4_p_a, model_4_p_u, model_4_a_u)))
rm(model_4_productive, model_4_unproductive, model_4_all, model_4_p_a, model_4_p_u, model_4_a_u)

model_combined_productive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_combined/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_combined_unproductive <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_combined/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_combined_all <- head(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_combined/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)), 100)
model_combined_p_a <- process_model(model_combined_productive, model_combined_all, 'productive - all', 'combined')
model_combined_p_u <- process_model(model_combined_productive, model_combined_unproductive, 'productive - unproductive', 'combined')
model_combined_a_u <- process_model(model_combined_all, model_combined_unproductive, 'all - unproductive', 'combined')
model_combined <- as.data.frame(do.call("rbind", list(model_combined_p_a, model_combined_p_u, model_combined_a_u)))
rm(model_combined_productive, model_combined_unproductive, model_combined_all, model_combined_p_a, model_combined_p_u, model_combined_a_u)

models <- do.call('rbind', list(model_0, model_1, model_2, model_3, model_4, model_combined))
rm(model_0, model_1, model_2, model_3, model_4, model_combined)

