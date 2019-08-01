# ----------
# GENERAL
# ----------
library(reshape2)
cc <- c(NA, rep("NULL", 6), rep(NA, 2))  # For AA and NT comparison

# ----------
# UNPRODUCTIVE
# ----------
model_0 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_0/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_1 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_1/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_2 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_2/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_3 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_3/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_4 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_4/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_combined <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_combined/pgen_estimate_unproductive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_default <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_igor/pgen_estimate_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
plot_title <- 'Normalized AA and NT sequence pgen compared between models'
plot_caption <- 'Subject-specific, combined and default models compared against each other. The models (except default) were trained using the\nunproductive sequences. All models evaluated the same combined data file that is also used for training the combined model.\nSpearman correlation method is used for the upper-right plots.'
output_filename <- '~/Downloads/evaluated_seqs_unproductive_rplot.png'

# ----------
# PRODUCTIVE
# ----------
model_0 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_0/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_1 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_1/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_2 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_2/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_3 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_3/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_4 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_4/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_combined <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_combined/pgen_estimate_productive_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_default <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_igor/pgen_estimate_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
plot_title <- 'Normalized AA and NT sequence pgen compared between models'
plot_caption <- 'Subject-specific, combined and default models compared against each other. The models (except default) were trained using the\nproductive sequences. All models evaluated the same combined data file that is also used for training the combined model.\nSpearman correlation method is used for the upper-right plots.'
output_filename <- '~/Downloads/evaluated_seqs_productive_rplot.png'

# ----------
# ALL
# ----------
model_0 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_0/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_1 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_1/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_2 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_2/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_3 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_3/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_4 <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_4/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_combined <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_combined/pgen_estimate_all_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
model_default <- melt(data.frame(read.table('~/Downloads/evaluated_seqs/evaluate_igor/pgen_estimate_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc)))
plot_title <- 'Normalized AA and NT sequence pgen compared between models'
plot_caption <- 'Subject-specific, combined and default models compared against each other. The models (except default) were trained using all\nsequences (productive and unproductive). All models evaluated the same combined data file that is also used for training the\ncombined model. Spearman correlation method is used for the upper-right plots.'
output_filename <- '~/Downloads/evaluated_seqs_all_rplot.png'
