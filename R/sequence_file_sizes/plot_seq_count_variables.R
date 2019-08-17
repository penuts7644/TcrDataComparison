# ----------
# PLOT VARIABLES
# ----------
plot_title <- 'Training sequence count for subject-specific models'
output_filename <- '~/Downloads/process_5_files/sequence_counts_rplot.png'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((0.75 - 0.15) * ((x - min(x)) / (max(x) - min(x))) + 0.15)
}

# ----------
# DATA PROCESSING
# ----------
counts <- data.frame(read.table('~/Downloads/process_5_files/sequence_counts.tsv', header=TRUE, sep='\t', check.names=FALSE))
counts <- counts[counts$type != 'evaluated', ]
counts <- counts[counts$model != 'combined', ]
counts$id_f = factor(counts$id, levels=c('500','1000','5000','10000', '50000', '100000', 'all'))

counts[counts$id == 'all' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == 'all' & counts$type == 'productive', 'count'])
counts[counts$id == 'all' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == 'all' & counts$type == 'unproductive', 'count'])
counts[counts$id == 'all' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == 'all' & counts$type == 'all', 'count'])

counts[counts$id == '100000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '100000' & counts$type == 'productive', 'count'])
counts[counts$id == '100000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '100000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '100000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '100000' & counts$type == 'all', 'count'])

counts[counts$id == '50000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '50000' & counts$type == 'productive', 'count'])
counts[counts$id == '50000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '50000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '50000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '50000' & counts$type == 'all', 'count'])

counts[counts$id == '10000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '10000' & counts$type == 'productive', 'count'])
counts[counts$id == '10000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '10000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '10000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '10000' & counts$type == 'all', 'count'])

counts[counts$id == '5000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '5000' & counts$type == 'productive', 'count'])
counts[counts$id == '5000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '5000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '5000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '5000' & counts$type == 'all', 'count'])

counts[counts$id == '1000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '1000' & counts$type == 'productive', 'count'])
counts[counts$id == '1000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '1000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '1000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '1000' & counts$type == 'all', 'count'])

counts[counts$id == '500' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '500' & counts$type == 'productive', 'count'])
counts[counts$id == '500' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '500' & counts$type == 'unproductive', 'count'])
counts[counts$id == '500' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '500' & counts$type == 'all', 'count'])
