# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)
library(scales)
library(seqinr)
library(gtools)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_v_reference <- function() {
  reference <- read.fasta(file = '~/Downloads/human_TRB/TRBV.fasta', seqtype = 'DNA', as.string = TRUE, strip.desc = TRUE, whole.header = TRUE)
  genes <- data.frame(v_gene_choice = character(0))
  genes$v_gene_choice <- as.character(genes$v_gene_choice)
  for (i in 1:length(reference)) {
    genes[i, ] <- c(strsplit(attr(x = reference[i], which = 'name'), '|', fixed = TRUE)[[1]][2])
  }
  genes$v_gene_choice <- as.factor(genes[mixedorder(gsub("[\\-]", '.', gsub("[a-zA-Z]", '', genes$v_gene_choice))), ])
  genes$id <- c(1:nrow(genes))
  genes <- genes[order(genes$v_gene_choice), ]
  rownames(genes) <- NULL
  return (genes)
}

process_model <- function(origin, y, x, combination, name) {
  origin <- as.data.frame(do.call("cbind", list(origin, y, x)))
  names(origin) <- c('V-gene', 'id', 'Y', 'X')
  origin$Y <- apply(origin, 1, function(i) {
    as.numeric(gsub('^%', '', i['Y']))
  })
  origin$X <- apply(origin, 1, function(i) {
    as.numeric(gsub('^%', '', i['X']))
  })
  origin$combination <- combination
  origin$name <- name
  origin <- origin[sample(nrow(origin)), ]
  return (origin)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'V-gene probability score'
plot_x <- 'V-gene probability score'
legend <- 'V-gene'
output_filename <- '~/Downloads/claim_1/v_gene_models_plot.png'

# ----------
# MODEL DATA
# ----------
ref_v_genes <- process_v_reference()

model_0_productive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_0/productive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_0_unproductive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_0/unproductive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_0_all <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_0/all_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_0_p_a <- process_model(ref_v_genes, model_0_productive, model_0_all, 'in (left) - all (bottom)', 'subject 1')
model_0_p_u <- process_model(ref_v_genes, model_0_productive, model_0_unproductive, 'in (left) - out (bottom)', 'subject 1')
model_0_a_u <- process_model(ref_v_genes, model_0_all, model_0_unproductive, 'all (left) - out (bottom)', 'subject 1')
model_0 <- as.data.frame(do.call("rbind", list(model_0_p_a, model_0_p_u, model_0_a_u)))
rm(model_0_productive, model_0_unproductive, model_0_all, model_0_p_a, model_0_p_u, model_0_a_u)

model_1_productive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_1/productive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_1_unproductive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_1/unproductive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_1_all <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_1/all_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_1_p_a <- process_model(ref_v_genes, model_1_productive, model_1_all, 'in (left) - all (bottom)', 'subject 2')
model_1_p_u <- process_model(ref_v_genes, model_1_productive, model_1_unproductive, 'in (left) - out (bottom)', 'subject 2')
model_1_a_u <- process_model(ref_v_genes, model_1_all, model_1_unproductive, 'all (left) - out (bottom)', 'subject 2')
model_1 <- as.data.frame(do.call("rbind", list(model_1_p_a, model_1_p_u, model_1_a_u)))
rm(model_1_productive, model_1_unproductive, model_1_all, model_1_p_a, model_1_p_u, model_1_a_u)

model_2_productive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_2/productive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_2_unproductive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_2/unproductive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_2_all <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_2/all_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_2_p_a <- process_model(ref_v_genes, model_2_productive, model_2_all, 'in (left) - all (bottom)', 'subject 3')
model_2_p_u <- process_model(ref_v_genes, model_2_productive, model_2_unproductive, 'in (left) - out (bottom)', 'subject 3')
model_2_a_u <- process_model(ref_v_genes, model_2_all, model_2_unproductive, 'all (left) - out (bottom)', 'subject 3')
model_2 <- as.data.frame(do.call("rbind", list(model_2_p_a, model_2_p_u, model_2_a_u)))
rm(model_2_productive, model_2_unproductive, model_2_all, model_2_p_a, model_2_p_u, model_2_a_u)

model_3_productive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_3/productive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_3_unproductive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_3/unproductive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_3_all <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_3/all_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_3_p_a <- process_model(ref_v_genes, model_3_productive, model_3_all, 'in (left) - all (bottom)', 'control 1')
model_3_p_u <- process_model(ref_v_genes, model_3_productive, model_3_unproductive, 'in (left) - out (bottom)', 'control 1')
model_3_a_u <- process_model(ref_v_genes, model_3_all, model_3_unproductive, 'all (left) - out (bottom)', 'control 1')
model_3 <- as.data.frame(do.call("rbind", list(model_3_p_a, model_3_p_u, model_3_a_u)))
rm(model_3_productive, model_3_unproductive, model_3_all, model_3_p_a, model_3_p_u, model_3_a_u)

model_4_productive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_4/productive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_4_unproductive <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_4/unproductive_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_4_all <- data.frame(t(read.table('~/Downloads/claim_1/models/subject_4/all_marginals.txt', sep=',', skip = 3, nrows = 1)), row.names = c(1:length(ref_v_genes$v_gene_choice)))
model_4_p_a <- process_model(ref_v_genes, model_4_productive, model_4_all, 'in (left) - all (bottom)', 'control 2')
model_4_p_u <- process_model(ref_v_genes, model_4_productive, model_4_unproductive, 'in (left) - out (bottom)', 'control 2')
model_4_a_u <- process_model(ref_v_genes, model_4_all, model_4_unproductive, 'all (left) - out (bottom)', 'control 2')
model_4 <- as.data.frame(do.call("rbind", list(model_4_p_a, model_4_p_u, model_4_a_u)))
rm(model_4_productive, model_4_unproductive, model_4_all, model_4_p_a, model_4_p_u, model_4_a_u)

models <- do.call('rbind', list(model_0, model_1, model_2, model_3, model_4))
rm(ref_v_genes, model_0, model_1, model_2, model_3, model_4)

# ----------
# MAKING THE PLOTS
# ----------
eval_compare <-
  ggplot(
    data = models,
    aes(
      x = X,
      y = Y,
      color = reorder(`V-gene`, id)
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = 'gray',
    size = 0.2
  ) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(
      size = 18
    ),
    legend.text = element_text(
      size = 15
    ),
    axis.title.x = element_text(
      size = 18
    ),
    axis.text.x = element_text(
      size = 13,
      angle = 60,
      hjust = 1
    ),
    axis.title.y = element_text(
      size = 18
    ),
    axis.text.y = element_text(
      size = 13
    ),
    strip.text.x = element_text(
      size = 15
    ),
    strip.text.y = element_text(
      size = 15
    ),
    legend.position = 'top',
    legend.direction = 'horizontal'
  ) +
  scale_x_sqrt(
    name = plot_x,
    limits = c(0, max(pmax(models$Y, models$X)))
  ) +
  scale_y_sqrt(
    name = plot_y,
    limits = c(0, max(pmax(models$Y, models$X)))
  ) +
  guides(
    colour = guide_legend(
      title = legend,
      override.aes = list(alpha = 1, size = 3)
    )
  ) +
  facet_grid(
    rows = vars(name),
    cols = vars(combination)
  )

jpeg(output_filename, width = 4000, height = 6000, res = 300)
eval_compare
dev.off()
