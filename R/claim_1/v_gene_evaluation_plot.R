# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)
library(scales)
library(gtools)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function(origin, y, x, combination, name) {
  v_genes <- unique(origin['v_gene_choice'])
  v_genes$Y <- 0
  v_genes$X <- 0
  v_genes$Y <- apply(v_genes, 1, function(i) {
    length(y[y$row_id %in% origin[origin$v_gene_choice == i['v_gene_choice'], 'row_id'], ])
  })
  v_genes$X <- apply(v_genes, 1, function(i) {
    length(x[x$row_id %in% origin[origin$v_gene_choice == i['v_gene_choice'], 'row_id'], ])
  })
  v_genes$Y <- v_genes$Y / sum(v_genes$Y)
  v_genes$X <- v_genes$X / sum(v_genes$X)
  v_genes <- as.data.frame(v_genes[mixedorder(gsub("[\\-]", '.', gsub("[a-zA-Z]", '', v_genes$v_gene_choice))), ])
  v_genes$id <- c(1:nrow(v_genes))
  v_genes <- v_genes[order(v_genes$v_gene_choice), ]
  rownames(v_genes) <- NULL
  v_genes$combination <- combination
  v_genes$name <- name
  names(v_genes)[1] <- 'V-gene'
  v_genes <- v_genes[sample(nrow(v_genes)), ]
  return (v_genes)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'V-gene occurrence (percent)'
plot_x <- 'V-gene occurrence (percent)'
legend <- 'V-gene'
output_filename <- '~/Downloads/claim_1/v_gene_evaluation_plot.png'

# ----------
# MODEL DATA
# ----------
cc1 <- c(NA, "NULL", NA, rep("NULL", 2), NA, "NULL")
cc2 <- c(NA, "NULL", NA, "NULL")

model_0_origin <- data.frame(read.table('~/Downloads/claim_1/models/subject_0/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
model_0_productive <- data.frame(read.table('~/Downloads/claim_1/models/subject_0/converted_full_length_productive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_0_unproductive <- data.frame(read.table('~/Downloads/claim_1/models/subject_0/converted_full_length_unproductive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_0_all <- data.frame(read.table('~/Downloads/claim_1/models/subject_0/converted_full_length.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_0_p_a <- process_model(model_0_origin, model_0_productive, model_0_all, 'in (left) - all (bottom)', 'subject 1')
model_0_p_u <- process_model(model_0_origin, model_0_productive, model_0_unproductive, 'in (left) - out (bottom)', 'subject 1')
model_0_a_u <- process_model(model_0_origin, model_0_all, model_0_unproductive, 'all (left) - out (bottom)', 'subject 1')
model_0 <- as.data.frame(do.call("rbind", list(model_0_p_a, model_0_p_u, model_0_a_u)))
rm(model_0_origin, model_0_productive, model_0_unproductive, model_0_all, model_0_p_a, model_0_p_u, model_0_a_u)

model_1_origin <- data.frame(read.table('~/Downloads/claim_1/models/subject_1/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
model_1_productive <- data.frame(read.table('~/Downloads/claim_1/models/subject_1/converted_full_length_productive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_1_unproductive <- data.frame(read.table('~/Downloads/claim_1/models/subject_1/converted_full_length_unproductive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_1_all <- data.frame(read.table('~/Downloads/claim_1/models/subject_1/converted_full_length.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_1_p_a <- process_model(model_1_origin, model_1_productive, model_1_all, 'in (left) - all (bottom)', 'subject 2')
model_1_p_u <- process_model(model_1_origin, model_1_productive, model_1_unproductive, 'in (left) - out (bottom)', 'subject 2')
model_1_a_u <- process_model(model_1_origin, model_1_all, model_1_unproductive, 'all (left) - out (bottom)', 'subject 2')
model_1 <- as.data.frame(do.call("rbind", list(model_1_p_a, model_1_p_u, model_1_a_u)))
rm(model_1_origin, model_1_productive, model_1_unproductive, model_1_all, model_1_p_a, model_1_p_u, model_1_a_u)

model_2_origin <- data.frame(read.table('~/Downloads/claim_1/models/subject_2/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
model_2_productive <- data.frame(read.table('~/Downloads/claim_1/models/subject_2/converted_full_length_productive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_2_unproductive <- data.frame(read.table('~/Downloads/claim_1/models/subject_2/converted_full_length_unproductive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_2_all <- data.frame(read.table('~/Downloads/claim_1/models/subject_2/converted_full_length.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_2_p_a <- process_model(model_2_origin, model_2_productive, model_2_all, 'in (left) - all (bottom)', 'subject 3')
model_2_p_u <- process_model(model_2_origin, model_2_productive, model_2_unproductive, 'in (left) - out (bottom)', 'subject 3')
model_2_a_u <- process_model(model_2_origin, model_2_all, model_2_unproductive, 'all (left) - out (bottom)', 'subject 3')
model_2 <- as.data.frame(do.call("rbind", list(model_2_p_a, model_2_p_u, model_2_a_u)))
rm(model_2_origin, model_2_productive, model_2_unproductive, model_2_all, model_2_p_a, model_2_p_u, model_2_a_u)

model_3_origin <- data.frame(read.table('~/Downloads/claim_1/models/subject_3/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
model_3_productive <- data.frame(read.table('~/Downloads/claim_1/models/subject_3/converted_full_length_productive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_3_unproductive <- data.frame(read.table('~/Downloads/claim_1/models/subject_3/converted_full_length_unproductive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_3_all <- data.frame(read.table('~/Downloads/claim_1/models/subject_3/converted_full_length.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_3_p_a <- process_model(model_3_origin, model_3_productive, model_3_all, 'in (left) - all (bottom)', 'control 1')
model_3_p_u <- process_model(model_3_origin, model_3_productive, model_3_unproductive, 'in (left) - out (bottom)', 'control 1')
model_3_a_u <- process_model(model_3_origin, model_3_all, model_3_unproductive, 'all (left) - out (bottom)', 'control 1')
model_3 <- as.data.frame(do.call("rbind", list(model_3_p_a, model_3_p_u, model_3_a_u)))
rm(model_3_origin, model_3_productive, model_3_unproductive, model_3_all, model_3_p_a, model_3_p_u, model_3_a_u)

model_4_origin <- data.frame(read.table('~/Downloads/claim_1/models/subject_4/converted_CDR3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc1))
model_4_productive <- data.frame(read.table('~/Downloads/claim_1/models/subject_4/converted_full_length_productive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_4_unproductive <- data.frame(read.table('~/Downloads/claim_1/models/subject_4/converted_full_length_unproductive.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_4_all <- data.frame(read.table('~/Downloads/claim_1/models/subject_4/converted_full_length.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc2))
model_4_p_a <- process_model(model_4_origin, model_4_productive, model_4_all, 'in (left) - all (bottom)', 'control 2')
model_4_p_u <- process_model(model_4_origin, model_4_productive, model_4_unproductive, 'in (left) - out (bottom)', 'control 2')
model_4_a_u <- process_model(model_4_origin, model_4_all, model_4_unproductive, 'all (left) - out (bottom)', 'control 2')
model_4 <- as.data.frame(do.call("rbind", list(model_4_p_a, model_4_p_u, model_4_a_u)))
rm(model_4_origin, model_4_productive, model_4_unproductive, model_4_all, model_4_p_a, model_4_p_u, model_4_a_u)

models <- do.call('rbind', list(model_0, model_1, model_2, model_3, model_4))
rm(cc1, cc2, model_0, model_1, model_2, model_3, model_4)

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
  scale_x_continuous(
    name = plot_x,
    trans = 'sqrt',
    limits = c(0, 0.12),
    labels = percent
  ) +
  scale_y_continuous(
    name = plot_y,
    trans = 'sqrt',
    limits = c(0, 0.12),
    labels = percent
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

jpeg(output_filename, width = 4000, height = 4000, res = 300)
eval_compare
dev.off()
