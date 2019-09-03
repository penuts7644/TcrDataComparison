# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_entropies <- function(data) {
  subject_list <- c('0', 'disease 1', '1', 'disease 2', '2', 'disease 3', '3', 'control 1', '4', 'control 2')
  frame_list <- c('productive', 'in', 'unproductive', 'out/stop', 'all', 'in/out/stop')
  tmp <- subset(data, grepl('0', Var1) & grepl('0', Var2) & Var1 != Var2)
  for (i in c('1', '2', '3', '4')) {
    tmp <- rbind(tmp, subset(data, grepl(i, Var1) & grepl(i, Var2) & Var1 != Var2))
  }
  for (i in 1:length(subject_list)) {
    tmp$L2[grepl(subject_list[i], tmp$Var1)] <- subject_list[i+1]
  }
  for (i in 1:length(frame_list)) {
    tmp$L3[grepl(frame_list[i], tmp$Var1)] <- frame_list[i+1]
    tmp$L4[grepl(frame_list[i], tmp$Var2)] <- frame_list[i+1]
  }
  names(tmp)[4] <- "Event level"
  names(tmp)[5] <- "Model ID"
  names(tmp)[6] <- "Frame ID"
  names(tmp)[7] <- "Frame type"
  tmp <- tmp[sample(nrow(tmp)), ]
  return (tmp)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'KL divergence (bits)'
plot_x <- 'Immune receptor component or event'
output_filename <- '~/Downloads/claim_1/model_entropies_plot.png'

# ----------
# MODEL DATA
# ----------
V <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_V.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
V_trim_3 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_V_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_D.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D_trim_3 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_D_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D_trim_5 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_D_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
J <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_J.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
J_trim_5 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_J_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
insert_length_VD <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_insert_length_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
insert_length_DJ <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_insert_length_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
dinuc_markov_VD <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_dinuc_markov_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
dinuc_markov_DJ <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_entropy_dinuc_markov_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
entropies <- process_entropies(
  melt(list(
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
)
rm(V, D, J, V_trim_3, D_trim_3, D_trim_5, J_trim_5, insert_length_VD, insert_length_DJ, dinuc_markov_VD, dinuc_markov_DJ)

# ----------
# MAKING THE PLOTS
# ----------
entr_compare <-
  ggplot(
    entropies,
    aes(
      x = `Event level`,
      y = value
    )
  ) +
  geom_boxplot(
    alpha = 0.5,
    width = 0.6,
    outlier.size = 0.2,
    outlier.shape = NA
  ) +
  geom_jitter(
    height = 0,
    width = .2,
    size = 4,
    aes(
      color = `Frame type`,
      shape = `Frame type`
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
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
    legend.title = element_text(
      size = 18
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  labs(
    y = plot_y,
    x = plot_x
  ) +
  scale_color_brewer(
    palette = 'Dark2'
  ) +
  facet_grid(
    rows = vars(`Model ID`),
    cols = vars(`Frame ID`)
  )

jpeg(output_filename, width = 3000, height = 4000, res = 300)
entr_compare
dev.off()
