# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)
library(scales)

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_entropies <- function(data, entropies) {
  subject_list <- c('0', 'subject 1', '1', 'subject 2', '2', 'subject 3', '3', 'control 1', '4', 'control 2')
  frame_list <- c('productive', 'in', 'unproductive', 'out', 'all', 'all')
  tmp <- subset(data, grepl('0', Var1) & grepl('0', Var2) & Var1 != Var2)
  for (i in c('1', '2', '3', '4')) {
    tmp <- rbind(tmp, subset(data, grepl(i, Var1) & grepl(i, Var2) & Var1 != Var2))
  }
  for (i in seq(1, length(subject_list), 2)) {
    tmp$L2[grepl(subject_list[i], tmp$Var1)] <- subject_list[i+1]
  }
  for (i in seq(1, length(frame_list), 2)) {
    tmp$L3[grepl(frame_list[i], tmp$Var1)] <- frame_list[i+1]
    tmp$L4[grepl(frame_list[i], tmp$Var2)] <- frame_list[i+1]
  }
  tmp <- transform(tmp, L5 = paste(L3, L4, sep=" - "))

  tmp2 <- merge(entropies, entropies, by = 'event', sort = TRUE)
  tmp2 <- merge(tmp2, entropies, by = 'event', sort = TRUE)
  tmp2 <- subset(tmp2, `id.x` != `id.y` & `id.x` != `id` & `id.y` != `id`)
  tmp2 <- subset(tmp2, `id.x` != `id.y` & `id.x` != `id` & `id.y` != `id`)
  tmp2 <- rbind(
    subset(tmp2, grepl('0', `id.x`) & grepl('0', `id.y`) & grepl('0', `id`)),
    subset(tmp2, grepl('1', `id.x`) & grepl('1', `id.y`) & grepl('1', `id`)),
    subset(tmp2, grepl('2', `id.x`) & grepl('2', `id.y`) & grepl('2', `id`)),
    subset(tmp2, grepl('3', `id.x`) & grepl('3', `id.y`) & grepl('3', `id`)),
    subset(tmp2, grepl('4', `id.x`) & grepl('4', `id.y`) & grepl('4', `id`))
  )
  tmp2 <- subset(tmp2, grepl('productive', `id.x`) & grepl('unproductive', `id.y`) & grepl('all', `id`))
  tmp2 <- subset(tmp2, `event` != 'total')
  tmp2$value <- apply(tmp2, 1, function(x) {
    return ((as.numeric(x['entropy.x']) + as.numeric(x['entropy.y']) + as.numeric(x['entropy'])) / 3)
  })
  for (i in seq(1, length(subject_list), 2)) {
    tmp2$L2[grepl(subject_list[i], tmp2$`id.x`)] <- subject_list[i+1]
  }
  tmp2 <- tmp2[, c(8, 1, 9)]
  tmp2$L5 <- 'average'
  names(tmp2) <- c('value', 'L1', 'L2', 'L5')
  tmp2$L1 <- as.character(tmp2$L1)
  tmp2$L1[tmp2$L1 == 'V'] <- 'V-gene'
  tmp2$L1[tmp2$L1 == 'D'] <- 'D-gene'
  tmp2$L1[tmp2$L1 == 'J'] <- 'J-gene'
  tmp2$L1[tmp2$L1 == 'V_trim_3'] <- 'V-gene trim (3)'
  tmp2$L1[tmp2$L1 == 'D_trim_3'] <- 'D-gene trim (3)'
  tmp2$L1[tmp2$L1 == 'D_trim_5'] <- 'D-gene trim (5)'
  tmp2$L1[tmp2$L1 == 'J_trim_5'] <- 'J-gene trim (5)'
  tmp2$L1[tmp2$L1 == 'insert_length_VD'] <- 'VD insert length'
  tmp2$L1[tmp2$L1 == 'insert_length_DJ'] <- 'DJ insert length'
  tmp2$L1[tmp2$L1 == 'dinuc_markov_VD'] <- 'DinucMarkov VD'
  tmp2$L1[tmp2$L1 == 'dinuc_markov_DJ'] <- 'DinucMarkov DJ'
  tmp2$L1 <- as.factor(tmp2$L1)
  tmp2 <- tmp2[order(-tmp2$value), ]

  tmp <- tmp[, c('value', 'L1', 'L2', 'L5')]
  tmp$L5[tmp$L5 == 'out - in'] <- 'in - out'
  tmp$L5[tmp$L5 == 'all - in'] <- 'in - all'
  tmp$L5[tmp$L5 == 'all - out'] <- 'out - all'
  tmp <- tmp[sample(nrow(tmp)), ]
  tmp <- rbind(tmp2, tmp)
  names(tmp)[2] <- "Event level"
  names(tmp)[3] <- "Model ID"
  names(tmp)[4] <- "Model comparison"
  return (tmp)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'KL divergence (bits)'
plot_x <- 'Immune receptor component/event'
output_filename <- '~/Downloads/claim_1/model_entropies_plot.png'

# ----------
# MODEL DATA
# ----------
V <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_V.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
V_trim_3 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_V_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_D.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D_trim_3 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_D_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D_trim_5 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_D_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
J <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_J.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
J_trim_5 <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_J_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
insert_length_VD <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_insert_length_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
insert_length_DJ <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_insert_length_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
dinuc_markov_VD <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_dinuc_markov_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
dinuc_markov_DJ <- data.matrix(read.table('~/Downloads/claim_1/entropies/calc_model_dkl_dinuc_markov_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
entropies <- data.frame(read.table('~/Downloads/claim_1/entropies/calc_model_entropy.tsv', header=TRUE, sep='\t', check.names=FALSE))

dkl <- process_entropies(
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
  )),
  entropies
)
rm(V, D, J, V_trim_3, D_trim_3, D_trim_5, J_trim_5, insert_length_VD, insert_length_DJ, dinuc_markov_VD, dinuc_markov_DJ, entropies)

# ----------
# MAKING THE PLOTS
# ----------
entr_compare <-
  ggplot(
    dkl,
    aes(
      x = factor(`Event level`, levels = unique(`Event level`)),
      y = value
    )
  ) +
  geom_jitter(
    height = 0,
    width = .3,
    size = 3,
    aes(
      color = `Model comparison`
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
  scale_y_continuous(
    name = plot_y,
    trans = 'log10',
    labels = trans_format("log10", math_format(1^.x))
  ) +
  labs(
    x = plot_x
  ) +
  scale_color_manual(
    values = c('#000000', '#1b9e77', '#d95f02', '#7570b3')
  ) +
  facet_wrap(
    vars(`Model ID`),
    nrow = 5
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
entr_compare
dev.off()
