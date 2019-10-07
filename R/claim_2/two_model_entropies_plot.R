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
  subject_list <- c('2', 'subject 3', '3', 'control 1')
  tmp <- subset(data, Var1 != Var2 & Var1 %in% c('2', '3') & Var2 %in% c('2', '3'))
  for (i in seq(1, length(subject_list), 2)) {
    tmp$Var1[subject_list[i] == tmp$Var1] <- subject_list[i+1]
    tmp$Var2[subject_list[i] == tmp$Var2] <- subject_list[i+1]
  }
  tmp <- transform(tmp, L2 = paste(Var1, Var2, sep=" - "))

  tmp2 <- subset(entropies, id %in% c('2', '3'))
  tmp2 <- subset(tmp2, `event` != 'total')
  for (i in seq(1, length(subject_list), 2)) {
    tmp2$id[subject_list[i] == tmp2$id] <- subject_list[i+1]
  }
  tmp2$event <- as.character(tmp2$event)
  tmp2$event[tmp2$event == 'V'] <- 'V-gene'
  tmp2$event[tmp2$event == 'D'] <- 'D-gene'
  tmp2$event[tmp2$event == 'J'] <- 'J-gene'
  tmp2$event[tmp2$event == 'V_trim_3'] <- 'V-gene trim (3)'
  tmp2$event[tmp2$event == 'D_trim_3'] <- 'D-gene trim (3)'
  tmp2$event[tmp2$event == 'D_trim_5'] <- 'D-gene trim (5)'
  tmp2$event[tmp2$event == 'J_trim_5'] <- 'J-gene trim (5)'
  tmp2$event[tmp2$event == 'insert_length_VD'] <- 'VD insert length'
  tmp2$event[tmp2$event == 'insert_length_DJ'] <- 'DJ insert length'
  tmp2$event[tmp2$event == 'dinuc_markov_VD'] <- 'DinucMarkov VD'
  tmp2$event[tmp2$event == 'dinuc_markov_DJ'] <- 'DinucMarkov DJ'
  tmp2$event <- as.factor(tmp2$event)
  tmp2 <- tmp2[order(-tmp2$entropy), ]

  # tmp$L2[tmp$L2 == 'control 1 - subject 3'] <- 'subject 3 - control 1'
  for (i in 1:nrow(tmp)) {
    tmp[i, 'value'] <- round((tmp[i, 'value'] / tmp2[tmp2$event == tmp[i, 'L1'] & tmp2$id == tmp[i, 'Var1'], 'entropy']), digits = 2)
  }
  tmp <- tmp[, c('value', 'L1', 'L2')]
  tmp <- tmp[order(match(tmp[, 'L1'], tmp2[, 'event'])), ]
  names(tmp)[2] <- "Event level"
  names(tmp)[3] <- "Compare"
  return (tmp)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'Entropy increase (percent)'
plot_x <- 'Immune receptor component/event'
output_filename <- '~/Downloads/claim_2/two_model_entropies_plot.png'

# ----------
# MODEL DATA
# ----------
V <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_V.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
V_trim_3 <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_V_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_D.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D_trim_3 <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_D_trim_3.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
D_trim_5 <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_D_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
J <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_J.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
J_trim_5 <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_J_trim_5.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
insert_length_VD <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_insert_length_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
insert_length_DJ <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_insert_length_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
dinuc_markov_VD <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_dinuc_markov_VD.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
dinuc_markov_DJ <- data.matrix(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_dkl_dinuc_markov_DJ.tsv', header=TRUE, row.names=1, sep='\t', check.names=FALSE))
entropies <- data.frame(read.table('~/Downloads/claim_2/entropies/emerson/calc_model_entropy.tsv', header=TRUE, sep='\t', check.names=FALSE))

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
      y = value,
      fill = `Compare`,
      label = percent(value)
    )
  ) +
  geom_col(
    position = position_dodge2(
      width = 0.8,
      preserve = 'single'
    )
  ) +
  geom_text(
    position = position_dodge2(
      width = 0.9,
      preserve = 'single'
    ),
    angle = 90,
    vjust = 0.5,
    hjust = -0.1,
    size = 4
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
    trans = 'sqrt',
    limits = c(0, 2),
    labels = percent
  ) +
  labs(
    x = plot_x
  ) +
  scale_fill_manual(
    values = c('#1b9e77', '#d95f02')
  )

jpeg(output_filename, width = 2000, height = 2000, res = 300)
entr_compare
dev.off()
