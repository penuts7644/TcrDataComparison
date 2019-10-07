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
process_model <- function(marginals, params) {
  data <- data.frame(gene = character(0), value = character(0))
  data$gene <- as.character(data$gene)
  data$value <- as.character(data$value)
  marginals <- strsplit(sub('%', '', readLines(marginals)[4]), ',')[[1]] # [4]
  params <- readLines(params)[3:149] # [3:170]
  for (i in 1:length(params)) {
    line <- strsplit(sub('%', '', params[i]), ';')[[1]][c(1, 3)]
    data[i, ] <- c(line[1], marginals[as.numeric(line[2]) + 1])
  }
  data <- data[mixedorder(gsub("[\\-]", '.', gsub("[a-zA-Z]", '', data$gene))), ]
  data$gene <- as.factor(data$gene)
  data$value <- as.numeric(data$value)
  data$id <- c(1:nrow(data))
  data <- data[sample(nrow(data)), ]
  return (data)
}

process_data <- function(y, x, combination, name) {
  data <- merge(y, x, by = c('gene', 'id'))
  names(data) <- c('gene', 'id', 'X', 'Y')
  data$combination <- combination
  data$name <- name
  data <- data[sample(nrow(data)), ]
  return (data)
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
model_0_productive <- process_model('~/Downloads/claim_1/models/subject_0/productive_marginals.txt', '~/Downloads/claim_1/models/subject_0/productive_params.txt')
model_0_unproductive <- process_model('~/Downloads/claim_1/models/subject_0/unproductive_marginals.txt', '~/Downloads/claim_1/models/subject_0/productive_params.txt')
model_0_all <- process_model('~/Downloads/claim_1/models/subject_0/all_marginals.txt', '~/Downloads/claim_1/models/subject_0/productive_params.txt')
model_0_p_a <- process_data(model_0_productive, model_0_all, 'in (left) - all (bottom)', 'subject 1')
model_0_p_u <- process_data(model_0_productive, model_0_unproductive, 'in (left) - out (bottom)', 'subject 1')
model_0_a_u <- process_data(model_0_all, model_0_unproductive, 'all (left) - out (bottom)', 'subject 1')
model_0 <- as.data.frame(do.call("rbind", list(model_0_p_a, model_0_p_u, model_0_a_u)))
rm(model_0_productive, model_0_unproductive, model_0_all, model_0_p_a, model_0_p_u, model_0_a_u)

model_1_productive <- process_model('~/Downloads/claim_1/models/subject_1/productive_marginals.txt', '~/Downloads/claim_1/models/subject_1/productive_params.txt')
model_1_unproductive <- process_model('~/Downloads/claim_1/models/subject_1/unproductive_marginals.txt', '~/Downloads/claim_1/models/subject_1/productive_params.txt')
model_1_all <- process_model('~/Downloads/claim_1/models/subject_1/all_marginals.txt', '~/Downloads/claim_1/models/subject_1/productive_params.txt')
model_1_p_a <- process_data(model_1_productive, model_1_all, 'in (left) - all (bottom)', 'subject 2')
model_1_p_u <- process_data(model_1_productive, model_1_unproductive, 'in (left) - out (bottom)', 'subject 2')
model_1_a_u <- process_data(model_1_all, model_1_unproductive, 'all (left) - out (bottom)', 'subject 2')
model_1 <- as.data.frame(do.call("rbind", list(model_1_p_a, model_1_p_u, model_1_a_u)))
rm(model_1_productive, model_1_unproductive, model_1_all, model_1_p_a, model_1_p_u, model_1_a_u)

model_2_productive <- process_model('~/Downloads/claim_1/models/subject_2/productive_marginals.txt', '~/Downloads/claim_1/models/subject_2/productive_params.txt')
model_2_unproductive <- process_model('~/Downloads/claim_1/models/subject_2/unproductive_marginals.txt', '~/Downloads/claim_1/models/subject_2/productive_params.txt')
model_2_all <- process_model('~/Downloads/claim_1/models/subject_2/all_marginals.txt', '~/Downloads/claim_1/models/subject_2/productive_params.txt')
model_2_p_a <- process_data(model_2_productive, model_2_all, 'in (left) - all (bottom)', 'subject 3')
model_2_p_u <- process_data(model_2_productive, model_2_unproductive, 'in (left) - out (bottom)', 'subject 3')
model_2_a_u <- process_data(model_2_all, model_2_unproductive, 'all (left) - out (bottom)', 'subject 3')
model_2 <- as.data.frame(do.call("rbind", list(model_2_p_a, model_2_p_u, model_2_a_u)))
rm(model_2_productive, model_2_unproductive, model_2_all, model_2_p_a, model_2_p_u, model_2_a_u)

model_3_productive <- process_model('~/Downloads/claim_1/models/subject_3/productive_marginals.txt', '~/Downloads/claim_1/models/subject_3/productive_params.txt')
model_3_unproductive <- process_model('~/Downloads/claim_1/models/subject_3/unproductive_marginals.txt', '~/Downloads/claim_1/models/subject_3/productive_params.txt')
model_3_all <- process_model('~/Downloads/claim_1/models/subject_3/all_marginals.txt', '~/Downloads/claim_1/models/subject_3/productive_params.txt')
model_3_p_a <- process_data(model_3_productive, model_3_all, 'in (left) - all (bottom)', 'control 1')
model_3_p_u <- process_data(model_3_productive, model_3_unproductive, 'in (left) - out (bottom)', 'control 1')
model_3_a_u <- process_data(model_3_all, model_3_unproductive, 'all (left) - out (bottom)', 'control 1')
model_3 <- as.data.frame(do.call("rbind", list(model_3_p_a, model_3_p_u, model_3_a_u)))
rm(model_3_productive, model_3_unproductive, model_3_all, model_3_p_a, model_3_p_u, model_3_a_u)

model_4_productive <- process_model('~/Downloads/claim_1/models/subject_4/productive_marginals.txt', '~/Downloads/claim_1/models/subject_4/productive_params.txt')
model_4_unproductive <- process_model('~/Downloads/claim_1/models/subject_4/unproductive_marginals.txt', '~/Downloads/claim_1/models/subject_4/productive_params.txt')
model_4_all <- process_model('~/Downloads/claim_1/models/subject_4/all_marginals.txt', '~/Downloads/claim_1/models/subject_4/productive_params.txt')
model_4_p_a <- process_data(model_4_productive, model_4_all, 'in (left) - all (bottom)', 'control 2')
model_4_p_u <- process_data(model_4_productive, model_4_unproductive, 'in (left) - out (bottom)', 'control 2')
model_4_a_u <- process_data(model_4_all, model_4_unproductive, 'all (left) - out (bottom)', 'control 2')
model_4 <- as.data.frame(do.call("rbind", list(model_4_p_a, model_4_p_u, model_4_a_u)))
rm(model_4_productive, model_4_unproductive, model_4_all, model_4_p_a, model_4_p_u, model_4_a_u)

models <- do.call('rbind', list(model_0, model_1, model_2, model_3, model_4))
rm(model_0, model_1, model_2, model_3, model_4)

for (j in unique(models$gene)) {
  if (sum(models[models$gene == j, c('Y', 'X')]) == 0) {
    models <- models[models$gene != j, ]
  }
}

# ----------
# MAKING THE PLOTS
# ----------
eval_compare <-
  ggplot(
    data = models,
    aes(
      x = X,
      y = Y,
      color = reorder(gene, id)
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

jpeg(output_filename, width = 4000, height = 4000, res = 300)
eval_compare
dev.off()
