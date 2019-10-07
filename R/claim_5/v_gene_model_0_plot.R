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
process_model <- function(marginals, params, subsample) {
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
  data$subsample <- subsample
  data <- data[sample(nrow(data)), ]
  return (data)
}

# ----------
# PLOT VARIABLES
# ----------
plot_y <- 'V-gene probability score'
plot_x <- 'V-gene'
legend <- 'Subsample'
output_filename <- '~/Downloads/claim_5/v_gene_model_0_plot.png'

# ----------
# MODEL DATA
# ----------
model_0_all <- process_model('~/Downloads/claim_5/models/subject_0/all/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/all/all_params.txt', 'all')
model_0_50000 <- process_model('~/Downloads/claim_5/models/subject_0/50000/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/50000/all_params.txt', '50000')
model_0_10000 <- process_model('~/Downloads/claim_5/models/subject_0/10000/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/10000/all_params.txt', '10000')
model_0_5000 <- process_model('~/Downloads/claim_5/models/subject_0/5000/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/5000/all_params.txt', '5000')
model_0_1000 <- process_model('~/Downloads/claim_5/models/subject_0/1000/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/1000/all_params.txt', '1000')
model_0_500 <- process_model('~/Downloads/claim_5/models/subject_0/500/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/500/all_params.txt', '500')
model_0_100 <- process_model('~/Downloads/claim_5/models/subject_0/100/all_marginals.txt', '~/Downloads/claim_5/models/subject_0/100/all_params.txt', '100')
models <- as.data.frame(do.call("rbind", list(model_0_all, model_0_50000, model_0_10000, model_0_5000, model_0_1000, model_0_500, model_0_100)))
rm(model_0_all, model_0_50000, model_0_10000, model_0_5000, model_0_1000, model_0_500, model_0_100)

for (j in unique(models$gene)) {
  if (sum(models[models$gene == j, 'value']) == 0) {
    models <- models[models$gene != j, ]
  }
}
models$subsample <- factor(as.factor(models$subsample), levels = c("100", "500", "1000", "5000", "10000", "50000", "all"))

# ----------
# MAKING THE PLOTS
# ----------
eval_compare <-
  ggplot(
    data = models,
    aes(
      x = reorder(gene, id),
      y = value,
      fill = subsample
    )
  ) +
  geom_bar(
    stat = 'identity'
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
      size = 11,
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
  labs(
    x = plot_x
  ) +
  scale_y_sqrt(
    name = plot_y,
    limits = c(0, max(models$value))
  ) +
  scale_fill_brewer(
    palette = 'Dark2'
  ) +
  guides(
    fill = guide_legend(
      title = legend,
      nrow = 1
    )
  )

jpeg(output_filename, width = 4000, height = 2000, res = 300)
eval_compare
dev.off()
