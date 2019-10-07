# ----------
# GENERAL
# ----------
library(reshape2)
library(ggplot2)
library(seqinr)
library(viridis)
library(gtools)
library(shadowtext)
MODELS <- c('subject 1', 'subject 2', 'subject 3', 'control 1', 'control 2')
COMPARE_METHOD <- '-' # or '/'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_v_reference <- function(params) {
  v_genes <- data.frame(gene = character(0), id = character(0))
  v_genes$gene <- as.character(v_genes$gene)
  v_genes$id <- as.character(v_genes$id)
  params_tmp <- readLines(params)[3:149] # [3:170]
  for (i in 1:length(params_tmp)) {
    line <- strsplit(sub('%', '', params_tmp[i]), ';')[[1]][c(1, 3)]
    v_genes[i, ] <- c(line[1], as.numeric(line[2]) + 1)
  }
  v_genes <- v_genes[mixedorder(gsub("[\\-]", '.', gsub("[a-zA-Z]", '', v_genes$gene))), ]
  v_genes$gene <- as.factor(v_genes$gene)
  v_genes$id <- as.numeric(v_genes$id)
  v_genes$order <- c(1:nrow(v_genes))
  v_genes <- v_genes[order(v_genes['id']), ]
  return (v_genes)
}

process_model <- function(marginals, params) {
  data <- data.frame()
  marginals_tmp <- readLines(marginals)[45:338] # [45:380]
  for (i in 1:length(marginals_tmp)) {
    if (grepl('%', marginals_tmp[i])) {
      data <- rbind(data, as.numeric(as.vector(strsplit(sub('%', '', marginals_tmp[i]), ','))[[1]]))
    }
  }
  params_tmp <- readLines(params)[194:214] # [193:213]
  tmp_params <- data.frame()
  for (i in 1:length(params_tmp)) {
    tmp_params <- rbind(tmp_params, as.numeric(as.vector(strsplit(sub('%', '', params_tmp[i]), ';'))[[1]]))
  }
  colnames(data) <- tmp_params[order(tmp_params[2]), 1]
  return (data)
}

process_data <- function(origin, y, x, combination) {
  if (COMPARE_METHOD == '/') {
    new_data <- log10(as.matrix(y) / as.matrix(x))
  } else {
    new_data <- as.matrix(y) - as.matrix(x)
  }
  new_data[is.na(new_data)] <- 0
  new_data[is.infinite(new_data)] <- 0
  origin <- melt(cbind(origin, new_data), id.vars = c('gene', 'id', 'order'))
  origin[2] <- NULL
  names(origin) <- c('V-gene', 'order', 'V-gene trim (3)', 'value')
  origin$combination <- combination
  origin <- origin[sample(nrow(origin)), ]
  return (origin)
}

for (i in 1:length(MODELS)) {

  # ----------
  # PLOT VARIABLES
  # ----------
  plot_y <- paste('V-gene (', MODELS[i], ')', sep = '')
  plot_x <- paste('V-gene trim (3) (', MODELS[i], ')', sep = '')
  if (COMPARE_METHOD == '/') {
    legend <- 'Probability score difference (Y-axis divided by X-axis), log scale'
    output_filename <- paste('~/Downloads/claim_1/v_gene_trim3_models_plot_divide_', i - 1,'.png', sep = '')
  } else {
    legend <- 'Probability score difference (Y-axis minus X-axis)'
    output_filename <- paste('~/Downloads/claim_1/v_gene_trim3_models_plot_substract_', i - 1,'.png', sep = '')
  }
  

  # ----------
  # MODEL DATA
  # ----------
  ref_v_genes <- process_v_reference(paste('~/Downloads/claim_1/models/subject_', i - 1,'/all_params.txt', sep = ''))

  model_productive <- process_model(paste('~/Downloads/claim_1/models/subject_', i - 1,'/productive_marginals.txt', sep = ''), paste('~/Downloads/claim_1/models/subject_', i - 1,'/productive_params.txt', sep = ''))
  model_unproductive <- process_model(paste('~/Downloads/claim_1/models/subject_', i - 1,'/unproductive_marginals.txt', sep = ''), paste('~/Downloads/claim_1/models/subject_', i - 1,'/unproductive_params.txt', sep = ''))
  model_all <- process_model(paste('~/Downloads/claim_1/models/subject_', i - 1,'/all_marginals.txt', sep = ''), paste('~/Downloads/claim_1/models/subject_', i - 1,'/all_params.txt', sep = ''))
  model_p_a <- process_data(ref_v_genes, model_productive, model_all, 'in - all')
  model_p_u <- process_data(ref_v_genes, model_productive, model_unproductive, 'in - out')
  model_a_u <- process_data(ref_v_genes, model_all, model_unproductive, 'all - out')
  models <- as.data.frame(do.call("rbind", list(model_p_a, model_p_u, model_a_u)))
  rm(ref_v_genes, model_productive, model_unproductive, model_all, model_p_a, model_p_u, model_a_u)

  for (j in unique(models$`V-gene`)) {
    if (sum(models[models$`V-gene` == j, 'value']) == 0) {
      models <- models[models$`V-gene` != j, ]
    }
  }

  # ----------
  # MAKING THE PLOTS
  # ----------
  heat_compare <-
    ggplot(
      models,
      aes(
        x = `V-gene trim (3)`,
        y = reorder(`V-gene`, -order)
      )
    ) +
    geom_raster(
      aes(
        fill = value
      )
    ) +
    geom_shadowtext(
      aes(
        label = round(value, digits = 2)
      ),
      size = 1.4
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      legend.title = element_text(
        size = 18
      ),
      legend.text = element_text(
        size = 15,
        angle = 60,
        hjust = 1
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
        size = 9
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
      y = plot_y,
      x = plot_x
    ) +
    scale_fill_viridis(
      name = legend,
      option = 'plasma'
    ) +
    facet_wrap(
      vars(combination)
    )

  jpeg(output_filename, width = 4000, height = 3000, res = 300)
  print(heat_compare)
  dev.off()
}
