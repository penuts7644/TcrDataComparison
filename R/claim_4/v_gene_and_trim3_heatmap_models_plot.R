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
  return (as.matrix(data))
}

process_data <- function(origin, x, name) {
  x[is.na(x)] <- 0
  x[is.infinite(x)] <- 0
  origin <- melt(cbind(origin, x), id.vars = c('gene', 'id', 'order'))
  origin[2] <- NULL
  names(origin) <- c('V-gene', 'order', 'V-gene trim (3)', 'value')
  origin$name <- name
  origin <- origin[sample(nrow(origin)), ]
  return (origin)
}

for (i in 1:length(MODELS)) {

  # ----------
  # PLOT VARIABLES
  # ----------
  plot_y <- paste('V-gene (', MODELS[i], ')', sep = '')
  plot_x <- paste('V-gene trim (3) (', MODELS[i], ')', sep = '')
  legend <- 'Probability score'
  output_filename <- paste('~/Downloads/claim_4/v_gene_trim3_models_plot_', i - 1,'.png', sep = '')
  
  # ----------
  # MODEL DATA
  # ----------
  ref_v_genes <- process_v_reference(paste('~/Downloads/claim_4/models/subject_', i - 1,'/all/all_params.txt', sep = ''))
  
  model_0_100 <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/100/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/100/all_params.txt', sep = '')), '100')
  model_0_500 <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/500/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/500/all_params.txt', sep = '')), '500')
  model_0_1000 <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/1000/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/1000/all_params.txt', sep = '')), '1000')
  model_0_5000 <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/5000/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/5000/all_params.txt', sep = '')), '5000')
  model_0_10000 <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/10000/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/10000/all_params.txt', sep = '')), '10000')
  model_0_50000 <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/50000/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/50000/all_params.txt', sep = '')), '50000')
  model_0_all <- process_data(ref_v_genes, process_model(paste('~/Downloads/claim_4/models/subject_', i - 1,'/all/all_marginals.txt', sep = ''), paste('~/Downloads/claim_4/models/subject_', i - 1,'/all/all_params.txt', sep = '')), 'all')
  
  models <- as.data.frame(do.call("rbind", list(model_0_100, model_0_500, model_0_1000, model_0_5000, model_0_10000, model_0_50000, model_0_all)))
  rm(ref_v_genes, model_0_100, model_0_500, model_0_1000, model_0_5000, model_0_10000, model_0_50000, model_0_all)
  
  sum_models <- data.frame()
  for (j in unique(models$`V-gene`)) {
    if (sum(models[models$`V-gene` == j, 'value']) == 0) {
      models <- models[models$`V-gene` != j, ]
    }
    sum_models <- rbind(sum_models, data.frame(j, sum(models[models$`V-gene` == j, 'value'])))
  }
  names(sum_models) <- c('gene', 'sum')
  
  models <- models[models$`V-gene` %in% unique(sum_models[order(-sum_models$sum), 'gene'])[1:10], ]
  models$name <- factor(as.factor(models$name), levels = c("100", "500", "1000", "5000", "10000", "50000", "all"))
  
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
      y = plot_y,
      x = plot_x
    ) +
    scale_fill_viridis(
      name = legend,
      option = 'plasma'
    ) +
    facet_wrap(
      vars(name),
      ncol = 1
    )
  
  jpeg(output_filename, width = 1600, height = 6000, res = 300)
  print(heat_compare)
  dev.off()
}
