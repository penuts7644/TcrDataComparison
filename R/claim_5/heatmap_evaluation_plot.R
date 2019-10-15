# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
library(gtools)
library(viridis)
library(shadowtext)
MODELS <- c('subject 1', 'subject 2', 'subject 3', 'control 1', 'control 2')
SUBSAMPLES <- c('all', '50000', '10000', '5000', '1000', '500', '100')
CORRELATION_METHOD <- 'pearson' # or 'spearman'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
process_model <- function() {
  cc <- c(NA, "NULL", NA, rep("NULL", 4), rep(NA, 2))
  sn <- list(`subject 1` = 'subject_0', `subject 2` = 'subject_1', `subject 3` = 'subject_2', `control 1` = 'subject_3', `control 2` = 'subject_4')
  comb <- expand.grid(MODELS, SUBSAMPLES, stringsAsFactors = FALSE)
  comb <- transform(comb, Var3 = paste(Var1, Var2, sep=" - "))[, 3]
  corr <- matrix(rep(NA, len=length(comb) * length(comb)), nrow = length(comb))
  diag(corr) <- 1
  colnames(corr) <- comb
  rownames(corr) <- comb
  for(i in 1:nrow(corr)) {
    row_info <- unlist(strsplit(rownames(corr)[i], " - ", fixed = TRUE))
    model_row <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', sn[row_info[1]], '/', row_info[2], '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
    model_row <- melt(model_row, id.vars = 'row_id')
    for(j in 1:ncol(corr)) {
      if (i == j) {
        next
      } else if (is.na(corr[j, i])) {
        col_info <- unlist(strsplit(colnames(corr)[j], " - ", fixed = TRUE))
        model_col <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', sn[col_info[1]], '/', col_info[2], '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
        model_col <- melt(model_col, id.vars = 'row_id')
        models <- merge(model_row, model_col, by=c('row_id', 'variable'))
        models <- models[!(models$row_id %in% models[is.na(models[, 3:4]), 'row_id']), ][, 2:4]
        corr[i, j] <- cor(models[, 2], models[, 3], method = CORRELATION_METHOD)
      } else {
        corr[i, j] <- corr[j, i]
      }
    }
  }
  corr <- melt(corr)
  return (corr)
}

# ----------
# EXTRACTING THE DATA
# ----------
correlations <- process_model()
correlations$Var1 <- factor(as.factor(correlations$Var1), levels = unique(mixedsort(as.vector(correlations$Var1))))
correlations$Var2 <- factor(as.factor(correlations$Var2), levels = unique(mixedsort(as.vector(correlations$Var2))))

plot_y <- 'Subject - subsample combination'
plot_x <- 'Subject - subsample combination'
legend <- paste('Correlation score (', CORRELATION_METHOD, ')', sep = '')
output_filename <- paste('~/Downloads/claim_5/subsample_heatmap_plot_', CORRELATION_METHOD, '.png', sep = '')

# ----------
# MAKING THE PLOTS
# ----------
heat_compare <-
  ggplot(
    correlations,
    aes(
      x = Var1,
      y = Var2
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
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
heat_compare
dev.off()
