# ----------
# GENERAL
# ----------
library(GGally)
library(reshape2)
library(gtools)
library(viridis)
MODELS <- c('subject 1', 'subject 2', 'subject 3', 'control 1', 'control 2')
SUBSAMPLES <- c('all', '50000', '10000', '5000', '1000', '500', '100')

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
    model_row <- melt(model_row[order(model_row$row_id), ][, 2:3])
    for(j in 1:ncol(corr)) {
      if (i == j) {
        next
      } else if (is.na(corr[j, i])) {
        col_info <- unlist(strsplit(colnames(corr)[j], " - ", fixed = TRUE))
        model_col <- data.frame(read.table(paste('~/Downloads/claim_5/evaluations/', sn[col_info[1]], '/', col_info[2], '/pgen_estimate_all_CDR3.tsv', sep = ''), header=TRUE, row.names=1, sep='\t', check.names=FALSE, colClasses=cc))
        model_col <- melt(model_col[order(model_col$row_id), ][, 2:3])
        models <- na.omit(cbind(model_row, model_col[2]))
        corr[i, j] <- cor(models[, 2], models[, 3], method = 'spearman')
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
plot_y <- 'Subject - subsample combination'
plot_x <- 'Subject - subsample combination'
legend <- 'Spearman correlation score'
output_filename <- '~/Downloads/claim_5/subsample_heatmap_plot.png'

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
  labs(
    y = plot_y,
    x = plot_x
  ) +
  scale_fill_viridis(
    name = legend,
    option = 'plasma'
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
print(heat_compare)
dev.off()
