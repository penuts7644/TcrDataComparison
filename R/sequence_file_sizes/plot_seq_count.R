# ----------
# GENERAL
# ----------
library(ggplot2)

# ----------
# PLOT VARIABLES
# ----------
plot_title <- 'Training sequence count for subject-specific models'
output_filename <- '~/Downloads/process_5_files/sequence_counts_rplot.png'

# ----------
# FUNCTIONS FOR PRE-PROCESSING
# ----------
normalize <- function(x) {
  return ((0.75 - 0.15) * ((x - min(x)) / (max(x) - min(x))) + 0.15)
}

# ----------
# DATA PROCESSING
# ----------
counts <- data.frame(read.table('~/Downloads/process_5_files/sequence_counts.tsv', header=TRUE, sep='\t', check.names=FALSE))
counts <- counts[counts$type != 'evaluated', ]
counts <- counts[counts$model != 'combined', ]
counts$id_f = factor(counts$id, levels=c('500','1000','5000','10000', '50000', '100000', 'all'))

counts[counts$id == 'all' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == 'all' & counts$type == 'productive', 'count'])
counts[counts$id == 'all' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == 'all' & counts$type == 'unproductive', 'count'])
counts[counts$id == 'all' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == 'all' & counts$type == 'all', 'count'])

counts[counts$id == '100000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '100000' & counts$type == 'productive', 'count'])
counts[counts$id == '100000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '100000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '100000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '100000' & counts$type == 'all', 'count'])

counts[counts$id == '50000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '50000' & counts$type == 'productive', 'count'])
counts[counts$id == '50000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '50000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '50000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '50000' & counts$type == 'all', 'count'])

counts[counts$id == '10000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '10000' & counts$type == 'productive', 'count'])
counts[counts$id == '10000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '10000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '10000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '10000' & counts$type == 'all', 'count'])

counts[counts$id == '5000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '5000' & counts$type == 'productive', 'count'])
counts[counts$id == '5000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '5000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '5000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '5000' & counts$type == 'all', 'count'])

counts[counts$id == '1000' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '1000' & counts$type == 'productive', 'count'])
counts[counts$id == '1000' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '1000' & counts$type == 'unproductive', 'count'])
counts[counts$id == '1000' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '1000' & counts$type == 'all', 'count'])

counts[counts$id == '500' & counts$type == 'productive', 'norm'] <- normalize(counts[counts$id == '500' & counts$type == 'productive', 'count'])
counts[counts$id == '500' & counts$type == 'unproductive', 'norm'] <- normalize(counts[counts$id == '500' & counts$type == 'unproductive', 'count'])
counts[counts$id == '500' & counts$type == 'all', 'norm'] <- normalize(counts[counts$id == '500' & counts$type == 'all', 'count'])

# ----------
# MAKING THE PLOT
# ----------
size_compare <-
  ggplot(
    counts,
    aes(
      x = model,
      y = norm,
      fill = model
    )
  ) +
  geom_bar(
    stat = 'identity',
    width = 0.8,
    position = 'dodge'
  ) +
  geom_text(
    aes(
      label = count
    ),
    vjust = -0.4,
    color = 'black',
    position = position_dodge(0.9),
    size = 5
  ) +
  theme_bw(
    plot.title = element_text(
      size = 20
    ),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(
      size = 15
    ),
    strip.text.x = element_text(
      size = 15
    ),
    strip.text.y = element_text(
      size = 15
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  ylim(0, 1) +
  labs(
    title = plot_title
  ) +
  scale_fill_brewer(
    palette = 'GnBu'
  ) +
  facet_grid(
    rows = vars(id_f),
    cols = vars(type)
  )

jpeg(output_filename, width = 4000, height = 4000, res = 300)
size_compare
dev.off()
