suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DropletUtils"))
suppressPackageStartupMessages(library("ggplot2"))

option_list = list(
  make_option(
    c("-r", "--rdsFile"),
    default=NA,
    type='character',
    help="File path to the Rds file containing a SingleCellExperiment object."
  ),
  make_option(
    c("-o", "--output"),
    dest="output",
    type='character',
    help="Output file path to save the SingleCellExperiment object as an Rds file.")
)

opt <- parse_args(OptionParser(option_list=option_list))
sce <- readRDS(opt$rdsFile)
bcrank <- barcodeRanks(counts(sce))
# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
bcrank.plot <- as.data.frame(bcrank[uniq,]@listData)
# Store information
sce@metadata$bcRankMetrics$knee <- bcrank@metadata$knee
sce@metadata$bcRankMetrics$inflection <- bcrank@metadata$inflection
sce@metadata$bcRankMetrics$rank <- bcrank.plot$rank
sce@metadata$bcRankMetrics$total <- bcrank.plot$total

#graph
sce@metadata$plots$bcRankMetrics <- ggplot(bcrank.plot,aes(rank,total)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  xlab("barcode rank") +
  ylab("count depth (total UMI count)") +
  ggtitle("Cells") +
  theme_bw() +
  geom_hline(yintercept=bcrank@metadata$inflection, linetype="dashed", 
             color = "darkgreen", size=1) +
  geom_hline(yintercept=bcrank@metadata$knee, linetype="dashed", 
             color = "dodgerblue", size=1)

# Data saving
saveRDS(sce, opt$output, compress=TRUE)