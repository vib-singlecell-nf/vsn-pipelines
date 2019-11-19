#!/usr/bin/env Rscript
args <- base::commandArgs(trailingOnly=TRUE)

base::suppressMessages(expr = { 
  library(DropletUtils)
  library(dplyr)
  library(ggplot2)
})

.reorder <- function(vals, lens, o) {
  out <- base::rep(vals, lens)
  out[o] <- out
  return(out)
}

# Source barcodeRanks (DropletUtils R package)
# Extend the function to take cell counts instead of SingleCellExperiment
BarcodeRanks <- function (cell.counts, lower = 100, fit.bounds = NULL, df = 20, ...) 
{
  o <- base::order(cell.counts$V1, decreasing = TRUE)
  stuff <- base::rle(cell.counts[[1]])
  run.rank <- base::cumsum(stuff$lengths) - (stuff$lengths - 1)/2
  run.totals <- stuff$values
  keep <- run.totals > lower
  if (sum(keep) < 3) {
    stop("insufficient unique points for computing knee/inflection points")
  }
  y <- base::log10(run.totals[keep])
  x <- base::log10(run.rank[keep])
  d1n <- base::diff(y)/base::diff(x)
  right.edge <- base::which.min(d1n)
  left.edge <- base::which.max(d1n[seq_len(right.edge)])
  if (is.null(fit.bounds)) {
    new.keep <- left.edge:right.edge
  }
  else {
    new.keep <- y > base::log10(fit.bounds[1]) & y < base::log10(fit.bounds[2])
  }
  fit <- stats::smooth.spline(x[new.keep], y[new.keep], df = df, 
                       ...)
  d1 <- stats::predict(fit, deriv = 1)$y
  d2 <- stats::predict(fit, deriv = 2)$y
  curvature <- d2/(1 + d1^2)^1.5
  knee <- 10^(y[which.min(curvature)])
  inflection <- 10^(y[right.edge])
  fitted.vals <- base::rep(NA_real_, length(keep))
  fitted.vals[keep][new.keep] <- 10^fitted(fit)
  out <- S4Vectors::DataFrame(rank = .reorder(run.rank, stuff$lengths, 
                                   o), total = .reorder(run.totals, stuff$lengths, o), 
                   fitted = .reorder(fitted.vals, stuff$lengths, o))
  rownames(out) <- cell.counts[[2]]
  metadata(out) <- list(knee = knee, inflection = inflection)
  out
}

# Calculate the Knee and Inflection point
print("Calculate the Knee and Inflection point")
cell.counts <- utils::read.table(file = args[1], header = F, sep = "\t", quote = '', stringsAsFactors = F)
br.out <- BarcodeRanks(cell.counts = cell.counts)
# I'll convert the list into a data frame
# and add the cell barcodes to the data frame
br.out.df <- as.data.frame(br.out)
br.out.df$barcode <- cell.counts$V2

x_knee <- br.out.df %>% 
  dplyr::filter(total > br.out@metadata$knee) %>% 
  dplyr::arrange(total) %>% 
  dplyr::select(rank) %>% 
  utils::head(1)
x_inflection <- br.out.df %>% 
  dplyr::filter(total > br.out@metadata$inflection) %>% 
  dplyr::arrange(total) %>% 
  dplyr::select(rank) %>% 
  utils::head(1)
padding <- length(br.out$rank) / 10

# Plot the Barcode Rank vs Total UMI
print("Plot the Barcode Rank vs Total UMI")
plot <- ggplot2::ggplot(br.out.df, aes(x = rank, y = total)) +
  ggplot2::geom_point() +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text = element_text(size = 10),
                 axis.title = element_text(size = 14),
                 title = element_text(size = 16)) +
  ggplot2::geom_hline(yintercept = br.out@metadata$knee, linetype = 2, colour = "dodgerblue") +
  ggplot2::geom_hline(yintercept = br.out@metadata$inflection, linetype = 2, colour = "forestgreen") +
  labs(x = "Rank", y = "Total", title = "Barcode Rank vs Total UMI") +
  ggplot2::annotate("text", label = paste0("Knee (", x_knee, ")"), x = x_knee$rank + padding, y = br.out@metadata$knee, size = 5) +
  ggplot2::annotate("text", label = paste0("Inflection (", x_inflection, ")"), x = x_inflection$rank + padding, y = br.out@metadata$inflection, size = 5)
ggplot2::ggsave(filename = paste0(args[2], ".barcode_rank_vs_total_umi_plot.png"), plot = plot, device = "png", width = 10, height = 10, type = "cairo")

# Saving selected cell barcodes
print("Save selected cell barcodes")
utils::write.table(x = cell.counts[1:x_knee$rank, "V2"], file = paste0(args[2],".selected_cell_barcodes_by_knee.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
utils::write.table(x = cell.counts[1:x_inflection$rank, "V2"], file = paste0(args[2],".selected_cell_barcodes_by_inflection.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
