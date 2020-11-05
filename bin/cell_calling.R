#!/usr/bin/env Rscript

print("##################################")
print("# ArchR: Arrow file -> QC plot   #")
print("##################################")

################################################################################

library("optparse")
parser <- OptionParser(
  prog = "run_archr_qc_plot",
  description = "Load an Arrow file, call cells, create QC plot."
)
parser <- add_option(
  parser,
  c("--arrow_file"),
  action = "store",
  default = NULL,
  help = "Input file [default]"
)
parser <- add_option(
  parser,
  c("--output_directory"),
  action = "store",
  default = NULL,
  help = "output directory prefix for plot file"
)
parser <- add_option(
  parser,
  c("--filter_tss"),
  action = "store",
  default = 4,
  help = "Filter on TSS Enrichment [%default]"
)
parser <- add_option(
  parser,
  c("--filter_frags"),
  action = "store",
  default = 1000,
  help = "Filter on fragments/cell [%default]"
)
parser <- add_option(
  parser,
  c("-g", "--genome"),
  action = "store",
  default = "hg38",
  help = "Genome to use, choose from (hg19, hg38, mm9, and mm10) [%default]"
)

parser <- add_option(
  parser,
  c("-t", "--threads"),
  action = "store",
  default = 1,
  help="Number of threads to use. [default %default]"
)
parser <- add_option(
  parser,
  c("-s", "--seed"),
  action = "store",
  default = 617,
  help="Seed. [default %default]"
)

args <- parse_args(parser)

cat("Parameters: \n")
print(args)


# ################################################################################
# args$arrow_file = "CNAG_1.arrow"
# args$genome = "hg38"
# args$filter_tss = 8
# args$filter_frags = 1000 
# ################################################################################

suppressMessages(library(ArchR, quietly=TRUE))

set.seed(args$seed)
addArchRThreads(threads=args$threads)
addArchRGenome(args$genome)



sample_id = gsub(".arrow$", "", basename(args$arrow_file))

nucLength = 147
allcells = ArchR:::.availableCells(args$arrow_file, passQC=FALSE)

fragSummary = ArchR:::.fastFragmentInfo(
    ArrowFile = args$arrow_file, 
    cellNames = allcells, 
    nucLength = nucLength,
    prefix = names(inputFile),
    logFile = NULL
)


ga = getGeneAnnotation()
tssSummary = ArchR:::.fastTSSEnrichment(
    TSS = ga$TSS,
    ArrowFile = args$arrow_file, 
    cellNames = allcells 
    )



# all(allcells == rownames(fragSummary$dfSummary))
# all(allcells == names(tssSummary$tssScores))

dat = data.frame(
                 barcode = sapply(strsplit(allcells,"#"),'[',2),
                 Keep = NA, #ifelse(allcells %in% filtcells,1,0),
                 nFrags = fragSummary$dfSummary$nFrags,
                 TSSEnrichment = tssSummary$tssScores,
                 stringsAsFactors=FALSE
                 )

### filtered cells from arrow file:
# filtcells = ArchR:::.availableCells(args$arrow_file, passQC=TRUE)

# determine based on the filters set here:
dat$Keep = as.integer(
                      dat$TSSEnrichment >= args$filter_tss &
                      dat$nFrags >= args$filter_frags
                      )

ggtitle <- sprintf("%s\n%s\n%s",
        paste0(sample_id, "\nnCells Pass Filter = ", sum(dat$Keep)),
        paste0("Median Fragments = ", median(dat$nFrags[dat$Keep==1])),
        paste0("Median TSS Enrichment = ", median(dat$TSSEnrichment[dat$Keep==1]))
      )


##################################################
library(viridis)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_TSS_nFrags = function(dat, filterTSS, filterFrags, ggtitle) {
    p2 = ggplot(data=dat, aes(x=log10(nFrags), y=TSSEnrichment,
                              color=get_density(log10(nFrags), TSSEnrichment, n=100))
                ) +
        geom_point(size=1) +
        labs(color="Density", title=ggtitle) +
        scale_color_viridis(option="cividis") +
        #scale_x_continuous(limits=c(1,NA), expand=c(0,0.03)) +
        scale_x_continuous(expand=expansion(mult=c(0.02,0.02))) +
        scale_y_continuous(expand=expansion(mult=c(0.02,0.02))) +
        xlab("log 10 (Unique fragments)") +
        ylab("TSS Enrichment") +
        geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
        geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25) +
        theme_minimal() +
        theme(panel.border=element_rect(color="black", fill=NA),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              aspect.ratio=1
              )
    return(p2)
}

p2 = plot_TSS_nFrags(dat, args$filter_tss, args$filter_frags, ggtitle)

pdf(paste0(args$output_directory,"/",sample_id,"-TSSEnrichment_vs_nFrags.pdf"),width=8,height=8)
p2
dev.off()

write.table(dat[dat$Keep==1,],
            paste0(args$output_directory,"/",sample_id,"-qc_stats.txt"), sep='\t', quote=FALSE)


