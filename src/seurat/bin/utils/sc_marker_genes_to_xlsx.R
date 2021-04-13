#!/usr/bin/env Rscript

library("optparse")
library("openxlsx")

option_list <- list(
    make_option(
        "--input",
        type = "character",
        dest = "input",
        help = "A Rds file containing results of Seurat::FindAllMarkers()"
    ),
    make_option(
        "--output",
        type = "character",
        dest = "output",
        help = "Output filename."
    )
)

args <- parse_args(OptionParser(option_list = option_list))
print(args)

markers <- tryCatch({
    readRDS(file = args$input)
}, error = function(e) {
    stop("VSN ERROR: Cannot read the given Rds file.")
})

markers.list <- list()

for (cluster in levels(markers$cluster)) {
    # Transform numeric cluster names to 'Cluster <number>'
    cluster.name <- cluster
    if (!grepl("[a-zA-Z]", cluster.name)) {
        cluster.name <- paste0("Cluster ", cluster.name)
    }
    markers.list[[cluster.name]] <- markers[markers$cluster == cluster,]
}

write.xlsx(markers.list, file = args$output)
