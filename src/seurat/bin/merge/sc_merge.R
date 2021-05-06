#!/usr/bin/env Rscript

library("optparse")
suppressMessages(library("Seurat", quietly = TRUE))

option_list <- list(
    make_option(
        "--output",
        type = "character",
        dest = "output",
        help = "Output filename."
    )
)

all.args <- parse_args(OptionParser(option_list = option_list), positional_arguments = T)
args <- all.args$options
positional.args <- all.args$args
print(args)

objects <- c()
for (file.path in positional.args) {
    object <- tryCatch({
        readRDS(file = file.path)
    }, error = function(e) {
        stop(paste0("VSN ERROR: Cannot read the given Rds file: ", file.path))
    })

    objects <- c(objects, object)
}

if (length(objects) <= 1) {
    stop(paste0("VSN ERROR: cannot merge ", length(objects), " objects. Needs at least 2 objects!"))
}

merged <- merge(objects[[1]], y = objects[-1])

saveRDS(
    object = merged,
    file = args$output
)
