suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing the seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  ),
  make_option(
    "--assay",
    default="RNA",
    type='character',
    help="Assay from the seurat object to put in the AnnData (default : RNA)"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

seuratObj <- readRDS(file = opt$seuratObj)

if(!is.na(opt$output)){
  h5Seurat_name <- gsub("\\.h5ad$",".h5Seurat",opt$output)
} else{
  h5Seurat_name <- paste0(seuratObj@project.name,".h5Seurat")
  opt$output <- paste0(seuratObj@project.name,"_",opt$assay,".h5ad")
}

if(length(seuratObj@misc)!=0 || length(seuratObj@tools)!=0){
  sink(paste0(seuratObj@project.name,"_diagnostics.txt"))
  print(seuratObj@tools)
  print(seuratObj@misc)
  sink()
  seuratObj@tools <- list()
  seuratObj@misc <- list()
}

SeuratDisk::SaveH5Seurat(seuratObj,h5Seurat_name)
  
if(opt$assay == "ALL"){
  for(i in names(seuratObj@assays)){
    h5ad_name <- gsub("ALL",i,opt$output)
    Convert(h5Seurat_name, dest = h5ad_name, assay=i)
  }
} else{
  Convert(h5Seurat_name, dest = opt$output, assay=opt$assay)
}
