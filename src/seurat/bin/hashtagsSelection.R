# Libraries loading
#suppressPackageStartupMessages(library(reticulate))
#use_python("/usr/bin/python3")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(patchwork))

# opt/arg parser definition
option_list = list(
  make_option("--inputSeuratRds",
              default=NA,
              type='character',
              help="Path to a RDS file containing a Seurat object"),
  make_option("--output",
              default=NA,
              type='character',
              help="Output name"),
  make_option("--hashtags",
              default = "1,2",
              type='character',
              help="Hashtag numbers, for Hashtags1, Hashtags4 and Hashtags7, --hashtags 1,4,7")
)

opt = parse_args(OptionParser(option_list=option_list))

# Input loading
seuratObj <- readRDS(opt$inputSeuratRds)

# Module code + report code

PresentHashingAntibodies <- paste0("Hashtag",strsplit(opt$hashtags,",")[[1]])

#COUNT DATA


### INSPECT THE HTO COUNTS IN THE WHITELIST OF HASHTAGS
HTOCounts <- Matrix::rowSums(seuratObj@assays$HTO@counts)
FalseDetectedAbs <- setdiff(names(HTOCounts), PresentHashingAntibodies)
TrueDetectedAbs <- intersect(PresentHashingAntibodies, names(HTOCounts))
seuratObj@tools$diagnostics[['UsedHashes']] <- TrueDetectedAbs

### FILTERED SEURAT OBJECT: only Hashtag Antibodies added to the experiment > Replace HTO assay
rawDataHTO <- as.matrix(seuratObj@assays$HTO@counts)
HashCounts <- rawDataHTO[c(TrueDetectedAbs),]
seuratObj[['HTO']] <- CreateAssayObject(counts = HashCounts,min.cells=1)

seuratObj@misc$subsamples <- paste0(seuratObj@project.name,".#",strsplit(opt$hashtags,",")[[1]])
seuratObj@misc$PresentHashingAntibodies <- PresentHashingAntibodies

qcfilename <- paste0(seuratObj@project.name,"_logQC.txt")
write("#QC1: ~~~~HASHTAGS SELECTION~~~~ < Intermediate QC",file =qcfilename)
write(paste("Hashtag names",paste(names(HTOCounts),collapse = " ")),file = qcfilename,append=T)
write(paste("HTO count",paste(HTOCounts,collapse = " ")),file = qcfilename,append = T)
write(paste("True Hashtags",paste(TrueDetectedAbs,collapse = " ")),file = qcfilename, append = T)
write(paste("False Hashtags",paste(FalseDetectedAbs,collapse = " ")),file = qcfilename, append = T)
write("#QC1: ~~~~HASHTAGS SELECTION~~~~ > Intermediate QC",file = qcfilename, append = T)

# Outputs writing into rds file
if(is.na(opt$output)){
  filename = paste0("SEURAT__HASHTAGS_SELECTION.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}
