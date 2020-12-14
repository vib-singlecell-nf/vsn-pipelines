#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--inputMatrix"),
              default=NA,
              type='character',
              help="Path to a .h5 count matrix or a folder containing the files : matrix.mtx(.gz), features.tsv(.gz or genes.tsv) and barcodes.tsv(.gz)"),
  make_option(c("-o", "--output"),
              default=NA,
              type='character',
              help="Output name"),
  make_option("--sample",
              default="sample1",
              type='character',
              help="Sample  name"),
  make_option(c("-f", "--minFeaturesGEX"),
              default=100,
              type='numeric',
              help="minimum number of features for a cell to be kept for the RNA assay"),
  make_option(c("-c", "--minCellsGEX"),
              default=2,
              type='numeric',
              help="minimum number of cells for a feature to be kept for the RNA assay")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)
# Functions
readSTARsolo <- function(sample) {
  if(is.null(sample.name)){
    sample.name <- sample
  }
  dir.data <- normalizePath(sample, mustWork = TRUE)
  m <- methods::as(Matrix::readMM(paste0(dir.data, "/","matrix.mtx")), "dgTMatrix")
  gene_info <- data.table::fread(paste0(dir.data, "/", "features.tsv"), header = FALSE,col.names=c("ID","Symbol"),stringsAsFactors = F)
  cell_info <- data.table::fread(paste0(dir.data, "/", "barcodes.tsv"), header = FALSE,stringsAsFactors = F,col.names = c("Barcode"))
  rownames(m) <- gene_info$Symbol
  colnames(m) <- cell_info$Barcode
  return(m)
}
getGoodName <- function(assay.name,features.names){
  if(assay.name == "Gene Expression"){
    new.name <- "RNA"
  } else if(assay.name == "Custom"){
    new.name <- "HTO"
  } else if(assay.name == "Antibody Capture"){
    if(length(grep("^Hashtag",features.names)) == length(features.names)){
      new.name <- "HTO"
    } else {
      new.name <- "ADT"
    }
  } else {
    print("The module allows only the assays of type : Gene Expression, Antibody Capture and Custom")
  }
  return(new.name)
}

createMultimodalSeurat <- function(assay.list,nbcellsRNA,nbfeaturesRNA, projName){

  if( is.null( names(assay.list) ) ) {
    # Only RNA in h5 file (so assay.list is no list, but a dgCMatrix)
    temp.seurat <- CreateSeuratObject(counts = assay.list,
                                      assay = 'RNA',
                                      min.features = nbfeaturesRNA,
                                      min.cells = nbcellsRNA,
                                      project = projName)
  }
  else {
    if("Gene Expression" %in% names(assay.list)){
      RNAindex <- which(names(assay.list) == "Gene Expression")
      assay.new.name <- getGoodName(names(assay.list)[RNAindex],
                                    assay.list[[RNAindex]]@Dimnames[[1]])
      temp.seurat <- CreateSeuratObject(counts = assay.list[[RNAindex]],
                                        assay = assay.new.name,
                                        min.features = nbfeaturesRNA,
                                        min.cells = nbcellsRNA,
                                        project = projName)

      assay.list <- assay.list[names(assay.list) != "Gene Expression"]
    }
    for(i in 1:length(assay.list)){
      assay.new.name <- getGoodName(names(assay.list)[i],
                                    assay.list[[i]]@Dimnames[[1]])
      if(exists("temp.seurat")){
        temp.seurat[[assay.new.name]] <- CreateAssayObject(counts = assay.list[[i]][,colnames(temp.seurat)])
      }else{
        temp.seurat <- CreateSeuratObject(counts = assay.list[[i]],
                                          assay = assay.new.name, 
                                          project = projName)
      }
    }
  }
  return(temp.seurat)
}

# Pipeline

if(grepl("\\.h5",as.character(opt$inputMatrix))){
  counts.matrix <- Read10X_h5(as.character(opt$inputMatrix))
  if(class(counts.matrix) == "list"){
    row.names(counts.matrix$`Gene Expression`) <- gsub("_{1,}","-",row.names(counts.matrix$`Gene Expression`))
    seuratobj.counts <- createMultimodalSeurat(assay.list= counts.matrix,
                                               nbfeaturesRNA= opt$minFeaturesGEX, 
                                               nbcellsRNA = opt$minCellsGEX,
                                               projName = opt$sample)
  }else{
    row.names(counts.matrix) <- gsub("_{1,}","-",row.names(counts.matrix))
    seuratobj.counts <- CreateSeuratObject(counts = counts.matrix,
                                           min.features= opt$minFeaturesGEX, 
                                           min.cells = opt$minCellsGEX,
                                           project = opt$sample)
  }
  
} else{
  file.names <- list.files(as.character(opt$inputMatrix))
  if("features.tsv.gz" %in% file.names){
    counts.matrix <- Read10X(as.character(opt$inputMatrix))
    if(class(counts.matrix) == "list"){
      row.names(counts.matrix$`Gene Expression`) <- gsub("_{1,}","-",row.names(counts.matrix$`Gene Expression`))
      seuratobj.counts <- createMultimodalSeurat(assay.list= counts.matrix,
                                                 nbfeaturesRNA= opt$minFeaturesGEX, 
                                                 nbcellsRNA = opt$minCellsGEX,
                                                 projName = opt$sample)
    }else{
      row.names(counts.matrix) <- gsub("_{1,}","-",row.names(counts.matrix))
      seuratobj.counts <- CreateSeuratObject(counts = counts.matrix,
                                             min.features= opt$minFeaturesGEX, 
                                             min.cells = opt$minCellsGEX,
                                             project = opt$sample)
    }
  } else if ("genes.tsv" %in% file.names){
    counts.matrix <- Read10X(as.character(opt$inputMatrix))
    row.names(counts.matrix) <- gsub("_{1,}","-",row.names(counts.matrix))
    seuratobj.counts <- CreateSeuratObject(counts = counts.matrix,
                                           min.features = opt$minFeaturesGEX,
                                           min.cells = opt$minCellsGEX,
                                           project = opt$sample)
  } else if ("features.tsv" %in% file.names){
    counts.matrix <- readSTARsolo(as.character(opt$inputMatrix))
    row.names(counts.matrix) <- gsub("_{1,}","-",row.names(counts.matrix))
    seuratobj.counts <- CreateSeuratObject(counts = counts.matrix,
                                           min.features = opt$minFeaturesGEX, 
                                           min.cells = opt$minCellsGEX,
                                           project = opt$sample)
  } else {
    stop("unrecognized input format")
  }
}

#diagnostics
diagnostics <- list()
for(i in names(seuratobj.counts@assays)){
  if(i == "RNA"){
    diagnostics[['dimRawDataRNA']]<-paste0(nrow(seuratobj.counts@assays$RNA)," genes - ",ncol(seuratobj.counts@assays$RNA)," cells")
    diagnostics[['nrGenes']]<-nrow(seuratobj.counts@assays$RNA)
    diagnostics[['nrCells']]<-ncol(seuratobj.counts@assays$RNA)
  }else if (i == "ADT"){
    diagnostics[['dimRawDataADT']]<-paste0(nrow(seuratobj.counts@assays$ADT)," antibodies - ",ncol(seuratobj.counts@assays$ADT)," cells")
  }else if(i == "HTO"){
    diagnostics[['dimRawDataHTO']]<-paste0(nrow(seuratobj.counts@assays$HTO)," hashes - ",ncol(seuratobj.counts@assays$HTO)," cells")
    diagnostics[['DetectedHashes']]<-nrow(seuratobj.counts@assays$HTO)
  }
}

diagnostics[['UsedCells']] <- ncol(seuratobj.counts)
seuratobj.counts@tools$diagnostics <- diagnostics

if(is.na(opt$output)){
  filename = paste0("SEURAT__SEURAT_OBJECT_BUILDER.rds")
  saveRDS(seuratobj.counts,file=filename, compress = T)
} else {
  saveRDS(seuratobj.counts,file=as.character(opt$output), compress = T)
}
