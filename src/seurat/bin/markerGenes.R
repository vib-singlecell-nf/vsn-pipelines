suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(limma))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="File path to the rds file containing a Seurat object."
  ),
  make_option(
    "--output",
    default=NA,
    type='character',
    help="output file name"
  ),
  make_option(
    "--markerGenes",
    default=NA,
    type='character',
    help="name of the marker genes"
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- list()
#opt$seuratObj <- "/home/quentinr/Documents/VIB_work/scRNA_seq/datasets/SAR13/work/3a/85fe7439ef422c4dea84f32451a41e/SAR13.SEURAT__FIND_VARIABLE_FEATURES_SCT.rds"
#opt$markerGenes <- "CCR10,CCR7,CD14,CD160,CD19,CD27,CD28,CD34,CD3D,CD4,CD79A,CD8A,CD8B,CST3,EBF1,ENTPD1,FCER1A,FCGR3A,FGR,FOXP3,GLIPR2,GNLY,ICOS,ID3,IFI30,IGJ,IL2RA,IL3RA,IL7R,ITGAX,KLRC1,KLRD1,LY86,LYZ,MS4A1,MS4A7,MZB1,NCAM1,NKG7,NRP1,PF4,PPBP,S100A4,SELL"

seuratObj <- readRDS(file = opt$seuratObj)
dir.create("Plots/RNA", recursive = T)
# the high variable features have been defined and a series of venndiagrams could be used:
print("# how many cell-cycle genes?")

genome_g2m.genes <- paste0("(^|-)",cc.genes.updated.2019$g2m.genes,"$")
genome_s.genes <- paste0("(^|-)",cc.genes.updated.2019$s.genes,"$")
HVG.g2m.genes <- unlist(lapply(genome_g2m.genes,function(x) grep(x,seuratObj@misc$HVG,value=T)))
HVG.s.genes <- unlist(lapply(genome_s.genes,function(x) grep(x,seuratObj@misc$HVG,value=T)))
perc.HVG.s <- round(100*length(HVG.s.genes)/length(cc.genes.updated.2019$s.genes))
perc.HVG.g2m <- round(100*length(HVG.g2m.genes)/length(cc.genes.updated.2019$g2m.genes))
seuratObj@tools$diagnostics[['percHVGs']]<-paste0(perc.HVG.s," % of S phase genes are detected as HVG")
seuratObj@tools$diagnostics[['percHVGg2m']]<-paste0(perc.HVG.g2m," % of G2M phase genes are detected as HVG")

if("rbc.genes" %in% names(seuratObj@misc)){
  print("# how many RBC genes?")
  genomeg_rbc.genes <- paste0("(^|-)",seuratObj@misc$rbc.genes,"$")
  HVG.RBC.genes <- unlist(lapply(genomeg_rbc.genes,function(x) grep(x,seuratObj@misc$HVG,value=T)))
  perc.HVG.rbc <- round(100*length(HVG.RBC.genes)/length(seuratObj@misc$rbc.genes))
  seuratObj@tools$diagnostics[['percHVGrbc']]<-paste0(perc.HVG.rbc," % of RBC genes are detected as HVG")

  nonHVG.RBC.genes.index <- unlist(lapply(1:length(genomeg_rbc.genes),
                                             function(x) ifelse(sum(grep(genomeg_rbc.genes[x],seuratObj@misc$HVG))==0,
                                                                x,
                                                                NA)))
  nonHVG.RBC.genes.index <- nonHVG.RBC.genes.index[!is.na(nonHVG.RBC.genes.index)]
  nonHVG.RBC.genes <- seuratObj@misc$rbc.genes[nonHVG.RBC.genes.index]
}
if("covid.genes" %in% names(seuratObj@misc)){
  print("# how many COVID genes?")
  genomeg_covid.genes <- paste0("(^|-)",seuratObj@misc$covid.genes,"$")
  HVG.COVID.genes <- unlist(lapply(genomeg_covid.genes,function(x) grep(x,seuratObj@misc$HVG,value=T)))
  perc.HVG.covid <- round(100*length(HVG.COVID.genes)/length(seuratObj@misc$covid.genes))
  seuratObj@tools$diagnostics[['percHVGcovid']]<-paste0(perc.HVG.covid," % of covid genes are detected as HVG")

  nonHVG.COVID.genes.index <- unlist(lapply(1:length(genomeg_covid.genes),
                                          function(x) ifelse(sum(grep(genomeg_covid.genes[x],seuratObj@misc$HVG))==0,
                                                             x,
                                                             NA)))
  nonHVG.COVID.genes.index <- nonHVG.COVID.genes.index[!is.na(nonHVG.COVID.genes.index)]
  nonHVG.COVID.genes <- seuratObj@misc$covid.genes[nonHVG.COVID.genes.index]
}

###################################
########## Building marker genes
###################################

# 4. how many marker genes or genes of interest?

seuratObj@misc$marker.genes <- strsplit(opt$markerGenes,",")[[1]]
genome_marker.genes <- paste0("(^|-)",seuratObj@misc$marker.genes,"$")
HVG.marker.genes <- unlist(lapply(genome_marker.genes,function(x) grep(x,seuratObj@misc$HVG,value = T)))

nonHVG.marker.genes.index <- unlist(lapply(1:length(genome_marker.genes),
                                     function(x) ifelse(sum(grep(genome_marker.genes[x],seuratObj@misc$HVG))==0,
                                                        x,
                                                        NA)))
nonHVG.marker.genes.index <- nonHVG.marker.genes.index[!is.na(nonHVG.marker.genes.index)]
nonHVG.marker.genes <- seuratObj@misc$marker.genes[nonHVG.marker.genes.index]

# marker genes that are not picked up in the HVG and thus are not used for clustering
perc.HVG.marker <- round(100*length(HVG.marker.genes)/length(seuratObj@misc$marker.genes))
seuratObj@tools$diagnostics[['percHVGmarker']] <- paste0(perc.HVG.marker," % of marker genes are detected as HVG")
seuratObj@tools$diagnostics[['missedHVGmarkers']] <- paste0(nonHVG.marker.genes,collapse = " ")

if("rbc.genes" %in% names(seuratObj@misc) && "covid.genes" %in% names(seuratObj@misc)){
  DF <- list(genes = c(seuratObj@misc$HVG,
                              nonHVG.RBC.genes,
                              nonHVG.COVID.genes))

  for(i in 1:length(DF$genes)){
    if(DF$genes[i] %in% seuratObj@misc$HVG){
      DF$HVG[i]<-TRUE
    } else{
      DF$HVG[i]<-FALSE
    }

    if(DF$genes[i] %in% HVG.RBC.genes){
      DF$RBC[i]<-TRUE
    } else{
      if(DF$genes[i] %in% nonHVG.RBC.genes){
        DF$RBC[i]<-TRUE
      }else {
        DF$RBC[i]<-FALSE
      }
    }

    if(DF$genes[i] %in% HVG.COVID.genes){
      DF$COVID[i]<-TRUE
    } else{
      if(DF$genes[i] %in% nonHVG.COVID.genes){
        DF$COVID[i]<-TRUE
      }else {
        DF$COVID[i]<-FALSE
      }
    }
  }

  DF <- as.data.frame(DF)
  rownames(DF)<-DF$genes
  DF<-DF[,-1]
  a <- vennCounts(DF)

  png(file="Plots/RNA/06_VennDiagram_HVG-RBC-COVID.png", width = 750, height = 750, type = 'cairo')
  vennDiagram(a, names = c("Highly variant genes","RBC genes","COVID genes"),
              circle.col = c("slateblue","darkorange","darkturquoise"))
  dev.off()

  pdf(file="Plots/RNA/06_VennDiagram_HVG-RBC-COVID.pdf", width = 10, height = 10)
  vennDiagram(a, names = c("Highly variant genes","RBC genes","COVID genes"),
              circle.col = c("slateblue","darkorange","darkturquoise"))
  dev.off()


  DFm <- list(genes = unique(c(seuratObj@misc$HVG,
                       nonHVG.RBC.genes,
                       nonHVG.COVID.genes,
                       nonHVG.marker.genes)))

  for(i in 1:length(DFm$genes)){
    if(DFm$genes[i] %in% seuratObj@misc$HVG){
      DFm$HVG[i]<-TRUE
    } else{
      DFm$HVG[i]<-FALSE
    }

    if(DFm$genes[i] %in% HVG.RBC.genes){
      DFm$RBC[i]<-TRUE
    } else{
      if(DFm$genes[i] %in% nonHVG.RBC.genes){
        DFm$RBC[i]<-TRUE
      }else {
        DFm$RBC[i]<-FALSE
      }
    }

    if(DFm$genes[i] %in% HVG.COVID.genes){
      DFm$COVID[i]<-TRUE
    } else{
      if(DFm$genes[i] %in% nonHVG.COVID.genes){
        DFm$COVID[i]<-TRUE
      }else {
        DFm$COVID[i]<-FALSE
      }
    }

    if(DFm$genes[i] %in% HVG.marker.genes){
      DFm$MARKER[i]<-TRUE
    } else{
      if(DFm$genes[i] %in% nonHVG.marker.genes){
        DFm$MARKER[i]<-TRUE
      }else {
        DFm$MARKER[i]<-FALSE
      }
    }
  }

  DFm<-as.data.frame(DFm)
  rownames(DFm)<-DFm$genes
  DFm<-DFm[,-1]
  a <- vennCounts(DFm)
  png(file="Plots/RNA/06_VennDiagram_HVG-RBC-COVID-Markers.png", width = 1000, height = 750, type = 'cairo')
  vennDiagram(a, names = c("Highly variant genes","RBC genes","COVID genes","Marker genes"),
              circle.col = c("slateblue","darkorange","darkturquoise","violet"))
  dev.off()

  pdf(file="Plots/RNA/06_VennDiagram_HVG-RBC-COVID-Markers.pdf", width = 12, height = 10)
  vennDiagram(a, names = c("Highly variant genes","RBC genes","COVID genes","Marker genes"),
              circle.col = c("slateblue","darkorange","darkturquoise","violet"))
  dev.off()
}else{

  DFn <- list(genes = unique(c(seuratObj@misc$HVG,
                               nonHVG.marker.genes)))

  for(i in 1:length(DFn$genes)){
    if(DFn$genes[i] %in% seuratObj@misc$HVG){
      DFn$HVG[i]<-TRUE
    } else{
      DFn$HVG[i]<-FALSE
    }

    if(DFn$genes[i] %in% HVG.marker.genes){
      DFn$MARKER[i]<-TRUE
    } else{
      if(DFn$genes[i] %in% nonHVG.marker.genes){
        DFn$MARKER[i]<-TRUE
      }else {
        DFn$MARKER[i]<-FALSE
      }
    }
  }

  DFn<-as.data.frame(DFn)
  rownames(DFn)<-DFn$genes
  DFn <-DFn[,-1]
  a <- vennCounts(DFn)
  pdf(file="Plots/RNA/06_VennDiagram_HVG-Markers.pdf", width = 12, height = 10)
  vennDiagram(a, names = c("Highly variant genes","Marker genes"),
              circle.col = c("slateblue","darkorange"))
  dev.off()
}


if(is.na(opt$output)){
  filename = paste0(seuratObj@project.name,"SEURAT___MARKER_GENES.rds")
  saveRDS(seuratObj,file=filename, compress = T)
} else {
  saveRDS(seuratObj,file=as.character(opt$output), compress = T)
}
