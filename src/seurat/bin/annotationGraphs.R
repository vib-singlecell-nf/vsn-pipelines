suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(MEM))
suppressPackageStartupMessages(library(future))

option_list = list(
  make_option(
    "--seuratObj",
    default=NA,
    type='character',
    help="Path to the rds file containing a Seurat object."
  ),
  make_option(
    "--markersFile",
    default=NA,
    type = "character",
    help= "Path to the PBMCmarkers.xlsx file"
  ),
  make_option(
    "--assay",
    default=NA,
    type = "character",
    help= "assay name"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

#################################################
### Cluster Annotation: provided marker genes
#################################################


seuratObj_hashing <- readRDS(opt$seuratObj)
dir.create(paste0("Plots/",opt$assay), recursive = T)

clusters.name <- paste0(opt$assay,"_clusters")

if(!is.na(opt$markersFile)){
  markers <- read.xlsx(opt$markersFile)
  colnames(markers) <- c('celltype', 'gene', 'marker')
  genes.cleaned <- c()
  for( gene in markers$gene ) {
    gene.full <- grep(paste0(gene, '$'), rownames(seuratObj_hashing), value = T)
    if( length(gene.full) != 1 ) {
      
      gene.full <- NA
    }
    genes.cleaned <- c(genes.cleaned, gene.full)
  }
  
  # Feature, Ridge and Violin plots
  # We will loop over the cell-types and create the plots in separate files
  markers$gene <- genes.cleaned
  det.markers <- markers %>% drop_na() %>% arrange(celltype)
  
  cell.ids <- unique(det.markers$celltype)
  for(i in cell.ids){
    pdf(file=paste0("Plots/",opt$assay,"/17_markers_",gsub(" ", "",i),".pdf")
        , width = 17*0.45, height = 12.4*0.45)
    
    print(RidgePlot(seuratObj_hashing, features = det.markers$gene[det.markers$celltype==i], ncol = 2))
    print(VlnPlot(seuratObj_hashing, features = det.markers$gene[det.markers$celltype==i], ncol = 2))
    print(FeaturePlot(seuratObj_hashing, features = det.markers$gene[det.markers$celltype==i], ncol = 2))
    
    dev.off()
  }
  
  # Seurat dotmap to help annotation
  marker.dot <- DotPlot(seuratObj_hashing, features = unique(det.markers$gene), group.by=clusters.name) +
    RotatedAxis() + theme(legend.position = "none")
  
  # Seurat heatmap to help annotation
  marker.heat <- DoHeatmap(subset(seuratObj_hashing,downsample = 100), group.by = clusters.name, features = unique(det.markers$gene), size = 3)
  
  pdf(file=paste0("Plots/",opt$assay,"/17_markermaps.pdf"), width = 17*0.45, height = 12.4*0.45)
  print(marker.dot)
  print(marker.heat)
  dev.off()
  
  # can we create a MEM heatmap that allows to see what is going on?
  # should work best on markers that are orthogonal (like flow panel)
  # could also be used with the top-marker genes
  # here with predefined marker genes
  
  MEM.input <- cbind(t(seuratObj_hashing@assays[[opt$assay]]@scale.data[unique(det.markers$gene),]),
                     cluster = seuratObj_hashing@meta.data[clusters.name][,1])
  
  # creates a series of tables ready for heatmap visualization.
  marker.MEM <- MEM(MEM.input)
  # the MEM is written out in a very inconvenient spot, can we do it better manually?
  # you can clearly see if it is informative or not
  
  marker.MEM.heatmap <- as.data.frame(marker.MEM$MEM_matrix) %>%
    tibble::rownames_to_column(var="cluster") %>%
    gather(key = "marker", value = "MEM.score",-one_of("cluster")) %>%
    mutate(cluster=as.factor(as.numeric(cluster)-1)) %>%
    ggplot(aes(cluster,marker,cluster)) +
    geom_tile(aes(fill=MEM.score),colour = "grey50") +
    theme_classic() +
    scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red" )
  
  
  pdf(file=paste0("Plots/",opt$assay,"/18_markerMEM.pdf"))
  print(marker.MEM.heatmap)
  dev.off()
}


#################################################
### Cluster Annotation: detection of marker genes
#################################################

Idents(seuratObj_hashing) <- clusters.name
exp.markers <- FindAllMarkers(seuratObj_hashing, assay = opt$assay)
top.exp.markers <- group_by(exp.markers,cluster) %>% top_n(n = 5, wt = avg_logFC)

# again we should try to build the heatmaps ourselves
exp.MEM.input <- cbind(t(seuratObj_hashing@assays[[opt$assay]]@scale.data[unique(top.exp.markers$gene),]),
                       cluster = seuratObj_hashing@meta.data[clusters.name][,1])
exp.MEM <- MEM(exp.MEM.input)

exp.MEM.heatmap <- as.data.frame(exp.MEM$MEM_matrix) %>%
  tibble::rownames_to_column(var="cluster") %>%
  gather(key = "marker", value = "MEM.score",-one_of("cluster")) %>%
  mutate(cluster=as.factor(as.numeric(cluster)-1)) %>%
  ggplot(aes(cluster,marker,cluster)) +
  geom_tile(aes(fill=MEM.score),colour = "grey50") +
  theme_classic() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red" )


pdf(file = paste0("Plots/",opt$assay,"/19_expMEM.pdf"))
print(exp.MEM.heatmap)
dev.off()

### Create list with markers
totalNrClusters <- seuratObj_hashing@meta.data[clusters.name][,1] %>% as.factor() %>% levels() %>% length()
markersList <- list()

for(i in 1:totalNrClusters){
  markersList[[i]] <- exp.markers[exp.markers$cluster == i - 1,]
}
names(markersList) <- paste0("Cluster ",0:(totalNrClusters-1))

### Write to Excel

write.xlsx(markersList, file = paste0("allClusters_",seuratObj_hashing@project.name,".xlsx"))