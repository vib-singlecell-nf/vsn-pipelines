##### Function drawVlnPlot
drawVlnPlot<-function(toPlot, fileName, colsToColor, png.device.type = NULL){
  toPlot<-toPlot[order(toPlot[,colsToColor[1]]),]
  p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()

  toPlot<-toPlot[order(toPlot[,colsToColor[2]]),]
  p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()

  toPlot<-toPlot[order(toPlot[,colsToColor[3]]),]
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()

  toPlot<-toPlot[order(toPlot[,colsToColor[4]]),]
  p_rbc <- ggplot(toPlot, aes(staticNr, percent.rbc)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, alpha=0.5, aes(color="percent.rbc")) +
    scale_color_manual(values=c("#00bfc4")) +
    theme_classic()

  toPlot<-toPlot[order(toPlot[,colsToColor[5]]),]
  p_COVID <- ggplot(toPlot, aes(staticNr, percent.COVID)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, alpha=0.5, aes(color="percent.COVID")) +
    scale_color_manual(values=c("#00bfc4")) +
    theme_classic()

  grid.arrange(p_nGene, p_nUMI,p_mito,p_rbc,p_COVID, ncol=5)
  if(fileName != ""){
    if(!is.null(png.device.type)){
      ggsave(grid.arrange(p_nGene, p_nUMI,p_mito,p_rbc,p_COVID, ncol=5),
             file=fileName,
             dpi=200,
             units = "mm",
             width = 500, height = 100,
             type = png.device.type)
    }else{
      ggsave(grid.arrange(p_nGene, p_nUMI,p_mito,p_rbc,p_COVID, ncol=5),
             file=fileName,
             dpi=200,
             units = "mm",
             width = 500, height = 100)
    }
  }
}

##### Function drawVlnPlot_out
drawVlnPlot_out<-function(toPlot, fileName, colsToColor, png.device.type = NULL){
  listplot <- list()

  toPlot<-toPlot[order(toPlot[,colsToColor[["nGene"]]]),]
  listplot[["nGene"]] <-ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[["nGene"]]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()

  toPlot<-toPlot[order(toPlot[,colsToColor[["nUMI"]]]),]
  listplot[["nUMI"]] <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="gray80") +
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[["nUMI"]]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()

  if("percent.mito" %in% colnames(toPlot)){
    toPlot<-toPlot[order(toPlot[,colsToColor[["percent.mito"]]]),]
    listplot[["percent.mito"]] <- ggplot(toPlot, aes(staticNr, percent.mito)) +
      geom_violin(fill="gray80") +
      geom_jitter(height = 0, width = 0.3, aes_string(col = colsToColor[["percent.mito"]]), alpha=0.5) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic()
  }

  if("percent.rbc" %in% colnames(toPlot)){
    toPlot<-toPlot[order(toPlot[,colsToColor[["percent.rbc"]]]),]
    listplot[["percent.rbc"]] <- ggplot(toPlot, aes(staticNr, percent.rbc)) +
      geom_violin(fill="gray80") +
      geom_jitter(height = 0, width = 0.3, aes_string(col = colsToColor[["percent.rbc"]]), alpha=0.5) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic()
  }

  if("percent.COVID" %in% colnames(toPlot)){
    toPlot<-toPlot[order(toPlot[,colsToColor[["percent.COVID"]]]),]
    listplot[["percent.COVID"]] <- ggplot(toPlot, aes(staticNr, percent.COVID)) +
      geom_violin(fill="gray80") +
      geom_jitter(height = 0, width = 0.3, aes_string(col = colsToColor[["percent.COVID"]]), alpha=0.5) +
      scale_color_manual(values=c("#00bfc4", "#F8766D")) +
      theme_classic()
  }

  if(fileName != ""){
      ggsave(do.call("grid.arrange", c(listplot, ncol=length(listplot))),
             file=fileName,
             dpi=200,
             units = "mm",
             width = 500, height = 100,
             type = png.device.type)
  }
}


##### Function inhouse HTO criteria: within estimation
# write function to estimate negative and within doublet
# just needs the maxval, fits a gaussian, returns label
# parameter that needs to be set is the quantile of the normal distribution
# q <- 0.999 # very stringent
within.crit <- function(x,q)
{
  y <- log2(x+1)
  fit <- fitdist(y,dist = "norm")
  lb <- fit$estimate[[1]]-qnorm(q)*fit$estimate[[2]]
  ub <- fit$estimate[[1]]+qnorm(q)*fit$estimate[[2]]
  z <- ifelse(y<lb,"Negative",
              ifelse(y>ub,"Doublet","Singlet"))
  return(z)
}

##### Function inhouse HTO criteria: between estimation
# B.1: how large should the difference be between the first and second hashtag for clear distinction?
# parameter to set is the percentage or fold change on the log2-scale, start with 20 %
#  p <- 0.8 # stringent enough?
dif.crit <- function(max,secmax,p)
{
  ifelse(log2(secmax+1)>p*log2(max+1),"too.close","ok")
}
# B.2: use a gaussian mixture on the second highest values and start from the negative distribution
# parameter to set is again the quantile...
# q <- 0.999 # again very stringent criterion
mix.crit <- function(secmax,q)
{
  mix <- tryCatch(mixtools::normalmixEM(log2(secmax+1),k=2),error = function(e) {secmax <<- rep(NA,length(secmax))})
  ifelse(is.na(secmax),
         "noNormalMixtures",
         ifelse(log2(secmax+1) > mix$mu[which.min(mix$mu)]+qnorm(q)*mix$sigma[which.min(mix$mu)],
                "Doublet",
                "Singlet")
  )
}

##### isOutlier function from scater (simplified) for Seurat pipeline uses.

is.Outlier <- function (metric, nmads = 3, type = c("both", "lower", "higher"),
                        log = FALSE, batch = NULL)
{
  if (log) {
    metric <- log2(metric)
  }
  N <- length(metric)
  if (nobatch <- is.null(batch)) {
    batch <- rep("1", N)
  }
  else {
    if (length(batch) != N) {
      stop("length of 'batch' must equal length of 'metric'")
    }
    batch <- as.character(batch)
  }
  all_batches <- sort(unique(batch))
  M <- metric
  B <- batch
  if (any(na.drop <- is.na(M))) {
    M <- M[!na.drop]
    B <- B[!na.drop]
    warning("missing values ignored during outlier detection")
  }
  by.batch <- split(M, B)
  empty <- rep(NA_real_, length(all_batches))
  names(empty) <- all_batches
  cur.med <- empty
  cur.med[names(by.batch)] <- unlist(lapply(by.batch, median,0))
  cur.mad <- empty
  cur.mad[names(by.batch)] <- unlist(mapply(mad,
                                            x = by.batch,
                                            center = cur.med[names(by.batch)],
                                            SIMPLIFY = FALSE))
  diff.val <- nmads * cur.mad
  upper.limit <- cur.med + diff.val
  lower.limit <- cur.med - diff.val
  type <- match.arg(type)
  if (type == "lower") {
    upper.limit[] <- Inf
  }
  else if (type == "higher") {
    lower.limit[] <- -Inf
  }
  collected <- (metric < lower.limit[batch] | upper.limit[batch] <
                  metric)
  names(collected) <- names(metric)
  all.threshold <- rbind(lower = lower.limit, higher = upper.limit)
  if (nobatch) {
    all.threshold <- drop(all.threshold)
  }
  if (log) {
    all.threshold <- 2^all.threshold
  }
  attr(collected, "thresholds") <- all.threshold
  collected
}


##### SHOW ALL GENE OUTLIERS = BLUE -- UMI OUTLIERS = PURPLE -- MITO OUTLIERS = ORANGE
drawVlnPlot_color<-function(toPlot,fileName, png.device.type = NULL){
  plotlist <- list()
  plotlist[["p_nGene"]] <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  plotlist[["p_nUMI"]] <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#ca2bd8"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

  titlefirstouts<-paste(sum(toPlot$final.drop),
                              "cells with low library sizes/expressed features are removed")

  if("percent.mito" %in% colnames(toPlot)){
    plotlist[["p_mito"]] <- ggplot(toPlot, aes(staticNr, percent.mito)) +
      geom_violin(fill="#e0e5e4") +
      geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
      scale_color_manual(values=c("#cecece", "#f26f18")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
    titlefirstouts<-paste(sum(toPlot$final.drop),
                                "cells with low library sizes/expressed features or high % mitochondrial genes are removed")
  }

  if(fileName != ""){
      ggsave( do.call("grid.arrange", c(plotlist, ncol=length(plotlist), top=titlefirstouts)),
              file=fileName,
              dpi=200,
              units = "mm",
              width = 500, height = 100,
              type = png.device.type)
  }
}

##### Do outlier cells with low Umi/nrGenes have increased mitochondrial genes?
##### Map the cells with low UMI or nGene (outliers) on the percent.mito charts > Mostly correspond
drawVlnPlot_color_nGene<-function(toPlot,fileName, png.device.type = NULL){
  toPlot<-toPlot[order(toPlot[,"nGene.drop"]),]
  p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

  if(fileName != ""){
    ggsave( grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3, top="Percentage mitochondrial genes of cells with low nGene and UMI count (blue)"),
            file=fileName,
            dpi=200,
            units = "mm",
            width = 500, height = 100,
            type = png.device.type)
  }
}
##### Map the cells with high proportional mito genes on the nGene/nUMI plots.
##### Is the cut-off ('nmads') for nGene reasonable?
drawVlnPlot_color_nGene_mito<-function(toPlot,fileName, png.device.type = NULL){
  p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nGeneMitodrop <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nUMIMitodrop <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  if(fileName != ""){
    ggsave( grid.arrange(p_nGene,p_nUMI, p_nGeneMitodrop, p_nUMIMitodrop,p_mito, ncol=5,top="nGene and UMI count of cells with high % mitochondrial genes (orange)"),
            file=fileName,
            dpi=200,
            units = "mm",
            width = 500, height = 100,
            type = png.device.type)
  }
}
##### Map the cells with high proportional mito genes on the other plots.
drawVlnPlot_color_mito_extra<-function(toPlot,fileName, png.device.type = NULL){
  toPlot<-toPlot[order(toPlot[,"mito.drop"]),]
  p_nGeneMitodrop <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nUMIMitodrop <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_RBCIMitodrop <- ggplot(toPlot, aes(staticNr, percent.rbc)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_COVIDMitodrop <- ggplot(toPlot, aes(staticNr, percent.COVID)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#f26f18"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  if(fileName != ""){
    ggsave( grid.arrange(p_nGeneMitodrop, p_nUMIMitodrop,p_mito,p_RBCIMitodrop,p_COVIDMitodrop, ncol=5, top=" high % mitochondrial genes (orange)"),
            file=fileName,
            dpi=200,
            units = "mm",
            width = 500, height = 100,
            type = png.device.type)
  }
}
##### Map the cells with outlier nGENE on the other plots.
drawVlnPlot_color_nGene_extra<-function(toPlot,fileName, png.device.type = NULL){
  toPlot<-toPlot[order(toPlot[,"nGene.drop"]),]
  p_nGeneMitodrop <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nUMIMitodrop <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_RBCIMitodrop <- ggplot(toPlot, aes(staticNr, percent.rbc)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_COVIDMitodrop <- ggplot(toPlot, aes(staticNr, percent.COVID)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#70bce5"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  grid.arrange(p_nGeneMitodrop, p_nUMIMitodrop,p_mito,p_RBCIMitodrop,p_COVIDMitodrop, ncol=5,
               top=" nGene outliers (blue)")
  if(fileName != ""){
    ggsave( grid.arrange(p_nGeneMitodrop, p_nUMIMitodrop,p_mito,p_RBCIMitodrop,p_COVIDMitodrop, ncol=5,top=" nGene outliers (blue)"),
            file=fileName,
            dpi=200,
            units = "mm",
            width = 500, height = 100,
            type = png.device.type)
  }
}
##### Map the cells with outlier nUMI on the other plots.
drawVlnPlot_color_nUMI_extra<-function(toPlot,fileName, png.device.type = NULL){
  toPlot<-toPlot[order(toPlot[,"nUMI.drop"]),]
  p_nGeneMitodrop <- ggplot(toPlot, aes(staticNr, nGene)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#ca2bd8"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_nUMIMitodrop <- ggplot(toPlot, aes(staticNr, nUMI)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#ca2bd8"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#ca2bd8")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_RBCIMitodrop <- ggplot(toPlot, aes(staticNr, percent.rbc)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#ca2bd8"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  p_COVIDMitodrop <- ggplot(toPlot, aes(staticNr, percent.COVID)) +
    geom_violin(fill="#e0e5e4") +
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop),size=0.8) +
    scale_color_manual(values=c("#cecece", "#ca2bd8"))+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  if(fileName != ""){

      ggsave( grid.arrange(p_nGeneMitodrop, p_nUMIMitodrop,p_mito,p_RBCIMitodrop,p_COVIDMitodrop, ncol=5,
                           top=" nUMI outliers (purple)"),file=fileName, dpi=200, units = "mm", width = 500, height = 100, type = png.device.type)
  }
}

##### Function drawUMI_mitoPlot
drawUMI_mitoPlot<-function(coordsTable, reductionType, clusterMatrix, columnName, titleInfo){
  columnNr<-which(colnames(clusterMatrix)==columnName)
  p <- ggplot()+
    geom_point(aes(x=sctTSNE_1,y=sctTSNE_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20)
  if(reductionType=="umap"){
    p <- ggplot()+
      geom_point(aes(x=sctUMAP_1,y=sctUMAP_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20)
  }
  p<-p +
    scale_colour_gradientn(colours = c("darkblue","cyan","green","yellow","orange","darkred")) +
    ggtitle(paste0(titleInfo," (",reductionType,")")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  return(p)
}
