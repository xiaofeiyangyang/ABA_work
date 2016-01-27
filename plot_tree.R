library(PhyloProfile)
library(ggtree)
library(ape)
data(fatp)
ATPphyloPlot <- PlotPhyloProfile(fatp$atpPhylo, speCol = fatp$specol, geneCol = fatp$genecol,
classCol = fatp$domain, legend.position = 'left')


ATPphyloPlot2 <- PlotPhyloProfile2(fatp$atpPhylo, speCol = fatp$specol, geneCol = fatp$genecol,
classCol = fatp$domain, legend.position = 'left')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##' Legend of species in phylogenetic profilings 
##'
##' Legend annotating the species categories, like the class or phylum.
##' @title Phylo legend
##' @param classCol A vector of colors with the names defining the species categories.
##' @param ... Parameters from geom_legend
##' @return ggplot2 object
##' @examples
##' data(fatp)
##' speLeg <- legend_spe(fatp$domain, legend.position = 'left')
##' \dontrun{
##' # plot legend
##' require(gridExtra)
##' grid.arrange(speLeg, ncol = 1)
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_manual
##' @export
##' 
legend_spe <- function(classCol, ...) {

  ## species object
  colMat <- data.frame(y = rep(0, length(classCol)),
                            x = 1:length(classCol),
                            fillCol = classCol)

  speBlockObj <- ggplot(colMat, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fillCol')) +
      scale_fill_manual(values = unname(classCol), name = 'Phylogeny', labels = names(classCol))

  speLegObj <- geom_legend(speBlockObj, ...)

  return(speLegObj)
}


##' Plot phylogenetic profiles with gene and species clusters.
##'
##' A combination plot of phylogenetic profiles with genes and species clusters.
##' 
##' @title Plot phylogenetic profiles
##' @param phyloData The phylogenetic profile data with 1 and 0 denoting the presence and absence of orthologous, respectively. The "phyloData" should be a numeric matrix, of which the row is gene and column is species.The "phyloData"has row names and column names which will be used for the dendrogram of row and column.
##' @param geneNameSize The size of gene names label, and the default value is 3.
##' @param geneNameCol The colour of gene names, and the default value is "grey55".
##' @param geneBetweenBlockCol The space color between gene blocks, and the default value is "NA" meaning no space color. If the number of genes is samll, for example less than 20, setting it as 'white' is fine. 
##' @param presentCol The color of present 1, the default value is "steelblue".
##' @param absentCol The color of present 0, the default value is "grey91".
##' @param speCol A vector of colors with names of species, which are the same as colnames of "phyloData" (may not in the same order). 
##' @param geneCol A vector of colors with names of genes, which are the same as rownames of "phyloData" (may not in the same order).
##' @param widthsShinkage The shinkage width vector.
##' @param heightsShinkage The shinkage width vector.
##' @inheritParams legend_spe
##' @return A plot object.
##' @examples
##' data(fatp)
##' ATPphyloPlot <- PlotPhyloProfile(fatp$atpPhylo, speCol = fatp$specol, geneCol = fatp$genecol,
##' classCol = fatp$domain, legend.position = 'left')
##' \dontrun{
##' # an example of saving output figures
##' cairo_pdf('FATPprofilePlot.pdf')
##' ATPphyloPlot <- PlotPhyloProfile(fatp$atpPhylo, speCol = fatp$specol, geneCol = fatp$genecol,
##' classCol = fatp$domain, legend.position = 'left')
##' dev.off()
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_text geom_tile geom_segment geom_point scale_fill_manual labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string coord_flip
##' @importFrom ggdendro dendro_data segment
##' @importFrom gridExtra grid.arrange
##' @importFrom reshape2 melt
##' @export
##' 


PlotPhyloProfile2 <- function(phyloData,
                             geneNameSize = 3,
                             geneNameCol = 'grey55',
                             geneBetweenBlockCol = NA,
                             presentCol = 'steelblue',
                             absentCol = 'grey91',
                             speCol,
                             geneCol,
                             classCol,
                             tipLabel,
                             widthsShinkage = c(0.3, 0.1, 9, 2.2),
                             heightsShinkage = c(0.3, 0.03, 0.05, 0.3),
                             ...) {

   require('grid')
   require('ggplot2')
   require('reshape2')
   require('ggdendro')
   require('gridExtra')
  #  phyloData <- as.matrix(proceMatrixUGTEcm)
  #  # phyloData <- fatp$atpPhylo[, 1:271]
  #  speCol <- tip_name
  #  geneCol <- colGene
  # legend.position = 'left'
  #  classCol = domain
   tip_label <- tipLabel


  # tmp1 <- ggtree(tree, branch.length="none", size = 0.5)
  # tmp2 <- tmp1$data
  # tmp3 <- tmp2[tmp2[, 7], ]
  # tip_label <- tmp3[order(tmp3[, 5], decreasing = TRUE), 6]
# cluster genes and species
  # hcGene <- hclust(dist(phyloData), method = 'average')
  # rowInd <- hcGene$order
  # hcSpe <- hclust(dist(t(phyloData)), method = 'average')
  # colInd <- hcSpe$order
  rowInd <-  nrow(phyloData):1
  colInd <- 1: ncol(phyloData)
  ## order 'phyloData'
  orderedPhyloData <- phyloData[rowInd,colInd]
  # orderedColNames <- colnames(orderedPhyloData)
  orderedColNames <- tip_label
  orderedRowNames <- rownames(orderedPhyloData)

  breaksRow <- 1:length(orderedRowNames)

  ## order 'geneCol'
  orderedGeneCol <- geneCol[match(rownames(orderedPhyloData), names(geneCol))]
  # orderedGeneCol <- rev(orderedGeneCol)
  ## order 'speCol'
  orderedSpeCol <- speCol[match(colnames(orderedPhyloData), names(speCol))]
  
  ## melt data for ggplot2
  colnames(orderedPhyloData) <- 1:ncol(orderedPhyloData)
  rownames(orderedPhyloData) <- 1:nrow(orderedPhyloData)
   orderedPhyloData <- melt(orderedPhyloData)
    # orderedPhyloData <- apply(orderedPhyloData, 2, rev)
  orderedPhyloData <- data.frame(geneNames = orderedPhyloData[, 1], speNames = orderedPhyloData[, 2], apData = factor(orderedPhyloData[, 3]))




  # ## plot spaces names
  # orderedRowNamesMat2 <- data.frame(
  #                                  x = seq(0.5, (length(tip_label)-0.5), 1),
  #                                  y = rep(0, length(tip_label)),
  #                                  fillName = tip_label)
  

  # geneNamesObj2 <- ggplot(orderedRowNamesMat2, aes_string('x', 'y', label = 'fillName')) +
  #   geom_text(size = 1, colour = geneNameCol, angle = 90) +
  #     labs(x = NULL, y = NULL) +

  #       scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  #         scale_x_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
  #           theme_pp(legend.position='none')


  orderedColNamesMat <- data.frame(y = rep(0, length(orderedColNames)),
                        x = seq(0.5, length(orderedColNames) - 0.5, 1),
                        fillName = orderedColNames)

  speciesNameSize <- 1


  speciesNamesObj <- ggplot(orderedColNamesMat, aes_string('x', 'y', label = 'fillName')) +
    geom_text(size = speciesNameSize, colour = geneNameCol, angle = 90) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, length(orderedColNames)), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
            theme_pp(legend.position='none')


  orderedRowNamesMat <- data.frame(x = rep(0, length(orderedRowNames)),
                        y = seq(0.5, length(orderedRowNames)-0.5, 1),
                        fillName = orderedRowNames)
  

  geneNamesObj <- ggplot(orderedRowNamesMat, aes_string('x', 'y', label = 'fillName')) +
    geom_text(size = 3, colour = geneNameCol) +
      labs(x = NULL, y = NULL) +
        scale_x_continuous(expand = c(0, 0), breaks = NULL) +
          scale_y_continuous(expand = c(0, 0), limits = c(0, length(orderedRowNames)), breaks = NULL) +
            theme_pp(legend.position='none')


  ## plot phylogenetic matrix
  phyloObj <- ggplot(orderedPhyloData, aes_string('speNames', 'geneNames')) +
    geom_tile(aes_string(fill = 'apData')) +
      scale_fill_manual(name = 'status', labels = c('absent', 'present'), values = c(absentCol, presentCol)) +
        labs(x = NULL, y = NULL) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
            scale_x_continuous(expand = c(0, 0), breaks = NULL) +
              theme_pp(legend.position='none')

  ## dendrogram plot for genes
  # hcTree <- as.hclust(tree)
  # ddata <- dendro_data(hcTree, type = 'rectangle')
  # segData <- segment(ddata)
  # segData[, c(1, 3)] <- segData[, c(1, 3)] - 0.5
  # geneDendroObj <- ggplot(segData) +
  #   geom_segment(aes_string(x = 'x', y = 'y', xend = 'xend', yend = 'yend')) +
  #     labs(x = NULL, y = NULL) +

      # coord_flip() +
        # scale_x_reverse(expand = c(0, 0), breaks = NULL) +
        # scale_y_reverse(expand = c(0, 0), breaks = NULL) +
          # scale_x_continuous(expand = c(0, 0), limits = c(0, ncol(phyloData)), breaks = NULL) +
              # theme_pp(legend.position='none')
   geneDendroObj <- ggplot(setDataTest) +
    geom_segment(aes_string(x = 'x', y = 'y', xend = 'xend', yend = 'yend'), size = 0.04) +
      labs(x = NULL, y = NULL) +
        scale_y_reverse(expand = c(0, 0), limits = c(length(treeNoLen$tip.label), 0), breaks = NULL) +
          scale_x_reverse(expand = c(0, 0), breaks = NULL) +
            coord_flip() +
              theme_pp(legend.position='none')      

  ## gene color block
  orderedGeneColMat <- data.frame(x = rep(0, length(orderedGeneCol)),
                                  y = 1:length(orderedGeneCol),
                                  fillCol = factor(orderedGeneCol))

  geneBlockObj <- ggplot(orderedGeneColMat, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fillCol'), color = geneBetweenBlockCol) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedGeneColMat$fillCol)) +
              theme_pp(legend.position='none')

  ## species color block
  orderedSpeColMat <- data.frame(y = rep(0, length(orderedSpeCol)),
                                 x = 1:length(orderedSpeCol),
                                 fillCol = factor(orderedSpeCol))

  speBlockObj <- ggplot(orderedSpeColMat, aes_string('x', 'y')) +
    geom_tile(aes_string(fill = 'fillCol')) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            scale_fill_manual(values = levels(orderedSpeColMat$fillCol)) +
              theme_pp(legend.position='none')

  ## species legend
  speLegObj <- legend_spe(classCol, ...)
 
  
  ## plot empty block
  emptyBlock <- geom_emptyblock()

  # plotRes <- marrangeGrob(
  #   list(empty, empty, empty, speBlockObj, geneDendroObj, geneNamesObj, geneBlockObj, phyloObj),
  #   ncol = 4,
  #   nrow = 2,
  #   widths = widthsShinkage,
  #   heights = heightsShinkage,
  #   top = NULL)

  plotRes <- grid.arrange(
 
    emptyBlock,
    emptyBlock,
    geneDendroObj,
    emptyBlock,


    emptyBlock,
    emptyBlock,
    speciesNamesObj,
    emptyBlock,

    emptyBlock,
    emptyBlock,
    speBlockObj,
    emptyBlock,
    
    geneNamesObj,
    geneBlockObj,
    phyloObj,
    speLegObj,
    ncol = 4,
    nrow = 4,
    widths = widthsShinkage,
    heights = heightsShinkage)

  return(plotRes)
} 
 cairo_pdf('figure4A.pdf', width = 20)
ATPphyloPlot2 <- PlotPhyloProfile2(phyloData, speCol = tip_name, geneCol = colGene,
classCol = domain, tipLabel = tip_label, legend.position = 'left')
 dev.off()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~process data~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/yangfang/mywork/eukaryotes/treeview/")
library("colorspace")
library(ggtree)
library('igraph')
tree <- read.tree("RAxML_bipartitions.T1")
 tmp1 <- ggtree(tree, branch.length="none", size = 0.5)
  tmp2 <- tmp1$data
  tmp3 <- tmp2[tmp2[, 7], ]
  tip_label <- tmp3[order(tmp3[, 5], decreasing = TRUE), 6]

UGTEcm <- read.table("UGT91family_0.txt", sep = "\t", header = TRUE, stringsAsFactor = FALSE)

proceUGTEcm <- UGTEcm[UGTEcm$ECM.ECM. == 'ECM',]
matrixUGTEcm <-proceUGTEcm[ , 13:length(proceUGTEcm[1,])-1]
proceMatrixUGTEcm <- matrixUGTEcm[ ,tip_label]
rownames(proceMatrixUGTEcm) <- proceUGTEcm$Gene.Symbol
ecmNum <- unique(proceUGTEcm$ECM.ID)
len_ecmNum <- length(ecmNum)
col <- rainbow_hcl(8)[1:len_ecmNum]
colGene <- NULL
for(i in 1:len_ecmNum){
	# i <- 1
	len <- length(proceUGTEcm$ECM.ID[proceUGTEcm$ECM.ID == i])
	colLen <- rep(col[i],len)
	colGene <- c(colGene,colLen)

}
#colGne 

UGTname <- proceUGTEcm$Gene.Symbol
 names(colGene) <- UGTname

phyloData <- as.matrix(proceMatrixUGTEcm)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



speccol <- rainbow_hcl(8)[6:8]

kingdom <- read.table("271euk_taxonomy4.txt", sep = "\t", header = FALSE, stringsAsFactor = FALSE)
specName <- kingdom$V1
phy <- kingdom$V9
for (i in 1: length(phy)){

	if ( is.na(phy[i])){
		phy[i] <- "#FFFFFF"
	} else if( phy[i] == "Fungi"){
		phy[i] <- speccol[1]
	} else if( phy[i] == "Viridiplantae"){
		phy[i] <- speccol[2]
	} else if( phy[i] == "Metazoa"){
		phy[i] <- speccol[3]
	} else{}
}
domain <- speccol[1:3]
names(domain) <- c("Fungi", "Viridiplantae", "Metazoa")
names(phy) <- specName

tip_name <- NULL
for(i in seq_along(tip_label)){
	for(j in seq_along(phy)){
		if(tip_label[i] == names(phy[j])){
			tip_name <- c(tip_name, phy[j])
		}
	}
}
# specols <- tip_label
# names(tip_name) <- specols









tree <- read.tree("RAxML_bipartitions.T1")
treeNoLen <- tree
treeNoLen$edge.length <- NULL
treeData <- ggplot(treeNoLen)$data
setData1 <- data.frame(x = treeData$x[treeData$parent],
	xend = treeData$x,
	y = treeData$y,
	yend = treeData$y)
setData2 <- data.frame(x = treeData$x[treeData$parent],
	xend = treeData$x[treeData$parent],
	y = treeData$y[treeData$parent],
	yend = treeData$y)
setData <- rbind(setData1, setData2)
setDataTest <- setData
setDataTest[, c(3, 4)] <- setData[, c(3, 4)] - 0.5

















theme_pp2 <- function(...) {
  theme_bw() %+replace%
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'mm'),
        axis.line = element_blank(),
        panel.margin = unit(0, 'mm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'line'),
        legend.margin = unit(0, 'mm'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        ...)
}


theme_pp <- function(...) {
  theme_bw() %+replace%
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0, 'mm'),
        axis.ticks.margin = unit(0, 'mm'),
        axis.line = element_blank(),
        panel.margin = unit(0, 'mm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), 'line'),
        legend.margin = unit(0, 'mm'),
        ...)
}


geom_emptyblock<- function(...) {
  emptyData <- data.frame(x = 1, y = 1)
  empty <- ggplot(emptyData) +
    geom_point(aes_string('x', 'y'), colour='white', ...) +
      labs(x = NULL, y = NULL) +
        scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          scale_x_continuous(expand = c(0, 0), breaks = NULL) +
            theme_pp(legend.position='none')

  return(empty)
}


geom_legend <- function(g, ...) {
  ## !!!Be aware of this dirty walkaround!!!
  pdf(file = NULL)
  g <- ggplotGrob(g + theme(legend.margin = unit(0, 'mm'), ...))$grobs
  dev.off()
  
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

  return(legend)
}
