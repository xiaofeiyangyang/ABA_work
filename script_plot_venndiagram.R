#############################################################################################
# venn diagram for prodiction by ABA in archaea
# read data
library(VennDiagram)
setwd("/home/yangfang/clime/ABA/archaea/output_ABA/")





cdpk <- read.table("ABA_CDPKs_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cipk <- read.table("ABA_CIPKs_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mapk <- read.table("ABA_MAPKs_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rlck <- read.table("ABA_RLCKs_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rlk <- read.table("ABA_RLKs_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cdpks <- as.character(na.omit(cdpk$Gene.Symbol[cdpk$LLR > 20]))
cipks <- as.character(na.omit(cipk$Gene.Symbol[cipk$LLR > 20]))
mapks <- as.character(na.omit(mapk$Gene.Symbol[mapk$LLR > 20]))
rlcks <- as.character(na.omit(rlck$Gene.Symbol[rlck$LLR > 20]))
rlks <- as.character(na.omit(rlk$Gene.Symbol[rlk$LLR > 20]))

list1 <- list(CDPKs = cdpks, CIPKs = cipks, MAPKs = mapks, RLCKs = rlcks, RLKs = rlks )
#venn for three
venn.diagram(list1, fill = rainbow(3), lty = "dotted", "out3.tiff")
#venn for five
venn.diagram(list1, fill = rainbow(5), lty = "dotted", "five_pathway.tiff", cat.cex = 1, lwd = 1,cat.dist = 0.1, margin = 0.1)

#lty = "blank",
#rotation.degree = 270,


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#venn diagram for prodiction by UGTs in eukaryotes
library(VennDiagram)

setwd("/home/yangfang/clime/ABA/eukaryotes/output_UGT/")
ugt71b <- read.table("UGT71Bfamily_0_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ugt71c <- read.table("UGT71Cfamily_0_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ugts4 <- read.table("4UGTS_ABA_0_0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ugt71bs <- as.character(na.omit(ugt71b$Gene.Symbol[ugt71b$LLR > 180]))
ugt71cs <- as.character(na.omit(ugt71c$Gene.Symbol[ugt71c$LLR > 180]))
ugts4s <- as.character(na.omit(ugts4$Gene.Symbol[ugts4$LLR > 180]))
list1 <- list(UGT71B = ugt71bs, UGT71C = ugt71cs, UGTS4 = ugts4s) 
venn.diagram(list1, fill = rainbow(3), "4UGTs.tiff")
