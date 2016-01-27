

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~annotation 610 bacteria tree with full names by phylum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggtree)
library(ggplot2)
# read data
setwd("/home/yangfang/mywork/archaea2/treeview/")
taxon <- read.csv("81archaea_taxonomy.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("RAxML_bipartitions.T1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/yangfang/mywork/eukaryotes/treeview/")
taxon <- read.csv("271euk_taxonomy.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("RAxML_bipartitions.T1")
#tree <- read.tree("RAxML_bipartitions.nwk")

phylum <- unique(taxon$V2)
uni_phylum <- phylum[!is.na(unique(phylum))]
#group1 <- taxon$V1[taxon$V2 %in% uni_phylum[1]]
cls <- sapply(uni_phylum, function(x){
	each_group <- taxon$V1[taxon$V2 %in% x]
})
cls2 <- cls
cls2 <- cls
tree <- groupOTU(tree,cls2)
library("colorspace")
#pdf(file="myplot3.pdf")
#test
ggtree(tree, aes(color=group),branch.length="none", layout = "fan",size = 0.5,ledderiza=FALSE)  + 
theme(legend.position="right")  + geom_text(aes(label=label, angle=angle+90),size = 3, vjust = 0.3, hjust = -.3)

#a test too
ggtree(tree, aes(color=group), layout = "circular",size= .5,right = TRUE, width =1)  + 
theme(legend.position="right") + geom_tiplab( aes(angle=angle+90),align = TRUE,hjust = -.1, vjust = 0.3, linetype="dotted", linesize = .5, size = 0.8)

#right to circular
cols <- rainbow(7,start=0.2)
ggtree(tree, aes(color=group),branch.length="none", layout = "circular",size = 0.3) + 
    geom_tiplab( aes(angle=angle+90),align = TRUE,hjust = -.3, vjust = 0.3, linetype="dotted", linesize = .3) +
         theme(legend.position="right") +
             scale_color_manual(values=cols, breaks=1:3,labels=c("Setosa", "Versicolor", "Virginica"))  




#dev.off()
p <- ggtree(tree, aes(color=group),branch.length="none", layout = "circular")  
p + geom_text(aes(label=label, angle=angle+90), vjust =-.2, na.rm = T)


ggtree(tree, aes(color=group, linetype=group)) + geom_text(aes(label=label),  hjust=-.25)  + theme(legend.position="right")
#################################################################################################################################\
################################################################################################################################

#~~~~~~~~~~~~~~~~~~~~annotation archaea tree with full names by phylum~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggtree)
library(ggplot2)
library(colorspace)
#read tree file fullnames file and taxonomy file
setwd("/home/yangfang/mywork/archaea2/treeview/")
taxon <- read.csv("81archaea_taxonomy.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("RAxML_bipartitions.T1")
names <- read.csv("merge_fullnames_archaea_update.csv", head=FALSE, stringsAsFactors=FALSE)
#process the data
re_names <- as.character(names$V1)
len <- length(names$V1)
#replacement the abbreviation names by full names
for (i in 1:len){
    # assemble regular expression
	abb_names <- paste0("\\b", re_names[i])
	abb_names <- paste0(abb_names, "\\b")
	# match the abbreviation get index
	index <- grep(abb_names, tree$tip.label)
	index2 <- grep(abb_names, taxon$V1)
	#replacement
	tree$tip.label[index] <- names$V2[i]
	taxon$V1[index2] <- names$V2[i]
}


#extract phylum
phylum <- unique(taxon$V2)
uni_phylum <- phylum[!is.na(unique(phylum))]
uni_phylum_num <- length(uni_phylum)
#group1 <- taxon$V1[taxon$V2 %in% uni_phylum[1]]
cls <- sapply(uni_phylum, function(x){
	each_group <- taxon$V1[taxon$V2 %in% x]
})
cls2 <- cls
#set group by groupOTU function
tree <- groupOTU(tree,cls2)

#pdf(file="81arc_circular.pdf")
#chose the color 
cols <- rainbow_hcl(uni_phylum_num+1)
#view the tree
ggtree(tree, aes(color=group),branch.length="none", layout = "circular", size = 0.5) + 
    geom_tiplab( aes(angle=angle+90),align = TRUE, hjust = -0.1, vjust = 0.3, linetype="dotted", linesize = 0.1, size = 2) +
         theme(legend.position="right") +
             scale_color_manual(values=cols, breaks=1:uni_phylum_num,labels=uni_phylum)  
#dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~annotation eukaryotes tree with full names by phylum~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggtree)
library(ggplot2)
library("colorspace")
#read tree file fullnames file and taxonomy file
setwd("/home/yangfang/mywork/eukaryotes/treeview/")
taxon <- read.csv("271euk_taxonomy2.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("RAxML_bipartitions.T1")
names <- read.csv("merge_fullnames_euk.csv", head=FALSE, stringsAsFactors=FALSE)
#process the data
re_names <- as.character(names$V1)
len <- length(names$V1)
#replacement the abbreviation names by full names
for (i in 1:len){
    # assemble regular expression
	abb_names <- paste0("\\b", re_names[i])
	abb_names <- paste0(abb_names, "\\b")
	# match the abbreviation get index
	index <- grep(abb_names, tree$tip.label)
	index2 <- grep(abb_names, taxon$V1)
	#replacement
	tree$tip.label[index] <- names$V2[i]
	taxon$V1[index2] <- names$V2[i]
}

for(i in seq_along(taxon$V3)){
	if(is.na(taxon$V3[i])){

	taxon$V3[i] <- taxon$V2[i]
	} 
else{}
}




#extract phylum
phylum <- unique(taxon$V3)
uni_phylum <- phylum[!is.na(unique(phylum))]
uni_phylum_num <- length(uni_phylum)
#group1 <- taxon$V1[taxon$V2 %in% uni_phylum[1]]
cls <- sapply(uni_phylum, function(x){
	each_group <- taxon$V1[taxon$V3 %in% x]
})
cls2 <- cls
#set group by groupOTU function
tree <- groupOTU(tree,cls2)

#pdf(file="270euk_circular3new.pdf", width= 15, height = 15)
#chose the color 
cols <- rainbow_hcl(uni_phylum_num+1)
#view the tree
ggtree(tree, aes(color=group), branch.length="none", layout = "circular", size = 0.3) + 
    geom_tiplab( aes(angle=angle+90),align = TRUE, hjust = 0, vjust = 0.3, linetype="dotted", linesize = 0.1, size = 1.5) +
         theme_tree(legend.position="right") +
             scale_color_manual(values=cols, breaks=1:uni_phylum_num+1,labels=uni_phylum)  
#dev.off()


ggtree(tree, branch.length="none", layout = "circular", size = 0.3) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  annotation bacteria tree with full names by phylums

library(ggtree)
library(ggplot2)
library("colorspace")
#read tree file fullnames file and taxonomy file
setwd("/home/yangfang/mywork/bactria2/viewtree/")
taxon <- read.csv("abb_610_names_taxonomy2.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
tree <- read.tree("RAxML_bipartitions.T1")
names <- read.csv("abb_600_fullnames.csv", head=TRUE, stringsAsFactors=FALSE)
#process the data
re_names <- as.character(names$V1)
len <- length(names$V1)
#replacement the abbreviation names by full names
for (i in 1:len){
    # assemble regular expression
	abb_names <- paste0("\\b", re_names[i])
	abb_names <- paste0(abb_names, "\\b")
	# match the abbreviation get index
	index <- grep(abb_names, tree$tip.label)
	index2 <- grep(abb_names, taxon$V1)
	#replacement
	tree$tip.label[index] <- names$V2[i]
	taxon$V1[index2] <- names$V2[i]
}


#extract phylum
phylum <- unique(taxon$V2)
uni_phylum <- phylum[!is.na(unique(phylum))]
uni_phylum_num <- length(uni_phylum)
#group1 <- taxon$V1[taxon$V3 %in% uni_phylum[1]]
cls <- sapply(uni_phylum, function(x){
	each_group <- taxon$V1[taxon$V2 %in% x]
})
cls2 <- cls
#set group by groupOTU function
tree <- groupOTU(tree,cls2)

#pdf(file="610bac_circular02.pdf")
#chose the color 
cols <- rainbow_hcl(uni_phylum_num+1,start=.3) 
#view the tree
ggtree(tree, aes(color=group),branch.length="none", layout = "circular", size = 0.5) + 
    geom_tiplab( aes(angle=angle+90),align = TRUE, hjust = -0.1, vjust = 0.3, size = 0.8) +
         theme(legend.position="right") +
             scale_color_manual(values=cols, breaks=1:uni_phylum_num,labels=uni_phylum)  
#dev.off()
ggtree(tree, branch.length = "none", layout = "circular", size = 0.5)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






pdf(file="138science_circular.pdf")
#ggtree(tree, branch.length = "none",layout = "circular")
ggtree(tree, aes(color=group),branch.length="none", layout = "circular")  + 
theme(legend.position="right")
# +  geom_text(aes(label=label, angle=angle), size=3, color="purple", vjust=-0.3)
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~circular plot tree by the way~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cls2 <- cls[1:15]
p <- ggtree(tree, branch.length = "none", layout = "circular") 
#pdf(file = "271eukcircular.pdf")
groupOTU(p, cls2) + aes(color=group) + geom_tiplab(size = 1, vjust = 1, hjust = 0.5) + scale_color_manual(values=c("black", rainbow_hcl(length(cls2)))) + theme(legend.position="right")
#dev.off()
 p <- ggtree(tree)
groupOTU(p, cls2) + aes(color=group) + geom_tiplab() 






Examples
file <- system.file("extdata/BEAST", "beast_mcc.tree", package="ggtree")
beast <- read.beast(file)
plot(beast, annotation="length_0.95_HPD", branch.length="none") + theme_tree()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~chose color~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(colorspace)
x = rainbow_hcl(64)
for(i in x ){
	write.table(i, file = "color_hcl2.txt",append =TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE )
}
