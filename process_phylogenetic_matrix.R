################################################################################


###################Get what we want organisms phylogentic matrix################

##~~~~~~~~~~~~~~~~~~~~~~~extract archaea phylogenetic matrix function~~~~~~~~~~~
extract_organisms <- function(pattern, all_organ){
	#extrcat all names
	name <- names(all_organ)
	#initialize
    indexsum <- vector(mode = "numeric", length = 0)
	for(i in pattern){
		gre <- grep(i,name)
		indexsum <- c(indexsum, gre) 
	}
    #update indexsum 
	indexsum <- c(1, indexsum)
	#get we pattern organisms
	organ <- all_organ[,indexsum]
    return(organ)
}
################################################################################
################################################################################



##~~~~~~~~~~~~~~~~~~~~~~write table archaea phylogentic matrix~~~~~~~~~~~~~~~~~~


#read the table and give some preprocessing requre function extract_organisms
ath834 <- read.table("834phy.txt", head = TRUE)
archaea <- read.table("600bac_abbnames.txt")
archaea <- as.character(archaea[,1])
#Regular expression to pattern names
#archaea <- paste("\\b", archaea, sep = "")
#archaea <- paste(archaea, "\\b", sep = "")

#extract organisms
#extract_archeae <- extract_organisms(archaea, ath834)
archaea <- c("symbol",archaea)
extract_archeae <- ath834[, names(ath834) %in% archaea]
#add Entrez
extract_archeae$Entrez <- c(1:nrow(extract_archeae))
#change Entrez to first col
len <- length(names(extract_archeae))
extract_archeae <- data.frame(extract_archeae[len],extract_archeae[1:len-1])
#remove "ath:"
symbol_names <-extract_archeae$symbol
symbol_names <- sub(pattern = "ath:", replacement = "", symbol_names)
extract_archeae$symbol <- symbol_names
#write the phylogenetic matrix
write.table(extract_archeae, file="600bac_matrix.txt", quote = FALSE, sep = "\t", row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####################################################################################
#a script change AT:XXXX names to we predict symbol
# read data
euk_mat <- read.table("600bac_matrix.txt", header = TRUE, stringsAsFactors = FALSE)
# read data
tag <- read.csv("112UGTs.csv", header = TRUE, stringsAsFactors = FALSE)
len <- length(tag$Locus.tag)
for(i in 1:len){
	index <- grep(tag$Locus.tag[i], euk_mat$symbol)
	euk_mat$symbol[index] <- tag$symbol[i]

}

write.table(euk_mat, file="600bac_matrix_for_UGTs.txt", quote = FALSE, sep = "\t", row.names = FALSE)





###################################################################################


#it's a script to conversion matrix names to what we want
ath834 <- read.table("834phy.txt", head = TRUE)
ath138 <- data.frame(ath834[1:139],ath834[338])
euk <- load("AddPhyloData0001.RData")

euk <- as.data.frame(AddPhyloData0001)

merge <- data.frame(ath138, euk[1:152])

names <- read.csv("138_fullnames.csv", head=FALSE, stringsAsFactors=FALSE)
re <- function(x){
	first_name <- strsplit(x, " ")[[1]][1]
	last_name <- strsplit(x," ")[[1]][2]
	re_names <- paste0(strtrim(first_name, 1), ".")
	re_names <- paste0(re_names, last_name)
	return(re_names)
}

re_names <- sapply(names$V2,re)	
re_names <- as.character(re_names)
match_names <- cbind(names$V1,re_names)
len <- length(names$V1)
for (i in 1:len){
	index <- grep(re_names[i], names(merge))
	names(merge)[index] <- names$V1[i]
}
names(merge)


write.table(merge, file="merge.txt", quote = FALSE, sep = "\t", row.names = FALSE)