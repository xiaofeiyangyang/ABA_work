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
archaea <- read.table("bacteria.txt")
archaea <- as.character(archaea[,1])
#Regular expression to pattern names
archaea <- paste("\\b", archaea, sep = "")
archaea <- paste(archaea, "\\b", sep = "")

#extract organisms
extract_archeae <- extract_organisms(archaea, ath834)
#add Entrez
extract_archeae$Entrez <- c(1:nrow(extract_archeae))
#change Entrez to first col
extract_archeae <- data.frame(extract_archeae[613],extract_archeae[1:612])
#remove "ath:"
symbol_names <-extract_archeae$symbol
symbol_names <- sub(pattern = "ath:", replacement = "", symbol_names)
extract_archeae$symbol <- symbol_names
#write the phylogenetic matrix
write.table(extract_archeae, file="610bacteria.txt", quote = FALSE, sep = "\t", row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




###################################################################################