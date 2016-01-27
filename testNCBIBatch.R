

mydata <- read.csv("30S12.txt", sep = "\t",head= FALSE)
num <- read.csv("archaea_taxonomy_num.txt",head =FALSE, sep = "\t")
taxon <- read.csv("81arc_taxonomy_number.txt",head =FALSE, sep = "\t")
num_in_mydata <- mydata$V2[mydata$V2 %in% taxon$V2]
num_in_mydata_un <- unique(num_in_mydata)
taxon_in_mydata <- mydata[mydata$V2 %in% num_in_mydata_un, 1:2]
taxon_in_mydata <- sapply(split(taxon_in_mydata[, 1], taxon_in_mydata[, 2]), '[[', 1)


num_in_mydata <- mydata$V2[mydata$V2 %in% num$V3]
un_num <- num$V3[!is.na(num$V3)]
num_in_mydata <- mydata$V2[mydata$V2 %in% un_num]

num_in_mydata_un <- unique(num_in_mydata)
mydata_uni <- mydata$V1[mydata$V2 %in% unique(mydata$V2)]
num_in_mydata <- mydata[as.character(mydata$V2) %in% num_in_mydata_un, 1:2]
num_in_mydata <- sapply(split(num_in_mydata[, 1], num_in_mydata[, 2]), '[[', 1)

library(stringr)

## transfer data.frame to character matrix
mydata <- apply(mydata, 1:2, function(x) {
	eachEle <- as.character(x)
	eachEle <- str_trim(x)

	return(eachEle)
})

## merge KEGG species ID with no_rank and species
taxList <- lapply(1:nrow(num), function(x){
	eachTax <- c(taxon[x, 2], num[x, 3])
	eachTax <- eachTax[!is.na(eachTax)]

	return(eachTax)
})
names(taxList) <- num[, 1]

## retrieve gi number
giList <- lapply(1:length(taxList), function(x) {
	eachData <- mydata[mydata[, 2] %in% taxList[[x]], ,drop = FALSE]

	return(eachData)
})
names(giList) <- num[, 1]

zeroTax <- taxList[which(sapply(giList, nrow) == 0)]
