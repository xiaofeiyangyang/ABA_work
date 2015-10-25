
library(ggtree)
library(ggplot2)
#read tree and  match names files
names <- read.csv("138_fullnames.csv", head=FALSE, stringsAsFactors=FALSE)
tree <- read.tree("species138.abbrev.manual_binary.nwk")
#a function to change full names to abbreviation example X.XXXX
re <- function(x){
	first_name <- strsplit(x, " ")[[1]][1]
	last_name <- strsplit(x," ")[[1]][2]
	re_names <- paste0(strtrim(first_name, 1), ".")
	re_names <- paste0(re_names, last_name)
	return(re_names)
}
#conversion to abb names
re_names <- sapply(names$V2,re)	
re_names <- as.character(re_names)
#cbind abb names
match_names <- cbind(names$V1,re_names)
len <- length(names$V1)
#match and replacement
for (i in 1:len){
	index <- grep(re_names[i], tree$tip.label)
	tree$tip.label[index] <- names$V1[i]
}
tree$tip.label