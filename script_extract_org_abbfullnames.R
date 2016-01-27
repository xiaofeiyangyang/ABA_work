############################################################################
##################extract organisms from save kegg html to local###########
#extract abbreviation and fullnames form kegg database
library(stringr)


##~~~~~~~~~~~~~A function use Regular expression to get organisms~~~~~~~~~~~~~~
extract_organ <- function(string){

    #Extract the organisms abbreviation 
	abb <- str_extract_all(string, "(?<=\\?org=)(.*?)(?='>)")

	#Extract the organisms fullname
	fullname <- str_extract_all(string, "(?<=www_bfind\\?T[\\d]{5}'>).*(?=</a>)")

	abb <- unlist(abb)

	fullname <- unlist(fullname)

	return(list(abb = abb, fullname= fullname))

}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##~~~~~~~~~~~~~~~~~~~~~read tabel and combine results~~~~~~~~~~~~~~~~~~~~~~~~~~


table <- read.csv("eukaryotes.csv", head = FALSE)

organ <- sapply(table, extract_organ)

results <- cbind(organ[[1]], organ[[2]])

write.csv(results, file = "output_eukaryotes.csv", row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##~~~~~~~~~~~~~~~~~~~~~~merge two table get the same organisms~~~~~~~~~~~~~~~~~~
x <- read.table("152eukaryotes.txt", head=FALSE)

y <- read.csv("output_eukaryotes.csv", head=TRUE)

mer <- merge(x, y,  all.x = T)

write.csv(mer, file="merge_fullnames_euk.csv", col.names = FALSE, row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



################################################################################
################################################################################