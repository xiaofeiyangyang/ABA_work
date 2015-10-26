library(stringr)
library(KEGGBioCycAPI)

require('Biostrings')
#read data and process
filenames <- dir()
len <- length(filenames)
#the abb names we want to remove
abb_names <- c("ert:", "gpa:", "ndl:", "rix:", "rum:", "tpq:", "euc:", "hci:", "ypx:","csk:")
# initiliza vector
lenabb <- length(abb_names)
index <- vector(mode = "numeric", length = 0)
k <- 1
for(i in 1:len){
	#read AA sequences
    seqs <- readAAStringSet(filenames[i], format="fasta")
    #found every  sequences and remove 
    for(j in 1:lenabb){

        dex <- grep(abb_names[j], names(seqs))
     
        if(length(dex) != 0){
        	index[j] <- dex
        }else{
        	index[j] <- NA
        }
        index <- index[!is.na(index)]
        #remove the match sequences
        deal_seqs <- seqs[-index]
    }
    #overwrite to filesp
    writeXStringSet(deal_seqs, filenames[i])
}
