
########################################################################################
##~~~~~~~~~~~~~~~~~~process download protein names match to organism abb names~~~~~~~~~~~
#read abbreviation bacteria and the protens sequence
files <- dir()
len <- length(files)
bac_names <- read.csv("merge_fullnames_archaea_update.csv", head = FALSE, stringsAsFactors = FALSE)
bac_names <- bac_names[['V2']]

for(k in 1:len){
	k <- 8
    pro_sequence <-readLines(files[k])
    #process the data
    len_pro <- length(pro_sequence)
    len_names <- length(bac_names)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~match to organism abb names~~~~~~~~~~~~~~~~~~~~~~~~
    #initilize
    update_seq <- 0
    j <- 1
    for(i in 1:len_pro){
    	#patter organism names
    	pattern <- length(grep("gi\\|.*", pro_sequence[i]))
    	if(pattern == 1){
    		#replace with the abbreviation organisms names
    		re <- sub(pattern = "gi\\|.*", replace = bac_names[j], pro_sequence[i])
            update_seq[i] <- re
            j <- j+1
    	}else{
    		#update which not pattern
    		update_seq[i] <- pro_sequence[i]
    	}
    }
    
    write.table(update_seq, file = paste0("process_",files[k]), row.names = F, col.names = F, quote = F)

}
