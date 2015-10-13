
########################################################################################
##~~~~~~~~~~~~~~~~~~process download protein names match to organism abb name~~~~~~~~~~~
#read abbreviation bacteria and the protens sequence
bac_names <- read.csv("merge_fullnames_bacteria_update.csv", head = TRUE, stringsAsFactors = FALSE)
pro_sequence <-readLines("polymerase_subunit_beta.txt")
#process the data
bac_names <- bac_names[['V1']]
len_pro <- length(pro_sequence)
len_names <- length(bac_names)
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

write.table(update_seq, file = "process_polymerase_subunit_beta.txt", row.names = F, col.names = F, quote = F)
