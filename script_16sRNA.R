## a script scaning the whole genome to match ko

files <- dir()

for(i in seq_along(files)){
	print(i)
	load(files[i])
	#initiliza
	name <- vector(mode = "character", length = 0)
	#grep KO
 	gr <- grep("K01977", eachSpeInfo)
 	if(length(gr) == 0){
 		#miss match

 	    name <- files[i]
 	} else{
		gr <- gr[1]
 		}
 	#write miss match name to file
 	write.table(name, file = "miss_name.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
}



################################################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~a script scaning the whole genome to match 16s ribosomal RNA 


files <- read.table("/home/yangfang/phylogenetic_tree/miss_name.txt")
files <- as.character(files$V1)
for(i in seq_along(files)){
	print(i)
	load(files[i])
	#initiliza
	name <- vector(mode = "character", length = 0)
	#grep KO
 	gr <- grep("16S ribosomal RNA", eachSpeInfo)
 	if(length(gr) == 0){
 		#miss match

 	    name <- files[i]
 	} else{
		gr <- gr[1]
 		}
 	#write miss match name to file
 	write.table(name, file = "/home/yangfang/phylogenetic_tree/miss_name2.txt", append = TRUE, row.names = FALSE, col.names = FALSE)

}







 eachSpeInfo[gr][[1]]$NTSEQ