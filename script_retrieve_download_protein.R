
source('../API_ncbi_prodatabase.R')
##~~~~~~~~~~~~~~a short script for split names to some group and download~~~~~~~
#read  data
bac_names <- read.csv("81_bac_fullnames.txt", head = FALSE)
#process the data
bac_names <- as.character(bac_names[[1]])


bac_names_split <- bac_names[1:10]
conserve_protein <- "30S Ribosomal protein S12"
# split organisms to group retrieve
#initilize vocter
i <- 1
t <- 0
#split  many group what you want
num <- 5
num_group <- ceiling(length(bac_names_split)/num)

for(j in 1:num_group){
	# set a time used
	time <- proc.time()
	a <- i +t
	b <- a+(num-1)
	split_names <- bac_names_split[a:b]
	split_names <- split_names[!is.na(split_names)]
	#download proteins
    proteins <- Get_protein(split_names, conserve_protein)
    
    #proteins <- Get_protein(bac_names_split, conserve_protein)

    hit_proteins <- proteins$proteins
    miss_retrieve <- proteins$miss_proteins
    #write the protein to the location

    
    write.table(hit_proteins, file="arc_ribosonmal_pro_s12.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)
    #write the protein to the location
    write.table(miss_retrieve, file="arc_ribosonmal_pro_s12_miss.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)
    elapsed_time <- proc.time() - time
	print(elapsed_time[3])
	#suspend 4 seconds because NCBI limit
    Sys.sleep(4)
	t <- b
	
}


#write.csv(group, file = paste0("names", as.character(j),".csv"),row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
