
source('../ABA_work/API_ncbi_prodatabase.R')
##~~~~~~~~~~~~~~a short script for split names to some group and download~~~~~~~
#read  data
bac_names <- read.csv("14eukaryotes_noabb.csv", head = FALSE)
#process the data
bac_names <- as.character(bac_names[[1]])


bac_names_split <- bac_names[1:14]
conserve_protein <- "Ribosomal protein L13"
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

    
    write.table(hit_proteins, file="euk14_Ribosomal_protein_s12.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)
    #write the protein to the location
    write.table(miss_retrieve, file="euk14_Ribosomal_protein_s12_miss.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)
    elapsed_time <- proc.time() - time
	print(elapsed_time[3])
	#suspend 4 seconds because NCBI limit
    Sys.sleep(4)
	t <- b
	
}

