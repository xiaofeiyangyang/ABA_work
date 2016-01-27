



#~~~~~~~~~~~~~~~~~download the ribosomal with  not split to group~~~~~~~~~~~~~~~~

source("../API_ncbi_prodatabase.R")
#read  data

bac_names <- read.csv("merge_fullnames_bacteria_update.csv", head = TRUE)
#process the data
bac_names <- bac_names[3]
bac_names <- apply(bac_names, 2 ,as.character)
bac_names <- bac_names[1]
conserve_protein <- "30S Ribosomal protein S12"
proteins <- Get_protein(bac_names, conserve_protein)
#write the protein to the location
write.table(proteins, file="proteins_s12.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~process download protein append miss protein~~~~~~~~~~

source("../API_ncbi_prodatabase.R")
#down load the miss protein 
pro_name <- read.csv("proteins_L1_miss.txt", head = FALSE, sep ="\t")
#pro_name_L3 <- read.csv("proteins_L3_miss.txt", head = FALSE, sep ="\t", stringsAsFactors = FALSE)
#pro_name_L3 <- pro_name_L3[['V1']]
pro_name <- as.character(pro_name$V2)

conserve_protein <- "30S Ribosomal protein s12"
proteins <- Get_protein(pro_name, conserve_protein)
hit_proteins <- proteins$proteins
miss_retrieve <- proteins$miss_proteins
write.table(hit_proteins, file="proteins_s12_add.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)
#write the protein to the location
write.table(miss_retrieve, file="proteins_s12_addmiss.txt",row.name = FALSE, col.name = FALSE, quote = FALSE, append = TRUE)




##############################################################################################
##############################################################################################


##~~~~~~~~~~~~~~~~~~~~~~~scrpit to check which proteins are not match by id~~~~~~~~~~~~~~~
  
  input_organisms <- bac_names[1:100]
  len_num <- length(input_organisms)
  proteins <- vector(mode = "character", length = 0)
  input_proteins <- "30S Ribosomal protein S12"
  base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # set the protein database
    db <- "protein"
    #query protein name of the orgaism
    organism_name <- input_organisms
    protein_name <- input_proteins
    for(i in 1:len_num){
    	time <-  proc.time()
    	query <- sprintf("%s[organism] AND %s[protein name]", organism_name[i], protein_name)
        #assemble the URL
        esearch_protein <- sprintf("esearch.fcgi?db=%s&term=%s", db, query)
        
        search_url <- paste0(base_url, esearch_protein)
        
        search_doc <- xmlParse(search_url)
         
        uids <- xpathSApply(search_doc, path = "//IdList/Id", fun = 'xmlValue')
        uid <- uids[1]
        #esearch_names <- xpathSApply(search_doc, path = "//Translation/From", fun = 'xmlValue')
        #esearch_names <- gsub(pattern = "\\[organism\\]", replacement ="", esearch_names)
        #if(length(esearch_names) == 0){
        #	esearch_names <- i
        #}else{
        #	esearch_names <- esearch_names
        #}
        proteins[i] <-uid
        elapsed_time <- proc.time() - time
	    print(elapsed_time[3])
	    print(i)
        Sys.sleep(1)    
    }  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




