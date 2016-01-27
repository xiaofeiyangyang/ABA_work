
#############a script retrieve  orgnism protein UID from NCBI################################

protein_name <- "30S Ribosomal protein S12"
organism <- read.csv("81_bac_fullnames.txt", head = FALSE)
#process the data
organism <- as.character(organism[[1]])
organism_name <- organism[34:81]
len_num <- length(organism_name)
ids <- vector(mod = "character", length=len_num)
for(i in 1:len_num){
	base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    db <- "protein"
    query <- sprintf("%s[organism] AND %s[protein name]", organism_name[i], protein_name)
    #assemble the URL
    esearch_protein <- sprintf("esearch.fcgi?db=%s&term=%s", db, query)
    
    search_url <- paste0(base_url, esearch_protein)
    
    search_doc <- xmlParse(search_url)
    
    # get the first UIDs 
    uids <- xpathSApply(search_doc, path = "//IdList/Id", fun = 'xmlValue')
    uid <- uids[1]
    print(uid)
    ids[i] <- rbind(uid)
    Sys.sleep(5)
}
ids <- cbind(ids)
write.table(ids, file= "81esearch_num.txt", quote = FALSE,append =TRUE)