#######################The API access for NCBI protein database#################
library(XML)
library(RCurl)


Get_protein <- function(input_organisms,input_proteins){
	#input_proteins <- "30S Ribosomal protein S12"
	#input_organisms <- "Aquifex aeolicus"
	
	len_num <- length(input_organisms)

	#initialize the vector
    proteins <- vector(mode = "character", length = len_num)
    miss_proteins <- vector(mode = "character", length = 0)
    #the base NCBI API EutilsURL
    base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # set the protein database
    db <- "protein"
    #query protein name of the orgaism
    organism_name <- input_organisms
    protein_name <- input_proteins
    for(i in 1:len_num){
    	
    	query <- sprintf("%s[organism] AND %s[protein name]", organism_name[i], protein_name)
        #assemble the URL
        esearch_protein <- sprintf("esearch.fcgi?db=%s&term=%s", db, query)
        
        search_url <- paste0(base_url, esearch_protein)
        
        search_doc <- xmlParse(search_url)
        
        # get the first UIDs 
        uids <- xpathSApply(search_doc, path = "//IdList/Id", fun = 'xmlValue')
        uid <- uids[1]
        
        #dowload  proteins
        # rtrieve the protein id
        if(uid != "NULL"){
        	#get the unique ids
         
            id <- uid
            #identify the type and mode
            rettype <- "fasta"
            retmod <- "text"
            
            #assemble the URL
            
            efetch_protein <- sprintf("efetch.fcgi?db=%s&id=%s&rettype=%s&retmod=%s", db, id, rettype, retmod)
            
            fetch_url <- paste(base_url, efetch_protein, sep = '')
            
            sequence <- getURL(fetch_url)
            

        }
        #not retrieve the protein id
        else if(uid == "NULL"){
        	sequence <- organism_name[i]
        	miss_proteins[i] <- organism_name[i]

        }
        #download the proteins 
        proteins[i] <- sequence
       
    }
    #miss_proteins <- miss_proteins[-which(is.na(miss_proteins))]
    my_list <- list(proteins = proteins, miss_proteins = miss_proteins)
  

    return(my_list)

}

#download.file(fetch_url, destfile = 'ribosomal12.txt', method = 'wget', mode = 'ab')


##########################################################################################