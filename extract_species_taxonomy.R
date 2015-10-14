#######################a script to extract abb orgnisms's taxonomy#################

#read abbreviation names and procession

org_names <- read.table("abb_610_names.txt",head = FALSE)
org_names <- as.character(org_names[[1]])
org_name_group <- org_names[1:610]
#initilize vactor
i <- 1
t <- 0
#split  many group what you want
num <- 5
num_group <- ceiling(length(org_name_group)/num)

for(j in 1:num_group){
	# set a time used
	time <- proc.time()
	a <- i +t
	b <- a+(num-1)
	split_names <- org_name_group[a:b]
	split_names <- split_names[!is.na(split_names)]
	print(split_names)
	#get names
	for(k in split_names){
        taxonomy_names <- Get_taxonomy(k)
        #rbind get names
        taxonomy_names<- rbind(taxonomy_names)
        #write to file
        write.table(taxonomy_names, file="610bacteria_taxonomy.txt", sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE, row.names = FALSE)  
	}
    
    elapsed_time <- proc.time() - time
	print(elapsed_time[3])
	#suspend 4 seconds because NCBI limit
    Sys.sleep(4)
	t <- b
	
}

#####################################################################################
#####################################################################################



##################function to extract species#######################################
library(XML)
library(RCurl)
Get_taxonomy <- function(abb_names){
    names_abb <- abb_names
    #assemble kegg url
    kegg_url <- sprintf("http://www.kegg.jp/kegg-bin/show_organism?org=%s", names_abb)
    #Parse kegg url
    kegg_doc <- htmlParse(kegg_url)
    #obtain organism taxonomy numbers
    taxon_num <- xpathSApply(kegg_doc, path = "//a", fun = 'xmlValue')[9]
    #use numbers to assemble ncbi url
    ncbi_url <- sprintf("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s", taxon_num)
    #Parse ncbi url
    ncbi_doc <- htmlParse(ncbi_url)
    #six taxonomy
    taxon_name <- c("phylum", "class", "order", "family", "genus","species")
    #initialize
    len <- length(taxon_name)
    taxon_org <- vector(mode = "character", length = len)
    for(i in 1:len){
    	#assemble the attribute
        att <- sprintf("//a[@alt='%s']",taxon_name[i])

        #get node set
        taxon<- getNodeSet(ncbi_doc, att)
         
        #extarct xmlvalue 
        location <- sapply(taxon, xmlValue)
        if(length(location) == 0){
         	location <- NA
        }else {
         	location <- location
        }
        taxon_org[i] <- location 

    }
    #aad organism names(abbreviation)
    taxon_six <- c(names_abb,taxon_org)
    #taxon_six <- taxon_org
    return(taxon_six)	
}


#########################################################################################
#########################################################################################



###############a script retrieve orgnisms's taxonomy number form ncbi###################
library(XML)
library(RCurl)

org_names <- read.table("archaea.txt",head = FALSE)
org_names <- as.character(org_names[[1]])
org_name_group <- org_names[1:81]
#initilize vactor
i <- 1
t <- 0
#split  many group what you want
num <- 5
num_group <- ceiling(length(org_name_group)/num)

for(j in 1:num_group){
	# set a time used
	time <- proc.time()
	a <- i +t
	b <- a+(num-1)
	split_names <- org_name_group[a:b]
	split_names <- split_names[!is.na(split_names)]
	print(split_names)
	#get names
	for(k in split_names){
        #assemble kegg url
        kegg_url <- sprintf("http://www.kegg.jp/kegg-bin/show_organism?org=%s", k)
        #Parse kegg url
        kegg_doc <- htmlParse(kegg_url)
        #obtain organism taxonomy numbers
        taxon_num <- xpathSApply(kegg_doc, path = "//a", fun = 'xmlValue')[9]
        #rbind get names
        taxonomy_names<- cbind(k,taxon_num)
        #write to file
        write.table(taxonomy_names, file="81arc_taxonomy_number.txt", sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE, row.names = FALSE)  
	}
    
    elapsed_time <- proc.time() - time
	print(elapsed_time[3])
	#suspend 4 seconds because NCBI limit
    Sys.sleep(3)
	t <- b
	
}


#########################################################################################
#########################################################################################

################################function to get species id###############################
library(XML)
library(RCurl)
library(stringr)
    
Get_species_num <- function(org_names){
	names_abb <- org_names
    #assemble kegg url
    kegg_url <- sprintf("http://www.kegg.jp/kegg-bin/show_organism?org=%s", names_abb)
    #Parse kegg url
    kegg_doc <- htmlParse(kegg_url)
    #obtain organism taxonomy numbers
    taxon_num <- xpathSApply(kegg_doc, path = "//a", fun = 'xmlValue')[9]
    #use numbers to assemble ncbi url
    ncbi_url <- sprintf("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s", taxon_num)
    #Parse ncbi url
    ncbi_doc <- htmlParse(ncbi_url)
    

    rank <- xpathSApply(ncbi_doc, path = "//em[3]", fun = 'getSibling')
    rank  <- sapply(rank, xmlValue)
    #getNodeSet(ncbi_doc, "//em")
    #six taxonomy
    taxon_name <- c("Superspecies","species","subspecies")
    #initialize
    len <- length(taxon_name)
    taxon_org <- vector(mode = "character", length = len)
    for(i in 1:len){
    	#assemble the attribute
        att <- sprintf("//a[@alt='%s']",taxon_name[i])

        #get node set
        taxon<- getNodeSet(ncbi_doc, att)
         
        #extarct xmlvalue 
        #location <- sapply(taxon, xmlValue)
        id <- sapply(taxon, xmlGetAttr,"href")
        uid <- str_extract_all(id, "(?<=&id=)(.*?)(?=&lvl)")

        if(length(uid) == 0){
         	uid <- NA
        }else {
         	uid <- uid
        }
        taxon_org[i] <- uid 

    }
    #aad organism names(abbreviation)
    taxon_spec <- c(names_abb,as.character(taxon_org), rank)
    taxon_spec <- rbind(taxon_spec)
    return(taxon_spec)


}



##~~~~~~~~~~~~~~~~~~~script to extract species id~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#read table and process
names <- read.csv("merge_fullnames_archaea_update.csv", header = FALSE)
names <- as.character(names[[2]])
names_spec <- names[1:81]
len <- length(names_spec)

for(i in 1:len){
	#used time
	time <- proc.time()
	#retrieve the sepcise
	taxon_spec <- Get_species_num(names_spec[i])
    #save to local disk
    write.table(taxon_spec, file = "archeae_taxonomy_num.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
    Sys.sleep(3)
    #retrieve used time 
    elapsed_time <- proc.time() - time
    print(elapsed_time[3])
}



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################################################################################
# read data
arc <- read.table("archeae_taxonomy_num.txt", header = FALSE, sep = "\t")
# read data
pro_data <- read.csv("30S12.txt", header = FALSE, sep = "\t")











