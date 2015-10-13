#####################API for NCBI Pubmed search paper#####################################
library(XML)
library(RCurl)
article_search <- function(art_name,pub_date){
	#the base NCBI API EutilsURL
    base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # set the protein database
    db <- "pubmed"
    retmax <- "200"
    #assemble the esearch URL
    query <- sprintf("%s [title] AND %s [pdat]", art_name, pub_date)
    esearch_url <- sprintf("esearch.fcgi?db=%s&term=%s&retmax=%s", db, query, retmax)
    esearch_url <- paste0(base_url, esearch_url)

    search_doc <- xmlParse(esearch_url)
    #retrieve the UIDs
    uids <- xpathSApply(search_doc, path = "//IdList/Id", fun = 'xmlValue')
    #id=4564646,646466,646466,........
    ids <- paste(uids, collapse = ",")
    #assemble the esummary URL
    esummary <- sprintf("esummary.fcgi?db=%s&id=%s",db, ids)

    esum_url <- paste0(base_url, esummary)
    #parse the URL
    esum_doc <- xmlParse(esum_url)
    #extract title date and source
    paper_names <- xpathSApply(esum_doc,"//Item[@Name='Title']", xmlValue)
    paper_date <- xpathSApply(esum_doc,"//Item[@Name='EPubDate']", xmlValue)
    paper_Source <- xpathSApply(esum_doc,"//Item[@Name='Source']", xmlValue)
    #get complete information of aritcles
    article <- paste0(paper_names, paper_date)
    article <- paste(article, paper_Source, sep = ",")
    return(article)    
}

#it's a test
art_name <- "epoxy resin"
pub_date <- "2015"
art <- article_search(art_name, pub_date)
write.csv(art, file = "epoxy_resin.csv")



#######################################################################################
