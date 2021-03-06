#######################################################################
## a script to retrieve one optimization KO number
library(KEGGBioCycAPI)
#read data
abb_names <- read.csv("276abb_names.csv",head = FALSE)
abb_names <- as.character(abb_names[[1]])
require('Biostrings')
len <- length(abb_names)
#get KO list
pro_kos <- ("K01409")
kegg_list <- getKEGGKO(pro_kos)
#intilization
ko_num <- vector(mode = "character", length = 0)
miss_names <- vector(mode = "character", length = 0)
k <- 1
#kegg_list <- strsplit(kegg_list, ":")
#kegg_list_abb <- sapply(kegg_list, function(v) return(v[1]))
for(i in 1: len){
    #assemble regular expression
    abb <- paste0(abb_names[i],":")
    abb <- paste0("^",abb )
    pro_ko <- kegg_list[grep(abb, kegg_list)]
    #print(pro_ko)
    if(length(pro_ko) == 0){
    	ko_num[i] <- NA
        miss_names[k] <- abb_names[i]
        k <- k+1
    }else{

    	ko_num[i] <- pro_ko[1]

    }
   ko_num <- ko_num[!is.na(ko_num)]
}
length(ko_num)
print(miss_names)




pro_seqs <- getKEGGGeneSeq(ko_num, n = 2)
#read to local disk
writeXStringSet(pro_seqs, "polymerase_subunit_alpha.fasta")
write.table(miss_names, file = "polymerase_subunit_alpha_miss.txt", quote = FALSE )




 ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###################################################################
###################################################################
#a function to download protein with two KO number(one Optimization 
#other is Alternative)
library(KEGGBioCycAPI)
Get_pro_kegg <- function(Optimization, Alternative , AbbNames){
    #initialize in put data
    pro_kos <- Optimization
    pro_ko_alternative <- Alternative
    abb_names <- AbbNames
  
    len <- length(abb_names)
    #get kegg list
    kegg_list <- getKEGGKO(pro_kos)
    kegg_list_alternative <- getKEGGKO(pro_ko_alternative)
    #initialize vector
    ko_num <- vector(mode = "character", length = 0)
    miss_names <- vector(mode = "character", length = 0)
    k <- 1
    #match kegg ko number
    for(i in 1: len){
        #assemble regular expression and match
        abb <- paste0(abb_names[i],":")
        abb <- paste0("^",abb )
        pro_ko <- kegg_list[grep(abb, kegg_list)]
        #dispose the miss match by alternative list
        if(length(pro_ko) == 0){
        	abb_alt <- paste0(abb_names[i],":")
            abb_alt <- paste0("^",abb_alt )
        	pro_ko_alt <- kegg_list_alternative[grep(abb_alt, kegg_list_alternative)]
        	# get also can't match organism
            if(length(pro_ko_alt) == 0){
            	ko_num[i] <- NA
            	miss_names[k] <- abb_names[i]
                k <- k+1
            }else{
               #get the hit organism by alternative kegg list
                ko_num[i] <- pro_ko_alt[1]
    
            }
        }else{
            #get the hit organism by optimization kegg list 
        	ko_num[i] <- pro_ko[1]
    
        }
        #dislodge the NA
        ko_num <- ko_num[!is.na(ko_num)]
    }
    #retrun a list content include hit KO number and the miss hit organism names
    return(list(miss_names = miss_names, ko_num = ko_num))
}

#read the data and processed
abb_names <- read.csv("276abb_names.csv",head = FALSE)
abb_names <- as.character(abb_names[[1]])
require('Biostrings')

kegg_num <- Get_pro_kegg("K02863 ","K02865",abb_names)
ko_split <- kegg_num$ko_num
AAseqs <- getKEGGGeneSeq(ko_split, n = 2)

writeXStringSet(AAseqs, "euk_ribosomal_protein_L1.fasta")
write.table(kegg_num$miss_names, file = "euk_ribosomal_protein_L1_miss.txt", quote = FALSE )
	
#####################################################################
#a script deal with too many mitochondrial protein by we get from kegg
library(stringr)
library(KEGGBioCycAPI)
#read data and process
require('Biostrings')
fullseq <- readBStringSet("euk_synthetase_ser.fasta")
seqs <- readLines("euk_synthetase_ser.fasta")
#extract abb names with mitochondrial
mito <- seqs[grep("mitochondrial",seqs)]
str_extract_all(mito[1], "(?<=\\>)(.*?)(?=:)")

mito_names <- sapply(mito,function(x){
	each_name <- str_extract_all(x, "(?<=\\>)(.*?)(?=:)")
	return(each_name)
})
mito_names <- as.character(mito_names)
#download mitochondrial abb proteins
abb_names <- mito_names
len <- length(abb_names)
#get KO list
pro_kos <- ("K01875")
kegg_list <- getKEGGKO(pro_kos)
#intilization
ko_num <- list()
miss_names <- vector(mode = "character", length = 0)
k <- 1
k <- 1
#kegg_list <- strsplit(kegg_list, ":")
#kegg_list_abb <- sapply(kegg_list, function(v) return(v[1]))
for(i in 1: len){
    #assemble regular expression
    abb <- paste0(abb_names[i],":")
    abb <- paste0("^",abb )
    pro_ko <- kegg_list[grep(abb, kegg_list)]
   
    if(length(pro_ko) == 0){
    	ko_num[[i]] <- NA
        miss_names[k] <- abb_names[i]
        k <- k+1
    }else{

    	ko_num[[i]] <- pro_ko

    }
   ko_num <- ko_num[!is.na(ko_num)]
}

ko_nums <- unlist(ko_num)
#remove the mitochondrial proteins frome download
pro_seqs <- getKEGGGeneSeq(ko_nums, n = 2)
gre <- grep("mitochondrial",names(pro_seqs))
#should manual check length  seqs
seqs <- pro_seqs[-gre]
seqs <- seqs[-c(2,4,22,35)]

#remove mitochondrial proteins forme local sequences
gre2 <- grep("mitochondrial",names(fullseq))

rm_mito_full <- fullseq[-gre2]
#merge to protein sequnces
mer <- c(seqs,rm_mito_full)

#write to file
writeXStringSet(mer, "euk_synthetase_ser_remito.fasta")
#write.table(miss_names, file = "euk_synthetase_leu_miss.txt", quote = FALSE )
