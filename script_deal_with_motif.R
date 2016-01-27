#~~~~~~~~~~~~~~~~~~~~~~~~~script get motif frome KEGG Database~~~~~~~~~~~~~~~~~~~~~~

library(KEGGAPI)
#read csv
ugt <- read.csv("112UGTs.csv", stringsAsFactors = FALSE)
#process
atName <- ugt$Locus.tag
atName <- paste0("ath:", atName)
utName <- ugt$UGT.symbol
len <- length(atName)
for(i in 1:len){
	print(i)
	#download all motif in UGTs
    motif <- getKEGGGeneMotif(atName[i])
    motif <- as.data.frame(motif)
    
    #index <- grep("pf:UDPGT", motif$Motif)
    
    #motifUDPGT <- motif[index,]
    #motifUDPGT$name <- utName[i]

    # add UGTs names
    motif$name <- utName[i]
    #write table to file
    write.table(motif, file = "112_all_motif.txt", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")	
}



#~~~~~~~~~~~~~~~~~~ process domain text to which domain in the same length ~~~~~~~~~~~~~~~~~~~~~~~~~
library(Biostrings)
# read data and precess
domain <- read.table("112_all_motif.txt", sep = '\t', stringsAsFactors = FALSE)
ugt <- read.csv("112UGTs.csv", stringsAsFactors = FALSE)
atName <- ugt$Locus.tag
atName <- paste0("ath:", atName)
utName <- ugt$UGT.symbol
#remove UGT81,82,83 because not have UDPGT motif
utName <- utName[- c(81,82,83)]

for(j in seq_along(utName)){
	#extract sub-family respectively
    eachDomain <- domain[domain$V7 %in% utName[j],]
    # get the UDPGT motif index
    index <- grep("pf:UDPGT", eachDomain$V1)
    #the number of motif
    dimRow <- 1:dim(eachDomain)[1]
    #a list contain the length of each domain
    domainLen <- sapply(dimRow, function(x){
    	eachseq <- seq(eachDomain$V2[x], eachDomain$V3[x])
    	})
    #remove the UDPGT domain
    otherRow <- dimRow[-index]
    k <- 1
    rem <- vector(mode = "numeric", length = 0)
    #check the whether other the domain with UDPGT are intersect 
    for(i in otherRow){
    
    	if(length(intersect(domainLen[[index]],domainLen[[i]])) != 0){
    		rem[k] <- i
    		k <- k + 1
    	} else{}
    	# if other domain intersect with UDPGT will be removed
    	if(length(rem) == 0){
    		eachDomain <- eachDomain

    	} else{

    			eachDomain <- eachDomain[-rem,]	
    		}
    		
    
    }
    #write the data to file
    write.table(eachDomain, file = "112_all_motif_rem.txt", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}




##########################################################################################
##########################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~deal data of domain to iTOL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#read data

library(Biostrings)

domain <- read.table("112_all_motif_rem_deal.txt", sep = '\t', stringsAsFactors = FALSE)
# get the AA numbers of each protein
aaLen <- fasta.seqlengths("process_AtUGTPROT.fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, seqtype="B", use.names=TRUE)
# get AA names
UGTName <- names(aaLen)

for(k in seq_along(UGTName)){
	# get the index of each protein
	nameAddD <- paste0(UGTName[k],"\\b")
	nameAddD <- paste0("\\b",nameAddD)
    ind <- grep(nameAddD,names(aaLen))
    ugtAALen <- aaLen[[ind]]
    # each domain of all family
    eachDomain <- domain[domain$V7 %in% UGTName[k],]
    
    dimRow <- 1:dim(eachDomain)[1]
    spr <- vector(mode = "character", length = 0)
    #assemble the Name and the AA length
    ugtLen <- paste(UGTName[k], ugtAALen,sep = ",")

    #only the one domain and more than one domain
    if(length(dimRow) == 1){
    	i <- 1
    	# get color and shape
    	att <- attr(eachDomain[i,]$V1)
    	shape <- att$Shape
    	col <- att$col
    	#get domain form and to
        form <- eachDomain[i,]$V2
        to <- eachDomain[i,]$V3
        #the motif name
        motif <- eachDomain[i,]$V1
        #assemble to result
        sp <- sprintf("%s|%s|%s|%s|%s", shape, form, to, col,motif)
    }else {
    	#more than one domain
    	for(i in dimRow){

    	    att <- attr(eachDomain[i,]$V1)
    	    shape <- att$Shape
    	    col <- att$col
            form <- eachDomain[i,]$V2
            to <- eachDomain[i,]$V3
            motif <- eachDomain[i,]$V1
            spr[i] <- sprintf("%s|%s|%s|%s|%s", shape, form, to, col,motif)
            
        }
        j <- 1
        if(j < length(spr)){
        	sp <- paste(spr[j],spr[j+1] ,sep = ",")
        	j <- j +1
        }
    }
    # add the names and length to result
    asembleSp <- paste(ugtLen, sp, sep = ",")
    #write the data
    write.table(asembleSp, file = "112_all_motif_finally.txt", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# a function get color and  shape
attr <- function(motifName){
	motif <- motifName 
	if(motif == "pf:UDPGT"){
		col <- "#FF0000"
		Shape <- "RE"
	} else if(motif == "pf:Glyco_transf_28"){
		col <- "#00FF00"
		Shape <- "HH"		

	} else if(motif == "pf:MatB"){
		col <- "#0000FF"
		Shape <- "HV"		
	} else if(motif == "pf:Glyco_trans_4_4"){
		col <- "#FFFF00"
		Shape <- "EL"

	} else if(motif == "pf:HisKA_2"){
		col <- "#00FFFF"
		Shape <- "DI"

	} else if(motif == "pf:Peptidase_M10_C"){
		col <- "#FF00FF"
		Shape <- "TR"

	} else if(motif == "pf:Glyco_trans_1_3"){
		col <- "#C0C0C0"
		Shape <- "TL"

	} else if(motif == "pf:SbmA_BacA"){
		col <- "#8A2BE2"
		Shape <- "PL"

	} else if(motif == "pf:Glyco_tran_28_C"){
		col <- "#A52A2A"
		Shape <- "PR"

	} else if(motif == "pf:YbgT_YccB"){
		col <- "#8E8E38"
		Shape <- "PU"

	} else if(motif == "pf:NSP2_assoc"){
		col <- "#7EC0EE"
		Shape <- "PD"

	} else if(motif == "pf:Tape_meas_lam_C"){
		col <- "#9F79EE"
		Shape <- "OC"

	} else if(motif == "pf:Arg_repressor_C"){
		col <- "#8B814C"
		Shape <- "GP"

	}else if(motif == "pf:MGDG_synth"){
		col <- "#008B00"
		Shape <- "OC"

	}
	return(list(col = col, Shape = Shape))
}


#~~~~~~~~~~~~~~~~~~~~~script to color label and clade~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("/home/yangfang/mywork/112UGTphylogentic/")
library("colorspace")
ugt <- read.csv("112UGTs.csv", stringsAsFactors = FALSE)
atName <- ugt$Locus.tag
utName <- ugt$UGT.symbol
num <- 71:92
num <- num[-7]
len <- length(num)
color <- rainbow_hcl(len)
for(i in 1:len){

	index <- grep(num[i], utName)
 	label <- paste(utName[index],"label", color[i], sep = "\t")
 	write.table(label, file = "label.txt", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#~~~~~~~~~~~~~~~~~~~script to Color strips~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/yangfang/mywork/112UGTphylogentic/")
library("colorspace")
ugt <- read.csv("112UGTs.csv", stringsAsFactors = FALSE)
atName <- ugt$Locus.tag
utName <- ugt$UGT.symbol
num <- 71:92
num <- num[-7]
len <- length(num)
color <- rainbow_hcl(len)
for(i in 1:len){

    index <- grep(num[i], utName)
    label <- paste(utName[index], color[i], sep = ",")
    write.table(label, file = "strips.txt", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
}