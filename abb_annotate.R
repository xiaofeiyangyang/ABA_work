setwd("/home/yangfang/mywork/abb_annotation/")

abb <- read.table("270euk_abb_names.txt",header = FALSE, sep ="\t",stringsAsFactors = FALSE)

full <- read.csv("wKEGGSpe.csv",header = TRUE, stringsAsFactors = FALSE)


logs <- full$abbName %in% abb$V1

merge <- full[logs,]
