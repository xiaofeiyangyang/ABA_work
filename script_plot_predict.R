library(ggplot2)
setwd("/home/yangfang/clime/ABA/eukaryotes3/output_UGT2/matrix")
num <- 1
files <- dir()
for(i in seq_along(files)){

	ugt71 <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    ugt71s <- as.character(ugt71$Gene.Symbol[ugt71$LLR > 20])
    num[i] <- length(ugt71s)
}

name <- paste0("UGT", 71:92)
name <- name[-7]
name <- paste0(name, "sub-family")
predict <- data.frame(name, num)


p <- ggplot(data=predict, aes(x=name, y=num))
p <- p + geom_bar( stat="identity", width=0.4, fill="cornflowerblue")
p <- p + geom_text(label=predict$num, colour="blue", vjust=-1, size=3)
p <- p + labs(x="Arabidopsis thaliana UGTs Predict LLR > 20 ", y="Number of predict genes")
p <- p + theme(axis.text.x=element_text(family="myFont2",face="bold",size=8,angle=45,))
pdf("predict genes")
p
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p <- ggplot() + 
geom_bar(aes(x=name, y=num, angle=45), stat="identity", fill="cornflowerblue", width=0.5) + 
theme(axis.text.x=element_text(family="myFont2",face="bold",size=8,angle=45,)) +

p <- p + geom_text(label=predict$num, colour="blue", vjust=-1)

library(ggplot2)
city<-c( '北京', '上海', '天津')
rate<-c(0.53, 0.29, 0.18)
data1<-data.frame(city,rate)
colnames(data1)<-c('city','rate')

p <- ggplot(data = data1, aes(x=city,y=rate))  
p <- p + geom_bar( stat="identity" , width = 0.4, fill = "cornflowerblue") 
p <- p + geom_text(label=paste(data1$rate * 100, "%", sep = "") ,colour = "blue", vjust=-1) 
p <- p + labs(x="",y="份额\n",title = "各省份额\n")  
p <- p + scale_y_continuous(limits=c(0, max(data1$rate)*1.1),labels = percent, breaks = seq(0, 2, 0.1)) 
p <- p + theme( plot.title = element_text(size = 16, face = "bold"))

p