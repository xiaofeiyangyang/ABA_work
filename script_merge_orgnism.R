

############################process bacteria data##################################

##~~~~~~~~~~~~~~~~~~~~~read tabel and combine results~~~~~~~~~~~~~~~~~~~~~~~~~~
table <- read.csv("bacteria.csv", head = FALSE)

organ <- sapply(table, extract_organ)

results <- cbind(organ[[1]], organ[[2]])

write.csv(results, file = "output_bacteria.csv", row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



##~~~~~~~~~~~~~~~~~~~~~~merge two table get the same organisms~~~~~~~~~~~~~~~~~~
x <- read.csv("full_name.csv", head=FALSE)

y <- read.csv("output_bacteria.csv", head=TRUE)

mer <- merge(x, y,  all.x = T)



















