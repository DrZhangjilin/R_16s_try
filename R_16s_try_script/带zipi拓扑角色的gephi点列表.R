setwd("G:/R/Rdata4_new/1zipi_new")
Cz <-read.delim("G:/R/Rdata4_new/1zipi_new/1zi_pi_resultC.txt", sep = '\t', check.names = FALSE)
#C <- read.table("G:/R/Rdata4_new/zipi_new/C_Gephi_allnode.csv",header = T,sep=",")
C <- read.table("G:/R/Rdata4_new/1zipi_new/node_C.csv",header = T,sep=",")
colnames(C)[1] <- "ID"
colnames(Cz)[1] <- "ID"
Cz1<-Cz[match(C$ID,Cz$ID),]
CB <- cbind(Cz1,C)
CB <- CB[,-6]
CB1 <- CB[!is.na(CB$ID),]
CB1$role <- CB1$ID
#CB1$role[which(CB1$within_module_connectivities < 2.5 & CB1$among_module_connectivities < 0.62)] <- ""
CB1$role[which(CB1$degree <= 6)] <- ""
CB1$role[which(CB1$clustering <=0.09)] <- ""
CB1$role[which(CB1$weighted.degree <= 6)] <- ""
CB1$role[which(CB1$colsenesscentrality<=0.14)] <- ""
CB1$role[which(CB1$betweenesscentrality>=0.05)] <- ""
CB1 <- CB1[-which(colnames(CB1) %in% c("timeset","modularity_class","within_module_connectivities","among_module_connectivities","degree.1"))]
write.table(CB1, 'C_node_new.csv', sep = ',', row.names = FALSE, quote = FALSE)

Cz <-read.delim("G:/R/Rdata4_new/1zipi_new/1zi_pi_resultDL.txt", sep = '\t', check.names = FALSE)
#C <- read.table("G:/R/Rdata4_new/zipi_new/DL_Gephi_allnode.csv",header = T,sep=",")
C <- read.table("G:/R/Rdata4_new/1zipi_new/node_DL.csv",header = T,sep=",")
colnames(C)[1] <- "ID"
colnames(Cz)[1] <- "ID"
Cz1<-Cz[match(C$ID,Cz$ID),]
CB <- cbind(Cz1,C)
CB <- CB[,-6]
CB1 <- CB[!is.na(CB$ID),]
CB1$role <- CB1$ID
#CB1$role[which(CB1$within_module_connectivities < 2.5 & CB1$among_module_connectivities < 0.62)] <- ""
CB1$role[which(CB1$degree <= 6)] <- ""
CB1$role[which(CB1$clustering <=0.09)] <- ""
CB1$role[which(CB1$weighted.degree <= 6)] <- ""
CB1$role[which(CB1$colsenesscentrality<=0.14)] <- ""
CB1$role[which(CB1$betweenesscentrality>=0.05)] <- ""
CB1 <- CB1[-which(colnames(CB1) %in% c("timeset","modularity_class","within_module_connectivities","among_module_connectivities","degree.1"))]
write.table(CB1, 'DL_node_new.csv', sep = ',', row.names = FALSE, quote = FALSE)

Cz <-read.delim("G:/R/Rdata4_new/1zipi_new/1zi_pi_resultNL.txt", sep = '\t', check.names = FALSE)
#C <- read.table("G:/R/Rdata4_new/zipi_new/NL_Gephi_allnode.csv",header = T,sep=",")
C <- read.table("G:/R/Rdata4_new/1zipi_new/node_NL.csv",header = T,sep=",")
colnames(C)[1] <- "ID"
colnames(Cz)[1] <- "ID"
Cz1<-Cz[match(C$ID,Cz$ID),]
CB <- cbind(Cz1,C)
CB <- CB[,-6]
CB1 <- CB[!is.na(CB$ID),]
CB1$role <- CB1$ID
#CB1$role[which(CB1$within_module_connectivities < 2.5 & CB1$among_module_connectivities < 0.62)] <- ""
CB1$role[which(CB1$degree <= 6)] <- ""
CB1$role[which(CB1$clustering <=0.09)] <- ""
CB1$role[which(CB1$weighted.degree <= 6)] <- ""
CB1$role[which(CB1$colsenesscentrality<=0.14)] <- ""
CB1$role[which(CB1$betweenesscentrality>=0.05)] <- ""
CB1 <- CB1[-which(colnames(CB1) %in% c("timeset","modularity_class","within_module_connectivities","among_module_connectivities","degree.1"))]
write.table(CB1, 'NL_node_new.csv', sep = ',', row.names = FALSE, quote = FALSE)

Cz <-read.delim("G:/R/Rdata4_new/1zipi_new/1zi_pi_resultNP.txt", sep = '\t', check.names = FALSE)
#C <- read.table("G:/R/Rdata4_new/zipi_new/NP_Gephi_allnode.csv",header = T,sep=",")
C <- read.table("G:/R/Rdata4_new/1zipi_new/node_NP.csv",header = T,sep=",")
colnames(C)[1] <- "ID"
colnames(Cz)[1] <- "ID"
Cz1<-Cz[match(C$ID,Cz$ID),]
CB <- cbind(Cz1,C)
CB <- CB[,-6]
CB1 <- CB[!is.na(CB$ID),]
CB1$role <- CB1$ID
#CB1$role[which(CB1$within_module_connectivities < 2.5 & CB1$among_module_connectivities < 0.62)] <- ""
CB1$role[which(CB1$degree <= 6)] <- ""
CB1$role[which(CB1$clustering <=0.09)] <- ""
CB1$role[which(CB1$weighted.degree <= 6)] <- ""
CB1$role[which(CB1$colsenesscentrality<=0.14)] <- ""
CB1$role[which(CB1$betweenesscentrality>=0.05)] <- ""
CB1 <- CB1[-which(colnames(CB1) %in% c("timeset","modularity_class","within_module_connectivities","among_module_connectivities","degree.1"))]
write.table(CB1, 'NP_node_new.csv', sep = ',', row.names = FALSE, quote = FALSE)

Cz <-read.delim("G:/R/Rdata4_new/1zipi_new/1zi_pi_resultNPDL.txt", sep = '\t', check.names = FALSE)
#C <- read.table("G:/R/Rdata4_new/zipi_new/NPDL_Gephi_allnode.csv",header = T,sep=",")
C <- read.table("G:/R/Rdata4_new/1zipi_new/node_NPDL.csv",header = T,sep=",")
colnames(C)[1] <- "ID"
colnames(Cz)[1] <- "ID"
Cz1<-Cz[match(C$ID,Cz$ID),]
CB <- cbind(Cz1,C)
CB <- CB[,-6]
CB1 <- CB[!is.na(CB$ID),]
CB1$role <- CB1$ID
#CB1$role[which(CB1$within_module_connectivities < 2.5 & CB1$among_module_connectivities < 0.62)] <- ""
CB1$role[which(CB1$degree <= 6)] <- ""
CB1$role[which(CB1$clustering <=0.09)] <- ""
CB1$role[which(CB1$weighted.degree <= 6)] <- ""
CB1$role[which(CB1$colsenesscentrality<=0.14)] <- ""
CB1$role[which(CB1$betweenesscentrality>=0.05)] <- ""
CB1 <- CB1[-which(colnames(CB1) %in% c("timeset","modularity_class","within_module_connectivities","among_module_connectivities","degree.1"))]
write.table(CB1, 'NPDL_node_new.csv', sep = ',', row.names = FALSE, quote = FALSE)

Cz <-read.delim("G:/R/Rdata4_new/1zipi_new/1zi_pi_resultNPNL.txt", sep = '\t', check.names = FALSE)
#C <- read.table("G:/R/Rdata4_new/zipi_new/NPNL_Gephi_allnode.csv",header = T,sep=",")
C <- read.table("G:/R/Rdata4_new/1zipi_new/node_NPNL.csv",header = T,sep=",")
colnames(C)[1] <- "ID"
colnames(Cz)[1] <- "ID"
Cz1<-Cz[match(C$ID,Cz$ID),]
CB <- cbind(Cz1,C)
CB <- CB[,-6]
CB1 <- CB[!is.na(CB$ID),]
CB1$role <- CB1$ID
#CB1$role[which(CB1$within_module_connectivities < 2.5 & CB1$among_module_connectivities < 0.62)] <- ""
CB1$role[which(CB1$degree <= 6)] <- ""
CB1$role[which(CB1$clustering <=0.09)] <- ""
CB1$role[which(CB1$weighted.degree <= 6)] <- ""
CB1$role[which(CB1$colsenesscentrality<=0.14)] <- ""
CB1$role[which(CB1$betweenesscentrality>=0.05)] <- ""
CB1 <- CB1[-which(colnames(CB1) %in% c("timeset","modularity_class","within_module_connectivities","among_module_connectivities","degree.1"))]
write.table(CB1, 'NPNL_node_new.csv', sep = ',', row.names = FALSE, quote = FALSE)
