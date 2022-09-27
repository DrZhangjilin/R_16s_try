setwd("G:/R/Rdata4_new/1zipi_new")
C1 <- read.table("G:/R/Rdata4_new/1zipi_new/C_node_tax.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/C_node_new.csv",header = T,row.names=1,sep=",")
dat <- merge(C1,C2,by='row.names')#按照行名合并文件
dat1 <- dat[,c(1,8:14,31)]
write.table(dat1, 'C_new_core.csv', sep = ',', row.names = FALSE, quote = FALSE)

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/DL_node_tax.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/DL_node_new.csv",header = T,row.names=1,sep=",")
dat <- merge(C1,C2,by='row.names')#按照行名合并文件
dat1 <- dat[,c(1,8:14,31)]
write.table(dat1, 'DL_new_core.csv', sep = ',', row.names = FALSE, quote = FALSE)

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/NL_node_tax.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/NL_node_new.csv",header = T,row.names=1,sep=",")
dat <- merge(C1,C2,by='row.names')#按照行名合并文件
dat1 <- dat[,c(1,8:14,31)]
write.table(dat1, 'NL_new_core.csv', sep = ',', row.names = FALSE, quote = FALSE)

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/NP_node_tax.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/NP_node_new.csv",header = T,row.names=1,sep=",")
dat <- merge(C1,C2,by='row.names')#按照行名合并文件
dat1 <- dat[,c(1,8:14,31)]
write.table(dat1, 'NP_new_core.csv', sep = ',', row.names = FALSE, quote = FALSE)

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/NPDL_node_tax.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/NPDL_node_new.csv",header = T,row.names=1,sep=",")
dat <- merge(C1,C2,by='row.names')#按照行名合并文件
dat1 <- dat[,c(1,8:14,31)]
write.table(dat1, 'NPDL_new_core.csv', sep = ',', row.names = FALSE, quote = FALSE)

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/NPNL_node_tax.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/NPNL_node_new.csv",header = T,row.names=1,sep=",")
dat <- merge(C1,C2,by='row.names')#按照行名合并文件
dat1 <- dat[,c(1,8:14,31)]
write.table(dat1, 'NPNL_new_core.csv', sep = ',', row.names = FALSE, quote = FALSE)