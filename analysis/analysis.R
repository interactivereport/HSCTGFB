library(reshape2)
library(FunPat)
library(outliers)

data = read.table("../Summary/Gene-fpkm-table.flt.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta = read.table("sample.annotation.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta = meta [, c("Tissue", "Time", "Compound")]

exp = as.matrix(t(data[,8:dim(data)[2]])) # dim(data)[2] - number of columns

#exp = as.matrix(t(data[1:3000,8:dim(data)[2]]))

exp1 = exp[!(row.names(exp) %in% c("5_S6", "7_S7")),]

exp1 = log(exp1+1e-4)

exp2 = merge(meta,exp1, by="row.names",all.x=TRUE)

meanNoOutliers <- function (x) {
	val = mean(rm.outlier(x, fill = FALSE, median = FALSE, opposite = FALSE))
	return(val)
}

blood_PFE = exp2[exp2$Tissue=='Blood' & exp2$Compound!="SDD",]
blood_PFE_avg <- aggregate(blood_PFE[5:dim(blood_PFE)[2]], list(Time=blood_PFE$Time), FUN = mean)
blood_PFE_t <- t(blood_PFE_avg[,2:ncol(blood_PFE_avg)])
colnames(blood_PFE_t) <- blood_PFE_avg[,1]


blood_SDD = exp2[exp2$Tissue=='Blood' & exp2$Compound!="PF-06826647",]
blood_SDD_avg <- aggregate(blood_SDD[5:dim(blood_SDD)[2]], list(Time=blood_SDD$Time), FUN = mean)
blood_SDD_t <- t(blood_SDD_avg[,2:ncol(blood_SDD_avg)])
colnames(blood_SDD_t) <- blood_SDD_avg[,1]


rep = transform(merge(blood_PFE_t[,2], blood_SDD_t[,2],by=0,all.x=TRUE), row.names=Row.names, Row.names=NULL)

rank.res<-SEL.TS.AREA(replicates=rep,data1=blood_PFE_t,data2=blood_SDD_t, htmlreport=TRUE,link_to_Entrez=NA,kplot=100)


# Spleen
spleen_PFE = exp2[exp2$Tissue=='Spleen' & exp2$Compound!="SDD",]
spleen_PFE_avg <- aggregate(spleen_PFE[5:dim(spleen_PFE)[2]], list(Time=spleen_PFE$Time), FUN = meanNoOutliers)
spleen_PFE_t <- t(spleen_PFE_avg[,2:ncol(spleen_PFE_avg)])
colnames(spleen_PFE_t) <- spleen_PFE_avg[,1]


spleen_SDD = exp2[exp2$Tissue=='Spleen' & exp2$Compound!="PF-06826647",]
spleen_SDD_avg <- aggregate(spleen_SDD[5:dim(spleen_SDD)[2]], list(Time=spleen_SDD$Time), FUN = meanNoOutliers)
spleen_SDD_t <- t(spleen_SDD_avg[,2:ncol(spleen_SDD_avg)])
colnames(spleen_SDD_t) <- spleen_SDD_avg[,1]


rep = transform(merge(spleen_PFE_t[,2], spleen_SDD_t[,2],by=0,all.x=TRUE), row.names=Row.names, Row.names=NULL)

rank.res<-SEL.TS.AREA(replicates=rep,data1=spleen_PFE_t,data2=spleen_SDD_t, htmlreport=TRUE,link_to_Entrez=NA,kplot=100)

