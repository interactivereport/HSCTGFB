library(edgeR)
library(plyr)

rawCount = read.table(file="Gene-count-table.flt.HSC.txt",header=TRUE,sep="\t",check.names=FALSE)
#rownames(rawCount) <- gsub("\\.\\d+", "",  rawCount$gene, perl=TRUE)
geneAnno = read.table(file="/ngs-pubdata/data/annotation/gene/Ensembl_hg38_gene.txt", header=TRUE,sep="\t",quote="\"")
rownames(geneAnno) <- geneAnno$gene_id
rawCount$gene = sub ("\\.[0-9]+","", rawCount$gene, perl=TRUE)


Tissue = "HSC"
meta = read.table("sample.annotation.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta = meta [meta$Tissue==Tissue, c("Tissue", "timepoint", "treatment")]
meta$sample_id = rownames(meta)


for (timepoint in c("24","48","72")) {

	print (timepoint)

	for (treatment in c("TGFb", "TGFb_and_JQ1", "Thrombin")) {
		targets = meta[meta$timepoint==timepoint&(meta$treatment==treatment|meta$treatment=="untreated"),]
	
	# levels setting is very important to have the contrast right
		mytreatment <- factor(targets$treatment, levels=c("untreated", treatment))
	
		design <- model.matrix(~mytreatment)
		
		count <- (rawCount[c(as.character(targets$sample_id))])
		rownames(count) <- rawCount$gene
		
		e <- DGEList(count=count)
		e <- calcNormFactors(e, "TMM")
		
		e <- estimateGLMCommonDisp(e, design)
		e <- estimateGLMTrendedDisp(e, design)
		e <- estimateGLMTagwiseDisp(e, design)
		
		fit <- glmFit(e, design)
		contrast <- paste("mytreatment",treatment, sep="");
		lrt <- glmLRT(fit, coef=contrast)
		
		topTable <- topTags(lrt,n=Inf,adjust.method="BH")$table
	
		mtable <-merge (topTable, geneAnno,  by="row.names", all.x=TRUE, sort=FALSE)
		mtable$gene_id <- mtable$Row.names
		mtable <- mtable[c("gene_id","gene_name","description","biotype","logFC", "logCPM", "LR", "PValue", "FDR")]
		
		write.table(format(mtable, digits=4), file=paste(treatment,"_t",timepoint,".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
	}
}
