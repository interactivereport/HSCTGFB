# remove sample "5" (wrong label), "7" (outlier)

# FunPat - prone to outliers
#Example C shows an experimental design where there are two replicates for all the time points. In this
#case we suggest using the mean expression in data1, data2 and replicates since the average provides more
#robust expression values to test, thus a better estimation of the p-values for the gene ranking. As seen for
#example A, for replicates the mean expression in T0 and C0 can be used. It is worth noting that, since
#replicates is used to describe the variability of the diㄦences observed between data1 and data2, whose
#expression values are averaged across the replicates, it is important that the error variance estimated from
#replicates reporting the averaged expression values of the replicates, otherwise the estimated variance tend
#to be over-estimated. The over-estimation of the variance does not act the gene ranking, but leads to
#higher p-values thus a lower number of selected seeds.
#
#
#replic<-Simdata$replicates
#nC<-Simdata$ctrl
#nD<-Simdata$treat
#rank.res<-SEL.TS.AREA(replicates=replic,data1=nC,data2=nD)
#
#source("http://bioconductor.org/biocLite.R")
#biocLite("ReportingTools")
#biocLite("hwriter")
#biocLite("Rgraphviz")
#library(ReportingTools)
#library(hwriter)
#library(Rgraphviz)


# filter out lowly expressed genes
./filterRPKM.pl ../Summary/Gene-fpkm-table.txt  1  9 Spleen > Gene-fpkm-table.flt.Spleen.txt

awk 'BEGIN{FS="\t";OFS="\t"} FNR == NR {hash[$1]=$1; next} $1 in hash {print $0}' Gene-fpkm-table.flt.Spleen.txt ../Summary/Gene-count-table.txt  > Gene-count-table.flt.Spleen.txt

./filterRPKM.pl ../Summary/Gene-fpkm-table.txt  1  9 Blood > Gene-fpkm-table.flt.Blood.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR == NR {hash[$1]=$1; next} $1 in hash {print $0}' Gene-fpkm-table.flt.Blood.txt  ../Summary/Gene-count-table.txt  > Gene-count-table.flt.Blood.txt

# Differential Gene Expression
Rscript edgeR.Blood.R
Rscript edgeR.Spleen.R

/ngs-pubdata/home/baohongz/bin/vennPlot.pl Spleen.comp.txt 0.05 0 /ngs-pubdata/data/annotation/gene/Ensembl_mm10_gene.txt > venn/Spleen.comp.csv


# gene filter
./filterRPKM.pl ../Summary/Gene-fpkm-table.txt 2  9 HSC | awk 'BEGIN{FS="\t";OFS="\t"} FNR == NR {hash[$1]=$1; next} $1 in hash {print $0}' - ../Summary/Gene-count-table.txt  > Gene-count-table.flt.HSC.txt

Rscript edgeR.section.R


/ngs-pubdata/home/baohongz/bin/vennPlot.pl t48_t72.comp.txt 0.05 1 /ngs-pubdata/data/annotation/gene/Ensembl_hg38_gene.txt > venn/t48_t72.comp.csv


# generate the DE gene list 
/ngs-pubdata/home/baohongz/QuickRNASeq1.2/comparativeList.pl --comparison_f comp1.txt --FDR_cutoff 0.0001 --logFC_cutoff 1



awk 'BEGIN{FS="\t";OFS="\t"} {if ($9<=0.0001 && ($5>= 1 || $5<=-1)) print $0}' TGFb_t24.txt | cut -f 1 > TGFb_t24.DEgene.txt
awk 'BEGIN{FS="\t";OFS="\t"} {if ($9>0.05) print $0}' TGFb_and_JQ1_t24.txt  | cut -f 1 > TGFb_and_JQ1_t24.nonDEgene.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR == NR {hash[$1]=1; next} $1 in hash {print $1}' TGFb_t24.DEgene.txt TGFb_and_JQ1_t24.nonDEgene.txt > TGFb_t24.rev.txt


awk 'BEGIN{FS="\t";OFS="\t"} {if ($9<=0.0001 && ($5>= 1 || $5<=-1)) print $0}' TGFb_t48.txt | cut -f 1 > TGFb_t48.DEgene.txt
awk 'BEGIN{FS="\t";OFS="\t"} {if ($9>0.05) print $0}' TGFb_and_JQ1_t48.txt  | cut -f 1 > TGFb_and_JQ1_t48.nonDEgene.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR == NR {hash[$1]=1; next} $1 in hash {print $1}' TGFb_t48.DEgene.txt TGFb_and_JQ1_t48.nonDEgene.txt > TGFb_t48.rev.txt


awk 'BEGIN{FS="\t";OFS="\t"} {if ($9<=0.0001 && ($5>= 1 || $5<=-1)) print $0}' TGFb_t72.txt | cut -f 1 > TGFb_t72.DEgene.txt
awk 'BEGIN{FS="\t";OFS="\t"} {if ($9>0.05) print $0}' TGFb_and_JQ1_t72.txt  | cut -f 1 > TGFb_and_JQ1_t72.nonDEgene.txt
awk 'BEGIN{FS="\t";OFS="\t"} FNR == NR {hash[$1]=1; next} $1 in hash {print $1}' TGFb_t72.DEgene.txt TGFb_and_JQ1_t72.nonDEgene.txt > TGFb_t72.rev.txt


cat *.rev.txt | sort | uniq > reverse.txt


/ngs-pubdata/home/baohongz/QuickRNASeq1.2/comparativeList.pl --comparison_f comp2.txt --geneset_f=reverse.txt 

# create Venn Diagram
/ngs-pubdata/home/baohongz/QuickRNASeq2/vennPlot.pl --comparison_f compall.txt --annotation_f /ngs-pubdata/data/annotation/gene/Ensembl_v81_GRCh38.p3_Gencode.v23_human.txt.gz



# collate the lists including reversals (not differentially expressed)
awk 'BEGIN{FS="\t";OFS="\t"} {if ($5>0) $5=9999} {if ($5<=0) $5=-9999} {if ($9>0.05) print $1,$2,$3,$4,$5,$6,$7,$8,0}' TGFb_and_JQ1_t24.txt > TGFb_and_JQ1_t24_nonDE.txt
awk 'BEGIN{FS="\t";OFS="\t"} {if ($5>0) $5=9999} {if ($5<=0) $5=-9999} {if ($9>0.05) print $1,$2,$3,$4,$5,$6,$7,$8,0}' TGFb_and_JQ1_t48.txt > TGFb_and_JQ1_t48_nonDE.txt
awk 'BEGIN{FS="\t";OFS="\t"} {if ($5>0) $5=9999} {if ($5<=0) $5=-9999} {if ($9>0.05) print $1,$2,$3,$4,$5,$6,$7,$8,0}' TGFb_and_JQ1_t72.txt > TGFb_and_JQ1_t72_nonDE.txt

/ngs-pubdata/home/baohongz/QuickRNASeq2/vennPlot.pl --comparison_f compall.txt --annotation_f /ngs-pubdata/data/annotation/gene/Ensembl_v81_GRCh38.p3_Gencode.v23_human.txt.gz
