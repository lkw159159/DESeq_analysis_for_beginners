#1-1
BiocManager::install("DESeq2")
BiocManager::install("tidyverse")
BiocManager::install("openxlsx")

library(DESeq2); library(tidyverse); library(openxlsx)

#2-1
## Data load
if (!require('zebrafishRNASeq', quietly = TRUE)){
  BiocManager::install("zebrafishRNASeq")
  library(zebrafishRNASeq)}else{
    library(zebrafishRNASeq)
  }
data(zfGenes)

dim(zfGenes);head(zfGenes)

#2-2
## Generation of metadata
metadata_mat = data.frame(Sample = colnames(zfGenes),
                    Condition = c(rep('affected',3), rep('unaffected',3)),
                    Batch = c('A','A','B','B','B','A'),
                    Seq.date = rep('22.11.15',6),
                    Seq.platform = rep('HiSeq4000',6))

#2-3
## Make a DeseqDataset from matrix
dds_mat = DESeqDataSetFromMatrix(countData = zfGenes,
                                 colData = metadata_mat,
                                 design = ~ Condition)

## Pre-filtering
dds_mat = dds_mat[apply(counts(dds_mat),1,function(x) sum(x == 0)) < ncol(dds_mat) / 3,]
dim(dds_mat)

## specify the reference level
dds_mat$Condition = relevel(dds_mat$Condition,ref = 'unaffected')
dds_mat$Condition

## Run DESeq2 
dds_mat = DESeq(dds_mat)

## Show the statistical characteristics
res_mat = results(dds_mat, alpha=0.05)
summary(res_mat)
vsd_mat=vst(dds_mat,blind=F)

#2-4
## Write the result 
out_mat = cbind(counts(dds_mat,normalized=T),res_mat)
write.xlsx(out_mat,'DESeq2 result from matrix.xlsx')


#3-1
## Make an transcriptome database from used GFF file
library(tximport); library(readr); library(GenomicFeatures)
txdb = makeTxDbFromGFF("Homo_sapiens.GRCh38.106.chr.gtf.gz")
k=keys(txdb, keytype = 'TXNAME')
tx2gene= select(txdb,k,"GENEID",'TXNAME')

#3-2
## load the metadata
metadata_txi=read.xlsx('metadata.xlsx')

#3-3
## Tximport with result files of SALMON (quant.sf)
txi=tximport(metadata_txi$quant_file,type="salmon",tx2gene=tx2gene,countsFromAbundance ="scaledTPM",ignoreTxVersion = T)

colnames(txi$counts)=metaData$Sample
head(txi$counts)

## Make a DeseqDataset from matrix
dds_mat = DESeqDataSetFromTximport(txi = txi,
                                   colData = metadata_txi,
                                   design = ~ Condition)




#4-1
## Make some plots to analyze
library(ggplot2); library(grid); library(gridBase); library(ggrepel); library(calibrate)

x11()             ## Open plot window
par(mfrow=c(2,2)) ## Divide the plot window by 4

### MA plot
plotMA(res_mat)

### Count comparison plot
plotCounts(dds_mat,gene=which.min(res_mat$padj),intgroup='Condition') 

### Volcano plot 
with(res_mat, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-10,10),col='gray'))
with(subset(res_mat, padj<0.05 ), points(log2FoldChange, -log10(padj), pch=20,))
with(subset(res_mat, padj<0.05 & log2FoldChange>=1.0),points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(res_mat, padj<0.05 & log2FoldChange<=-1.0),points(log2FoldChange, -log10(padj), pch=20, col="blue"))
legend("bottomleft",legend=c(paste0("Up-regulated: ",length(rownames(subset(res_mat, padj<0.05 & log2FoldChange>=1.0)))),
                             paste0("Down-regulated: ",length(rownames(subset(res_mat, padj<0.05 & log2FoldChange<=-1.0))))),
       col=c("red","blue"),pch=19,border="white",box.lty=0,cex=1)

### PCA plot
p1 = plotPCA(vsd_mat,intgroup=c('Condition','Batch'))

####do not change below lines, just run them
plot.new()
vps = baseViewports()
pushViewport(vps$figure)
vp1 = plotViewport(c(1.8,1,0,1))
print(p1,vp=vp1) 
#################################



