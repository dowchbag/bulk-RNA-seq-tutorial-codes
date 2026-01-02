#load libraires 
library(DESeq2)
library(tidyverse)
library(airway)
library(ggrepel)
library(tximport)
library(readr)
library(pheatmap)
library(colorspace)
library(biomaRt)
library(tximportData)
library(CompDTUReg)
library(txdbmaker)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(EnhancedVolcano)

srr_ids <- c(
  "SRR24860170", "SRR24860169", "SRR24860168", "SRR24860167",
  "SRR24860166", "SRR24860165", "SRR24860164", "SRR24860163",
  "SRR24860162", "SRR24860161", "SRR24860160", "SRR24860159", 
  "SRR24860158", "SRR24860157", "SRR24860156", "SRR24860155", 
  "SRR24860154", "SRR24860153", "SRR24860152"
)

quants <- paste0('data/', srr_ids, '_quant/quant.sf')


tx2gene <- read.table('C:/Users/dowpa/OneDrive - Nexus365/nim-lab/references/ensembl.G3Cm39.tx2genename.txt', sep = '\t', header = TRUE)
tx2gene <- tx2gene[, c(2, 1)] 
tx2gene <- tx2gene[tx2gene$Gene.name != '', ]

#create txi object
txi.salmon <- tximport(quants, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
colnames(txi.salmon$counts) <- srr_ids

#create colData 
phag <- rep(
  c("phagocytic", "non-phagocytic", "phagocytic", "non-phagocytic"),
  times = c(4, 6, 4, 5)
)
apoe <- rep(
  c("APOE3", 'APOE4'),
  times = c(10, 9)
)

colData <- data.frame(
  phag = phag,
  apoe = apoe,
  row.names = srr_ids
)

#create DeSeqDataSet
dds <- DESeqDataSetFromTximport(txi.salmon,
                                colData = colData, 
                                design = ~ phag)

#filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,] 

#set factor level 
dds$condition <- relevel(dds$phag, ref = 'non-phagocytic')

#Run DESeq 
dds <- DESeq(dds)
res <- results(dds)
res

#data transformation
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)




###FIGURE 1F
plotPCA(vsd, intgroup = c('phag', 'apoe'))
