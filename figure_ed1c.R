phag_srr_ids <- c("SRR24860170", "SRR24860169", "SRR24860168", "SRR24860167",
                  "SRR24860160", "SRR24860159", "SRR24860158", "SRR24860157"
)

phag_quants <- paste0('data/', phag_srr_ids, '_quant/quant.sf')

#create txi object
txi.salmon.phag <- tximport(phag_quants, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
colnames(txi.salmon.phag$counts) <- phag_srr_ids

phag_colData <- data.frame(
  apoe = rep(c("APOE3", "APOE4"), times = c(4,4)),
  row.names = phag_srr_ids
)

#create DeSeqDataSet
phag_dds <- DESeqDataSetFromTximport(txi.salmon.phag,
                                colData = phag_colData, 
                                design = ~ apoe)

#filter
keep <- rowSums(counts(phag_dds)) >= 10
phag_dds <- phag_dds[keep,] 

#set factor level 
phag_dds$condition <- relevel(phag_dds$apoe, ref = 'APOE3')

#Run DESeq 
phag_dds <- DESeq(phag_dds)
phag_res <- results(phag_dds, alpha = 0.05)
phag_res <- lfcShrink(phag_dds, coef="apoe_APOE4_vs_APOE3", type="normal", res=phag_res)



####MANUALLY GENERATING VOLCANO PLOT USING GGPLOT2 

de <- as.data.frame(phag_res)
de <- de[complete.cases(de), ] 
de$minuslog10P <- -log10(de$padj)
View(de)

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.5 & de$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.5 & de$pvalue < 0.05] <- "DOWN"

#define labelled genes 
label_genes <- c(
  "Mid1", "Npy", "Ch25h", "Il1rn", "H2-Q7", "Ccrl2",
  "Lgals3", "Lilrb4b", "Atf3", "Ccl5", "Cxcl10", "Clec7a",
  "Grn", "Irf7", "Ccl9", "B2m",
  "Rassf4", "Asns", "Depp1", "Lgi4", "Ccnd2"
)

de$delabel <- NA 
de$delabel[rownames(de) %in% label_genes] <- rownames(de)[rownames(de) %in% label_genes]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  scale_color_manual(values=c("blue", "black", "red"))





#### AUTOMATICALLY GENERATING A PLOT USING ENHANCEDVOLCANO

EnhancedVolcano(phag_res,
                lab = rownames(phag_res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                pCutoff = 0.05, 
                FCcutoff = 0.5, 
                selectLab = label_genes, 
                drawConnectors = TRUE)





