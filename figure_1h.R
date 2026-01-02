###FIGURE 1H 
#get expression matrix
expr <- assay(vsd)
expr <- expr[ , colData$phag == 'phagocytic']

#define genes for each phenotype
phag_genes <- c("Fcgr1","Thbs1","Rab7b","Sod1","Washc5","Snx3","Itgb2")
autophagy_genes <- c("Gramd1a","Smcr8","Ubqln2","Tollip","Trim12a","Lamp2","Becn1","Irgm1")
ifny_genes <- c("Kif16b","Acod1","Ccl9","Ccl4","Ccl3","Capg","H2-Q7","Gbp2","Irgm1","Jak2")
antigen_genes <- c("Fcer1g","H2-Q7","H2-T22","B2m","H2-D1","Fcgr1","Psmb8","Thbs1","H2-Q4","H2-T24")

#combine into a named list for iteration
pathways <- list(
  Phagocytosis = phag_genes,
  Autophagosome = autophagy_genes,
  IFNg_signaling = ifny_genes,
  Antigen_presentation = antigen_genes
)

# 3) --- subset matrix to genes present, warn about missing genes ---
pathway_mats <- list()
for(p in names(pathways)) {
  genes <- pathways[[p]]
  present <- genes[genes %in% rownames(expr)]
  mat <- expr[present, , drop=FALSE]
  # preserve the order in your gene list (like the figure)
  mat <- mat[match(genes, rownames(mat), nomatch = 0), , drop = FALSE]
  pathway_mats[[p]] <- mat
  missing <- setdiff(genes, rownames(mat))
  if(length(missing)>0) message(p, ": missing genes -> ", paste(missing, collapse = ", "))
}



#z-score
pathway_z <- lapply(pathway_mats, function(m) {
  t(scale(t(m)))
})
View(pathway_z$Phagocytosis)

#define colour mapping
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


#create heatmap 
annotation_col <- data.frame(
  APOE = rep(c("APOE3", "APOE4"), times = c(4,4))
)
rownames(annotation_col) <- colnames(m)
pheatmap(
  pathway_z[["Phagocytosis"]],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(50),
  show_colnames = FALSE,
  main = "Phagocytosis"
)


pheatmap(
  pathway_z[["Autophagosome"]],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(50),
  show_colnames = FALSE,
  main = "Autophagosome"
)

pheatmap(
  pathway_z[["IFNg_signaling"]],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(50),
  show_colnames = FALSE,
  main = "IFNg_signaling"
)

pheatmap(
  pathway_z[["Antigen_presentation"]],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(50),
  show_colnames = FALSE,
  main = "Antigen_presentation"
)
