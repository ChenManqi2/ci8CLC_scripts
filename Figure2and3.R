########################################################################################################
########>>>>>                        Figure2 downstream analysis                           <<<<<########


# - 1st. Configuration
# - 2nd. Gene: 8CmCherry pos vs 8CmCherry and ci8CLCs vs prEpiSC --> Fig2O-P, Fig3G-I
# - 3rd. Repeat(Unique-mapping): 8CmCherry pos vs 8CmCherry and ci8CLCs vs prEpiSC --> Fig2Q, Fig3K
# - 4th. Repeat(Random-mapping): ci8CLCs vs prEpiSC --> Fig3J
# - 5th. Comparison of early human embryos, prEpiSCs, ci8CLCs, hESCs, DUX4OE hESCs, and FSHD muscle cells --> Fig3L


### ==================
### 1st. Configuration
### ==================

### >>> 1. save data
setwd("/home/data/wanglab/YMZ_hTPSCs_1st_batch/analysis/rnaseq/scripts")
#save.image("Figure2and3.RData")
#load("Figure2and3.RData")

### >>> 2. load package
library(dplyr)
library(tibble)
library(tidyr)
library(Seurat)
library(SeuratObject)
library(stringr)
library(scater)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(sva)
library(factoextra)
library(FactoMineR)

### >>> 3. load function
CountToTpm <- function(count, length){
  for (i in seq(1, ncol(count))){
    numer <- log(count[, i]) - log(length)
    denom <- log(sum(exp(numer)))
    count[, i] <- exp(numer - denom + log(1e6))
  }
  return(count)
}
CountToCpm <- function(geneCounts){
  for (i in seq(1, ncol(geneCounts))){
    temp.numer <- log(geneCounts[, i])
    temp.denom <- log(sum(exp(temp.numer)))
    geneCounts[, i] <- exp(temp.numer + log(1e6) - temp.denom)
  }
  # Return the resutls
  return(geneCounts)
}
DbiIDtrans <- function(genelist, intype, outtype, species){
  # Load the required packages
  # load the packages
  library("org.Hs.eg.db")
  library("org.Mm.eg.db")
  library("org.Rn.eg.db")
  library("AnnotationDbi")
  # Judge the species
  if (species == "human"){
    database <- "org.Hs.eg.db"
    database_short <- "hsa"
  } else if (species == "mouse"){
    database <- "org.Mm.eg.db"
    database_short <- "mmu"
  } else if (species == "rat"){
    database <-  "org.Rn.eg.db"
    database_short <- "rno"
  } else {
    print("You may need to define a new species dataset in the function")
  }
  # Transforming(keytypes list: keytypes(org.Hs.eg.db))
  if (species == "human"){
    newgenelist <- mapIds(x = org.Hs.eg.db, keys = unique(as.character(genelist)),
                          column = outtype, keytype = intype, multiVals = "first")
  } else if (species == "mouse"){
    newgenelist <- mapIds(x = org.Mm.eg.db, keys = unique(as.character(genelist)),
                          column = outtype, keytype = intype, multiVals = "first")
  } else if (species == "rat"){
    newgenelist <- mapIds(x = org.Rn.eg.db, keys = unique(as.character(genelist)),
                          column = outtype, keytype = intype, multiVals = "first")
  } else {
    print("You may need to define a new species dataset in the function")
  }
  # Retuen the results
  return(newgenelist)
}
edgeR.Ge.ET.LRT <- function(count, meta, g1, g2, lfc, sig, nor, with.rep){
  library("edgeR")
  library("dplyr")
  library("ggplot2")
  library("cowplot")
  con.index <- grep(paste("^", g2, "$", sep = ""), as.character(meta$SampleGroup))
  exp.index <- grep(paste("^", g1, "$", sep = ""), as.character(meta$SampleGroup))
  count <- count[, c(con.index, exp.index)]
  meta <- meta[c(con.index, exp.index), ]
  nor <- nor[, c(con.index, exp.index)]
  nor$SYMBOL <- rownames(nor)
  if (isTRUE(with.rep)){
    contrast <- factor(c(rep(1, length(con.index)), rep(2, length(exp.index))))
    deg <- DGEList(counts = count, group = contrast)
  } else if (isFALSE(with.rep)){
    deg <- DGEList(counts = count, samples = meta)
    deg$samples$group <- meta$SampleGroup
  }
  keep <- filterByExpr(deg, min.count = 10)
  deg <- deg[keep, , keep.lib.sizes = FALSE]
  deg <- calcNormFactors(deg, method = "TMM")
  lcpm <- as.data.frame(cpm(deg, log = F))
  lcpm$SYMBOL <- rownames(lcpm)
  if (isTRUE(with.rep)){
    design <- model.matrix(~contrast)
    deg <- estimateDisp(deg, design)
    deg.fit <- glmFit(deg, design)
    deg.test <- glmLRT(deg.fit, coef = 2)
  } else if (isFALSE(with.rep)){
    design <- model.matrix(~meta$SampleGroup)
    deg <- estimateGLMCommonDisp(deg, method = "deviance", robust = TRUE, subset = NULL)
    deg.test <- exactTest(deg, pair = c(g2, g1), rejection.region = "doubletail")
  }
  test.res <- topTags(deg.test, sort.by = "logFC", adjust.method = "fdr", n = nrow(deg))
  test.res$table$SYMBOL <- rownames(test.res$table)
  compare.info <- data.frame(Name = c("Comparison", "Testing for DE genes", "Method to adjust p.value"),
                             info = c(paste(g1, "vs", g2, sep = " "), test.res$test, test.res$adjust.method))
  res <- merge(test.res$table, lcpm, by = "SYMBOL")
  res <- merge(res, nor, by = "SYMBOL")
  colnames(res) <- gsub("\\.x", "\\.TMM", colnames(res))
  colnames(res) <- gsub("\\.y", "\\.TPM", colnames(res))
  res <- res %>% mutate(Log2MeanTMM = log2(rowMeans(res[, grep("TMM", colnames(res))])),
                        Log2MeanTPM = log2(rowMeans(res[, grep("TPM", colnames(res))])),
                        color = case_when((logFC >= lfc & FDR <= sig) ~ "Up.regulated", 
                                          (logFC <= -lfc & FDR <= sig) ~ "Down.regulated",
                                          (abs(logFC) < lfc | FDR > sig) ~ "Non.changed"),
                        trans = case_when((logFC >= lfc & FDR <= sig) ~ 1, 
                                          (logFC <= -lfc & FDR <= sig) ~ 1,
                                          (abs(logFC) < lfc | FDR > sig) ~ 0.5))
  all.sig <- subset(res, FDR < sig & abs(logFC) > lfc)
  up.sig <- subset(res, FDR < sig & logFC > lfc)
  down.sig <- subset(res, FDR < sig & logFC < -(lfc))
  all.res <- list(compare.info, lcpm, res, all.sig, up.sig, down.sig)
  names(all.res) <- c("compare.info", "lcpm", "res", "all.sig", "up.sig", "down.sig")
  return(all.res)
}
Deseq2.Ge <- function(count, meta, g1, g2, lfc, sig, nor){
  library("DESeq2")
  library("apeglm")
  con.index <- grep(paste("^", g2, "$", sep = ""), as.character(meta$SampleGroup))
  exp.index <- grep(paste("^", g1, "$", sep = ""), as.character(meta$SampleGroup))
  count <- count[, c(con.index, exp.index)]
  meta <- meta[c(con.index, exp.index), ]
  nor <- nor[, c(con.index, exp.index)]
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~SampleGroup)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds$SampleGroup <- relevel(dds$SampleGroup, g2)
  dds <- estimateSizeFactors(dds)
  count.nor <- counts(dds, normalized = TRUE)
  count.nor <- as.data.frame(count.nor)
  count.nor$SYMBOL <- rownames(count.nor)
  nor$SYMBOL <- rownames(nor)
  dds <- DESeq(dds)
  res <- results(dds, pAdjustMethod = "fdr", contrast = c("SampleGroup", g1, g2), parallel = T)
  summary(res)
  mcols(res, use.names = T)
  dea.info <- data.frame(mcols(res, use.names = T))
  resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")
  res <- as.data.frame(res)
  res$SYMBOL <- rownames(res)
  res <- merge(res, count.nor, by = "SYMBOL")
  res <- merge(res, nor, by = "SYMBOL")
  colnames(res) <- gsub("\\.x", "\\.RLE", colnames(res))
  colnames(res) <- gsub("\\.y", "\\.TPM", colnames(res))
  all.sig <- subset(res, abs(log2FoldChange) >= lfc & padj <= sig)
  up.sig <- subset(res, log2FoldChange >= lfc & padj <= sig)
  down.sig <- subset(res, log2FoldChange <= -(lfc) & padj <= sig)
  all.res <- list(dea.info, count.nor, res, all.sig, up.sig, down.sig, resLFC)
  names(all.res) <- c("dea.info", "count.nor", "res", "all.sig", "up.sig", "down.sig", "resLFC")
  return(all.res)
}
Deseq2.Te <- function(count, meta, g1, g2, lfc, sig, nor){
  library("DESeq2")
  library("apeglm")
  con.index <- grep(paste("^", g2, "$", sep = ""), as.character(meta$SampleGroup))
  exp.index <- grep(paste("^", g1, "$", sep = ""), as.character(meta$SampleGroup))
  count <- count[, c(con.index, exp.index)]
  meta <- meta[c(con.index, exp.index), ]
  nor <- nor[, c(con.index, exp.index)]
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~SampleGroup)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep, ]
  dds$SampleGroup <- relevel(dds$SampleGroup, g2)
  dds <- estimateSizeFactors(dds)
  count.nor <- counts(dds, normalized = TRUE)
  count.nor <- as.data.frame(count.nor)
  count.nor$Name <- rownames(count.nor)
  nor$Name <- rownames(nor)
  dds <- DESeq(dds)
  res <- results(dds, pAdjustMethod = "fdr", contrast = c("SampleGroup", g1, g2), parallel = F)
  summary(res)
  dea.info <- data.frame(mcols(res, use.names = T))
  resLFC <- lfcShrink(dds, coef = 2, type = "apeglm")
  resLFC <- as.data.frame(resLFC)
  res <- as.data.frame(res)
  res$Name <- rownames(res)
  res <- merge(res, count.nor, by = "Name")
  res <- merge(res, nor, by = "Name")
  colnames(res) <- gsub("\\.x", "\\.RLE", colnames(res))
  colnames(res) <- gsub("\\.y", "\\.CPM", colnames(res))
  all.sig <- subset(res, abs(log2FoldChange) >= lfc & padj <= sig)
  up.sig <- subset(res, log2FoldChange >= lfc & padj <= sig)
  down.sig <- subset(res, log2FoldChange <= -(lfc) & padj <= sig)
  all.res <- list(dea.info, count.nor, res, all.sig, up.sig, down.sig, resLFC)
  names(all.res) <- c("dea.info", "count.nor", "res", "all.sig", "up.sig", "down.sig", "resLFC")
  return(all.res)
}
pca.plot <- function(data, co.vec){
  fviz_pca_ind(PCA(t(data), scale.unit = TRUE, ncp = 5, graph = FALSE), axes = c(1, 2),
               geom.ind = c("point", "text"),
               pointshape = 16, pointsize = 6, alpha.ind = 1,
               col.ind = co.vec, 
               invisible =  c("quali"), 
               repel = T, max.overlaps = 100, 
               label = "all", labelsize = 3) + 
    ggtitle("2D PCA-plot") + theme_1
}
theme_1 <- theme_bw() + 
  theme(axis.title.x = element_text(face="plain", colour = "#000000", size = 12, angle = 0),
        axis.title.y = element_text(face="plain", colour = "#000000", size = 12, angle = 90),
        axis.text.x  = element_text(face="plain", colour = "#000000", size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(face="plain", colour = "#000000", size = 12, angle = 0, hjust = 1))

### =================================================================================
### 2nd. Gene: 8CmCherry pos vs 8CmCherry and ci8CLCs vs prEpiSC --> Fig2O-P, Fig3G-I
### =================================================================================

### >>> 1. load gene set
ge.set <- list()
ge.set$e3.gmt <- read.gmt("/home/cmq/bioinfo/EmboSC/E-MATB-3929/tables/combined_with_GSE36552/hs_embryo_8cell_specific_gene.gmt")

### >>> 2. load count table
# hTPSC data
# batch 1
tpsc.ge <- list()
tpsc.ge$co <- read.table("/home/data/wanglab/YMZ_hTPSCs_1st_batch/analysis/rnaseq/results/featurecounts/top10_gene/all_samples_gene_count_name_matrix.txt", header = T, row.names = 1, stringsAsFactors = F)
colnames(tpsc.ge$co)[-1:-5] <- c("Ctrl", "hTPSC_1", "hTPSC_2", "hTPSC_3") 
tpsc.ge$tpm <- CountToTpm(tpsc.ge$co[,-1:-5], tpsc.ge$co$Length)
# batch 2
tpsc.b2.ge <- list()
tpsc.b2.ge$co <- read.table("/home/data/wanglab/YMZ_hTPSCs_2nd_batch/bulk_rna_seq/analysis/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt", header = T, row.names = 1, stringsAsFactors = F)
colnames(tpsc.b2.ge$co)[-1:-5] <- c("hTPSC_b2_1", "Ctrl_b2_1", "hTPSC_b2_2", "Ctrl_b2_2")
tpsc.b2.ge$co <- tpsc.b2.ge$co[,c(colnames(tpsc.b2.ge$co)[1:5], "Ctrl_b2_1", "Ctrl_b2_2", "hTPSC_b2_1", "hTPSC_b2_2")]
tpsc.b2.ge$tpm <- CountToTpm(tpsc.b2.ge$co[,-1:-5], tpsc.b2.ge$co$Length)
# batch 1 and batch2
tpsc.b1b2.ge <- list()
tpsc.b1b2.ge$co <- cbind(tpsc.ge$co, tpsc.b2.ge$co)
tpsc.b1b2.ge$co <- tpsc.b1b2.ge$co[,c(1:6,15:16,8:9,17:18)]
tpsc.b1b2.ge$tpm <- CountToTpm(tpsc.b1b2.ge$co[,-1:-5], tpsc.b1b2.ge$co$Length)

### >>> 3. create metadata
dea.meta <- list()
# tpsc1/2/3 vs ctrl separately
dea.meta$tpsc2 <- data.frame(SampleID = colnames(tpsc.ge$co)[c(-1:-5)], SampleGroup = colnames(tpsc.ge$co)[c(-1:-5)], 
                             Batch = rep("1", 4), row.names = colnames(tpsc.ge$co)[c(-1:-5)])
# batch 1 and batch 2 
dea.meta$tpsc.b1b2 <- data.frame(SampleID = colnames(tpsc.b1b2.ge$tpm), 
                                 SampleGroup = str_split_fixed(colnames(tpsc.b1b2.ge$tpm), "_", "2")[,1], 
                                 Replicates = c(rep(1:3,2), 4), 
                                 row.names = colnames(tpsc.b1b2.ge$tpm))

### >>> 4. run DEA and visualization
# tpsc1(8C mCherry positive) vs ctrl 
tpsc1.vs.ctrl.ge <- edgeR.Ge.ET.LRT(tpsc.ge$co[,c(-1:-5)], dea.meta$tpsc2, "hTPSC_1", "Ctrl", 1, 0.05, tpsc.ge$tpm, F)
# visualization of MA plot
pd <- subset(tpsc1.vs.ctrl.ge$res, Ctrl.TPM >= 1 | hTPSC_1.TPM >= 1)
pd[,c(1:2,4:5)] %>% 
  mutate(logFDR=(-log10(FDR))) %>% 
  mutate(color = case_when((logFC >= 1 & FDR <= 0.05) ~ "Up-regulated",
                           (abs(logFC) < 1 | FDR > 0.05) ~ "No-changed",
                           (logFC <= -1 & FDR <= 0.05) ~ "Down-regulated")) -> pd
pd <- pd[c(grep("No-changed", pd$color), grep("Down-regulated", pd$color), grep("Up-regulated", pd$color)),]
pd.anno <- subset(pd, (SYMBOL %in% c("LEUTX","ZSCAN4","TPRX1","ZSCAN5B","TRIM43","TRIM43B","MBD3L2", "MBD3L3", "H3.Y","UBTFL7",
                                     "PRAMEF1","PRAMEF12", "PRAMEF2","PRAMEF4", "RFPL2","RFPL4A","SLC34A2","CCNA1", "KDM4E")) | logFC>=5)
#pdf("results/R/Graphs/v2/hTPSC_1_vs_Ctrl_differential_expressed_genes_plot.pdf", height = 6, width = 7)
pd %>% 
  ggplot(aes(x = logFC, y = logFDR)) + 
  geom_point(aes(color = color), size = 1.5, alpha = 0.5) + 
  geom_label_repel(aes(label = SYMBOL), data = pd.anno, 
                   box.padding = 0.1, segment.color = "#130f40", 
                   label.size = 0.1, 
                   size = 2.5, color = "#000000") +
  scale_color_manual(values = c(`Up-regulated` = "#d31d00", `Down-regulated` = "#0043d3", `No-changed` = "#b2bec3")) + 
  geom_hline(yintercept = c(-log10(0.05)), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) + 
  geom_vline(xintercept = c(-1,1), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) + 
  xlab("log2(Fold of Change)") + ylab("-log10(FDR)") + 
  theme_classic() +
  labs(color = "Color") +
  theme(legend.position = "bottom", legend.box="vertical")
#dev.off()

# tpsc vs ctrl
tpsc.b1b2.vs.ctrl.ge <- Deseq2.Ge(tpsc.b1b2.ge$co[,-1:-5], dea.meta$tpsc.b1b2, "hTPSC", "Ctrl", 1, 0.05, tpsc.b1b2.ge$tpm)
pd <- tpsc.b1b2.vs.ctrl.ge$res
pd[,c(1:3,6:7)] %>% 
  mutate(logFDR=(-log10(padj))) %>% 
  mutate(color = case_when((log2FoldChange >= 2 & padj <= 0.05) ~ "Up-regulated",
                           (abs(log2FoldChange) < 2 | padj > 0.05) ~ "No-changed",
                           (log2FoldChange <= -2 & padj <= 0.05) ~ "Down-regulated")) -> pd
pd <- pd[c(grep("No-changed", pd$color), grep("Down-regulated", pd$color), grep("Up-regulated", pd$color)),]
pd.anno <- subset(pd, (SYMBOL %in% c("LEUTX","ZSCAN4","TPRX1","ZSCAN5B","TRIM43","TRIM43B","MBD3L2","DUXA","UBTFL1",
                                     "KDM4E","PRAMEF1","PRAMEF12","RFPL2","RFPL4A","SLC34A2","CCNA1") | (log2FoldChange>=5 & logFDR >= 10)))
#pdf("results/R/Graphs/v2/batch1_2_hTPSC_vs_Ctrl_differential_expressed_genes_plot.pdf", height = 6, width = 7)
pd %>% 
  ggplot(aes(x = log2FoldChange, y = logFDR)) + 
  geom_point(aes(color = color), size = 1.5, alpha = 0.5) + 
  geom_label_repel(aes(label = SYMBOL), data = pd.anno, 
                   box.padding = 0.1, segment.color = "#130f40", 
                   label.size = 0.1, 
                   size = 2.5, color = "#000000", max.overlaps = 100) +
  scale_color_manual(values = c(`Up-regulated` = "#d31d00", `Down-regulated` = "#0043d3", `No-changed` = "#b2bec3")) + 
  geom_hline(yintercept = c(-log10(0.05)), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) + 
  geom_vline(xintercept = c(-2,2), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) + 
  xlab("log2(Fold of Change)") + ylab("-log10(FDR)") + 
  theme_classic() +
  labs(color = "Color") +
  theme(legend.position = "bottom", legend.box="vertical")
#dev.off()


### >>> 5. GSEA on 8Cell specific genes
deg.list <- list()
deg.rank <- list()
gsea.result <- list()
# htpsc1(8C mCherry positive) vs Ctrl
deg.list$tpsc1.vs.ctrl.ge <- tpsc1.vs.ctrl.ge$all.sig$SYMBOL
deg.rank$tpsc1.vs.ctrl.ge <- tpsc1.vs.ctrl.ge$all.sig$logFC
names(deg.rank$tpsc1.vs.ctrl.ge) <- as.vector(deg.list$tpsc1.vs.ctrl.ge)
deg.rank$tpsc1.vs.ctrl.ge <- sort(deg.rank$tpsc1.vs.ctrl.ge, decreasing = T)
gsea.result$tpsc1.vs.ctrl.ge <- GSEA(deg.rank$tpsc1.vs.ctrl.ge, TERM2GENE = ge.set$e3.gmt, verbose=FALSE, pvalueCutoff = 0.5)
# htpsc vs ctrl
deg.list$tpsc.b1b2 <- tpsc.b1b2.vs.ctrl.ge$all.sig$SYMBOL
deg.rank$tpsc.b1b2 <- tpsc.b1b2.vs.ctrl.ge$all.sig$log2FoldChange
names(deg.rank$tpsc.b1b2) <- as.vector(deg.list$tpsc.b1b2)
deg.rank$tpsc.b1b2 <- sort(deg.rank$tpsc.b1b2, decreasing = T)
# run GSEA
gsea.result <- list()
for (i in seq(1, length(deg.rank))) {
  gsea.result[[i]] <- GSEA(deg.rank[[i]], TERM2GENE = ge.set$e3.gmt, verbose=FALSE, pvalueCutoff = 0.5)
  names(gsea.result)[i] <- names(deg.rank)[i]
  #pdf(paste("results/R/Graphs/v2/", names(deg.rank)[i], "_vs_ctrl_dea_gsea_on_8cell_specific_genes.pdf", sep = ""), 
  #    height = 4, width = 4)
  print(gseaplot2(gsea.result[[i]], geneSetID = "human 8cell specific genes", pvalue_table = T))
  #dev.off()
}; rm(i)


### >>> 6. 
# load pluripotent genes
ge.set$pluri <- c("UTF1", "KLF4", "DNMT3L", "ZFP42", "DNMT3A", "BMP2", "DNMT3B", "SALL4", "PRDM14", "NANOG", "TFCP2L1", "POU5F1", "ESRG", "IL6R", "SOX2", 
                  "TDGF1", "FGF4", "IFI16", "IFITM1", "IFITM2", "IFITM3", "GDF3")
# load embryo scRNAseq data
emb.ge$tpm.nt <- readRDS("/home/cmq/bioinfo/project-cmq/gse36552/results/R/RDS/hs.sc.ge.tpm.rds")
cell.meta$nt <- readRDS("/home/cmq/bioinfo/project-cmq/gse36552/results/R/RDS/hs.sc.meta.full.rds")
cell.meta$nt <- subset(cell.meta$nt, !(cell_id %in% c("EPI_4", "EPI_5")))
emb.ge$tpm.nt <- emb.ge$tpm.nt[,grep("ave|med|var|EPI_4|EPI_5", colnames(emb.ge$tpm.nt), invert = T)]
emb.ge$tpm.nt %>% mutate(Oocyte_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="Oocyte",]$cell_id]), 
                         Zygote_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="Zygote",]$cell_id]), 
                         `2Cell_ave`=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="2-cell",]$cell_id]),
                         `4Cell_ave`=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="4-cell",]$cell_id]),
                         `8Cell_ave`=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="8-cell",]$cell_id]),
                         Morula_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="Morula",]$cell_id]),
                         EPI_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="EPI",]$cell_id]),
                         TE_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$cell.type=="TE",]$cell_id]), 
                         hESC_p0_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$embryo.passage=="hESC P0",]$cell_id]), 
                         hESC_P10_ave=rowMeans(emb.ge$tpm.nt[,cell.meta$nt[cell.meta$nt$embryo.passage=="hESC P10",]$cell_id])) -> emb.ge$tpm.nt
emb.ge$tpm.nt$ENSEMBL <- rownames(emb.ge$tpm.nt)
emb.ge$tpm.nt$SYMBOL <- DbiIDtrans(emb.ge$tpm.nt$ENSEMBL, "ENSEMBL", "SYMBOL", "human")
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000254764",]$SYMBOL <- "TRIM53CP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000237194",]$SYMBOL <- "SNAI1P1"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000269118",]$SYMBOL <- "FAM90A28P"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000267908",]$SYMBOL <- "ZSCAN5DP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000166013",]$SYMBOL <- "TRIM53BP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000144188",]$SYMBOL <- "TRIM43CP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000226185",]$SYMBOL <- "TRIM64FP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000251360",]$SYMBOL <- "KHDC1P1"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000254751",]$SYMBOL <- "TRIM64DP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000284306",]$SYMBOL <- "SSU72P2"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000204455",]$SYMBOL <- "TRIM51BP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000188831",]$SYMBOL <- "DPPA3P2"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000249910",]$SYMBOL <- "TRIM51CP"
emb.ge$tpm.nt[emb.ge$tpm.nt$ENSEMBL=="ENSG00000229437",]$SYMBOL <- "LINC00366"
emb.ge$tpm.nt <- emb.ge$tpm.nt[!emb.ge$tpm.nt$SYMBOL%in%NA,]
emb.ge$tpm.nt <- emb.ge$tpm.nt[!emb.ge$tpm.nt$SYMBOL %in% emb.ge$tpm.nt[duplicated(emb.ge$tpm.nt$SYMBOL),]$SYMBOL,]
rownames(emb.ge$tpm.nt) <- emb.ge$tpm.nt$SYMBOL
# combine embryo and hTPSC tpm table
emb.tpsc.ge <- list()
# raw TPM
emb.tpsc.ge$tpm <- cbind(tpsc.ge$tpm, tpsc.b2.ge$tpm)
emb.tpsc.ge$tpm[intersect(rownames(emb.tpsc.ge$tpm), rownames(emb.ge$tpm.nt)),] %>% rownames_to_column("SYMBOL") %>% 
  merge(emb.ge$tpm.nt[intersect(rownames(emb.tpsc.ge$tpm), rownames(emb.ge$tpm.nt)),grep("ave|SYMBOL", colnames(emb.ge$tpm.nt))]) -> emb.tpsc.ge$tpm
rownames(emb.tpsc.ge$tpm) <- emb.tpsc.ge$tpm$SYMBOL
emb.tpsc.ge$tpm <- emb.tpsc.ge$tpm[,c("Ctrl", "Ctrl_b2_1", "Ctrl_b2_2", "hTPSC_2", "hTPSC_3", "hTPSC_b2_1", "hTPSC_b2_2", colnames(emb.tpsc.ge$tpm)[10:19])]
# zscore of raw TPM
emb.tpsc.ge$tpm.zs <- as.data.frame(t(apply(emb.tpsc.ge$tpm, 1, function(x){(x - mean(x))/(sd(x))})))
# batch effect corrected TPM
sva.meta <- data.frame(SampleID=colnames(emb.tpsc.ge$tpm), 
                       SampleInfo=c("hTPSC", rep("hTPSC_b2",2), rep("hTPSC",2), rep("hTPSC_b2",2), rep("Nature",10)), 
                       Batch=c("hTPSC", rep("hTPSC_b2",2), rep("hTPSC",2), rep("hTPSC_b2",2), rep("Nature",10)))
sva.modcombat <- model.matrix(~1, data = sva.meta)
emb.tpsc.ge$tpm.rb <- ComBat(dat = as.matrix(emb.tpsc.ge$tpm), batch = sva.meta$Batch,
                             mod = sva.modcombat, par.prior = TRUE, prior.plots = FALSE) %>% as.data.frame() 
emb.tpsc.ge$tpm.rb[emb.tpsc.ge$tpm.rb<0] <- 0
# zscore of batch effect corrected TPM
emb.tpsc.ge$tpm.rb.zs <- as.data.frame(t(apply(emb.tpsc.ge$tpm.rb, 1, function(x){(x - mean(x))/(sd(x))})))
# plotting
pd <- emb.tpsc.ge$tpm.rb.zs[ge.set$pluri,]
#pdf("results/R/Graphs/v2/hTPSC_and_embryo_expression_level_of_pluripotent_genes_zscore_of_batch_corrected_tpm_heatmap.pdf", width = 10, height = 10)
p1 <- Heatmap(as.matrix(pd[,grep("hTPSC|Ctrl", colnames(pd))]),
              # colour and legend name setting
              col = colorRamp2(c(-2, 0, 2), c("#1078d4", "#ffffff", "#f03039")), 
              name = "Cell Line: Zscore", use_raster = F, 
              # clustering setting
              cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
              clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
              cluster_columns = F, show_column_dend = F, column_dend_height = unit(1, "cm"), column_dend_side = "top",
              clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
              # names, rotation and position of columns and row names
              show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
              show_row_names = F, row_names_rot = 0, row_names_side = "right", 
              # column split and annotation
              column_order = colnames(pd[,grep("hTPSC|Ctrl", colnames(pd))]), 
              # size of heatmap, also including heatmap_height and heatmap_width
              width = unit(1.5, "cm"), height = unit(10, "cm")) 
p1 + Heatmap(as.matrix(pd[,grep("ave", colnames(pd))]),
             # colour and legend name setting
             col = colorRamp2(c(-2, 0, 2), c("#1078d4", "#ffffff", "#f03039")), 
             name = "Embryo: Zscore", use_raster = F, 
             # clustering setting
             cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
             clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
             cluster_columns = F, show_column_dend = F, column_dend_height = unit(1, "cm"), column_dend_side = "top",
             clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
             # names, rotation and position of columns and row names
             show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
             show_row_names = T, row_names_rot = 0, row_names_side = "right", 
             # column split and annotation
             column_order = colnames(pd[,grep("ave", colnames(pd))]), 
             # size of heatmap, also including heatmap_height and heatmap_width
             width = unit(3, "cm"), height = unit(10, "cm")) 
dev.off()


### ================================================================================================
### 3rd. Repeat(Unique-mapping): 8CmCherry pos vs 8CmCherry and ci8CLCs vs prEpiSC --> Fig2Q, Fig3K
### ================================================================================================

### >>> 1. load embryo repeat count table
emb.re <- list()
emb.re$co.gse101571 <- read.table("/home2/data/publicdata/GSE101571/rnaseq/analysis/results/featurecounts/top10_repeat/all_samples_repeat_count_matrix.txt", row.names = 1, stringsAsFactors = F, header = T)
cell.meta <- list()
cell.meta$gse101571 <- read.table("/home2/data/publicdata/GSE101571/rnaseq/analysis/sample_meta.txt", header = F)
colnames(emb.re$co.gse101571)[-1:-5] <- str_split_fixed(str_split_fixed(colnames(emb.re$co.gse101571)[-1:-5], "repeat.", 2)[,2], "_1_val", 2)[,1]
if (all(colnames(emb.re$co.gse101571)[-1:-5]==cell.meta$gse101571$V1)) {
  colnames(emb.re$co.gse101571)[-1:-5] <- cell.meta$gse101571$V2
}
emb.re$tpm.gse101571 <- CountToTpm(emb.re$co.gse101571[,-1:-5], emb.re$co.gse101571$Length)

### >>> 2. load hTPSC repeat count table
tpsc.re <- list()
tpsc.re$raw <- read.table("/home/data/wanglab/YMZ_hTPSCs_1st_batch/analysis/rnaseq/results/featurecounts/top10_repeat/all_samples_repeat_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(tpsc.re$raw)[-1:-6] <- c("Ctrl", "hTPSC_1", "hTPSC_2", "hTPSC_3")
tpsc.re$tpm <- CountToTpm(tpsc.re$raw[,-1:-6], tpsc.re$raw$Length)
rownames(tpsc.re$tpm) <- tpsc.re$raw$Geneid
tpsc.b2.re <- list()
tpsc.b2.re$raw <- read.table("/home/data/wanglab/YMZ_hTPSCs_2nd_batch/bulk_rna_seq/analysis/results/featurecounts/repeat_uni/all_samples_repeat_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(tpsc.b2.re$raw)[-1:-6] <- c("hTPSC_b2_1", "Ctrl_b2_1", "hTPSC_b2_2", "Ctrl_b2_2")
tpsc.b2.re$tpm <- CountToTpm(tpsc.b2.re$raw[,-1:-6], tpsc.b2.re$raw$Length)
rownames(tpsc.b2.re$tpm) <- tpsc.b2.re$raw$Geneid

### >>> 3. merge embryo and hTPSC repeat tpm
emb.tpsc.re <- list()
emb.tpsc.re$tpm <- cbind(emb.re$tpm.gse101571, tpsc.re$tpm, tpsc.b2.re$tpm)
emb.tpsc.re$tpm %>% rownames_to_column("TeID") %>% 
  mutate(`8Cell_ave`=rowMeans(emb.tpsc.re$tpm[,grep("8Cell", colnames(emb.tpsc.re$tpm))]),
         ICM_ave=rowMeans(emb.tpsc.re$tpm[,grep("ICM", colnames(emb.tpsc.re$tpm))]), 
         Ctrl_ave=rowMeans(emb.tpsc.re$tpm[,grep("Ctrl", colnames(emb.tpsc.re$tpm))]), 
         hTPSC_ave=rowMeans(emb.tpsc.re$tpm[,grep("hTPSC_[2-3]|hTPSC_b", colnames(emb.tpsc.re$tpm))])) -> emb.tpsc.re$tpm

### >>> 4. create repeat annotation files
re.anno <- data.frame(TeID=emb.tpsc.re$tpm$TeID)
re.anno <- separate(re.anno, "TeID", into = c("Class", "Family", "Name", "Locus"), sep = ":", remove = F)
# calculate number of repeat locus
tmp <- table(re.anno$Name) %>% as.data.frame()
colnames(tmp) <- c("Name", "Loci.Number")
# merge raw annotation file and loci number table
re.anno <- merge(re.anno, tmp, "Name")
rownames(re.anno) <- re.anno$TeID

### >>> 5. merge repeat annotation file and repeat tpm table
emb.tpsc.re$tpm <- merge(emb.tpsc.re$tpm, re.anno, "TeID")

### >>> 6. calculate activate locus number and ratio
emb.tpsc.re$tpm_a01 <- emb.tpsc.re$tpm %>% group_by(Name) %>% 
  mutate(`8Cell_an`=sum(`8Cell_ave`>=0.1), `8Cell_af`=sum(`8Cell_ave`>=0.1)/Loci.Number, 
         ICM_an=sum(ICM_ave>=0.1), ICM_af=sum(ICM_ave>=0.1)/Loci.Number, 
         Ctrl_an=sum(Ctrl_ave>=0.1), Ctrl_af=sum(Ctrl_ave>=0.1)/Loci.Number, 
         hTPSC_an=sum(hTPSC_ave>=0.1), hTPSC_af=sum(hTPSC_ave>=0.1)/Loci.Number,
         Ctrl_1_an=sum(Ctrl>=0.1), Ctrl_1_af=sum(Ctrl>=0.1)/Loci.Number, 
         hTPSC_1_an=sum(hTPSC_1>=0.1), hTPSC_1_af=sum(hTPSC_1>=0.1)/Loci.Number) %>% as.data.frame()
emb.tpsc.re$tpm_a01 <- emb.tpsc.re$tpm_a01[,c(27:29,32:43)] %>% unique()
emb.tpsc.re$tpm_a01 <- emb.tpsc.re$tpm_a01 %>% filter(!Class %in% c("DNA?", "LTR?", "RC?", "SINE?"))

### >>> 7. barplot hTPSC_1 vs Ctrl_1 --> Fig2Q
# data preparation
pd <- subset(emb.tpsc.re$tpm_a01, Name %in% c("HERVH-int", "LTR7", "LTR7B", "MLT2A1", "HERVL-int", "LTR12C"))
pd <- data.frame(Name=c(pd$Name, pd$Name), fraction=c(pd$hTPSC_1_af, pd$Ctrl_1_af), 
                 number=c(pd$hTPSC_1_an, pd$Ctrl_1_an), sample=c(rep("hTPSC_1", length(pd$Name)), rep("Ctrl_1", length(pd$Name))))
pd %>% gather(key = "type", value = "y", -Name, -sample) -> pd
pd %>% mutate(y=case_when(type=="fraction" ~ y*1000, 
                          type=="number" ~ y)) %>% 
  mutate(type=factor(type, levels = c("number", "fraction"))) -> pd
#pdf("results/R/Graphs/v2/active_repeat_unique/hTPSC_1_and_Ctrl_1_active_repeats_barplot_TPM0.1.pdf", height = 4.5, width = 8)
pd %>% 
  ggplot(aes(x=type, y=y, fill=sample)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  scale_y_continuous(name="number", sec.axis = sec_axis(~ 0.001*., name="fraction")) + 
  facet_wrap(.~Name, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
        legend.text = element_text(face = "plain", colour = "#000000", size = 10, angle = 0), 
        strip.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"))
#dev.off()

### >>> 8. dotplot htpsc vs ctrl (ci8CLC vs prEpiSC) --> Fig3K
# load prEpiSC and primed repeat raw count table
prEpiSC.re <- list()
prEpiSC.re$raw <- read.table("/home/data/wanglab/ICM_like/analysis/rnaseq/hs38/results/featurecounts/repeat_uni/all_samples_repeat_count_matrix.txt", row.names = 1, header = T, stringsAsFactors = F)
colnames(prEpiSC.re$raw)[-1:-5] <- str_split_fixed(str_split_fixed(colnames(prEpiSC.re$raw)[-1:-5], "repeat.", 2)[,2], "_1", 2)[,1]
prEpiSC.re$tpm <- CountToTpm(prEpiSC.re$raw[,-1:-5], prEpiSC.re$raw$Length)
# create table for plot
emb.tpsc.re$tpm2 <- prEpiSC.re$tpm %>% rownames_to_column("TeID") %>% merge(emb.tpsc.re$tpm, by = "TeID")
pd.repeat <- c("LTR12C", "MLT2A1", "LTR7", "LTR7Y", "HERVL-int")
# plotting
subset(emb.tpsc.re$tpm2, Name %in% pd.repeat) %>% 
  subset(hTPSC_ave >= 0.1 | Ctrl_ave >= 0.1) %>%  
  mutate(color=case_when(hTPSC_ave>=Ctrl_ave ~ "Up", 
                         Ctrl_ave>hTPSC_ave ~ "Down")) %>% 
  group_by(color) %>% as.data.frame() -> pd
print(pd %>% 
        ggplot(aes(x=log2(Ctrl_ave+0.1), y=log2(hTPSC_ave+0.1), color=color)) + 
        geom_point() + 
        geom_abline(intercept = 0, slope = 1, colour = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
        facet_wrap(.~Name) + 
        theme_bw() + 
        theme(axis.title.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 90),
              axis.title.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              axis.text.x = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              axis.text.y = element_text(face = "plain", colour = "#000000", size = 12, angle = 0),
              strip.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid")))

### =========================================================
### 4th. Repeat(Random-mapping): ci8CLCs vs prEpiSC --> Fig3J
### =========================================================
### >>> 1. read raw count table
tpsc.re <- list()
tpsc.re$raw <- read.table("/home/data/wanglab/YMZ_hTPSCs_1st_batch/analysis/rnaseq/results/featurecounts/random_repeat/all_samples_repeat_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(tpsc.re$raw)[-1:-6] <- c("Ctrl", "hTPSC_1", "hTPSC_2", "hTPSC_3")
tpsc.b2.re <- list()
tpsc.b2.re$raw <- read.table("/home/data/wanglab/YMZ_hTPSCs_2nd_batch/bulk_rna_seq/analysis/results/featurecounts/repeat/all_samples_repeat_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(tpsc.b2.re$raw)[-1:-6] <- c("hTPSC_b2_1", "Ctrl_b2_1", "hTPSC_b2_2", "Ctrl_b2_2")

### >>> 2. data pre-processing
# merge batch1 and batch2 raw count table
tpsc.b1b2.re <- list()
tpsc.b1b2.re$raw <- merge(tpsc.re$raw, tpsc.b2.re$raw, c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
tpsc.b1b2.re$raw <- tpsc.b1b2.re$raw[,c(1:7,12,14,9:11,13)]
colnames(tpsc.b1b2.re$raw)
# merge count of all locus into subfamily 
tpsc.b1b2.re$raw %>% separate("Geneid", into = c("Class", "Family", "Name", "Locus"), sep = ":", remove = T) -> tpsc.b1b2.re$name
tpsc.b1b2.re$name %>% group_by(Name) %>% mutate(Ctrl=sum(Ctrl), Ctrl_b2_1=sum(Ctrl_b2_1), Ctrl_b2_2=sum(Ctrl_b2_2), hTPSC_2=sum(hTPSC_2), 
                                                hTPSC_3=sum(hTPSC_3), hTPSC_b2_1=sum(hTPSC_b2_1), hTPSC_b2_2=sum(hTPSC_b2_2)) %>% as.data.frame() -> tpsc.b1b2.re$name
colnames(tpsc.b1b2.re$name)
tpsc.b1b2.re$name <- tpsc.b1b2.re$name[,c(3,10:16)] %>% unique()
rownames(tpsc.b1b2.re$name) <- tpsc.b1b2.re$name$Name
tpsc.b1b2.re$name <- tpsc.b1b2.re$name[,-1]
# calculate cpm
re.cpm <- list()
re.cpm$tpsc.b1b2 <- CountToCpm(tpsc.b1b2.re$name)
all(dea.meta$tpsc.b1b2$SampleID==colnames(re.cpm$tpsc.b1b2))
# run DEseq2
tpsc.b1b2.vs.ctrl.re <- Deseq2.Te(tpsc.b1b2.re$name, dea.meta$tpsc.b1b2, "hTPSC", "Ctrl", 1, 0.05, re.cpm$tpsc.b1b2)

### >>> 3. visualization
pd <- tpsc.b1b2.vs.ctrl.re$res
pd[,c(1:3,6:7)] %>% 
  mutate(logPvalue=(-log10(pvalue))) %>% 
  mutate(color = case_when((log2FoldChange >= 1 & pvalue <= 0.05) ~ "Up-regulated",
                           (abs(log2FoldChange) < 1 | pvalue > 0.05) ~ "No-changed",
                           (log2FoldChange <= -1 & pvalue <= 0.05) ~ "Down-regulated")) -> pd
pd <- pd[c(grep("No-changed", pd$color), grep("Down-regulated", pd$color), grep("Up-regulated", pd$color)),]
pd.anno <- subset(pd, (Name %in% c("MLT2A1", "ERVL-int", "MLT2A2", "HSAT2", "ACRO1", "LTR7", "HERVH-int", "LTR5_Hs", 
                                   "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F") | (log2FoldChange>=1 & logPvalue>=5)))
#pdf("results/R/Graphs/v2/batch1_2_hTPSC_vs_Ctrl_differential_expressed_repeats_name_plot_random.pdf", height = 6, width = 7)
pd %>% 
  ggplot(aes(x = log2FoldChange, y = logPvalue)) + 
  geom_point(aes(color = color), size = 1.5, alpha = 0.5) + 
  geom_label_repel(aes(label = Name), data = pd.anno, 
                   box.padding = 0.1, segment.color = "#130f40", 
                   label.size = 0.1, 
                   size = 2.5, color = "#000000", max.overlaps = 100) +
  scale_color_manual(values = c(`Up-regulated` = "#d31d00", `Down-regulated` = "#0043d3", `No-changed` = "#b2bec3")) + 
  geom_hline(yintercept = c(-log10(0.05)), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) + 
  geom_vline(xintercept = c(-1,1), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) + 
  xlab("log2(Fold of Change)") + ylab("-log10(Pvalue)") + 
  xlim(c(-6,6)) + 
  theme_classic() +
  labs(color = "Color") +
  theme(legend.position = "bottom", legend.box="vertical")
#dev.off()


### ===============================================================================================================
### 5th. Comparison of early human embryos, prEpiSCs, ci8CLCs, hESCs, DUX4OE hESCs, and FSHD muscle cells --> Fig3L
### ===============================================================================================================
### >>> 1. load raw count
# human embryo data
emb.ge <- list()
emb.ge$co.gse101571 <- read.table("/home2/data/publicdata/GSE101571/rnaseq/analysis/results/featurecounts/top10_gene/all_samples_gene_count_name_matrix.txt",
                                  header = T, sep = "\t", row.names = 1)
colnames(emb.ge$co.gse101571)[-1:-5] <- str_split_fixed(str_split_fixed(colnames(emb.ge$co.gse101571)[-1:-5], "gene.", 2)[,2], "_1_val", 2)[,1]
if (all(colnames(emb.ge$co.gse101571)[-1:-5]==cell.meta$gse101571$V1)) {
  colnames(emb.ge$co.gse101571)[-1:-5] <- cell.meta$gse101571$V2
}
emb.ge$co.gse101571 <- emb.ge$co.gse101571[,grep("2Cell_3PN", colnames(emb.ge$co.gse101571), invert = T)]
emb.ge$tpm.gse101571 <- CountToTpm(emb.ge$co.gse101571[,-1:-5], emb.ge$co.gse101571$Length)
# DUX4 oe data
dux4.oe.ge <- list()
dux4.oe.ge$gse85632.co <- read.table("/home2/data/publicdata/GSE85632/analysis/hipsc_rnaseq/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt", header = T, row.names = 1, stringsAsFactors = F)
colnames(dux4.oe.ge$gse85632.co)[-1:-5] <- str_split_fixed(str_split_fixed(colnames(dux4.oe.ge$gse85632.co)[-1:-5], "gene.", 2)[,2], "_merge", 2)[,1]
dux4.oe.ge$gse153301.co <- read.table("/home2/data/publicdata/GSE153301/analysis/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt", header = T, row.names = 1, stringsAsFactors = F)
colnames(dux4.oe.ge$gse153301.co)[-1:-5] <- c("WT_S1", "WT_S2", "WT_S3", "FSHD_S4", "FSHD_S5", "FSHD_S6")
# calculate tpm value
for (i in seq(1, length(dux4.oe.ge))) {
  dux4.oe.ge[[i+2]] <- CountToTpm(dux4.oe.ge[[i]][,-1:-5], dux4.oe.ge[[i]]$Length)
}; rm(i)
names(dux4.oe.ge)[3:4] <- c("gse85632.tpm", "gse153301.tpm")

### >>> 2. combine embryo, hTPSC and dux4 oe data
# merge tpm value
mix.ge <- list()
mix.ge$tpm <- cbind(tpsc.ge$tpm, tpsc.b2.ge$tpm, emb.ge$tpm.gse101571, dux4.oe.ge$gse85632.tpm, dux4.oe.ge$gse153301.tpm)
# create sample meta
mix.meta <- data.frame(Sample=colnames(mix.ge$tpm), Source=c(rep("WangLab1", 4), rep("WangLab2", 4), rep("Embryo", 12), rep("GSE85632", 8), rep("GSE153301", 6)))
mix.meta %>% mutate(Type=case_when(str_detect(Sample, "Ctrl") ~ "Ctrl", str_detect(Sample, "hTPSC") ~ "hTPSC",
                                   str_detect(Sample, "GV") ~ "GV_oocyte", str_detect(Sample, "M2") ~ "M2_oocyte",
                                   str_detect(Sample, "2Cell") ~ "2Cell", str_detect(Sample, "4Cell") ~ "4Cell",
                                   str_detect(Sample, "8Cell") ~ "8Cell", str_detect(Sample, "ICM") ~ "ICM",
                                   str_detect(Sample, "DUX4CA_14h") ~ "DUX4CA_14h", str_detect(Sample, "DUX4CA_24h") ~ "DUX4CA_24h",
                                   str_detect(Sample, "Luciferase_14h") ~ "Luciferase_14h", str_detect(Sample, "Luciferase_24h") ~ "Luciferase_24h",
                                   str_detect(Sample, "WT_") ~ "WT", str_detect(Sample, "FSHD") ~ "FSHD")) -> mix.meta
rownames(mix.meta) <- mix.meta$Sample
# create metadata
mix.ge$tpm.1 <- mix.ge$tpm[,grep("hTPSC_1", colnames(mix.ge$tpm), invert = T)]
sva.meta <- data.frame(SampleID=mix.meta[colnames(mix.ge$tpm.1),]$Sample,
                       SampleInfo=mix.meta[colnames(mix.ge$tpm.1),]$Source,
                       Batch=mix.meta[colnames(mix.ge$tpm.1),]$Source)
sva.modcombat <- model.matrix(~1, data = sva.meta)
mix.ge$tpm.1.rb <- ComBat(dat = as.matrix(mix.ge$tpm.1), batch = sva.meta$Batch,
                          mod = sva.modcombat, par.prior = TRUE, prior.plots = FALSE) %>% as.data.frame()
mix.ge$tpm.1.rb[mix.ge$tpm.1.rb<0] <- 0
# run pca
pd <- mix.ge$tpm.rb[,grep("hTPSC_1", colnames(mix.ge$tpm.rb), invert = T)]
pd <- subset(pd, rowSums(pd)>=10)
pca.plot(pd, mix.meta[colnames(pd),]$Type)


