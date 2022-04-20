########################################################################################################
########>>>>>                        Figure1 downstream analysis                           <<<<<########


# - 1st. Configuration
# - 2nd. Analysis of human embryo scRNA-seq and define 8Cell specific genes --> Fig1A
# - 3rd. Analysis of 10x prEpiSC scRNA-seq --> Fig1E-I
# - 4th. Comparison between prEpiSC and naive or conventional hESC --> Fig1B-D


### ==================
### 1st. Configuration
### ==================

### >>> 1. save data
setwd("/home/data/wanglab/YMZ_hTPSCs_1st_batch/analysis/rnaseq/scripts")
#save.image("Figure1.RData")
#load("Figure1.RData")

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

### >>> 3. load function
CountToTpm <- function(count, length){
  for (i in seq(1, ncol(count))){
    numer <- log(count[, i]) - log(length)
    denom <- log(sum(exp(numer)))
    count[, i] <- exp(numer - denom + log(1e6))
  }
  return(count)
}
edgeR <- function(count, meta, g1, g2, lfc, sig){
  library("edgeR")
  # obtain samples index
  g1.cols <- grep(paste("^", g1, "$", sep = ""), meta$Stage)
  g2.cols <- grep(paste("^", g2, "$", sep = ""), meta$Stage)
  # re-organize count data and meta data
  tar.meta <- meta[c(g1.cols, g2.cols), ]
  print(table(tar.meta$Stage))
  tar.count <- count[, c(g1.cols, g2.cols)]
  tar.count <- subset(tar.count, rowSums(tar.count) >= 10)
  # differential expression analysis
  contrast <- factor(c(rep(1, length(g2.cols)), rep(2, length(g1.cols))))
  dge <- DGEList(counts = tar.count, samples = tar.meta, group = contrast)
  dge <- calcNormFactors(dge)
  cdr <- scale(colMeans(dge$counts > 0))
  design <- model.matrix(~ cdr + dge$samples$Stage)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf)
  degs.up <- subset(tt$table, logFC >= lfc & PValue <= sig)
  degs.down  <- subset(tt$table, logFC <= -lfc & PValue <= sig)
  degs.all <- subset(tt$table, abs(logFC) >= lfc & PValue <= sig)
  degs.non <- subset(tt$table, abs(logFC) < lfc | PValue > sig)
  degs.all.res <- list(tt$table, degs.all, degs.up, degs.down, degs.non)
  names(degs.all.res) <- c("all", "degs.all", "degs.up", "degs.down", "degs.non")
  return(degs.all.res)
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
  all.sig <- subset(res, abs(log2FoldChange) >= lfc & pvalue <= sig)
  up.sig <- subset(res, log2FoldChange >= lfc & pvalue <= sig)
  down.sig <- subset(res, log2FoldChange <= -(lfc) & pvalue <= sig)
  return.all <- list(dea.info, count.nor, res, all.sig, up.sig, down.sig, resLFC)
  names(return.all) <- c("dea.info", "count.nor", "res", "all.sig", "up.sig", "down.sig", "resLFC")
  return(return.all)
}

### =================================================================================
### 2nd. Analysis of human embryo scRNA-seq and define 8Cell specific genes --> Fig1A
### =================================================================================
### >>> 1. gene annotation file
ge.anno <- readRDS("/home/cx/work/scRNA-seq/ensembl_hg38_v34_gene_anno.rds")
ge.anno$ENSEMBL <- rownames(ge.anno)
### >>> 2. sample annotation
s.anno <- list()
# 2013 Nature
s.anno$`2013nature` <- read.csv("/home/cmq/bioinfo/EmboSC/GSE36552/upstream_from_cx/cell_SRR_cellname.txt", stringsAsFactors = F, sep = "\t", header = F)
colnames(s.anno$`2013nature`) <- c("ID", "Cell.raw", "Cell")
s.anno$`2013nature` %>% 
  mutate(Stage = case_when(str_detect(Cell, regex("Oocyte")) ~ "Oocyte", str_detect(Cell, regex("Zygote")) ~ "Zygote", 
                           str_detect(Cell, regex("2-cell")) ~ "2Cell", str_detect(Cell, regex("4-cell")) ~ "4Cell", 
                           str_detect(Cell, regex("8-cell")) ~ "E3", str_detect(Cell, regex("Morula")) ~ "E4", 
                           str_detect(Cell, regex("LB")) ~ "LB", str_detect(Cell, regex("hESC")) ~ "hESC")) %>% 
  mutate(Embryo = str_split_fixed(str_split_fixed(Cell, " ", 2)[,2], "\\.", 2)[,1]) -> s.anno$`2013nature`
s.anno$`2013nature` <- s.anno$`2013nature`[,c(3, 5, 1, 4)]
s.anno$`2013nature` <- subset(s.anno$`2013nature`, !Stage%in%c("hESC"))

# 2016 Cell
s.anno$`2016cell` <- read.table("/home/cmq/bioinfo/EmboSC/E-MATB-3929/upstream_from_cx/sample_anno.txt", header = T, stringsAsFactors = F, sep = "\t")
colnames(s.anno$`2016cell`) <- c("Cell", "Embryo", "ID")
s.anno$`2016cell` %>% mutate(Stage = str_split_fixed(s.anno$`2016cell`$Embryo, "\\.", 2)[,1]) -> s.anno$`2016cell`

### >>> 3. raw count table
emb.co <- list()
# 2013 Nature
emb.co$`2013nature` <- read.table("/home/cmq/bioinfo/EmboSC/GSE36552/featurecounts/all_samples_gene_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(emb.co$`2013nature`)[-1:-6] <- str_split_fixed(str_split_fixed(colnames(emb.co$`2013nature`)[-1:-6], "star.", 2)[,2], "_", 2)[,1]
emb.co$`2013nature` <- emb.co$`2013nature`[,c(colnames(emb.co$`2013nature`)[1:6], s.anno$`2013nature`$ID)]
# 2016 Cell
emb.co$`2016cell` <- read.table("/home/cmq/bioinfo/EmboSC/E-MATB-3929/featurecounts/all_samples_gene_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(emb.co$`2016cell`)[-1:-6] <- str_split_fixed(str_split_fixed(colnames(emb.co$`2016cell`)[-1:-6], "\\.", 13)[,12], "Align", 2)[,1]

### >>> 4. transform geneid into gene symbol
for (i in seq(1, length(emb.co))) {
  emb.co[[i]] <- merge(emb.co[[i]], ge.anno[,5:6])
  emb.co[[i]] <- emb.co[[i]][emb.co[[i]]$Geneid %in% ge.anno[!ge.anno$Symbol %in% ge.anno[duplicated(ge.anno$Symbol),]$Symbol,]$Geneid,]
  rownames(emb.co[[i]]) <- emb.co[[i]]$Symbol
  emb.co[[i]] <- emb.co[[i]][,-grep("Geneid|Symbol", colnames(emb.co[[i]]))]
}; rm(i)
emb.tpm$cb[c("DPRX", "AQP3", "BAK1"),grep("ave|med", colnames(emb.tpm$cb))]

### >>> 5. create seurat object
emb.sr <- list()
for (i in seq(1, length(s.anno))) {
  rownames(s.anno[[i]]) <- s.anno[[i]]$ID
  emb.sr[[i]] <- CreateSeuratObject(counts = emb.co[[i]][,-1:-5], project = names(s.anno)[i], meta.data = s.anno[[i]]) %>% 
    PercentageFeatureSet(pattern = "^MT-", col.name = "pect.mt")
  emb.sr[[i]]$Dataset <- names(s.anno)[i]
  # raw quality control 
  #pdf(paste("graphs/combined_with_GSE36552/", names(s.anno)[i], "_raw_vlnplot.pdf", sep = ""), height = 5, width = 15)
  print(VlnPlot(emb.sr[[i]], features = c("nFeature_RNA", "nCount_RNA", "pect.mt"), ncol = 3, group.by = "Dataset"))
  #dev.off()
  # filter cells
  qc.ncount <- isOutlier(emb.sr[[i]]$nCount_RNA, log=TRUE, type="lower", nmads = 4)
  qc.nfeature <- isOutlier(emb.sr[[i]]$nFeature_RNA, log=TRUE, type="lower", nmads = 4)
  qc.mito <- isOutlier(emb.sr[[i]]$pect.mt, type="higher", nmads = 4)
  emb.sr[[i]] <- emb.sr[[i]][,!(qc.ncount | qc.nfeature | qc.mito)]
  # quality control plot after filtering cells
  #pdf(paste("graphs/combined_with_GSE36552/", names(s.anno)[i], "_filtered_vlnplot.pdf", sep = ""), height = 5, width = 15)
  print(VlnPlot(emb.sr[[i]], features = c("nFeature_RNA", "nCount_RNA", "pect.mt"), ncol = 3, group.by = "Dataset"))
  #dev.off()
  names(emb.sr)[i] <- names(s.anno)[i]
}; rm(i)

### >>> 6. merge two datasets
emb.sr.m <- merge(emb.sr$`2013nature`, emb.sr$`2016cell`)
emb.sr.m$Source <- paste(emb.sr.m$Stage, emb.sr.m$Dataset, sep = ".")
# normalization, find variable genes and scale data
emb.sr.m %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = rownames(emb.sr.m)) -> emb.sr.m
# pca and define dimensions
emb.sr.m <- RunPCA(emb.sr.m, features = VariableFeatures(object = emb.sr.m))
#pdf("graphs/combined_with_GSE36552/combeined_pca_group_by_dataset.pdf", height = 5, width = 6.5)
DimPlot(emb.sr.m, reduction = "pca", group.by = "Dataset")
#dev.off()
#pdf("graphs/combined_with_GSE36552/combeined_pca_group_by_dataset_and_stage.pdf", height = 5, width = 6.5)
DimPlot(emb.sr.m, reduction = "pca", group.by = "Source")
#dev.off()
# output merged seurat object
#saveRDS(emb.sr.m, "/home/cmq/bioinfo/EmboSC/E-MATB-3929/rds/E-MATB-3929_GSE36552_merged_seurat_object.rds")

### >>> 7. create metadata for edger
emb.meta <- rbind(s.anno$`2013nature`, s.anno$`2016cell`)
emb.meta <- emb.meta[emb.meta$ID %in% emb.sr.m@meta.data$ID,]
emb.meta <- emb.meta[!emb.meta$Stage %in% c("LB", "Oocyte"),] %>% arrange(Stage)
emb.meta <- emb.meta[c(1514:1516, 1:1513),]
emb.meta$CellName <- ave(emb.meta$Stage, emb.meta$Stage, FUN = function(i) paste0(i, '_', seq_along(i)))
rownames(emb.meta) <- emb.meta$CellName
emb.meta %>% mutate(Dataset=case_when(str_detect(ID, regex("^SRR")) ~ "nature", str_detect(ID, regex("^ERR")) ~ "cell")) %>% 
  mutate(Source=paste(Dataset, Stage, sep = ".")) -> emb.meta

### >>> 8. create count table for edger
emb.co$cb <- merge(emb.co$`2016cell`[,colnames(emb.co$`2016cell`) %in% c(colnames(emb.co$`2016cell`)[1:5], emb.meta$ID)] %>% rownames_to_column("Geneid"), 
                   emb.co$`2013nature`[,colnames(emb.co$`2013nature`) %in% c(colnames(emb.co$`2013nature`)[1:5], emb.meta$ID)] %>% rownames_to_column("Geneid"))
rownames(emb.co$cb) <- emb.co$cb$Geneid
emb.co$cb <- emb.co$cb[,emb.meta$ID]
if (all(colnames(emb.co$cb)==emb.meta$ID)) {
  colnames(emb.co$cb) <- emb.meta$CellName
}

### >>> 9. create tpm table
emb.tpm <- list()
emb.tpm$cb <- CountToTpm(emb.co$cb, emb.co$`2016cell`[rownames(emb.co$cb),]$Length)
# zscore of tpm
emb.tpm$cb.zs <- as.data.frame(t(apply(emb.tpm$cb, 1, function(x){(x - mean(x))/(sd(x))})))
# average of tpm
emb.tpm$cb %>% mutate(Zygote_ave = rowMeans(emb.tpm$cb[, grep("Zygote", colnames(emb.tpm$cb))]), 
                      `2Cell_ave` = rowMeans(emb.tpm$cb[, grep("2Cell", colnames(emb.tpm$cb))]), 
                      `4Cell_ave` = rowMeans(emb.tpm$cb[, grep("4Cell", colnames(emb.tpm$cb))]), 
                      E3_ave = rowMeans(emb.tpm$cb[, grep("E3", colnames(emb.tpm$cb))]), 
                      E4_ave = rowMeans(emb.tpm$cb[, grep("E4", colnames(emb.tpm$cb))]), 
                      E5_ave = rowMeans(emb.tpm$cb[, grep("E5", colnames(emb.tpm$cb))]), 
                      E6_ave = rowMeans(emb.tpm$cb[, grep("E6", colnames(emb.tpm$cb))]), 
                      E7_ave = rowMeans(emb.tpm$cb[, grep("E7", colnames(emb.tpm$cb))]),
                      Zygote_med = rowMedians(as.matrix(emb.tpm$cb[, grep("Zygote", colnames(emb.tpm$cb))])), 
                      `2Cell_med` = rowMedians(as.matrix(emb.tpm$cb[, grep("2Cell", colnames(emb.tpm$cb))])), 
                      `4Cell_med` = rowMedians(as.matrix(emb.tpm$cb[, grep("4Cell", colnames(emb.tpm$cb))])), 
                      E3_med = rowMedians(as.matrix(emb.tpm$cb[, grep("E3", colnames(emb.tpm$cb))])), 
                      E4_med = rowMedians(as.matrix(emb.tpm$cb[, grep("E4", colnames(emb.tpm$cb))])), 
                      E5_med = rowMedians(as.matrix(emb.tpm$cb[, grep("E5", colnames(emb.tpm$cb))])), 
                      E6_med = rowMedians(as.matrix(emb.tpm$cb[, grep("E6", colnames(emb.tpm$cb))])), 
                      E7_med = rowMedians(as.matrix(emb.tpm$cb[, grep("E7", colnames(emb.tpm$cb))]))) -> emb.tpm$cb

### >>> 10. run edger for differential expression analysis
# zygote vs E3
zy.vs.e3.ge <- edgeR(emb.co$cb, emb.meta, "Zygote", "E3", 1, 0.05)
# 2cell vs E3
e3.vs.c2.ge <- edgeR(emb.co$cb, emb.meta, "2Cell", "E3", 1, 0.05)
# 4cell vs E3
e3.vs.c4.ge <- edgeR(emb.co$cb, emb.meta, "4Cell", "E3", 1, 0.05)
# 8cell vs E4
e4.vs.e3.ge <- edgeR(emb.co$cb, emb.meta, "E3", "E4", 1, 0.05)
# 8cell vs E5
e5.vs.e3.ge <- edgeR(emb.co$cb, emb.meta, "E3", "E5", 1, 0.05)
# 8cell vs E6
e6.vs.e3.ge <- edgeR(emb.co$cb, emb.meta, "E3", "E6", 1, 0.05)
# 8cell vs E7
e7.vs.e3.ge <- edgeR(emb.co$cb, emb.meta, "E3", "E7", 1, 0.05)

### >>> 11. 8Cell specific genes
# define
c8.spe.ge <- Reduce(intersect, list(rownames(zy.vs.e3.ge$degs.down), rownames(e3.vs.c2.ge$degs.up), 
                                    rownames(e3.vs.c4.ge$degs.up), rownames(e4.vs.e3.ge$degs.down), 
                                    rownames(e5.vs.e3.ge$degs.down), rownames(e6.vs.e3.ge$degs.down), 
                                    rownames(e7.vs.e3.ge$degs.down)))
c8.spe.ge <- emb.tpm$cb[c8.spe.ge,]
c8.spe.ge <- subset(c8.spe.ge, E3_ave >= 5)
# visualize
pd <- emb.tpm$cb.zs[rownames(c8.spe.ge),]
#pdf("graphs/combined_with_GSE36552/8Cell_specific_genes_zscore_of_TPM.pdf", height = 8, width = 16)
Heatmap(as.matrix(pd),
        # colour and legend name setting
        col = colorRamp2(c(-1, 0, 2), c("#1078d4", "#ffffff", "#f03039")), 
        name = "ZScore", use_raster = T, 
        # clustering setting
        cluster_rows = T, show_row_dend = F, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        cluster_columns = F, show_column_dend = T, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        # names, rotation and position of columns and row names
        show_column_names = F, column_names_rot = 45, column_names_side = "bottom",
        show_row_names = F, row_names_rot = 0, row_names_side = "right", 
        # column split and annotation
        column_order = colnames(pd), 
        column_split = factor(str_split_fixed(colnames(pd), "_", 2)[,1], 
                              levels = c("Zygote", "2Cell", "4Cell", "E3", "E4", "E5", "E6", "E7")), 
        top_annotation = HeatmapAnnotation(Stage = c(rep("Zygote", length(grep("Zygote", colnames(pd)))), 
                                                     rep("2Cell", length(grep("2Cell", colnames(pd)))), 
                                                     rep("4Cell", length(grep("4Cell", colnames(pd)))), 
                                                     rep("E3", length(grep("E3_[0-9]", colnames(pd)))), 
                                                     rep("E4", length(grep("E4_[0-9]", colnames(pd)))),
                                                     rep("E5", length(grep("E5_[0-9]", colnames(pd)))), 
                                                     rep("E6", length(grep("E6_[0-9]", colnames(pd)))),
                                                     rep("E7", length(grep("E7_[0-9]", colnames(pd))))), 
                                           col = list(Stage = c("Zygote"="#f0b60e", "2Cell"="#c5c4c1", "4Cell"="#39d104",  
                                                                "E3"="#A6CEE3", "E4"="#1F78B4", "E5"="#CAB2D6", 
                                                                "E6"="#6A3D9A", "E7"="#FB9A99")), 
                                           simple_anno_size = unit(0.25, "cm")), 
        # size of heatmap, also including heatmap_height and heatmap_width
        width = unit(15, "cm"), height = unit(10, "cm"))
#dev.off
# output 
write.table(c8.spe.ge %>% rownames_to_column("SYMBOL"), "tables/combined_with_GSE36552/8Cell_specific_genes.txt", 
            col.names = T, row.names = F, sep = "\t", quote = F)


### ==================================================
### 3rd. Analysis of 10x prEpiSC scRNA-seq --> Fig1E-I
### ==================================================

### >>> 1. annotate human embryo cells with cell annotation from CSC2021 
# load cell annotation from CSC2021
cell.anno <- read.table("/home/cmq/bioinfo/EmboSC/E-MATB-3929/tables/cell_annotation_fromCSC_2016P.txt", header = T, stringsAsFactors = F, sep = "\t")
cell.anno <- subset(cell.anno, Dataset %in% c("Petropoulos2016", "Yan2013"))
cell.anno <- cell.anno[,c(15,21)]
colnames(cell.anno)[1] <- "ID"
# add cell annotation to seurat object
emb.sr.m.an <- emb.sr.m
tmp <- merge(emb.sr.m.an@meta.data, cell.anno, all.x = T)
tmp %>% 
  mutate(UMAP.cluster=case_when(Stage=="Oocyte" ~ "Oocyte", 
                                Stage=="Zygote" ~ "Zygote",
                                Stage=="2Cell" ~ "2Cell", 
                                Stage=="4Cell" ~ "4Cell", 
                                TRUE ~ UMAP.cluster)) -> tmp
rownames(tmp) <- tmp$ID
tmp <- tmp[colnames(emb.sr.m.an),]
if(all(colnames(emb.sr.m.an)==rownames(tmp))){
  emb.sr.m.an@meta.data <- tmp
}
emb.sr.m.an <- subset(emb.sr.m.an, ID %in% tmp[!tmp$UMAP.cluster%in%c(NA, "B1.EPI", "EPI.early_TE", "EPI.PrE", "EPI.PrE.TE", "PrE.TE") & !tmp$Source%in%c("E3.2013nature", "LB.2013nature", "E4.2013nature"),]$ID)
# dimension reduction
emb.sr.m.an %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = rownames(emb.sr.m.an)) -> emb.sr.m.an
emb.sr.m.an <- RunPCA(emb.sr.m.an)
ElbowPlot(emb.sr.m.an, ndims = 50)
emb.sr.m.an <- RunUMAP(emb.sr.m.an, reduction = "pca", dims = 1:11) %>% RunTSNE(reduction = "pca", dims = 1:11) %>% 
  FindNeighbors(reduction = "pca", dims = 1:11) %>% FindClusters(resolution = 0.8)
DimPlot(emb.sr.m.an, reduction = "umap", group.by = c("UMAP.cluster", "Source"))
#saveRDS(emb.sr.m.an, "rds/E-MATB-3929_GSE36552_merged_and_annotated_with_2021CSC_seurat_object.rds")

### >>> 2. gene set preparation
ge.set <- list()
ge.set$e3.spe <- read.table("/home/cmq/bioinfo/EmboSC/E-MATB-3929/tables/combined_with_GSE36552/8Cell_specific_genes.txt", header = T, stringsAsFactors = F)$SYMBOL
ge.set$ddr <- read.table("/home/cmq/document/gene_list/dna_repair/Human_DNA_Repair_Genes.txt")
colnames(ge.set$ddr) <- "SYMBOL"
ge.set$ddr %>% mutate(SYMBOL_alias=case_when(SYMBOL=="H2AX" ~ "H2AFX", SYMBOL=="TFIIH" ~ "GTF2H1", SYMBOL=="MRE11A" ~ "MRE11", TRUE ~ SYMBOL)) -> ge.set$ddr

### >>> 3. create seurat object
nc.sr <- Read10X("/home/data/wanglab/YMZ_hTPSCs_2nd_batch/naivecell_scrna/rawdata/filtered_feature_bc_matrix/") %>% as.data.frame()
rownames(nc.sr)[rownames(nc.sr)=="H3.Y"] <- "H3Y1"
nc.sr <- CreateSeuratObject(counts = nc.sr, project = "Naive Cell") %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "pect.mt")

### >>> 4. filter cells
# raw cells qc
#pdf(paste(seurat.out, "raw_cells_quality_control.pdf", sep = ""), height = 4, width = 8)
VlnPlot(nc.sr, features = c("nFeature_RNA", "nCount_RNA", "pect.mt"), ncol = 3)
#dev.off()
# filter cells
nc.sr <- subset(nc.sr, nFeature_RNA>=2000 & nCount_RNA>=5000 & pect.mt<20)
# filtered cells qc
#pdf(paste(seurat.out, "filtered_cells_quality_control.pdf", sep = ""), height = 4, width = 8)
VlnPlot(nc.sr, features = c("nFeature_RNA", "nCount_RNA", "pect.mt"), ncol = 3)
#dev.off()

### >>> 5. normalization, find hvg, run pca and dimension reduction
nc.sr %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(nc.sr)) -> nc.sr
nc.sr <- RunPCA(nc.sr, features = VariableFeatures(object = nc.sr))
ElbowPlot(nc.sr, ndims = 50)
# dimension reduction and clustering
nc.sr <- RunUMAP(nc.sr, reduction = "pca", dims = 1:20) %>% RunTSNE(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>% FindClusters(resolution = 0.1)

### >>> 6. cell type annotation
# find markers
nc.sr.marker <- FindAllMarkers(nc.sr, only.pos = T)
# visualize pluripotent and totipotent marker genes
VlnPlot(nc.sr, features = c("POU5F1", "SOX2", "NANOG", "TFCP2L1", "KLF4", "KLF17", "ZSCAN4", "LEUTX", "TPRX1", "H3Y1"), stack = T, flip = T)
# add initiative cell type information to metadata
DimPlot(nc.sr, reduction = "umap", label = TRUE, repel = TRUE, group.by = c("seurat_clusters"))
DimPlot(nc.sr, reduction = "tsne", label = TRUE, repel = TRUE, group.by = c("seurat_clusters"))
nc.sr@meta.data$Cell <- rownames(nc.sr@meta.data)
nc.sr@meta.data$Data <- "CellLine"
nc.sr@meta.data %>% mutate(RawType=case_when(seurat_clusters=="0" ~ "Naive Cells",
                                             seurat_clusters=="1" ~ "TELCs",
                                             seurat_clusters=="2" ~ "8CLCs")) -> nc.sr@meta.data
nc.sr@meta.data$Group <- paste("Cluster", nc.sr@meta.data$seurat_clusters, sep = " ")

### >>> 6. define 8CLCs and TELCs via clustering imputed 8CLCs, TELCs and human embryo
# load embryo seurat object
emb.sr <- readRDS("/home/cmq/bioinfo/EmboSC/E-MATB-3929/rds/E-MATB-3929_GSE36552_merged_and_annotated_with_2021CSC_seurat_object.rds")
emb.sr@meta.data$Data <- "Embryo"
emb.sr@meta.data %>%
  mutate(RawType=case_when(UMAP.cluster=="early_TE" ~ "Early TE",
                           UMAP.cluster=="medium_TE" ~ "Medium TE",
                           UMAP.cluster=="late_TE" ~ "Late TE",
                           UMAP.cluster=="EightCells" ~ "8Cell",
                           UMAP.cluster=="B1_B2" ~ "B1 B2",
                           TRUE ~ UMAP.cluster)) -> emb.sr@meta.data
emb.sr@meta.data$Group <- emb.sr@meta.data$RawType
emb.sr@meta.data$CellType <- emb.sr@meta.data$RawType
# subset 8CLCs and TELCs
nc.sr.lc <- subset(nc.sr, RawType %in% c("8CLCs", "TELCs"))
nc.sr.lc <- NormalizeData(nc.sr.lc) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# combine raw TELCs, 8CLCs and embryo
ol.hvg.int <- SelectIntegrationFeatures(object.list = list(nc.sr.lc, emb.sr), nfeatures = 2000)
lc.emb.int.an <- FindIntegrationAnchors(object.list = list(nc.sr.lc, emb.sr), anchor.features = ol.hvg.int, k.filter = NA)
lc.emb.int.sr <- IntegrateData(anchorset = lc.emb.int.an, k.weight = 65)
DefaultAssay(lc.emb.int.sr) <- "integrated"
lc.emb.int.sr <- ScaleData(lc.emb.int.sr, verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
lc.emb.int.sr <- RunUMAP(lc.emb.int.sr, reduction = "pca", dims = 1:15) %>%
  RunTSNE(reduction = "pca", dims = 1:15) %>% FindNeighbors(reduction = "pca", dims = 1:15) %>% FindClusters(resolution = 0.8)
lc.emb.int.sr@meta.data %>%
  mutate(Group=factor(Group, levels = c("Oocyte", "Zygote", "2Cell", "4Cell", "8Cell", "Morula",
                                        "B1 B2", "EPI", "PrE", "Early TE", "Medium TE", "Late TE", 
                                        "Cluster 1", "Cluster 2"))) -> lc.emb.int.sr@meta.data
# define TELCs and 8CLCs
table(lc.emb.int.sr@meta.data[,c(9,6)])
lc.emb.int.sr@meta.data %>%
  mutate(CellType=case_when(RawType=="TELCs" & seurat_clusters%in%c("0", "1", "3", "6", "7") ~ "TELCs",
                            RawType=="8CLCs" & seurat_clusters%in%c("9") ~ "8CLCs",
                            RawType=="TELCs" & seurat_clusters%in%c("4", "10") ~ "Intermediate",
                            RawType=="8CLCs" & seurat_clusters%in%c("4", "5", "10") ~ "Intermediate",
                            TRUE ~ RawType)) -> lc.emb.int.sr@meta.data
lc.emb.int.sr@meta.data %>%
  mutate(CellType=factor(CellType, levels = c("Oocyte", "Zygote", "2Cell", "4Cell", "8Cell", "Morula", "B1 B2", "EPI",
                                              "PrE", "Early TE", "Medium TE", "Late TE", "Intermediate", "8CLCs", "TELCs"))) -> lc.emb.int.sr@meta.data
Idents(lc.emb.int.sr) <- lc.emb.int.sr@meta.data$CellType
# add cell type to metadata
nc.sr@meta.data %>%
  mutate(CellType=case_when(Cell %in% lc.emb.int.sr@meta.data[lc.emb.int.sr@meta.data$CellType=="Intermediate",]$Cell ~ "Intermediate",
                            Cell %in% lc.emb.int.sr@meta.data[lc.emb.int.sr@meta.data$CellType=="8CLCs",]$Cell ~ "8CLCs",
                            Cell %in% lc.emb.int.sr@meta.data[lc.emb.int.sr@meta.data$CellType=="TELCs",]$Cell ~ "TELCs",
                            TRUE ~ "Naive Cells")) -> nc.sr@meta.data
# rename idents
Idents(nc.sr) <- nc.sr$CellType
# visualization of combined cell line and embryo
# tsne dimension plot
#pdf(paste(seurat.out, "combined_cell_line_and_embryo_tsne_plot.pdf", sep = ""), height = 5, width = 18.5)
print(DimPlot(lc.emb.int.sr, reduction = "tsne", label = TRUE, repel = TRUE, group.by = c("Group", "seurat_clusters", "CellType"),
              cols = c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C", "#B2DF8A", "#33A02C",
                       "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#09E9E2", "#04918D",
                       "#F1F144", "#979704", "#9C9C9B", "#525252", "#848587", "#bec0c2", "#266c3f", "#70a081")))
#dev.off()
# visualization pluripotent and totipotent genes
#pdf(paste(seurat.out, "combined_cell_line_and_embryo_expression_level_of_pluripotent_and_totipotent_markers_tsne_featureplot_V1.pdf", sep = ""), height = 10, width = 30)
FeaturePlot(lc.emb.int.sr, features = c("NANOG", "KLF4", "KLF17", "ZSCAN4", "LEUTX", "TPRX1"), ncol = 3, order = T, reduction = "tsne")
#dev.off()
# visualization of marker genes in three groups via dotplot
pd.ge <- c("POU5F1","SOX2","NANOG","TFCP2L1","KLF4","DPPA5","IFITM3","ESRG","PRDM14","UTF1","NODAL",
           "IFI16","KRT8","KRT18","FABP5","CD24","HAND1","GATA2","GATA3","CDX2","DLX3","CCNA1","H3Y1",
           "MBD3L2","MBD3L3","KLF17","SLC34A2","GDF15","ZSCAN4","LEUTX","TPRX1","TRIM43","PRAMEF12")
pd.sr <- merge(lc.emb.int.sr, subset(nc.sr, Group=="Cluster 0"))
pd.sr@meta.data %>%
  mutate(Type.pd=case_when(Group=="Cluster 0" ~ "Naive Cells",
                           CellType=="Intermediate" & Group%in%c("Cluster 1", "Cluster 2") ~ "Intermediate",
                           TRUE ~ CellType)) %>%
  mutate(Type.pd=factor(Type.pd, levels=c("Oocyte", "Zygote", "2Cell", "4Cell", "8Cell", "Morula", "B1 B2", "EPI",
                                          "PrE", "Early TE", "Medium TE", "Late TE", "Naive Cells", "Intermediate", "8CLCs", "TELCs")))-> pd.sr@meta.data
#pdf(paste(seurat.out, "combined_cell_line_and_embryo_expression_level_of_marker_genes_dotplot_v2.pdf", sep = ""), height = 5, width = 12)
DotPlot(pd.sr, features = pd.ge, assay = "RNA", group.by = "Type.pd",
        cols = c("#eaeaea", "#e74c3c"), scale = T, col.min = -1) + RotatedAxis() +
  scale_x_discrete(labels = pd.ge) +
  labs(x = "Marker genes", y = "Cluster") +
  theme(axis.title.x = element_text(face="plain", color="#000000", size=12, angle=0),
        axis.title.y = element_text(face="plain", color="#000000", size=12, angle=90),
        axis.text.x = element_text(face="plain", color="#000000", size=11, angle=45),
        axis.text.y = element_text(face="plain", color="#000000", size=11, angle=0),
        legend.box = "horizontal")
#dev.off()


### ==========================================================================
### 4th.Comparison between prEpiSC and naive or conventional hESC --> Fig1B-D
### ==========================================================================

### >>> 1. load gene count
ge.co <- list()
# prEpiSC
ge.co$prepi <- read.table("/home/cx/work/hPSC/wang-PSC/bulk/results/featurecounts/all_samples_gene_count_matrix.txt", header = T, stringsAsFactors = F)
colnames(ge.co$prepi)[-1:-6] <- str_split_fixed(str_split_fixed(colnames(ge.co$prepi)[-1:-6], "star.", 2)[,2], "_1", 2)[,1]
# ctrl samples in hTPSC bulk RNA-seq data
ge.co$ctrl.b1 <- read.table("/home/data/wanglab/YMZ_hTPSCs_1st_batch/analysis/rnaseq/results/featurecounts/top10_gene/all_samples_gene_count_name_matrix.txt", header = T, stringsAsFactors = F, row.names = 1)
colnames(ge.co$ctrl.b1)[-1:-5] <- c("Ctrl", "hTPSC_1", "hTPSC_2", "hTPSC_3")
ge.co$ctrl.b2 <- read.table("/home/data/wanglab/YMZ_hTPSCs_2nd_batch/bulk_rna_seq/analysis/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt", header = T, stringsAsFactors = F, row.names = 1)
colnames(ge.co$ctrl.b2)[-1:-5] <- c("hTPSC_b2_1", "Ctrl_b2_1", "hTPSC_b2_2", "Ctrl_b2_2")
# psc
geo <- c("GSE135989", "GSE150772", "GSE147839", "GSE153212", "GSE93241", "GSE87452", "E-MTAB-4461", "E-MTAB-2857")
ge.co$psc <- list()
for (i in seq(1, length(geo))) {
  if (geo[i]=="GSE93241") {
    ge.co$psc[[i]] <- read.table(paste("/home2/data/publicdata/", geo[i], "/analysis/results/featurecounts/top10_gene/all_samples_gene_count_name_matrix.txt", sep = ""),
                                 header = T, stringsAsFactors = F, row.names = 1)
    colnames(ge.co$psc[[i]])[-1:-5] <- str_split_fixed(str_split_fixed(colnames(ge.co$psc[[i]])[-1:-5], "gene.", 2)[,2], "_", 2)[,1]
    names(ge.co$psc)[i] <- geo[i]
  }else {
    ge.co$psc[[i]] <- read.table(paste("/home2/data/publicdata/", geo[i], "/analysis/results/featurecounts/gene/all_samples_gene_count_name_matrix.txt", sep = ""),
                                 header = T, stringsAsFactors = F, row.names = 1)
    colnames(ge.co$psc[[i]])[-1:-5] <- str_split_fixed(str_split_fixed(colnames(ge.co$psc[[i]])[-1:-5], "gene.", 2)[,2], "_", 2)[,1]
    names(ge.co$psc)[i] <- geo[i]
  }
}; rm(i)
ge.co$psc <- Reduce(cbind, ge.co$psc)
ge.co$psc <- ge.co$psc[,c(1:5, grep("Chr|Start|End|Strand|Length", colnames(ge.co$psc), invert = T))]
# transform geneid to gene symbol
ge.anno <- readRDS("/home/cmq/bioinfo/project-yangmz/chimera/icm_like/emb_psc/rawdata/metadata/ENSEMBL_hg38_v34_Gene_anno.rds")
ge.anno %>% rownames_to_column("ENSEMBL") -> ge.anno
for (i in 1) {
  ge.co[[i]] <- merge(ge.co[[i]], ge.anno[,c(6,7)])
  ge.co[[i]] <- ge.co[[i]][!ge.co[[i]]$Symbol %in% ge.anno[duplicated(ge.anno$Symbol),]$Symbol,]
  rownames(ge.co[[i]]) <- ge.co[[i]]$Symbol
  ge.co[[i]] <- ge.co[[i]][,-grep("Geneid|Symbol", colnames(ge.co[[i]]))]
}; rm(i)

### >>> 2. all psc sample metadata
sp.meta <- list()
for (i in seq(1, length(geo))) {
  sp.meta[[i]] <- read.table(paste("/home2/data/publicdata/", geo[i], "/analysis/sample_meta.txt", sep = ""))
  colnames(sp.meta[[i]]) <- c("Sample", "SRR", "Author")
}; rm(i)
sp.meta <- Reduce(rbind, sp.meta)
# change psc sample name in count table
ge.co$psc <- ge.co$psc[, c("Chr", "Start", "End", "Strand", "Length", sp.meta$SRR)]
if(all(colnames(ge.co$psc)[-1:-5]==sp.meta$SRR)){
  colnames(ge.co$psc)[-1:-5] <- sp.meta$Sample
}

### >>> 3. calculate tpm
ge.tpm <- list()
for (i in seq(1, length(ge.co))) {
  ge.tpm[[i]] <- CountToTpm(ge.co[[i]][,-1:-5], ge.co[[i]]$Length)
  names(ge.tpm)[i] <- names(ge.co)[i]
}; rm(i)

### >>> 4.  DEA of genes between prEpiSC and primed H9 hESC
# count table
dea.co <- list()
dea.co$prEpiSC.vs.primedH9 <- cbind(ge.co$prepi[intersect(rownames(ge.co$prepi), rownames(ge.co$psc)),-8:-9],
                                    ge.co$psc[intersect(rownames(ge.co$prepi), rownames(ge.co$psc)),c(1:5, grep("^Primed_H9", colnames(ge.co$psc)))])
dea.co$prEpiSC.vs.primedH9 <- dea.co$prEpiSC.vs.primedH9[,-8:-12]
# metadata
dea.meta <- list()
dea.meta$prEpiSC.vs.primedH9 <- data.frame(SampleID = colnames(dea.co$prEpiSC.vs.primedH9)[-1:-5],
                                           SampleGroup = c(rep("WC", 2), rep("Primed_H9", 3)),
                                           Replicates = c("R1", "R2", "R1", "R2", "R3"),
                                           row.names = colnames(dea.co$prEpiSC.vs.primedH9)[-1:-5])
prEpiSC.vs.primedH9.ge <- Deseq2.Ge(dea.co$prEpiSC.vs.primedH9[,-1:-5], dea.meta$prEpiSC.vs.primedH9, "WC", "Primed_H9",
                                    1, 0.05, "results/R/tables/deseq2/prEpiSC_vs_primedH9_ge", dea.tpm$prEpiSC.vs.primedH9)
# tpm table
dea.tpm <- list()
dea.tpm$prEpiSC.vs.primedH9 <- CountToTpm(dea.co$prEpiSC.vs.primedH9[,-1:-5], dea.co$prEpiSC.vs.primedH9$Length)
# run Deseq2
prEpiSC.vs.primedH9.ge <- Deseq2.Ge(dea.co$prEpiSC.vs.primedH9[,-1:-5], dea.meta$prEpiSC.vs.primedH9, "WC", "Primed_H9",
                                    1, 0.05, dea.tpm$prEpiSC.vs.primedH9)

### >>> 5. create differential expressed genes table for visualization 
pd <- prEpiSC.vs.primedH9.ge$res[,c(1:3,7,13:17)]
pd %>% mutate(Log10MeanTPM=log10(rowMeans(pd[,grep("TPM", colnames(pd))]))) %>%
  mutate(color = case_when((log2FoldChange >= 1 & padj <= 0.05) ~ "Up-regulated",
                           (abs(log2FoldChange) < 1 | padj > 0.05) ~ "No-changed",
                           (log2FoldChange <= -1 & padj <= 0.05) ~ "Down-regulated")) -> pd
pd <- pd[c(grep("No-changed", pd$color), grep("Down-regulated", pd$color), grep("Up-regulated", pd$color)),]
pd.anno <- pd[pd$SYMBOL %in% c("LEUTX","ZSCAN4","TPRX1","ZSCAN5B","TRIM43","TRIM43B","MBD3L2","DUXA","UBTFL1",
                               "KDM4E","PRAMEF1","PRAMEF12","RFPL2","RFPL4A","SLC34A2","CCNA1"),]
# visualization of MA plot
#pdf("results/R/graphs/prEpiSC_vs_PrimedH9/prEpiSC_vs_Primed_H9_hESC_MA_plot.pdf", height = 10, width = 10)
pd %>%
  ggplot(aes(x = Log10MeanTPM, y = log2FoldChange)) +
  geom_point(aes(color = color), size = 1.5, alpha = 0.5) +
  geom_label_repel(aes(label = SYMBOL), data = pd.anno,
                   box.padding = 0.1, segment.color = "#130f40",
                   label.size = 0.1,
                   size = 2.5, color = "#000000") +
  scale_color_manual(values = c(`Up-regulated` = "#d31d00", `Down-regulated` = "#0043d3", `No-changed` = "#b2bec3")) +
  geom_hline(yintercept = c(-1,1), color = "#000000", linetype = "longdash", size = 1, alpha = 0.5) +
  xlab("log10(Mean TPM)") + ylab("log2(Fold of Change)") +
  theme_classic() +
  labs(color = "Color") +
  theme(legend.position = "bottom", legend.box="vertical")
#dev.off()

### >>> 6. visualize expression level of DUX4 related genes via heatmap
# load dux4 genes
dux4.ge <- read.table("/home/cmq/bioinfo/project-yangmz/chimera/icm_like/emb_psc/rawdata/dux4_related_gene.txt", header = F, stringsAsFactors = F)
# load embryo gene TPM table
ge.tpm$emb2 <- read.table("/home/cmq/bioinfo/EmboSC/E-MATB-3929/tables/combined_with_GSE36552/all_gene_tpm.txt", header = T, stringsAsFactors = F)
rownames(ge.tpm$emb2) <- ge.tpm$emb2$SYMBOL
# create table for plot
dea.tpm$prEpiSC.vs.primedH9 <- cbind(dea.tpm$prEpiSC.vs.primedH9, ge.tpm$emb2[rownames(dea.tpm$prEpiSC.vs.primedH9),grep("ave", colnames(ge.tpm$emb2))])
dea.tpm$prEpiSC.vs.primedH9.zs <- as.data.frame(t(apply(dea.tpm$prEpiSC.vs.primedH9[,-grep("ave", colnames(dea.tpm$prEpiSC.vs.primedH9))],
                                                        1, function(x){(x - mean(x))/(sd(x))})))
dea.tpm$emb.zs <- as.data.frame(t(apply(dea.tpm$prEpiSC.vs.primedH9[,grep("ave", colnames(dea.tpm$prEpiSC.vs.primedH9))],
                                        1, function(x){(x - mean(x))/(sd(x))})))
# visualization
pd <- dea.tpm$prEpiSC.vs.primedH9.zs[rownames(dea.tpm$prEpiSC.vs.primedH9.zs) %in% dux4.ge$V1,]
pd <- cbind(pd, dea.tpm$emb.zs[rownames(pd),])
pd <- subset(pd, rowSums(pd)!="NaN")
p1 <- Heatmap(as.matrix(pd[,-grep("ave", colnames(pd))]),
              # colour and legend name setting
              col = colorRamp2(c(-2, 0, 2), c("#1078d4", "#ffffff", "#f03039")),
              name = "CellLine: Zscore", use_raster = F,
              # clustering setting
              cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"),
              clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
              cluster_columns = F, show_column_dend = F, column_dend_height = unit(1, "cm"), column_dend_side = "top",
              clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
              # names, rotation and position of columns and row names
              show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
              show_row_names = F, row_names_rot = 0, row_names_side = "right",
              # annotation
              column_order = colnames(pd[,-grep("ave", colnames(pd))]),
              # size of heatmap, also including heatmap_height and heatmap_width
              width = unit(3, "cm"), height = unit(8, "cm"))
#pdf("results/R/graphs/prEpiSC_vs_PrimedH9/prEpiSC_and_Primed_H9_hESC_expression_level_of_DUX4_target_top100_genes_zscore_heatmap.pdf", height = 10, width = 10)
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
             show_row_names = F, row_names_rot = 0, row_names_side = "right",
             # annotation
             column_order = colnames(pd[,grep("ave", colnames(pd))]),
             right_annotation = rowAnnotation(mark = anno_mark(at = match(c("LEUTX","ZSCAN4","TPRX1","ZSCAN5B","TRIM43","TRIM43B","MBD3L2","DUXA","UBTFL1",
                                                                            "KDM4E","PRAMEF1","PRAMEF12","RFPL2","RFPL4A","SLC34A2","CCNA1"), rownames(pd)),
                                                               labels = c("LEUTX","ZSCAN4","TPRX1","ZSCAN5B","TRIM43","TRIM43B","MBD3L2","DUXA","UBTFL1",
                                                                          "KDM4E","PRAMEF1","PRAMEF12","RFPL2","RFPL4A","SLC34A2","CCNA1"))),
             # size of heatmap, also including heatmap_height and heatmap_width
             width = unit(4, "cm"), height = unit(8, "cm"))
#dev.off()

### >>> 7. visualize 8Cell marker genes in different hPSC
ge.tpm.cb <- cbind(ge.tpm$prepi[intersect(rownames(ge.tpm$prepi), rownames(ge.tpm$psc)),], ge.tpm$psc[intersect(rownames(ge.tpm$prepi), rownames(ge.tpm$psc)),])
pd.marker <- ge.tpm.cb[,grep("WC_001|WC_002|HNES|5iLA|naive|NK2|Primed_NK2|FTW|EPS", colnames(ge.tpm.cb))]
pd.marker <- pd.marker[,grep("iPS", colnames(pd.marker), invert = T)]
pd.marker.zs <- as.data.frame(t(apply(pd.marker, 1, function(x){(x - mean(x))/(sd(x))})))
pd <- pd.marker.zs[rownames(pd.marker.zs) %in% c("RFPL4A", "PRAMEF12", "PRAMEF2", "TRIM43", "RFPL4B", "TRIM43B", "CCNA1", "RP11-554D14.4", "MBD3L2", "PRAMEF1", "RFPL2", "LEUTX", "DUXA", "TPRX1", "SNAI1",
                                                 "RP11-432MB.17", "RP11-321E2.13", "ZSCAN4", "KDM4E", "SLC34A2", "FAM90A27P", "KLF17"),]
#pdf("results/R/graphs/expression_level_of_8Cell_marker_in_PSCs_zscore_of_raw_tpm.pdf", height = 10, width = 10)
Heatmap(as.matrix(pd),
        # colour and legend name setting
        col = colorRamp2(c(-2, 0, 2), c("#1078d4", "#ffffff", "#f03039")),
        name = "Zscore", use_raster = F,
        # clustering setting
        cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"),
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        cluster_columns = F, show_column_dend = F, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        # names, rotation and position of columns and row names
        show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
        show_row_names = T, row_names_rot = 0, row_names_side = "right",
        # annotation
        column_order = colnames(pd),
        # size of heatmap, also including heatmap_height and heatmap_width
        width = unit(12, "cm"), height = unit(9, "cm"))
#dev.off()

pd.marker <- cbind(ge.tpm$prepi[intersect(rownames(ge.tpm$prepi), rownames(ge.tpm$psc)),], ge.tpm$psc[intersect(rownames(ge.tpm$prepi), rownames(ge.tpm$psc)),])
pd.marker <- pd.marker[,grep("WC_001|WC_002|HNES|5iLA|naive|NK2|Primed_NK2|FTW", colnames(pd.marker))]
pd.marker.zs <- as.data.frame(t(apply(pd.marker, 1, function(x){(x - mean(x))/(sd(x))})))

### >>> 2. visualization
# zscore of raw tpm
pd <- pd.marker.zs[rownames(pd.marker.zs) %in% c("RFPL4A", "PRAMEF12", "PRAMEF2", "TRIM43", "RFPL4B", "TRIM43B", "CCNA1", "RP11-554D14.4", "MBD3L2", "PRAMEF1", "RFPL2", "LEUTX", "DUXA", "TPRX1", "SNAI1", 
                                                 "RP11-432MB.17", "RP11-321E2.13", "ZSCAN4", "KDM4E", "SLC34A2", "FAM90A27P", "KLF17"),]
#pdf("results/R/graphs/expression_level_of_8Cell_marker_in_PSCs_zscore_of_raw_tpm.pdf", height = 10, width = 10)
Heatmap(as.matrix(pd),
        # colour and legend name setting
        col = colorRamp2(c(-2, 0, 2), c("#1078d4", "#ffffff", "#f03039")), 
        name = "Zscore", use_raster = F, 
        # clustering setting
        cluster_rows = T, show_row_dend = T, row_dend_side = "left", row_dend_width = unit(1, "cm"), 
        clustering_distance_rows = "euclidean", clustering_method_rows = "complete",
        cluster_columns = F, show_column_dend = F, column_dend_height = unit(1, "cm"), column_dend_side = "top",
        clustering_distance_columns = "euclidean", clustering_method_columns = "complete",
        # names, rotation and position of columns and row names
        show_column_names = T, column_names_rot = 45, column_names_side = "bottom",
        show_row_names = T, row_names_rot = 0, row_names_side = "right", 
        # annotation
        column_order = colnames(pd), 
        # size of heatmap, also including heatmap_height and heatmap_width
        width = unit(12, "cm"), height = unit(9, "cm")) 
#dev.off()
