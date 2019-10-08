# generate results files for app.R
library(DESeq2)
library(hypeR)
library(qusage)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

HALLMARK <- msigdb_download_one(species="Homo sapiens", category="H")
laurenti <- qusage::read.gmt("Laurenti_Muschen.ensembl.gmt")

data(sample_table)
rownames(sample_table) <- sample_table$Name
sample_table$Name <- NULL
sample_table$SampleID <- NULL
sample_table$batch <- c(rep("nov_18", 27), rep("mar_19",32))
sample_table$stage <- c("w12","w12","w10","w10","w10","CS19",
                        "CS19","CS19","CB","CB","CB","CS16",
                        "CS16","CB","CB","CB","w13","w13",
                        "CS20","CS20","CS20","w8","w8","w8",
                        "BM","BM","BM","BM","BM","BM","CS19"
                        ,"CS19","BM","BM","BM","BM","BM","BM",
                        "CB","CB","CB","CS22","CS22","CS22","CS20"
                        ,"CS20","CS20","w13","w13","w13","CS23",
                        "CS23","CS23","CS22","CS22","CS22",
                        "w13","w13","w13")


sample_table$pools <- c("W12_13","W12_13","W8_10","W8_10","W8_10","CS19_23","CS19_23","CS19_23","CB",
                        "CB","CB","CS16","CS16","CB","CB","CB","W12_13","W12_13","CS19_23","CS19_23",
                        "CS19_23","W8_10","W8_10","W8_10","BM","BM","BM","BM","BM","BM","CS19_23",
                        "CS19_23","BM","BM","BM","BM","BM","BM","CB","CB","CB","CS19_23","CS19_23",
                        "CS19_23","CS19_23","CS19_23","CS19_23","W12_13","W12_13","W12_13","CS19_23",
                        "CS19_23","CS19_23","CS19_23","CS19_23","CS19_23","W12_13","W12_13","W12_13")
sample_table$age <- case_when(str_detect(sample_table$pools,"W|CS") ~ "fetal", TRUE ~ "adult"  )

# remove CB and CS
sample_table <- sample_table[!str_detect(sample_table$pools,"CB|CS16"),]

# set diff groups to be adult and fetal
sample_table$batch <- str_replace(sample_table$batch, "mar", "jun")
# previously identified outliers

sample_table <- sample_table[!rownames(sample_table) %in%
               c("CB2_E4_8W_P12","CB2_E6_BM_P12","CB2_B2_CS19_P12",
                                "CB2_A7_10w_P12","CB2_D9_CS20_P12"),]

data(big_counts)
big_counts <- big_counts[, which(colnames(big_counts) %in% rownames(sample_table))]

# HSC CS vs Adult
## DGE
hsc_cs_v_adult_coldata <- rownames_to_column(sample_table,"sample_id") %>%
  dplyr::filter(Phenotype == "CD34+CD38-CD45RA-") %>%
  dplyr::filter(pools == "CS19_23" | pools == "BM")

rownames(hsc_cs_v_adult_coldata) <- hsc_cs_v_adult_coldata$sample_id

counts <- big_counts[, colnames(big_counts) %in% rownames(hsc_cs_v_adult_coldata)]                           

dds <- DESeqDataSetFromMatrix(counts, colData = hsc_cs_v_adult_coldata, design = ~ pools)
dds$pools <- relevel(dds$pools, ref = "BM")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

ens_to_external <- readRDS('data/ens_to_external.RDS')

res <- data.frame(res) %>% rownames_to_column("ensembl_gene_id") 
my_colnames <- colnames(res)
res <- left_join(res,ens_to_external,by = "ensembl_gene_id") %>% dplyr::select("ensembl_gene_id","external_gene_name",
                                                                               my_colnames)
res <- res[!is.na(res$padj),]
res <- distinct(res)
res$external_gene_name <- case_when(duplicated(res$external_gene_name) ~ res$ensembl_gene_id,
                                    is.na(res$external_gene_name) ~ res$ensembl_gene_id,
                                    TRUE ~ as.character(res$external_gene_name))
rownames(res) <- res$external_gene_name
res$external_gene_name <- NULL

saveRDS(res, 'shiny/data/hsc_cs_v_adult.rds')

## Hallmark GSEA
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
hallmark_hyp_obj <- hypeR(geneList, HALLMARK, test = "kstest", fdr_cutoff = 0.01, do_plots = TRUE)
saveRDS(hallmark_hyp_obj, 'shiny/data/hsc_cs_v_adult_hallmark.rds')
## Laurenti GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$ensembl_gene_id
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
laurenti_hyp_obj <- hypeR(geneList, laurenti, test="kstest", fdr_cutoff=0.01, do_plots = TRUE)
saveRDS(laurenti_hyp_obj, 'shiny/data/hsc_cs_v_adult_laurenti.rds')
## GO
gene <- names(geneList)[abs(geneList) > 2]
my_ego2 <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

saveRDS(my_ego2, file = 'shiny/data/hsc_cs_v_adult_ego.rds')

# IL7R CS vs Adult

il7r_cs_v_adult_coldata <- rownames_to_column(sample_table,"sample_id") %>%
  dplyr::filter(Phenotype == "IL7R+Kit+") %>%
  dplyr::filter(pools == "CS19_23" | pools == "BM")

rownames(il7r_cs_v_adult_coldata) <- il7r_cs_v_adult_coldata$sample_id

counts <- big_counts[, colnames(big_counts) %in% rownames(il7r_cs_v_adult_coldata)]                           

dds <- DESeqDataSetFromMatrix(counts, colData = il7r_cs_v_adult_coldata, design = ~ pools)
dds$pools <- relevel(dds$pools, ref = "BM")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

# ens_to_external <- readRDS('data/ens_to_external.RDS')

res <- data.frame(res) %>% rownames_to_column("ensembl_gene_id") 
my_colnames <- colnames(res)
res <- left_join(res,ens_to_external,by = "ensembl_gene_id") %>% dplyr::select("ensembl_gene_id","external_gene_name",
                                                                               my_colnames)
res <- res[!is.na(res$padj),]
res <- distinct(res)
res$external_gene_name <- case_when(duplicated(res$external_gene_name) ~ res$ensembl_gene_id,
                                    is.na(res$external_gene_name) ~ res$ensembl_gene_id,
                                    TRUE ~ as.character(res$external_gene_name))
rownames(res) <- res$external_gene_name
res$external_gene_name <- NULL

saveRDS(res, 'shiny/data/il7r_cs_v_adult.rds')

## Hallmark GSEA
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
hallmark_hyp_obj <- hypeR(geneList, HALLMARK, test = "kstest", fdr_cutoff = 0.01, do_plots = TRUE)
saveRDS(hallmark_hyp_obj, 'shiny/data/il7r_cs_v_adult_hallmark.rds')
## Laurenti GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$ensembl_gene_id
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
laurenti_hyp_obj <- hypeR(geneList, laurenti, test="kstest", fdr_cutoff=0.01, do_plots = TRUE)
saveRDS(laurenti_hyp_obj, 'shiny/data/il7r_cs_v_adult_laurenti.rds')
## GO
gene <- names(geneList)[abs(geneList) > 2]
my_ego2 <- enrichGO(gene          = gene,
                    universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

saveRDS(my_ego2, file = 'shiny/data/il7r_cs_v_adult_ego.rds')

# proB CS vs Adult

proB_cs_v_adult_coldata <- rownames_to_column(sample_table,"sample_id") %>%
  dplyr::filter(Phenotype == "CD34+CD19+") %>%
  dplyr::filter(pools == "CS19_23" | pools == "BM")

rownames(proB_cs_v_adult_coldata) <- proB_cs_v_adult_coldata$sample_id

counts <- big_counts[, colnames(big_counts) %in% rownames(proB_cs_v_adult_coldata)]                           

dds <- DESeqDataSetFromMatrix(counts, colData = proB_cs_v_adult_coldata, design = ~ pools)
dds$pools <- relevel(dds$pools, ref = "BM")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

# ens_to_external <- readRDS('data/ens_to_external.RDS')

res <- data.frame(res) %>% rownames_to_column("ensembl_gene_id") 
my_colnames <- colnames(res)
res <- left_join(res,ens_to_external,by = "ensembl_gene_id") %>% dplyr::select("ensembl_gene_id","external_gene_name",
                                                                               my_colnames)
res <- res[!is.na(res$padj),]
res <- distinct(res)
res$external_gene_name <- case_when(duplicated(res$external_gene_name) ~ res$ensembl_gene_id,
                                    is.na(res$external_gene_name) ~ res$ensembl_gene_id,
                                    TRUE ~ as.character(res$external_gene_name))
rownames(res) <- res$external_gene_name
res$external_gene_name <- NULL

saveRDS(res, 'shiny/data/proB_cs_v_adult.rds')

## Hallmark GSEA
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
hallmark_hyp_obj <- hypeR(geneList, HALLMARK, test = "kstest", fdr_cutoff = 0.01, do_plots = TRUE)
saveRDS(hallmark_hyp_obj, 'shiny/data/proB_cs_v_adult_hallmark.rds')
## Laurenti GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$ensembl_gene_id
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
laurenti_hyp_obj <- hypeR(geneList, laurenti, test="kstest", fdr_cutoff=0.01, do_plots = TRUE)
saveRDS(laurenti_hyp_obj, 'shiny/data/proB_cs_v_adult_laurenti.rds')
## GO
gene <- names(geneList)[abs(geneList) > 2]
my_ego2 <- enrichGO(gene          = gene,
                    universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

saveRDS(my_ego2, file = 'shiny/data/proB_cs_v_adult_ego.rds')


# HSC W12-13 vs Adult

hsc_w12_13_v_adult_coldata <- rownames_to_column(sample_table,"sample_id") %>%
  dplyr::filter(Phenotype == "CD34+CD38-CD45RA-") %>%
  dplyr::filter(pools == "W12_13" | pools == "BM")

rownames(hsc_w12_13_v_adult_coldata) <- hsc_w12_13_v_adult_coldata$sample_id

counts <- big_counts[, colnames(big_counts) %in% rownames(hsc_w12_13_v_adult_coldata)]                           

dds <- DESeqDataSetFromMatrix(counts, colData = hsc_w12_13_v_adult_coldata, design = ~ pools)
dds$pools <- relevel(dds$pools, ref = "BM")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

# ens_to_external <- readRDS('data/ens_to_external.RDS')

res <- data.frame(res) %>% rownames_to_column("ensembl_gene_id") 
my_colnames <- colnames(res)
res <- left_join(res,ens_to_external,by = "ensembl_gene_id") %>% dplyr::select("ensembl_gene_id","external_gene_name",
                                                                               my_colnames)
res <- res[!is.na(res$padj),]
res <- distinct(res)
res$external_gene_name <- case_when(duplicated(res$external_gene_name) ~ res$ensembl_gene_id,
                                    is.na(res$external_gene_name) ~ res$ensembl_gene_id,
                                    TRUE ~ as.character(res$external_gene_name))
rownames(res) <- res$external_gene_name
res$external_gene_name <- NULL

saveRDS(res, 'shiny/data/hsc_w12_13_v_adult.rds')

## Hallmark GSEA
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
hallmark_hyp_obj <- hypeR(geneList, HALLMARK, test = "kstest", fdr_cutoff = 0.01, do_plots = TRUE)
saveRDS(hallmark_hyp_obj, 'shiny/data/hsc_w12_13_v_adult_hallmark.rds')
## Laurenti GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$ensembl_gene_id
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
laurenti_hyp_obj <- hypeR(geneList, laurenti, test="kstest", fdr_cutoff=0.01, do_plots = TRUE)
saveRDS(laurenti_hyp_obj, 'shiny/data/hsc_w12_13_v_adult_laurenti.rds')
## GO
gene <- names(geneList)[abs(geneList) > 2]
my_ego2 <- enrichGO(gene          = gene,
                    universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

saveRDS(my_ego2, file = 'shiny/data/hsc_w12_13_v_adult_ego.rds')

# IL7R W12-13 vs Adult

il7r_w12_13_v_adult_coldata <- rownames_to_column(sample_table,"sample_id") %>%
  dplyr::filter(Phenotype == "IL7R+Kit+") %>%
  dplyr::filter(pools == "W12_13" | pools == "BM")

rownames(il7r_w12_13_v_adult_coldata) <- il7r_w12_13_v_adult_coldata$sample_id

counts <- big_counts[, colnames(big_counts) %in% rownames(il7r_w12_13_v_adult_coldata)]                           

dds <- DESeqDataSetFromMatrix(counts, colData = il7r_w12_13_v_adult_coldata, design = ~ pools)
dds$pools <- relevel(dds$pools, ref = "BM")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

# ens_to_external <- readRDS('data/ens_to_external.RDS')

res <- data.frame(res) %>% rownames_to_column("ensembl_gene_id") 
my_colnames <- colnames(res)
res <- left_join(res,ens_to_external,by = "ensembl_gene_id") %>% dplyr::select("ensembl_gene_id","external_gene_name",
                                                                               my_colnames)
res <- res[!is.na(res$padj),]
res <- distinct(res)
res$external_gene_name <- case_when(duplicated(res$external_gene_name) ~ res$ensembl_gene_id,
                                    is.na(res$external_gene_name) ~ res$ensembl_gene_id,
                                    TRUE ~ as.character(res$external_gene_name))
rownames(res) <- res$external_gene_name
res$external_gene_name <- NULL

saveRDS(res, 'shiny/data/il7r_w12_13_v_adult.rds')

## Hallmark GSEA
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
hallmark_hyp_obj <- hypeR(geneList, HALLMARK, test = "kstest", fdr_cutoff = 0.01, do_plots = TRUE)
saveRDS(hallmark_hyp_obj, 'shiny/data/il7r_w12_13_v_adult_hallmark.rds')
## Laurenti GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$ensembl_gene_id
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
laurenti_hyp_obj <- hypeR(geneList, laurenti, test="kstest", fdr_cutoff=0.01, do_plots = TRUE)
saveRDS(laurenti_hyp_obj, 'shiny/data/il7r_w12_13_v_adult_laurenti.rds')
## GO
gene <- names(geneList)[abs(geneList) > 2]
my_ego2 <- enrichGO(gene          = gene,
                    universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

saveRDS(my_ego2, file = 'shiny/data/il7r_w12_13_v_adult_ego.rds')

# proB W12-13 vs Adult

proB_w12_13_v_adult_coldata <- rownames_to_column(sample_table,"sample_id") %>%
  dplyr::filter(Phenotype == "CD34+CD19+") %>%
  dplyr::filter(pools == "W12_13" | pools == "BM")

rownames(proB_w12_13_v_adult_coldata) <- proB_w12_13_v_adult_coldata$sample_id

counts <- big_counts[, colnames(big_counts) %in% rownames(proB_w12_13_v_adult_coldata)]                           

dds <- DESeqDataSetFromMatrix(counts, colData = proB_w12_13_v_adult_coldata, design = ~ pools)
dds$pools <- relevel(dds$pools, ref = "BM")
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

# ens_to_external <- readRDS('data/ens_to_external.RDS')

res <- data.frame(res) %>% rownames_to_column("ensembl_gene_id") 
my_colnames <- colnames(res)
res <- left_join(res,ens_to_external,by = "ensembl_gene_id") %>% dplyr::select("ensembl_gene_id","external_gene_name",
                                                                               my_colnames)
res <- res[!is.na(res$padj),]
res <- distinct(res)
res$external_gene_name <- case_when(duplicated(res$external_gene_name) ~ res$ensembl_gene_id,
                                    is.na(res$external_gene_name) ~ res$ensembl_gene_id,
                                    TRUE ~ as.character(res$external_gene_name))
rownames(res) <- res$external_gene_name
res$external_gene_name <- NULL

saveRDS(res, 'shiny/data/proB_w12_13_v_adult.rds')

## Hallmark GSEA
geneList <- res$log2FoldChange
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
hallmark_hyp_obj <- hypeR(geneList, HALLMARK, test = "kstest", fdr_cutoff = 0.01, do_plots = TRUE)
saveRDS(hallmark_hyp_obj, 'shiny/data/proB_w12_13_v_adult_hallmark.rds')
## Laurenti GSEA
geneList <- res$log2FoldChange
names(geneList) <- res$ensembl_gene_id
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
laurenti_hyp_obj <- hypeR(geneList, laurenti, test="kstest", fdr_cutoff=0.01, do_plots = TRUE)
saveRDS(laurenti_hyp_obj, 'shiny/data/proB_w12_13_v_adult_laurenti.rds')
## GO
gene <- names(geneList)[abs(geneList) > 2]
my_ego2 <- enrichGO(gene          = gene,
                    universe      = names(geneList),
                    OrgDb         = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

saveRDS(my_ego2, file = 'shiny/data/proB_w12_13_v_adult_ego.rds')


# PCA
dds <- DESeqDataSetFromMatrix(big_counts, colData = sample_table, design = ~ pools)
rld <- vst(dds)

# calculate the variance for each gene
rv <- rowVars(assay(rkd))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

intgroup.df <- as.data.frame(colData(rld)[, c("cell_type", "pools"), drop=FALSE])

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],PC4=pca$x[,4],PC5=pca$x[,5],
                intgroup.df, name=colnames(rld))

p1 <- ggplot(d, aes(PC1,PC2, shape=pools, colour = cell_type)) + 
  geom_point(size = 3 ) + theme_classic()
p2 <- ggplot(d, aes(PC2,PC3, shape=pools, colour = cell_type)) + 
  geom_point(size = 3 ) + theme_classic()
p3 <- ggplot(d, aes(PC3,PC4, shape=pools, colour = cell_type)) + 
  geom_point(size = 3 ) + theme_classic()
p4 <- ggplot(d, aes(PC4,PC5, shape=pools, colour = cell_type)) + 
  geom_point(size = 3 ) + theme_classic()

pca_plot <- grid.arrange(p1,p2,p3,p4)
saveRDS(pca_plot, 'shiny/data/pca_plot.rds')
