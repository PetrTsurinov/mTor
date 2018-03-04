setwd('./Autism/GSE65106_RAW/')

library(dplyr)
library(affy)
library(tkWidgets)
library(tcltk)
library(oligo)
library(limma)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(affyQCReport)
library(biomaRt)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)


celFiles <- list.celfiles()
celFiles <- celFiles[1:33]
affyRaw <- read.celfiles(celFiles)
eset <- rma(affyRaw)
write.exprs(eset,file="data.txt")

my_frame <- data.frame(exprs(eset))
#my_frame <- write.table('./data.txt', sep = )
#Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
#all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
#write.table(all,file="data.ann.txt",sep="\t")

###QC
plotMA3by2(as.matrix(my_frame))
boxplot(log(my_frame))
hist(my_frame)

###analisys by limma

old_name = colnames(my_frame)
old_name[grepl('hESC', old_name)] <- 'hESC'
old_name[grepl('Fibroblast', old_name)] <- 'Fibroblast'
old_name[grepl('iPSC_', old_name)] <- 'iPSC'

colnames(my_frame) <- old_name

old_name

design <- cbind(CTRL=1, 
                iPSCvsCTRL=c(rep(0,20),rep(1,13)), 
                hESCvsCTRL=c(rep(0,16),rep(1,4),rep(0,13)))
fit <- lmFit(eset, design)
fit <- eBayes(fit)
res_hESC = as.data.frame(topTable(fit, coef="hESCvsCTRL", number = 20000))
res_iPSC = as.data.frame(topTable(fit, coef="iPSCvsCTRL", number = 20000))
res_iPSC$probe_id <- rownames(res_iPSC)
res_hESC$probe_id <- rownames(res_hESC)
png(filename = '~/Dropbox/CardioCenter/projects/2018/Philip/Results/Hist.png', width = 800, height = 600, res = 150)
hist(res$logFC, 
     main = 'Histogram of logFC', 
     xlab = 'LogGC', 
     ylab = 'Frequency of genes')
dev.off()
res <- res %>% filter(abs(logFC) > 0.5)
res_hESC <- res_hESC %>% filter(abs(logFC) > 0.5)
res_iPSC <- res_iPSC %>% filter(abs(logFC) > 0.5)
png(filename = '~/Dropbox/CardioCenter/projects/2018/Philip/Results/Hist2.png', width = 800, height = 600, res = 150)
hist(res_hESC$logFC, 
     main = 'Histogram of logFC', 
     xlab = 'LogGC', 
     ylab = 'Frequency of genes')
dev.off()
min(abs(res_hESC$logFC))
####Annotation
Ann <- read.csv('../HuGene-1_0-st-v1-na36-hg19-transcript-csv/Annotaion.csv', header = F, sep = ';', stringsAsFactors = F)
head(Ann)
colnames(Ann) <- c('probe_id', 'ENST')
Ann$probe_id <- as.character(Ann$probe_id)
dt_hESC <- merge(res_hESC, Ann, by = 'probe_id')
dt_iPSC <- merge(res_iPSC, Ann, by = 'probe_id')

mart <- useDataset(dataset = "hsapiens_gene_ensembl", 
                   mart    = useMart("ENSEMBL_MART_ENSEMBL",
                                     host    = "www.ensembl.org"))
resultTable <- getBM(attributes = c("ensembl_transcript_id","external_gene_name"),
                     mart       = mart)
resultTable <- resultTable[!resultTable$ensembl_transcript_id == "",]

resultTable[resultTable$ensembl_transcript_id == "ENST00000375078",]


dt_hESC_pval <- dt_hESC %>% filter(P.Value < 0.05)
dt_iPSC_pval <- dt_iPSC %>% filter(P.Value < 0.05)


dt <- dt[order(dt$P.Value, decreasing = F),]
dt$Gene_Symbol <- NA
dt <- dt[!is.na(dt$ENST),]
dt <- dt %>% filter(P.Value < 0.05)
png(filename = '~/Dropbox/CardioCenter/projects/2018/Philip/Results/Hist3.png', width = 800, height = 600, res = 150)
hist(dt$logFC, 
     main = 'Histogram of logFC', 
     xlab = 'LogGC', 
     ylab = 'Frequency of genes')
dev.off()

dt <- dt_hESC_pval
dt$Gene_Symbol <- NA
for(i in 1:nrow(dt)){
  geneNameList = character()
  dtENST <- unlist(strsplit(dt[i,]$ENST, split = ','))
  for(j in 1:length(dtENST)){
    geneNameList <- c(geneNameList,resultTable[resultTable$ensembl_transcript_id == dtENST[j],]$external_gene_name)
  }
  if(length(unique(geneNameList)) > 0){
    dt$Gene_Symbol[i] <- paste(unique(geneNameList), collapse = ',')}
}

my_dt_sub = data.frame()
for(i in 1:nrow(dt)){
  for(gene in unlist(strsplit(dt$Gene_Symbol[i], ','))){
    print(gene)
    if(gene %in% mTOR_path){
      print(gene)
      my_dt_sub <- rbind(my_diff_sub, my_diff[i,])
    }
  }
  
  }


my_dt_sub[my_dt_sub$PValue < 0.05,]
















write.csv(file = '../Results/gene_144.csv', x = dt[!is.na(dt$Gene_Symbol),c(9,2:7)], row.names = F)

CADgenes <- read_file("../../../../Papers/2017/microRNA_network/Final_version/CADgensEXP.txt")
CADgenes <- unlist(strsplit(CADgenes, " ", fixed = TRUE))

# the CADgene taking from Filip data
filipCAD = 0
filipCADgene = character()
for(gene in dt$Gene_Symbol){
  if(gene %in% CADgenes){
    filipCAD = filipCAD + 1
    filipCADgene = c(filipCADgene, gene)
  }
}

cat(filipCADgene)

miR21gene <- read_file("../../../../Papers/2017/microRNA_network/miR_Targets/hsa-miR-21-.txt")
miR21gene <- unlist(strsplit(miR21gene, "\r\n", fixed = TRUE))

# the check mir-21 targets in filipCADgene

for(gene in filipCADgene){
  if(gene %in% miR21gene){
    print(gene)
  }
}

# the check mir-21 targets in filipGene

for(gene in dt$Gene_Symbol){
  if(gene %in% miR21gene){
    print(gene)
  }
}

# the same analysis for mir-223

miR223gene <- read_file("./hsa-miR-223.txt")
miR223gene <- unlist(strsplit(miR223gene, "\r\n", fixed = TRUE))

for(gene in filipCADgene){
  if(gene %in% miR223gene){
    print(gene)
  }
}

for(gene in dt$Gene_Symbol){
  if(gene %in% miR223gene){
    print(gene)
  }
}

####heatmap

library(pheatmap)

new_frame = data.frame()
rbind(new_frame, my_frame[4,])
for(i in 1:nrow(my_frame)){
  if(rownames(my_frame)[i] %in% dt$probe_id){
    new_frame = rbind(new_frame, my_frame[i,])
  }
}

pheatmap(as.matrix(new_frame))


data_log <- log2(new_frame)

# Center and scale data by row
scale_rows <- function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_norm <- scale_rows(data_log)

png('../Results/heatmap.png', width = 800, height = 600, res = 150)
pheatmap(data_norm, cluster_rows=T, cluster_cols=T, show_rownames = F, legend = F)
dev.off()
plotMDS(data_norm, col = c(rep('red',6), rep('green',6)))

#geting cluster

aaa <- pheatmap(data_norm, cluster_rows=T, cluster_cols=T)
cluster <- cutree(aaa$tree_row, k = 2)

cluster <- data.frame(probe_id = names(cluster),
                      cluster = cluster)


dt_clus <- merge.data.frame(dt, cluster, by = 'probe_id', all.x = T)
dt_clus <- dt_clus[!is.na(dt_clus$Gene_Symbol),]


cat(dt_clus$Gene_Symbol[dt_clus$cluster == 2], sep = ', ')



######The Enrichment games

GeneCl_one <- dt_clus$Gene_Symbol[dt_clus$cluster == 1]
GeneCl_two <- dt_clus$Gene_Symbol[dt_clus$cluster == 2]


GeneCl_oneID <- (bitr(GeneCl_one, fromType = 'SYMBOL',
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db))$ENTREZID

GeneCl_twoID <- (bitr(GeneCl_two, fromType = 'SYMBOL',
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db))$ENTREZID


ego1 <- enrichGO(gene          = GeneCl_oneID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
write.csv(x = head(ego1, 40), file = '../Results/enGO_clust1.csv')

dotplot(ego1)
enrichMap(ego1, vertex.label.font = 0.01)
cnetplot(ego1)
plotGOgraph(ego1)

aa <- enrichMap(ego1)
aa_edge_list <- as_edgelist(aa)
write.csv(x = data.frame(E1=aa_edge_list[,1],
                         E2=aa_edge_list[,2]), file = '../Results/graph_ego1.csv', row.names = F)

ego2 <- enrichGO(gene          = GeneCl_twoID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

dotplot(ego2)
enrichMap(ego2)

head(ego2)


kk2 <- enrichKEGG(gene         = GeneCl_twoID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
head(kk2, 20)
browseKEGG(kk2, 'hsa04650')

kk1 <- enrichKEGG(gene         = GeneCl_oneID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
head(kk2, 20)
browseKEGG(kk2, 'hsa04650')

write.csv(x = head(ego1, 40), file = '../Results/enGO_clust1.csv')


kk2_data <- as.data.frame(kk2)
kk2_data$geneSymbol <- NA
for(i in 1:nrow(kk2_data)){
  GeneIDs <- unlist(strsplit(kk2_data$geneID[i], '/'))
  GeneSymbols <- (bitr(GeneIDs, fromType = "ENTREZID",
                       toType = c('SYMBOL'),
                       OrgDb = org.Hs.eg.db))$SYMBOL
  GeneSymbols <- paste(GeneSymbols, collapse = '/')
  kk2_data$geneSymbol[i] <- GeneSymbols
}
kk2_data$geneID <- kk2_data$geneSymbol
kk2_data$geneSymbol <- NULL

write.csv(x = kk2_data, file = '../Results/kk2_clust1.csv')

##common expression of CTRL
png('../Results/dotplotkk.png', width = 800, height = 600, res = 100)
dotplot(kk2)
dev.off()
enrichMap(kk2)
expCTRL <- apply(my_frame[7:ncol(my_frame)],1,median)
sum(expCTRL < 3)

####for network analysis
diffGenList = character()
for(gene in dt_clus$Gene_Symbol){
  diffGenList = c(diffGenList, unlist(strsplit(gene, ',')))
}
write.table(x = diffGenList, file = './diffGenLis_144.txt', eol = ' ', row.names = F, col.names = F, quote = F)



