setwd('./Autism/GSE65106_RAW/')

source("https://bioconductor.org/biocLite.R")
biocLite("hgu133plus2.db")
library("hgu133plus2.db")
x <- hgu133plus2ALIAS2PROBE
as.list(x)[[100000]]
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


#celFiles <- list.celfiles()
#celFiles <- celFiles[1:33]
#affyRaw <- read.celfiles(celFiles)
#eset <- rma(affyRaw)
write.exprs(eset,file="data.txt")

my_frame <- read.table(file = './big_data/embryonic_data 2.txt', header = T, sep = '\t')
my_frame_fibro <- read.table(file = './big_data/fibroblast_data.txt', header = T, sep = '\t')

my_frame_merded <- cbind.data.frame(my_frame, my_frame_fibro)

colnames(my_frame_merded)[2:104] <- paste0(colnames(my_frame_merded)[2:104], '_ESC')
colnames(my_frame_merded)[106:548] <- paste0(colnames(my_frame_merded)[106:548], '_fibrob')

my_frame_merded <- my_frame_merded[-105]

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

design <- model.matrix(~colnames(my_frame)) #col_names должно быть маркером группы

fit <- lmFit(eset, design)
fit <- eBayes(fit)

my_coef = 2

dataf = as.data.frame(topTable(fit, coef=my_coef, number = nrow(fit[[1]])))
colnames(design)[my_coef]

####Контрасты#####
fit2 <- contrasts.fit(fit, c(0,-1,1))
fit <- eBayes(fit2)
dataf = as.data.frame(topTable(fit, number = nrow(fit[[1]])))

###Аннотация####
Ann <- read.csv('./big_data/annotation.csv.csv', header = F, sep = ',', stringsAsFactors = F)
colnames(Ann) <- c('probe_id', 'GeneID')
Ann$probe_id <- as.character(Ann$probe_id)
#########

get_results <- function(res, Ann){
  res$probe_id <- rownames(res)
  res <- res %>% filter(abs(logFC) > 0.5)
  res <- res %>% filter(P.Value < 0.05)
  ####Annotation
  dt <- merge(res, Ann, by = 'probe_id')
  print(head(dt))
  mart <- useDataset(dataset = "hsapiens_gene_ensembl", 
                     mart    = useMart("ENSEMBL_MART_ENSEMBL",
                                       host    = "www.ensembl.org"))
  resultTable <- getBM(attributes = c("ensembl_transcript_id","external_gene_name"),
                       mart       = mart)
  resultTable <- resultTable[!resultTable$ensembl_transcript_id == "",]
  
  dt$Gene_Symbol <- NA
  dt <- dt[!is.na(dt$ENST),]

  for(i in 1:nrow(dt)){
    geneNameList = character()
    dtENST <- unlist(strsplit(dt[i,]$ENST, split = ','))
    for(j in 1:length(dtENST)){
      geneNameList <- c(geneNameList,resultTable[resultTable$ensembl_transcript_id == dtENST[j],]$external_gene_name)
    }
    if(length(unique(geneNameList)) > 0){
      dt$Gene_Symbol[i] <- paste(unique(geneNameList), collapse = ',')}
  }
  
  return(dt)
  } #пришиватель аннотации
get_mTOR <- function(dt){
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
}

dt_iPSCvsESC_my_function = get_results(res = dataf, Ann = Ann)


dt <- dt_iPSCvsESC_my_function
dt <- dt_iPSC_my_function
dt <- dt_hESC_my_function

my_dt_sub = data.frame(matrix(ncol = ncol(dt)))
colnames(my_dt_sub) = colnames(dt)
for(i in 1:nrow(dt)){
  for(gene in unlist(strsplit(dt$Gene_Symbol[i], ','))){
    print(gene)
    if(gene %in% mTOR_path){
      print('1')
      my_dt_sub <- rbind(my_dt_sub, dt[i,])
    }
  }

}

my_dt_sub <- my_dt_sub[-1,]
mTORsub_iPSC <- my_dt_sub
mTORsub_hESC <- my_dt_sub


head(dt_hESC_my_function$Gene_Symbol)
head(dt_hESC$Gene_Symbol)  

#####to ML

colnames(my_frame_merded)[1] <- 'probe_id'
to_ML <- merge.data.frame(my_frame_merded, Ann, by = 'probe_id', all.x = T)


#get_results_ML <- function(res, Ann){
  #res$probe_id <- rownames(res)
  ####Annotation
  dt <- merge(res, Ann, by = 'probe_id')
  print(head(dt))
  mart <- useDataset(dataset = "hsapiens_gene_ensembl", 
                     mart    = useMart("ENSEMBL_MART_ENSEMBL",
                                       host    = "www.ensembl.org"))
  resultTable <- getBM(attributes = c("ensembl_transcript_id","external_gene_name"),
                       mart       = mart)
  resultTable <- resultTable[!resultTable$ensembl_transcript_id == "",]
  
  dt$Gene_Symbol <- NA
  dt <- dt[!is.na(dt$ENST),]
  
  for(i in 1:nrow(dt)){
    geneNameList = character()
    dtENST <- unlist(strsplit(dt[i,]$ENST, split = ','))
    for(j in 1:length(dtENST)){
      geneNameList <- c(geneNameList,resultTable[resultTable$ensembl_transcript_id == dtENST[j],]$external_gene_name)
    }
    if(length(unique(geneNameList)) > 0){
      dt$Gene_Symbol[i] <- paste(unique(geneNameList), collapse = ',')}
  }
  
  return(dt)
} 

my_frame_merded$probe_id

#colnames(my_frame) <- paste0(colnames(my_frame),1:33)
to_ML <- get_results_ML(my_frame_merded, Ann)
write_csv(x = to_ML, path = './to_ML.csv')


###correlation

library(Hmisc)

my_TF = c('POU5F1', 'SOX2', 'NANOG', 'KLF4', 'MYC')
mTOR_path

res_matrix <- matrix(ncol = length(my_TF), nrow = length(mTOR_path))
colnames(res_matrix) = my_TF
rownames(res_matrix) = mTOR_path


random_corr_distribution = function(tests, data){
  res = vector()
  for(i in 1:tests){
    random_set <- sample(1:nrow(data), size = 2, replace = F)
    
    our_gene1 <- to_ML[random_set[1],][-c(1,numcol, numcol-1)]
    our_gene2 <- to_ML[random_set[2],][-c(1,numcol, numcol-1)]
    
    our_gene1 <- apply(our_gene1, 2, sum)[apply(design,1,sum) == 1]
    our_gene2 <- apply(our_gene2, 2, sum)[apply(design,1,sum) == 1]
    
    res <- c(res,rcorr(our_gene1, our_gene2, type = 'spearman')$r[1,2])
    
  }
  return(res)
}
get_p.value = function(gene1,gene2,to_ML){
  numcol = ncol(to_ML)
  our_gene1 <- to_ML[grep(paste0('\\',gene1,'\\b'),to_ML$Gene_Symbol),][-c(1,numcol, numcol-1)]
  our_gene2 <- to_ML[grep(paste0('\\',gene2,'\\b'),to_ML$Gene_Symbol),][-c(1,numcol, numcol-1)]
  
  our_gene1 <- apply(our_gene1, 2, mean)[apply(design,1,sum) == 1]
  our_gene2 <- apply(our_gene2, 2, mean)[apply(design,1,sum) == 1]
  
  test_ourGenes <- rcorr(our_gene1, our_gene2, type = 'spearman')$r[1,2]
  my_random <- random_corr_distribution(2500,to_ML)
  
  z=(test_ourGenes - mean(my_random))/sd(my_random)
  results = pnorm(-abs(z)) #p-value
  return(results)
}
aaa <- random_corr_distribution(2500,to_ML)
sink('../../set_corr.csv')
for(TF in my_TF){
  for(gene in mTOR_path){
    my_p.value <- get_p.value(TF, gene, to_ML, aaa)
    res_matrix[grep(gene, rownames(res_matrix)), grep(TF, colnames(res_matrix))] <- my_p.value
      cat(TF,gene,my_p.value, sep = ',', end = '\n')
  }
}
sink()


data <- read.csv('../../set_corr.csv', header = F, stringsAsFactors = F)
data$V4 <- stats::p.adjust(data$V3 , method = 'BH')
colnames(data) = c('TF', 'gene_mTOR', 'p_value', 'p_adj_BH')
res_rnaseq_inter <- data

write.csv(file = './set_corr.adj.csv', x = data, row.names = F)

TF = 'POU5F1'
gene = 'FBXW11'


get_p.value(gene1, gene2, to_ML)


