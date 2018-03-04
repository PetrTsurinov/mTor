library(dplyr)
library(edgeR)
library(ggrepel)
library("MASS")
library(biomaRt)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

dt <- read.csv('~/Dropbox/BioHack/new_table2.csv')
colnames(dt)[96] <- 'GeneID'
dt2 <- dt
data <- dt2
rownames(data) = data$GeneID
data <- data[-96]
data <- as.matrix(data)
dez <- factor(ifelse(grepl("OSKM", colnames(data)),"IM","C"))
col <- ifelse(dez=="IM","red","green")

  edger <- DGEList(data, group=dez)  
  
  ## Filtering
  
  edger <- edger[rowSums(edger$counts) > 1,] #оставляет только строки сумма которых больше 1 (удаляем 0 и 1 значения)
  o <- order(rowSums(edger$counts), decreasing=TRUE) # перестаовка в порядке убывания (Зачем? Затем что написано ниже)
  edger <- edger[o,]
  d <- duplicated(rownames(edger)) #определеяет дупликаты, при этом дупликаты те элементы у которых индекс выше, т.е. они ниже по списку 
  edger <- edger[!d,]
  
  
  ## Нормализация
  
  edger <- calcNormFactors(edger, method = "RLE")
  
  # nrow(edger)
  # plotMDS(edger, col = col)

  
dez <- relevel(dez, ref="C") #перетасовывает переменные  

  
  design <- model.matrix(~dez)
  rownames(design) <- colnames(edger)
  edger <- estimateDisp(edger, design)
  print(design)
  
  ####Fitting
    fit <-glmFit(edger, design)
    lrt <- glmLRT(fit)
  
  res <- topTags(lrt, n = nrow(lrt$table))
  res <- as.data.frame(res)

diff_data <- as.data.frame(cpm(edger))
colnames(diff_data) <- paste0(dez,1:96)
GeneCl_two <- rownames(diff_data)
GeneCl_twoID <- bitr(GeneCl_two, fromType = 'ENSEMBL',
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db)

colnames(GeneCl_twoID)[1] <- 'GeneID'
diff_data$GeneID <- rownames(diff_data)
my_diff <- merge.data.frame(x = diff_data, y = GeneCl_twoID, by = 'GeneID')


diff_data$GeneID <- rownames(diff_data)
res$GeneID <- rownames(res)
ebdata <- merge(diff_data, res[c(4,5,6)], by = 'GeneID')
ebdata <- ebdata[ebdata$FDR < 0.05,]




ss

###Какая-то херня######


mTOR_path <- as.character(read.table('../../geneset.txt',header = F)$V1)
my_diff_sub = as.data.frame(matrix(nrow = 1, ncol = ncol(my_diff)))
colnames(my_diff_sub) = colnames(my_diff)
for(i in 1:nrow(my_diff)){
  if(my_diff$SYMBOL[i] %in% mTOR_path){
    my_diff_sub <- rbind(my_diff_sub, my_diff[i,])
  }
}

my_diff_sub %>% filter(PValue < 0.05) %>% select_('SYMBOL')

res.data = data.frame()
edger_data <- edger$counts
for(i in my_diff_sub$GeneID){
  res.data = rbind(res.data, edger_data[rownames(edger_data) == i,])
}

sum(apply(res.data, 1, median) > 100)


############

library(pheatmap)
edata = as.matrix(ebdata[c(2:97)])
edata = edata/mean(edata)

data_log <- log2(edata)

# Center and scale data by row
scale_rows <- function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

data_norm <- scale_rows(edata)

pheatmap(data_norm, cluster_rows=T, cluster_cols=T, show_rownames = F, legend = F)




