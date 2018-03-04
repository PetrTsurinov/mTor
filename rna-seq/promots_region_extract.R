all_promotors= read.csv('../../../BioHack/promotors.csv')
head(all_promotors)
mTOR_path

GeneCl_two <- mTOR_path
GeneCl_twoID <- bitr(GeneCl_two, fromType = 'SYMBOL',
                     toType = c("ENTREZID"),
                     OrgDb = org.Hs.eg.db)

mTORpromots = data.frame(matrix(ncol = ncol(all_promotors)))
colnames(mTORpromots) = colnames(all_promotors)
for(i in GeneCl_twoID$ENTREZID){
  mTORpromots = rbind(mTORpromots,all_promotors[all_promotors$geneID == i,])
}
mTORpromots=mTORpromots[-1,-1]

GeneCl_two <- mTORpromots$geneID
GeneCl_twoID <- bitr(GeneCl_two, fromType = 'ENTREZID',
                     toType = c("SYMBOL"),
                     OrgDb = org.Hs.eg.db)
mTORpromots$Gene_SYMBOL = GeneCl_twoID$SYMBOL

write.table(mTORpromots, file = './mTORpromots.bed', sep = '\t', row.names = F, quote = F)


by_TF <- lapply(unique(res_rnaseq_inter$TF), function(x){res_rnaseq_inter$gene_mTOR[res_rnaseq_inter$TF == x]})


####to_.bed
setwd('../../')
num_TF = 1
for(i in by_TF){
  mTORpromots %>%  
    dplyr::filter(grepl(paste0(i, collapse = '|'), mTORpromots$Gene_SYMBOL)) %>% 
    write.table(file = paste0('./',unique(res_rnaseq_inter$TF)[num_TF],'promotors','.bed'), 
                sep = '\t', 
                row.names = F, 
                quote = F)
  num_TF = num_TF + 1
}


