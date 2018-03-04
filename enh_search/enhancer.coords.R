library(stringr)

mouse.pos <- fread("../data/mouse_positive.bed") %>% setNames(c("chr", "start", "end"))
mouse.pos$start <- as.character(mouse.pos$start)
mouse.pos$end <- as.character(mouse.pos$end)

mouse.all <- fread("../data/mouse.enh.genes.csv", header = T)

x <- str_split_fixed(mouse.all$V2, ":", 2)
y <- as.data.frame(str_split_fixed(x[,2], "-", 2)) %>% mutate_all(funs(gsub(",", "", .)))
mouse.all <- cbind(x[,1], y, mouse.all %>% select(3)) %>% setNames(c("chr", "start", "end", "genes")) %>% 
  mutate(gene1 = sub("-.*$", "", genes), gene2 = sub("^[^-]*-", "", genes)) %>% 
  filter(!grepl("\\(intragenic\\)", genes))

write.table(mouse.all %>% select(-genes), "mouse.pos.enh.txt", sep = "\t", col.names = T, row.names = F, quote = F)

homo.all <- fread("../data/human.enh.genes.csv", header = F)

x <- str_split_fixed(homo.all$V2, ":", 2)
y <- as.data.frame(str_split_fixed(x[,2], "-", 2)) %>% mutate_all(funs(gsub(",", "", .)))
homo.all <- cbind(x[,1], y, homo.all %>% select(3)) %>% setNames(c("chr", "start", "end", "genes")) %>% 
  mutate(gene1 = sub("-.*$", "", genes), gene2 = sub("^[^-]*-", "", genes)) %>%
  filter(!grepl("\\(intragenic\\)", genes))

write.table(homo.all %>% select(-genes), "homo.pos.enh.txt", sep = "\t", col.names = T, row.names = F, quote = F)


mm9.genes <- fread("../data/mm9.genes.gtf")
xx <- str_split_fixed(mm9.genes$V9, ";", 10)
gid <- xx[,5]
gid <- sub('.*"([^"]*)".*', "\\1", gid)
mm9.genes$V9 <- gid
names(mm9.genes)[9] <- "gene"
write.table(mm9.genes %>% 
              select(1,4,5,7,9) %>% 
              setNames(c("chr", "start", "end", "strand", "gene_id")), "mm9.genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

genes1 <- merge(mouse.all, mm9.genes, by.x = "gene1", by.y = "gene") %>% filter(chr == V1)
genes1$start <- as.integer(genes1$start)
genes1$dist <- genes1$start - genes1$V4
plot(genes1$dist, type = "h")
mean(abs(genes1$dist))
median(abs(genes1$dist))
boxplot(genes1$dist, outline = F)


genes2 <- merge(mouse.all, mm9.genes, by.x = "gene2", by.y = "gene") %>% filter(chr == V1)
genes2$start <- as.integer(genes2$start)
genes2$end <- as.integer(genes2$end)
genes2$dist <- ifelse(abs(genes2$start - genes2$V4) > abs(genes2$end - genes2$V5), 
                      abs(genes2$start - genes2$V4),
                      abs(genes2$end - genes2$V5))
plot(genes2$dist, type = "h")
mean(abs(genes2$dist))
median(abs(genes2$dist))
boxplot(genes2$dist, outline = F)

enh.dist <- c(genes1$dist, genes2$dist)[c(genes1$dist, genes2$dist) > 0]
median(enh.dist)
mean(enh.dist)

pdf("mm9.enh.dist.pdf")
options(scipen=99)
boxplot(enh.dist, outline = F, main = "mm9 enhancer distance")
dev.off()






mm10.genes <- fread("~/IMG/data/mus/mm10.gene.gtf")
xx <- str_split_fixed(mm10.genes$V9, ";", 10)
gid <- xx[,3]
gid <- sub('.*"([^"]*)".*', "\\1", gid)
mm10.genes$V9 <- gid
names(mm10.genes)[9] <- "gene"
write.table(mm10.genes %>% 
              select(1,4,5,7,9) %>% 
              setNames(c("chr", "start", "end", "strand", "gene_id")), "mm10.genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)





hg19.genes <- fread("../data/hg19.genes.gtf")
yy <- str_split_fixed(hg19.genes$V9, ";", 10)
gid <- yy[,5]
gid <- sub('.*"([^"]*)".*', "\\1", gid)

hg19.genes$V9 <- gid
names(hg19.genes)[9] <- "gene"
write.table(hg19.genes %>% 
              select(1,4,5,7,9) %>% 
              setNames(c("chr", "start", "end", "strand", "gene_id")), "hg19.genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)



genes1 <- merge(homo.all, hg19.genes, by.x = "gene1", by.y = "gene") %>% dplyr::slice(-c(426, 429))
genes1$start <- as.integer(genes1$start)
genes1$dist <- genes1$start - genes1$V4
plot(genes1$dist, type = "h")
mean(abs(genes1$dist))
boxplot(genes1$dist, outline = F)


genes2 <- merge(homo.all, hg19.genes, by.x = "gene2", by.y = "gene") %>% filter(chr == V1)
genes2$end <- as.integer(genes2$end)
genes2$dist <- genes2$V4 - genes2$end
plot(genes2$dist, type = "h")
mean(abs(genes1$dist))
boxplot(genes1$dist, outline = F)

enh.dist <- c(genes1$dist, genes2$dist)[c(genes1$dist, genes2$dist) > 0]
median(enh.dist)
mean(enh.dist)

pdf("hg19.enh.dist.pdf")
options(scipen=99)
boxplot(enh.dist, outline = F, main = "hg19 enhancer distance")
dev.off()

