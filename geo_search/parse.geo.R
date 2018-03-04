library(dplyr)
library(GenomicRanges)
library(data.table)
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# biocLite("GEOmetadb")
library(edgeR)
library(GEOmetadb)

mtorc1 <- c(
  "Mtor",
  "Rptor",
  "Akt1s1",
  "Mlst8",
  "Deptor",
  "Telo2",
  "Tti1"
)

mtorc2 <- c(
  "Mtor",
  "Mapkap1",
  "Rictor",
  "Prr5",
  "Deptor",
  "Mlst8",
  "Telo2",
  "Tti1"
)

getSQLiteFile()

if(file.exists('GEOmetadb.sqlite')) {
  a <- columnDescriptions()
  b <- geoConvert("Mus",'GSM')
} else {
  print("use getSQLiteFile() to get a copy of the GEOmetadb SQLite file
and then rerun the example")
}

con <- dbConnect(SQLite(),'GEOmetadb.sqlite')

sql <- paste("SELECT DISTINCT gse.title,gse.gse",
             "FROM",
             "  gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
             "  JOIN gse ON gse_gsm.gse=gse.gse",
             "  JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             "  gse.title LIKE '%Chip-seq%' AND",
             "  gse.title LIKE '%Oct4%' AND",
             "  gpl.organism LIKE '%Homo sapiens%' OR",
             "  gpl.organism LIKE '%Mus musculus%' OR",
             "  gpl.organism LIKE '%Drosophila%' OR",
             "  gpl.organism LIKE '%Rattus norvegicus%'",sep=" ")

sql.hum <- paste("SELECT DISTINCT gse.title,gse.gse",
                 "FROM",
                 "  gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
                 "  JOIN gse ON gse_gsm.gse=gse.gse",
                 "  JOIN gse_gpl ON gse_gpl.gse=gse.gse",
                 "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
                 "WHERE",
                 "  gse.title LIKE '%Chip-seq%' AND",
                 "  gse.title LIKE '%Sox2%' AND",
                 "  gpl.organism LIKE '%Homo sapiens%'",sep=" ")

query.chip.human <- dbGetQuery(con,sql.hum)

sql.chip <- paste("SELECT DISTINCT gse.title, gse.gse, gsm.title,gsm.gsm",
                 "FROM",
                 "  gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
                 "  JOIN gse ON gse_gsm.gse=gse.gse",
                 "  JOIN gse_gpl ON gse_gpl.gse=gse.gse",
                 "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
                 "WHERE",
                 "  gse.title LIKE '%Chip-seq%' AND",
                 "  (gsm.title LIKE '%Sox2%' OR gsm.title LIKE '%Oct4%' OR gsm.title LIKE '%Nanog%'",
                 "  OR gsm.title LIKE '%Klf4%' OR gsm.title LIKE '%C-myc%') AND",
                 "  (gpl.organism LIKE '%Mus musculus%' OR gpl.organism LIKE '%Homo sapiens%' OR",
                 "  gpl.organism LIKE '%Rattus norvegicus%' OR gpl.organism LIKE '%melanogaster%')", sep=" ")

query.chip <- dbGetQuery(con,sql.chip)


sql.exp <- paste("SELECT DISTINCT gse.title, gse.gse, gsm.title,gsm.gsm",
             "FROM",
             "  gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
             "  JOIN gse ON gse_gsm.gse=gse.gse",
             "  JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             "  gsm.title LIKE '%RNA-seq%' AND",
             "  (gsm.title LIKE '%ESC%' OR gsm.title LIKE '%MEF%' OR gsm.title LIKE '%IPSc%'",
             "  OR gsm.title LIKE '%Klf4%' OR gsm.title LIKE '%C-myc%') AND",
             "  (gpl.organism LIKE '%Mus musculus%' OR gpl.organism LIKE '%Homo sapiens%' OR",
             "  gpl.organism LIKE '%Rattus norvegicus%')", sep=" ")

query.exp <- dbGetQuery(con,sql.exp)


q.c.title <- query.chip %>% group_by(gse, s.title) %>% summarize( n = n())

q.e.title <- query.exp %>%setNames(c("gse.title", "gse", "gsm.title", "gsm")) %>%  group_by(gse, gse.title) %>% summarize( n = n())

write.table(q.c.title[1:20,1:2], "chip.seq.table.1.csv", sep = "\t", row.names = F, col.names = T, quote = F)

write.table(q.c.title[21:45,1:2], "chip.seq.table.2.csv", sep = "\t", row.names = F, col.names = T, quote = F)

write.table(q.e.title, "rna.seq.table.csv", sep = "\t", row.names = F, col.names = T, quote = F)


