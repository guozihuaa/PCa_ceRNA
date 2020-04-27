rm(list = ls())
library(readr)
library(tidyverse)
library(readxl)
setwd("F:/WORK/Lab Work/indivial/GUOZIHU/")

mirwalk5utr <- read_tsv(file = "hsa_miRWalk_5UTR/hsa_miRWalk_5UTR.txt" )

mirwalk5utr2 <- dplyr::select(mirwalk5utr, miRNA, Genesymbol) %>% 
                distinct()

rm(mirwalk5utr)
gc()

mirwalkcds <- read_tsv(file = "hsa_miRWalk_CDS/hsa_miRWalk_CDS.txt")

mirwalkcds2 <- dplyr::select(mirwalkcds, miRNA, Genesymbol) %>% 
  distinct()

rm(mirwalkcds)
gc()


mir3utr <- read_tsv(file = "hsa_miRWalk_3UTR/hsa_miRWalk_3UTR.txt")

mirwalk3utr2 <- dplyr::select(mir3utr, miRNA, Genesymbol) %>% 
  distinct()

rm(mir3utr)
gc()


mirna_mran <- bind_rows(mirwalk3utr2, mirwalk5utr2, mirwalkcds2) %>% 
              distinct()

save(mirna_mran, file = "2/rawdata/miRNA_mRNA_inte.RData")


############## lncbase

lncbaseraw <- read_csv(file = "lncBaseV2_predicted_data (2)/lncBaseV2_predicted_human_data.csv")
names(lncbaseraw) <- c("trans", "Gene_ID", "miRNA", "score")

lncbaseraw2 <- lncbaseraw[grep("hsa", lncbaseraw$miRNA), ]

lncbaseraw2 <- distinct(lncbaseraw2)

mirnas <- str_extract(lncbaseraw2$Gene_ID, "\\([^()]+\\)")

k <- substring(mirnas, 2, nchar(mirnas)-1)

genes <- str_remove(lncbaseraw2$miRNA, "\\([^()]+\\)")


lncbase_inter <- data.frame(mirna = genes, target = k, stringsAsFactors = F) %>% 
                 distinct()


################# starbase 数据

starbase <- read_xlsx(path = "starBase_Human_Pan_Cancer_miRNA_LncRNA.xlsx") %>% 
            dplyr::select(mirna = name, target = geneName)

############### miRcode

mircodeH <- read_tsv("mircode_highconsfamilies.txt/mircode_highconsfamilies.txt") %>% 
            dplyr::select(microrna, gene_symbol)

mircodeM <- read_tsv("mircode_medconsfamilies.txt/mircode_medconsfamilies.txt") %>% 
  dplyr::select(microrna, gene_symbol)

mircode <- bind_rows(mircodeH, mircodeM) %>% 
           distinct()

allmir <- str_remove(mircode$microrna, "miR-") %>% 
          str_split("/")

alltarget <- rep(mircode$gene_symbol, sapply(allmir, length))

mircodeDF <- data.frame(mirna = unlist(allmir), target = alltarget, stringsAsFactors = F) %>% 
              distinct()



mirname_con <- function(mirname){
  
  tmps <- str_extract_all(str_extract(mirname, "[a-f]+"), "[a-f]")[[1]]
  
  return(str_replace(mirname, "[a-f]+", tmps))
  
}

containab_re <- filter(mircodeDF, grepl("ab",mirna))

splitrname <- lapply(unique(containab_re$mirna), function(x){mirname_con(x)})
names(splitrname) <-  unique(containab_re$mirna)

newmirname <- lapply(containab_re$mirna, function(x){splitrname[[x]]})

newmircodeDF2 <- data.frame(mirna  = unlist(newmirname),
                            target = rep(containab_re$target, sapply(newmirname, length)),
                            stringsAsFactors = F)


newmircodeDF3 <- filter(mircodeDF, !grepl("ab",mirna))


newmircodeDF <- bind_rows(newmircodeDF2, newmircodeDF3) %>% 
               distinct() %>% 
               mutate(mirna = str_c("hsa-miR-", mirna, sep = ""))


##################################### lncSNP
library(annotables)
data(grch37)

lincpida <- read_tsv("lncRNASNP2/human_lncrna_alias.txt")

aliass <- str_split(lincpida$alias, ",|;")

lincpidanew <- data.frame(id1 = rep(lincpida$transcript, sapply(aliass, length)),
                          lncIDs = unlist(aliass),
                          stringsAsFactors = F) %>%
              dplyr::filter(grepl("ENSG", lncIDs)) %>% 
              mutate(lncIDs = tools::file_path_sans_ext(lncIDs)) %>% 
              distinct()

lncsnpH <- read_tsv("lncRNASNP2/mirnas_lncrnas_conserved.txt", col_names = F)

lncsnp <- read_tsv("lncRNASNP2/mirnas_lncrnas_validated.txt", col_names = F) %>% 
           dplyr::select(X1, X2) %>% 
          bind_rows(lncsnpH) %>% 
         distinct() %>% 
         inner_join(lincpidanew, by = c("X1" = "id1")) %>% 
         inner_join(grch37, by = c("lncIDs" = "ensgene")) %>% 
         dplyr::select(mirna = X2, target = symbol) %>% 
        distinct()

#################

all_mir_lncRNA <- bind_rows(lncsnp, newmircodeDF, starbase, lncbase_inter) %>% 
                  distinct()

save(all_mir_lncRNA, file = "2/rawdata/miRNA_lncRNA_inter.RData")

load("rawdata/miRNA_lncRNA_inter.RData")
