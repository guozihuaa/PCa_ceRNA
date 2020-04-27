####### 找到合适的miRNA 数据集
dbDEMC <- read_csv("rawdata/miRExpAll.csv", col_names = F) %>% 
  group_split(X3)


names(dbDEMC) <- sapply(dbDEMC, function(x){x$X4[1]})

mirnas <- lapply(dbDEMC, function(x){x$X2})[c("GSE76260",  "GSE21036","GSE31568",	"TCGA_PRAD"	)]



result <- sapply(mirnas, function(x) sapply(mirnas, function(y) length(intersect(x,y))))


############ GEO mRNA 表达数据处理



########## TCGA 前列腺癌
##### 临床数据

#clincila data
clinical <- read_csv("rawdata/clinical_information.csv")

filterclinical <- clinical

############


#### 从TCGA直接下载下载表达数据
query <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts"
)
GDCdownload(query)
data <- GDCprepare(query, 
                   save = TRUE, 
                   save.filename = "HTSeq_coutn_exp.rds",
                   summarizedExperiment = TRUE)

load("HTSeq_coutn_exp.rds")

lihc_ranexp2 <- assay(data)


clinical2 <- colData(data) %>% 
  as.data.frame() %>% 
  dplyr::select(sample, patient, barcode, shortLetterCode)

save(lihc_ranexp2, clinical2, file = "PRAD_RNA_and_clinical.RData")

load("PRAD_RNA_and_clinical.RData")

########## 下载miRNA数据
library(GDCRNATools)
gdcRNADownload(project.id     = 'TCGA-PRAD', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               directory      = 'data/',
               method = "gdc-client")


metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-PRAD',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

mirMatrix <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = 'data',
                         data.type = 'miRNAs',
                         organized = F)

mirExpr <- gdcVoomNormalization(counts = mirMatrix, filter = FALSE)

save(mirMatrix, mirExpr, file = "rawdata/mirexpr.RData")

############################### GEO 数据

####### GEO 数据预处理 GSE104131
library(GEOquery)
library(viper)
library(annotables)

data("grch37")

annotiondata <- dplyr::select(grch37, ensgene, symbol)

fpkm1 <- read_tsv("rawdata/mRNA/GSE104131/GSE104131_batch1.cufflinks.fpkms_tracking.TableFormatted.txt/GSE104131_batch1.cufflinks.fpkms_tracking.TableFormatted.txt")
fpkm2 <- read_tsv("rawdata/mRNA/GSE104131/GSE104131_batch2.cufflinks.fpkms_tracking.TableFormatted.txt/GSE104131_batch2.cufflinks.fpkms_tracking.TableFormatted.txt")

genenames <- fpkm1$GID

fpkm1 <- dplyr::select(fpkm1, -GENESYMBOL, -Locus) %>% 
         group_by(GID) %>% 
         summarise_all(mean)


rownames(fpkm1) <- fpkm1$GID

fpkm1 <- dplyr::select(fpkm1, -GID)
names(fpkm1) <- str_c("Sample_", names(fpkm1), sep = "")

gse1 <- getGEO(filename="rawdata/mRNA/GSE104131/GSE104131_series_matrix.txt.gz", getGPL = F)

temp <- gse1@phenoData@data %>% 
        dplyr::select(title, geo_accession, characteristics_ch1) %>% 
        data.frame(stringsAsFactors = F)

temp$title <- as.character(temp$title)

temp$title[11] <- "Sample_881_B19"
temp$title[12] <- "Sample_881_B26"

intersect(temp$title, names(fpkm1))

sampleldescript <- mutate(temp, 
                          type = ifelse(grepl("normal",characteristics_ch1), "Normal", "Tumor")) %>% 
  dplyr::select(title, type)

save(fpkm1,sampleldescript, file = "rawdata/mRNA_GSE104131_clinin.RData")

load("rawdata/mRNA_GSE104131_clinin.RData")


####### GEO 数据预处理GSE89223

filenames <- list.files("rawdata/mRNA/GSE89223/GSE89223_RAW (2)/")

allpaths <- paste("rawdata/mRNA/GSE89223/GSE89223_RAW (2)/", filenames, sep = "")

explist <- list()
for (i in 1:length(filenames )) {
  explist[[i]] <- read_tsv(allpaths[i],  col_names = F)
}

neexplist <- lapply(explist, function(x){x[2]})

exprdata <- as.matrix(bind_cols(neexplist))

colnames(exprdata) <- str_remove(filenames, "_.*")
rownames(exprdata) <- explist[[1]]$X1

gsepath <- paste("rawdata/mRNA/GSE89223/GSE89223_series_matrix.txt.gz")

gse1 <- getGEO(filename = gsepath, getGPL = F)

temp <- gse1@phenoData@data %>% 
  dplyr::select(geo_accession, characteristics_ch1.1) %>% 
  data.frame(stringsAsFactors = F)

sampleldescript <- mutate(temp, 
                          type = ifelse(grepl("normal",characteristics_ch1.1), "Normal", "Tumor")) %>% 
  dplyr::select(title = geo_accession, type)

save(exprdata, sampleldescript, file = "rawdata/mRNA_GSE89223_clinin.RData")


load("rawdata/mRNA_GSE89223_clinin.RData")

####### GEO 数据miRNA处理GSE21036

gsepath <- paste("rawdata/miRNA/GSE21036/GSE21036_series_matrix.txt.gz")
#gplpath <- paste("rawdata/miRNA/GSE21036/")
gse1 <- getGEO(filename = gsepath, getGPL = F)
#gpldata1 <- getGEO(filename=gplpath, getGPL = F)

exprdata <-gse1@assayData[["exprs"]] %>% 
           as.data.frame() %>% 
          dplyr::select(-GSM526532)

sampleldescript <- gse1@phenoData@data %>% 
                   filter(geo_accession %in% names(exprdata)) %>% 
           mutate( type = ifelse(grepl("normal",title), "Normal", "Tumor")) %>% 
  dplyr::select(title = geo_accession, type)


save(exprdata, sampleldescript, file = "rawdata/miRNA_GSE21036_cline.RData")


####### GEO 数据miRNA处理GSE76260
gsepath <- paste("rawdata/miRNA/GSE76260/GSE76260_series_matrix.txt.gz")
gplpath <- read_tsv("rawdata/miRNA/GSE76260/GPL8179_humanMI_V2_R0_XS0000124-MAP.txt/GPL8179_humanMI_V2_R0_XS0000124-MAP.txt")
gse1 <- getGEO(filename = gsepath, getGPL = F)

exprdata <-gse1@assayData[["exprs"]] 

sampleldescript <- gse1@phenoData@data %>%  
  mutate( type = ifelse(grepl("normal",title), "Normal", "Tumor")) %>% 
  dplyr::select(title = geo_accession, type)


save(exprdata, sampleldescript, file = "rawdata/miRNA_GSE76260_clinic.RData")

load("rawdata/geo_lich_featureANDexp.RData")


