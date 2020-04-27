rm(list = ls())
library(annotables)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)

clinical <- read_csv("rawdata/clinical_information.csv")
filterclinical <- clinical %>% 
                  na.omit()
load( "PRAD_RNA_and_clinical.RData")
rm(lihc_ranexp2)
allnodes <- read.csv("result/network/nodes.csv", stringsAsFactors = F)

load("result/exprofile.RData")

exprofile2 <- column_to_rownames(exprofile, "sampleID")[, allnodes$symbol] %>% 
              rownames_to_column("sampleID")

siguregene <- colnames(exprofile2)[-1]
################### AUC
clinical2$sampleID <- str_sub(clinical2$sample,  1, -2)
joindata <- inner_join(exprofile2, clinical2, by = c("sampleID" = "sampleID"))

library(ROCR)

joindata <- mutate(joindata, 
                   meth = ifelse(shortLetterCode == "TP", 1, 0))

auclist <- c()

direct <- allnodes$dir 
names(direct) <- allnodes$symbol

for (j in 1:length(siguregene)) {
  if (direct[siguregene[j]] == "up") {
    pred <- prediction(joindata[[siguregene[j]]], joindata$meth)
  } else 
  {pred <- prediction(-joindata[[siguregene[j]]], joindata$meth)}
  perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
  
  auc <- performance(pred,'auc')
  auclist[j]=unlist(slot(auc,"y.values"))
  
}

aucdf <- data.frame(name = siguregene,
                    auc = auclist) %>%
  inner_join(allnodes, by = c("name" = "symbol")) 

###################### COXPH 回归

#### 临床数据预处理
tumorsample <- dplyr::filter(joindata, shortLetterCode == "TP")

joindata <- dplyr::filter(exprofile2, sampleID %in% tumorsample$sampleID) %>% 
  inner_join(clinical, by = c("sampleID" = "Sample ID"))

## OS.time
joindata$DFDtime <- joindata$`Disease Free (Months)`

## status
joindata$DFSstatus <- ifelse(joindata$`Disease Free Status` == "DiseaseFree", 0, 1)

## age

joindata$age <- joindata$`Diagnosis Age`

## T_stage
library(stringr) 

joindata$Tstage <- str_remove(joindata$`American Joint Committee on Cancer Tumor Stage Code`, "[abc]")

### Nstage
joindata$Nstage <- joindata$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`


## bind your clinical data
pheno <- joindata[,c("sampleID",
                     "age",
                  "DFDtime",
                  "DFSstatus",
                  "Tstage",
                  "Nstage")]
t_exp <- joindata[, allnodes$symbol]

save(pheno, t_exp, joindata, file = "rawdata/for_figure_plot.RData")

multi_results <- apply(t_exp , 2, function(gene){
  ## gene <- t_exp[1, ]
  group <- unlist(gene)
  survival_dat <- data.frame(expr = group,
                             age = pheno$age,
                             DFDtime = pheno$DFDtime,
                             DFSstatus = pheno$DFSstatus, 
                               Tstage = pheno$Tstage,
                             Nstage = pheno$Nstage,
                             stringsAsFactors = F)
  
  res.cox <- coxph(Surv(DFDtime, DFSstatus) ~ age + Tstage + Nstage + expr, 
                   data =  survival_dat)
  ## summary(res.cox)
  beta <- coef(res.cox)
  se <- sqrt(diag(vcov(res.cox)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  res <- as.data.frame(round(cbind(coef = beta,
                                   se = se,
                                   z = beta/se,
                                   p.value = 1 - pchisq((beta/se)^2, 1),
                                   HR = HR,
                                   HRse = HRse,
                                   HRz = (HR - 1) / HRse,
                                   HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                   HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                                   HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 5))

})

multi_results2 <- do.call(rbind, multi_results) %>% 
                  rownames_to_column("IDs") %>% 
                  dplyr::filter(grepl(".expr", IDs))
  

multi_results2$IDs <- str_remove(multi_results2$IDs, ".expr")

aucdf2 <- inner_join(aucdf, multi_results2, by = c("name" = "IDs"))



################ KM 

####### 

#raw_clinic_data$FreeSurvival <- as.numeric(raw_clinic_data$FreeSurvival)
getcutoff <- function(x, 
                      vv = c(0,0.25,0.5,0.75,1), 
                      down=2, 
                      up=4){
  upcut <- quantile(x, vv)[up]
  downcut <- quantile(x, vv)[down]
  
  y <- x[(x>downcut)&(x<upcut)]
  return(y)
}

siguregene <- names(joindata)[8:67]
##### 确定最优阈值

siguregene <- colnames(t_exp)

survival_dat <- data.frame(DFDtime = pheno$DFDtime,
                           DFSstatus = pheno$DFSstatus, 
                           stringsAsFactors = F)

pvaluelist <- c()
cut_remian <- c()

for (i in 1:length(siguregene)){
  
  tmpsigure <- bind_cols(survival_dat, t_exp[siguregene[i]])
  
  all_cutoff <- getcutoff(tmpsigure[[3]])
  
  tmpsurpvalue <- c()
  
  for (j in 1:length(all_cutoff)) {
    
    tmpsigure$group <- cut(unlist(tmpsigure[, 3]), 
                           c(min(tmpsigure[, 3]), all_cutoff[j], max(tmpsigure[, 3])),
                           include.lowest = TRUE)
    
    tmpsigure$group <- factor(tmpsigure$group,
                              labels=c("Low", 
                                       "High"))
    
    
    fittedSurv <- survfit(Surv(time=DFDtime, event=DFSstatus) ~ group, 
                          data=tmpsigure)
    
    diffSurv <- survdiff(Surv(time=DFDtime, event=DFSstatus) ~ group,           
                         data=tmpsigure)
    
    tmpsurpvalue[j] <- signif(1 - pchisq(diffSurv$chisq, length(diffSurv$n) - 1), 2)
    
    
  }
  
  pvaluelist[i] <- min(tmpsurpvalue)
  cut_remian[i] <- all_cutoff[which.min(tmpsurpvalue)]
  
}

pvaludf <- data.frame(name = siguregene, 
                      pvalue = pvaluelist,
                      cut = cut_remian) %>%
  inner_join(aucdf2, by = "name")



resultnode <- dplyr::filter(pvaludf, pvalue < 0.05, auc > 0.8, p.value < 0.05)

####
all_edges <- readxl::read_xlsx("result/network/network.xlsx") %>% 
             filter(fromNode %in% resultnode$name, toNode %in% resultnode$name)

write.csv(pvaludf, "result/network/pvalueDF.csv", quote = F, row.names = F)

write.csv(all_edges, "result/network/corenetwork.csv", quote = F, row.names = F)
write.csv(resultnode, "result/network/cornet_node.csv", quote = F, row.names = F)
