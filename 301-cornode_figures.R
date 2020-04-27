
library(Vennerable)
library(survminer)
library(survival)
rm(list = ls())
mydata <- read_csv("result/network/pvalueDF.csv")
load("rawdata/for_figure_plot.RData")

resultnode <- dplyr::filter(mydata, pvalue < 0.05, auc > 0.8, p.value < 0.05)

AUCgene <- dplyr::filter(mydata, auc > 0.8)

KMgne <- dplyr::filter(mydata, pvalue < 0.05)
COXgne <- dplyr::filter(mydata, p.value < 0.05)

######## 韦恩图
genelist <- list(AUC = AUCgene$name,
                 'COX regression' = COXgne$name,
                 'Log Rank-test' = KMgne$name
                 )

Vsteam <- Venn(genelist)

tiff("result/Figures/Venn.tiff", width = 2000, height = 2000, units = 'px', res = 300)
plot(Vsteam, doWeights = FALSE)
dev.off()



####
all_edges <- readxl::read_xlsx("result/network/network.xlsx") %>% 
  filter(fromNode %in% resultnode$name, toNode %in% resultnode$name)

finalnoudes <- c(all_edges$fromNode, all_edges$toNode) %>% 
               unique()

forplot <- dplyr::filter(mydata, name %in% finalnoudes)

joindata <- mutate(joindata, 
                   meth = ifelse(shortLetterCode == "TP", 1, 0))

##################### KM plot

####### 生存期分析

raw_clinic_data <- pheno



for (i in 1:nrow(forplot )){
  
  tmpsigure <- bind_cols(raw_clinic_data[, 3:4], t_exp[forplot$name[i]]) %>% 
               na.omit()
  
  
  tmpsigure$group <- cut(unlist(tmpsigure[, 3]), 
                         c(min(tmpsigure[, 3]), forplot$cut[i], max(tmpsigure[, 3])),
                         include.lowest = TRUE)
  
  tmpsigure$group <- factor(tmpsigure$group,
                            labels=c("Low", 
                                     "High"))
  
  
  fit <- survfit(Surv(DFDtime, DFSstatus) ~ group, 
                 data = tmpsigure)
  pfit1 <- ggsurvplot(fit,
                      palette = c("blue", "red"),
                      data = tmpsigure, 
                      #color = c("blue", "red"),
                      censor.shape="|", 
                      pval = TRUE,
                      censor.size = 3,
                      break.time.by = 30,
                      legend.title = paste(forplot$name[i], ":", sep = ""),
                      
                      #legend = "left",
                      #legend.labs = c("Advanced", "Early"),
                      xlab = "Time in months", 
                      
                      legend.labs = c("Low", "High"),
                      pval.size = 7,
                      pval.coord = c(0,0.1),
                      #ggtheme = theme_bw(),
                      
                      font.legend = 18,
                      font.x =  c(18, "bold", "black"),
                      font.y = c(18, "bold", "black"),
                      font.tickslab = c(15, "bold", "black")
                      # risk.table = TRUE,
                      # risk.table.col = "strata",
                      # risk.table.height = 0.25
                      #fontsize = 10
  )
  
  
  # ggsave(plot = print(pfit1), filename = paste(forplot$name[i], ".tiff", sep = "") )
  
  ggpubr::ggexport(plot = pfit1, 
                   filename = paste(forplot$name[i], ".tiff", sep = ""), 
                   device = "tiff", 
                   res = 300,
                   width = 2200, 
                   height = 2200)
  
}


######### 差异比较
library(ggpubr)

finaldata <- dplyr::filter(joindata, shortLetterCode %in% c("NT", "TP"))

finaldata$State <- factor(finaldata$shortLetterCode, 
                          levels = c("NT", "TP"), 
                          labels = c("Normal Tissue","Primary Tumor"))

raw_clinic_data <- finaldata

for (i in 1:nrow(forplot )){
  
  tmpsigure <- data.frame(State = raw_clinic_data[, "State"], 
                          rawexp = raw_clinic_data[,forplot$name[i]])
  
  tmpsigure$exp <- log10(tmpsigure$rawexp + abs(min(tmpsigure$rawexp)) + 1)
  
  my_comparisons <- list( c("Normal Tissue","Primary Tumor") )
  
  ylables <- max(tmpsigure$exp)
  
  p1 <- ggboxplot(tmpsigure, 
                 x = "State", 
                 y = "exp",
                 shape="State",
                 color = "State",
                 bxp.errorbar = T,
                 bxp.errorbar.width = 0.2,
                 size = 0.5,
                 #fill = "State",
                 palette = "aaas",
                 #font.label = list(size = 50, face = "bold"),
                 add = c("jitter"), 
                 font.main = c(14,"bold"),
                 font.x = c(16, "bold"),
                 font.y = c(16, "bold"),
                 font.legend = c(14, "bold"),
                 font.tickslab = c(14,"bold"),
                 ylab = paste("Log10  of ", forplot$name[i], " Expression") )+
    # Add significance levels
    stat_compare_means(method = "t.test",
                       label.y = ylables + 0.05,
                       size = 6) 
  
  # p <- ggboxplot(df1, x="dose", y="len", color = "dose",
  #                palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  #                add = "jitter", shape="dose")
  

  
  
  ggpubr::ggexport(plot = p1, 
                   filename = paste(forplot$name[i], "comp.tiff", sep = ""), 
                   device = "tiff", 
                   res = 300,
                   width = 2200, 
                   height = 2200)
  
}

# finaldata$PIK3CG <- log10(finaldata$PIK3CG + 1)
# finaldata$CXCR4 <- log10(finaldata$CXCR4 + 1)
# finaldata$IL7R <- log10(finaldata$IL7R + 1)

my_comparisons <- list( c("NT", "TP") )

p1 <- ggviolin(finaldata, 
               x = "State", 
               y = forplot$name[i],
               fill = "State",
               palette = c("#00AFBB", "#E7B800"),
               add = "boxplot", 
               add.params = list(fill = "white"),
               ylab = paste("Log10  of ", forplot$name[i], " Expression") )+
  stat_compare_means(comparisons = my_comparisons,
                     
                     label = "p.signif")+ # Add significance levels
  stat_compare_means(method = "t.test", label.y = 4)  


ggpubr::ggexport(plot = p1, 
                 filename = paste(forplot$name[i], "comp.tiff", sep = ""), 
                 device = "tiff", 
                 res = 300,
                 width = 2000, 
                 height = 2000)


#### 画AUC图
library(ROCR)
library(pROC)

# 总图
tiff("PlotAUC.tiff", width = 1800, height = 1500, units = 'px', res = 300)
roc1 <- roc(finaldata$meth, finaldata$`MAGI2-AS3`)   
roc2 <- roc(finaldata$meth, finaldata$SNHG3) 
roc3 <- roc(finaldata$meth, finaldata$`RP11-166D19.1`)
plot(roc1, col="blue")  
plot.roc(roc2, add=TRUE, col="red")
plot.roc(roc3, add=TRUE, col="green")


legend(0.65, 0.3, 
       c("MAGI2-AS3 (AUC=0.91)", "SNHG3 (AUC=0.92)", "RP11-166D19.1 (AUC=0.83)"), 
       col = c("blue", "red", "green"),
       text.col = c("blue", "red", "green"), 
       lty = c(1, 1, 1), 
       lwd = 2,
       pch = NA,
       merge = TRUE, bg = "gray90",
       cex = 1 )

dev.off()

finaldata2 <- finaldata[, c(forplot$name, "sample", "meth")]

tiff("result/Figures/PlotAUCmirna.tiff", width = 1800, height = 1500, units = 'px', res = 300)
roc1 <- roc(finaldata$meth, finaldata2$`hsa-miR-221-5p`)   
roc2 <- roc(finaldata$meth, finaldata2$`hsa-miR-222-3p`) 
plot(roc1, col="blue")  
plot.roc(roc2, add=TRUE, col="red")

legend(0.65, 0.3, 
       c("hsa-miR-221-5p (AUC=0.83)", "hsa-miR-222-3p (AUC=0.89)"), 
       col = c("blue", "red"),
       text.col = c("blue", "red"), 
       lty = c(1, 1, 1), 
       lwd = 2,
       pch = NA,
       merge = TRUE, bg = "gray90",
       cex = 1 )

dev.off()



##### mRNA

tiff("result/Figures/PlotAUCmRNA.tiff", width = 1800, height = 1500, units = 'px', res = 300)
roc1 <- roc(finaldata$meth, finaldata2$AOX1)   
roc2 <- roc(finaldata$meth, finaldata2$ALDH1A2) 
roc3 <- roc(finaldata$meth, finaldata2$ATP2B4)   
roc4 <- roc(finaldata$meth, finaldata2$TSPAN1) 
roc5 <- roc(finaldata$meth, finaldata2$SLC43A1)   
plot(roc1, col= rainbow(6)[1])  
plot.roc(roc2, add=TRUE, col=rainbow(6)[6])
plot.roc(roc3, add=TRUE, col=rainbow(6)[3])
plot.roc(roc4, add=TRUE, col=rainbow(6)[4])
plot.roc(roc5, add=TRUE, col=rainbow(6)[5])

legend(0.5, 0.4, 
       c("AOX1 (AUC=0.93)", "ALDH1A2 (AUC=0.86)", "ATP2B4 (AUC=0.87)", "TSPAN1 (AUC=0.82)", "SLC43A1 (AUC=0.86)"), 
       col = rainbow(6)[c(1,6,3,4,5)],
       text.col = rainbow(6)[c(1,6,3,4,5)], 
       lty = c(1, 1, 1, 1, 1), 
       lwd = 2,
       pch = NA,
       merge = TRUE, bg = "gray90",
       cex = 1 )

dev.off()




######### 不同分期
library(ggpubr)

finaldata$tumor_stage2 <- str_remove(finaldata$tumor_stage, "stage ") %>% 
  str_remove("a")%>% 
  str_remove("b")%>% 
  str_remove("c")

finaldata$tumor_stage2[finaldata$tumor_stage2 == "i/ii nos"] <- "i"

finaldata <- dplyr::filter(finaldata, tumor_stage2 != "not reported") %>%               dplyr::filter( tumor_stage2 != "0")


finaldata <- dplyr::filter(finaldata, shortLetterCode %in% c("TM", "TP"))

finaldata$Stage <- factor(finaldata$tumor_stage2, 
                          levels = c("i", "ii", "iii", "iv"), 
                          labels = c("Stage I","Stage II","Stage III","Stage IV"))



siguregene <- c(names(joindata)[8:67], names(joindata)[70:120])

pvalulist <- c()

for (j in 1:length(siguregene)){
  
  
  
  tmpsigure <- data.frame(Stage = raw_clinic_data[, "Stage"], 
                          rawexp = raw_clinic_data[,siguregene[j]])
  
  tmpsigure$exp <- log10(tmpsigure$rawexp + abs(min(tmpsigure$rawexp)) + 1)
  
  aaa <- summary(aov(exp~Stage, data = tmpsigure))
  
  pvalulist[j] <- aaa[[1]]$`Pr(>F)`[1]
  
}

newdf <- data.frame(name = siguregene, anovP= pvalulist) %>% 
  inner_join(selectnodes, by = "name") %>% 
  dplyr::filter(pvalue.x < 0.05, pvalue.y < 0.05, pvalulist < 0.05)

cenet <- read.csv("script/cerna/ceRNAedges.csv", stringsAsFactors = F) %>% 
  filter(fromNode %in% newdf$name, toNode %in% newdf$name)


aaa <- summary(aov(hsa-miR-141~Stage, data = finaldata))

aaa[[1]]$`Pr(>F)`[1]

raw_clinic_data <- finaldata

for (i in 1:nrow(forplot )){
  
  tmpsigure <- data.frame(Stage = raw_clinic_data[, "Stage"], 
                          rawexp = raw_clinic_data[,forplot$name[i]])
  
  tmpsigure$exp <- log10(tmpsigure$rawexp + abs(min(tmpsigure$rawexp)) + 1)
  
  
  my_comparisons <- list( c("Stage I","Stage II") ,c("Stage I","Stage III"), c("Stage I","Stage IV"))
  
  
  
  p1 <- ggviolin(tmpsigure, 
                 x = "Stage", 
                 y = "exp",
                 fill = "Stage",
                 palette = get_palette("npg", 4),
                 add = "boxplot", 
                 add.params = list(fill = "white"),
                 ylab = paste("Log10  of ", forplot$name[i], " Expression") )+
    stat_compare_means(comparisons = my_comparisons, label = "p.format")+ # Add significance levels
    stat_compare_means( method = "anova",  label.y = 0.6 )  
  
  
  
  
  ggpubr::ggexport(plot = p1, 
                   filename = paste(forplot$name[i], "comp.tiff", sep = ""), 
                   device = "tiff", 
                   res = 300,
                   width = 2000, 
                   height = 2000)
}




#####################################
### 森林图
#################
library(ggplot2)
forplot <- mutate(forplot, group = ifelse(HR>1, "high", "Low")) 
  

p<-ggplot(forplot,aes(x=name,y=HR,group=group, color=group))+
  #设置图形的数据集和xy轴变量，分组变量
  geom_errorbar(aes(ymin=HRCILL, ymax=HRCIUL), width=.4,size=0.9) +
  #设置误差图，误差图宽度，间距，粗细
  geom_point(size=4)+
  #设置误差图中心的点估计
  geom_hline(aes(yintercept=1), colour="gray50", linetype="dashed",size=1)+
  scale_y_continuous("Hazard ratios",limits=c(0.25,1.75))+
  xlab("Names")+
  coord_flip()+
  theme_bw()+
  #guides(fill=FALSE)+
  scale_color_manual(values=c('red','blue'))+
  theme(panel.border = element_rect(size = 1.5,fill = NA),
        legend.position="none",
        axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 12,colour = "black"),
        axis.title = element_text(face = "bold", size = 15)
        )
  
  ggsave("result/Figures/forestplot.tiff", p, width = 5, height = 10, dpi = 300)


#######################
  
  lncrnas <- filter(forplot, type == "lnc")

compete_relation <- read_csv("result/gdcceRNAout.csv") %>% 
                    dplyr::filter(lncRNAs %in% lncrnas$name)



count <- as.matrix(table(filterclinical$Nstage, filterclinical$Tstage))

matdata <- as.matrix(unclass(count))

fisher.test(matdata, simulate.p.value = T)





