library(metafor)
library(ggplot2)
dat = get(data(dat.raudenbush1985))
random = rma(yi, sei=se, data=dat)


library(ggplot2)
p<-ggplot(dat,aes(x=group,y=or,group=group))+
  #设置图形的数据集和xy轴变量，分组变量
  geom_errorbarh(aes(ymin=low, ymax=up), width=.2,size=0.5) +
  #设置误差图，误差图宽度，间距，粗细
  geom_point(size=3)+
  #设置误差图中心的点估计
  geom_hline(aes(yintercept=1), colour="black", linetype="dashed",size=1)+
  scale_y_continuous("Adjusted Odds Ratios and 95% CIs",limits=c(-1,3))
#设置x和y轴
p
