setwd("~/projects/ACE2/FPKM-UQ_processTCGA2/ACE2_expr_cor/")
# input all gene expression data (FPKM)
expr <- readRDS("Data_TCGALungCancerSampleExpression.rds")
head(expr)
names(expr) <- c("gene","name","gene_id","expr","type","case","sample","tissue")
# calculate the median expression of each gene in tumor samples and normal samples
mean_expr <- as.data.frame(expr %>% group_by(name,tissue,type) %>% summarise(mean=median(expr),na.rm=T))
write.table(mean_expr,"median_expr_all_gene_lung_cancer.txt")
# input median expression of all genes
mean_expr <- read.table("~/projects/ACE2/FPKM-UQ_processTCGA2/ACE2_expr_cor/median_expr_all_gene_lung_cancer.txt")
head(mean_expr)
# conver data frame from long to wide
library(reshape2)
mean_con <- dcast(mean_expr, name+type  ~ tissue, value.var="mean")
head(mean_con)
# plot median expression normal vs tumor in LUAD
p1 <- ggplot()+geom_point(data=subset(mean_con,type=="TCGA-LUAD"),aes(y=log2(`Primary Tumor`+1),
                                      x=log2(`Solid Tissue Normal`+1)),alpha=0.2,size=0.2,colour="#104E8B")+mytheme_violin+
  xlab("Expression level in normal tissue (log2(FPKM+1))")+ylab("Expression level in tumor (log2(FPKM+1))")+
  xlim(0,30)+ylim(0,30)+ geom_abline(intercept=0, slope=1)+ggtitle("LUAD RHOV and ACE2")+
  geom_point(data=subset(mean_con,type=="TCGA-LUAD" & name=="ACE2"),aes(y=log2(`Primary Tumor`+1),
                                                                        x=log2(`Solid Tissue Normal`+1)),colour="#CD0000",size=0.3)+
  geom_point(data=subset(mean_con,type=="TCGA-LUAD" & name=="RHOV"),aes(y=log2(`Primary Tumor`+1),
                                                                            x=log2(`Solid Tissue Normal`+1)),size=0.3,colour="#CD0000")
p1
# plot median expression normal vs tumor in LUSC
p2 <- ggplot()+geom_point(data=subset(mean_con,type=="TCGA-LUSC"),aes(y=log2(`Primary Tumor`+1),
                                                                x=log2(`Solid Tissue Normal`+1)),alpha=0.2,size=0.2,colour="#104E8B")+mytheme_violin+
  xlab("Expression level in normal tissue (log2(FPKM+1))")+ylab("Expression level in tumor (log2(FPKM+1))")+
  xlim(0,30)+ylim(0,30)+ geom_abline(intercept=0, slope=1)+ggtitle("LUSC RHOV and ACE2")+
  geom_point(data=subset(mean_con,type=="TCGA-LUSC" & name=="ACE2"),aes(y=log2(`Primary Tumor`+1),
                                                         x=log2(`Solid Tissue Normal`+1)),colour="#CD0000",size=0.3)+
  geom_point(data=subset(mean_con,type=="TCGA-LUSC" & name=="RHOV"),aes(y=log2(`Primary Tumor`+1),
                                                                          x=log2(`Solid Tissue Normal`+1)),size=0.3,colour="#CD0000")
p2

library("grid")
pdf("gene_expr_cor_lung_cancer.pdf",width=20,height=20)
print(p1,vp=viewport(.09,.095,x=.2,y=.9))
print(p2,vp=viewport(.09,.095,x=.33,y=.9))
dev.off()
