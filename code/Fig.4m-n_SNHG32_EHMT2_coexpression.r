setwd("/data/home/yuanmengqin/AD_GWAS/paper_eQTL_data/2019_Brain/RNA-seq/eQTL_lncRNA_coexpression")
exp=read.table("HIPPO_control_exp_rpkm.txt",sep="\t",header=T)
exp1=exp[-grep("_PAR_Y",rownames(exp)),]
exp_name=unlist(lapply(strsplit(rownames(exp1), ".", fixed = TRUE), "[", 1))
rownames(exp1)=exp_name
lnc=read.table("./SNHG32_SRG_pairs.txt",sep="\t",header=T)
out=c()
pdf("HIPPO_cor_ggplot.pdf",width=3,height=3.5)
for(i in 1:dim(lnc)[1]){
	df=as.data.frame(cbind(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],])))
	re=cor.test(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],]))
	re_sp=cor.test(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],]),method="spearman")
	out=rbind(out,c(lnc[i,1],lnc[i,2],re$estimate,re$p.value,re_sp$estimate,re_sp$p.value))
	colnames(df)=c("lnc","gene")
	p=ggplot(df,aes(x=lnc,y=gene))+
	geom_smooth(method="lm",  # 拟合曲线
              se=T,         # 是否添加置信区间
              color="#AEC5B5")+
  geom_point(colour="#0D5122",fill="#AEC5B5",shape=21,size=3)+             # 散点图
  #scale_color_manual(values=c("red"))+scale_fill_manual(values=c("#85A484"))+
  annotate(
    "text", x = min(df[,1],na.rm = T)+0.1, y = max(df[,2],na.rm = T)-0.1,                   # 坐标位置，左上角
    label = paste0("R=",round(re$estimate,2),                                   # 注释标记内容相关性+P值
                   ",P=",format(re$p.value, scientific = TRUE, digits = 4)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
	annotate(
    "text", x = min(df[,1],na.rm = T)+0.1, y = max(df[,2],na.rm = T)/2,                   # 坐标位置，左上角
    label = paste0("R=",round(re_sp$estimate,2),                                   # 注释标记内容相关性+P值
                   ",P=",format(re_sp$p.value, scientific = TRUE, digits = 4)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
	labs(x=lnc[i,1],y=lnc[i,2])+
  theme_classic()
  print(p)
}
dev.off()
write.table(out,"SNHG32_SRG_pairs_HIPPO_coexp.txt",sep="\t",quote=F,row.names=F)
####DLPFC
setwd("/data/home/yuanmengqin/AD_GWAS/paper_eQTL_data/2019_Brain/RNA-seq/eQTL_lncRNA_coexpression")
exp=read.table("DLPFC_control_exp_rpkm.txt",sep="\t",header=T)
exp1=exp[-grep("_PAR_Y",rownames(exp)),]
exp_name=unlist(lapply(strsplit(rownames(exp1), ".", fixed = TRUE), "[", 1))
rownames(exp1)=exp_name
lnc=read.table("./SNHG32_SRG_pairs.txt",sep="\t",header=T)
out=c()
pdf("DLPFC_cor_ggplot.pdf",width=3,height=3.5)
for(i in 1:dim(lnc)[1]){
	df=as.data.frame(cbind(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],])))
	re=cor.test(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],]))
	re_sp=cor.test(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],]),method="spearman")
	out=rbind(out,c(lnc[i,1],lnc[i,2],re$estimate,re$p.value,re_sp$estimate,re_sp$p.value))
	colnames(df)=c("lnc","gene")
	p=ggplot(df,aes(x=lnc,y=gene))+
	geom_smooth(method="lm",  # 拟合曲线
              se=T,         # 是否添加置信区间
              color="#AEC5B5")+
  geom_point(colour="#0D5122",fill="#AEC5B5",shape=21,size=3)+             # 散点图
  #scale_color_manual(values=c("red"))+scale_fill_manual(values=c("#85A484"))+
  annotate(
    "text", x = min(df[,1],na.rm = T)+0.1, y = max(df[,2],na.rm = T)-0.1,                   # 坐标位置，左上角
    label = paste0("R=",round(re$estimate,2),                                   # 注释标记内容相关性+P值
                   ",P=",format(re$p.value, scientific = TRUE, digits = 4)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
	annotate(
    "text", x = min(df[,1],na.rm = T)+0.1, y = max(df[,2],na.rm = T)/2,                   # 坐标位置，左上角
    label = paste0("R=",round(re_sp$estimate,2),                                   # 注释标记内容相关性+P值
                   ",P=",format(re_sp$p.value, scientific = TRUE, digits = 4)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
	labs(x=lnc[i,1],y=lnc[i,2])+
  theme_classic()
  print(p)
}
dev.off()
write.table(out,"SNHG32_SRG_pairs_DLPFC_coexp.txt",sep="\t",quote=F,row.names=F)
####caudate

setwd("/data/home/yuanmengqin/AD_GWAS/paper_eQTL_data/2022_NG/RNA-seq/coexpression")
data=read.table("../caudate_brainseq_phase3_hg38_rseGene_merged_n464_FPKM.txt",sep="\t")
con_sample=read.table("../caudate_Control_sample.txt",sep="\t",header=F)
data1=data[-grep("_PAR_Y",rownames(data)),]
exp_name=unlist(lapply(strsplit(rownames(data1), ".", fixed = TRUE), "[", 1))
rownames(data1)=exp_name
exp1=data1[,con_sample[,1]]
lnc=read.table("./SNHG32_SRG_pairs.txt",sep="\t",header=T)
out=c()
pdf("Caudate_cor_ggplot.pdf",width=3,height=3.5)
for(i in 1:dim(lnc)[1]){
	df=as.data.frame(cbind(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],])))
	re=cor.test(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],]))
	re_sp=cor.test(as.numeric(exp1[lnc[i,1],]),as.numeric(exp1[lnc[i,2],]),method="spearman")
	out=rbind(out,c(lnc[i,1],lnc[i,2],re$estimate,re$p.value,re_sp$estimate,re_sp$p.value))
	colnames(df)=c("lnc","gene")
	p=ggplot(df,aes(x=lnc,y=gene))+
	geom_smooth(method="lm",  # 拟合曲线
              se=T,         # 是否添加置信区间
              color="#AEC5B5")+
  geom_point(colour="#0D5122",fill="#AEC5B5",shape=21,size=3)+             # 散点图
  #scale_color_manual(values=c("red"))+scale_fill_manual(values=c("#85A484"))+
  annotate(
    "text", x = min(df[,1],na.rm = T)+0.1, y = max(df[,2],na.rm = T)-0.1,                   # 坐标位置，左上角
    label = paste0("R=",round(re$estimate,2),                                   # 注释标记内容相关性+P值
                   ",P=",format(re$p.value, scientific = TRUE, digits = 4)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
	annotate(
    "text", x = min(df[,1],na.rm = T)+0.1, y = max(df[,2],na.rm = T)/2,                   # 坐标位置，左上角
    label = paste0("R=",round(re_sp$estimate,2),                                   # 注释标记内容相关性+P值
                   ",P=",format(re_sp$p.value, scientific = TRUE, digits = 4)),   # P值采用科学计数法，保留2位小数
    size = 5,hjust=0,vjust=1)+
	labs(x=lnc[i,1],y=lnc[i,2])+
  theme_classic()
  print(p)
}
dev.off()
write.table(out,"SNHG32_SRG_pairs_caudate_coexp.txt",sep="\t",quote=F,row.names=F)