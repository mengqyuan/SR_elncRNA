#module ava
#module load r/4.0

library(ggplot2)

# 读取数据
lncRNA <- read.table("ATAC_seq_CAUD.txt",
                     head = TRUE, sep = "\t", quote = "", as.is = TRUE)

# 转换为因子变量
lncRNA$Group <- factor(lncRNA$Group, levels = c("WBG", "LBG", "ES_lncRNA","Brain_elncRNA", "SR_elncRNA"))

# 筛选相关组
temp <- lncRNA$Group %in% c("WBG", "LBG", "ES_lncRNA","Brain_elncRNA", "SR_elncRNA")
use <- lncRNA[temp, ]

# 绘制箱线图
p <- ggplot(use, aes(x = Group, y = ATAC_signial, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# 计算箱线图统计数据以确定Y轴范围
stats <- boxplot.stats(use$ATAC_signial)$stats
# 将Y轴范围设置为从最小值到最大值的2倍
ylim1 <- range(stats[1], stats[5]* 2)

# 添加Y轴范围限制
p1 <- p + coord_cartesian(ylim = ylim1)

# 输出为PDF文件
ggsave("ATAC_seq_CAUD.pdf", plot = p1, width = 5, height = 5)


# 读取数据
lncRNA <- read.table("/data/home/liuhaizhou/AD_GWAS/ATAC-seq/GSE14672/plot/ATAC_seq_CAUD.txt",
                     head = TRUE, sep = "\t", quote = "", as.is = TRUE)

lncRNA$Group<- factor(lncRNA$Group, level = c("SR_elncRNA", "WBG"))
res1 <- wilcox.test(ATAC_signial ~ Group, data = lncRNA, alternative = "greater") #p-value=8.565298e-72

# 读取数据
lncRNA <- read.table("/data/home/liuhaizhou/AD_GWAS/ATAC-seq/GSE14672/plot/ATAC_seq_CAUD.txt",
                     head = TRUE, sep = "\t", quote = "", as.is = TRUE)

lncRNA$Group<- factor(lncRNA$Group, level = c("SR_elncRNA", "LBG"))
res2 <- wilcox.test(ATAC_signial ~ Group, data = lncRNA, alternative = "greater") #p-value=3.632009e-43

# 读取数据
lncRNA <- read.table("/data/home/liuhaizhou/AD_GWAS/ATAC-seq/GSE14672/plot/ATAC_seq_CAUD.txt",
                     head = TRUE, sep = "\t", quote = "", as.is = TRUE)

lncRNA$Group<- factor(lncRNA$Group, level = c("SR_elncRNA", "ES_lncRNA"))
res3 <- wilcox.test(ATAC_signial ~ Group, data = lncRNA, alternative = "greater") #p-value=8.548155e-33

# 读取数据
lncRNA <- read.table("/data/home/liuhaizhou/AD_GWAS/ATAC-seq/GSE14672/plot/ATAC_seq_CAUD.txt",
                     head = TRUE, sep = "\t", quote = "", as.is = TRUE)

lncRNA$Group<- factor(lncRNA$Group, level = c("SR_elncRNA", "Brain_elncRNA"))
res4 <- wilcox.test(ATAC_signial ~ Group, data = lncRNA, alternative = "greater") #p-value=7.06251e-30

