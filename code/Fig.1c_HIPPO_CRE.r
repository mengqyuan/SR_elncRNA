module ava
module load r/4.0

library(ggplot2)

# 读取数据
lncRNA <- read.table("multi_omics_elncRNA_6_groups.txt",
                     head = TRUE, sep = "\t", quote = "", as.is = TRUE)

# 转换为因子变量
lncRNA$Group <- factor(lncRNA$Group, levels = c("WBG lncRNA", "LBG lncRNA", "eQTL lncRNA","Brain elncRNA", "SCZ elncRNA"))

# 筛选相关组
temp <- lncRNA$Group %in% c("WBG lncRNA", "LBG lncRNA", "eQTL lncRNA", "Brain elncRNA", "SCZ elncRNA")
use <- lncRNA[temp, ]


# 绘制箱线图
p <- ggplot(use, aes(x = Group, y = tss_HIPPO_pcHIC_po, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

# 计算箱线图统计数据以确定Y轴范围
stats <- boxplot.stats(use$tss_HIPPO_pcHIC_po)$stats
# 将Y轴范围设置为从最小值到最大值的2倍
ylim1 <- range(stats[1], stats[5]* 3)

# 添加Y轴范围限制
p1 <- p + coord_cartesian(ylim = ylim1)

# 输出为PDF文件
ggsave("SCZ elncRNA_tss_HIPPO_pcHIC_po.pdf", plot = p1, width = 4, height = 5)
---------------------------------------------------------------------------------------------------------------------
lncRNA$Group<- factor(lncRNA$Group, level = c("SCZ elncRNA", "WBG lncRNA"))
res1 <- wilcox.test(ACTAC_HIPPO1 c, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.002954
res2 <- wilcox.test(ACTAC_HIPPO2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.00149
res3 <- wilcox.test(ACTAC_HIPPO3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.001005
res4 <- wilcox.test(ACTAC_Caudate1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0004262
res5 <- wilcox.test(ACTAC_Caudate2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0006941
res6 <- wilcox.test(ACTAC_Caudate3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0005202
res7 <- wilcox.test(tss_CP_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 1.969e-07
res8 <- wilcox.test(tss_GZ_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 8.431e-08
res9 <- wilcox.test(tss_capHIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.02408
res10 <- wilcox.test(tss_Fantom5 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.9164
res11 <- wilcox.test(tss_DLPFC_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 9.826e-05
res12 <- wilcox.test(tss_HIPPO_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 6.505e-05
res13 <- wilcox.test(HIPPO_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.02066
res14 <- wilcox.test(DLPFC_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0006254
res15 <- wilcox.test(CAUDATE_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.002221
---------------------------------------------------------------------------------------------------------------------
lncRNA$Group<- factor(lncRNA$Group, level = c("SCZ elncRNA", "LBG lncRNA"))
res1 <- wilcox.test(ACTAC_HIPPO1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.02004
res2 <- wilcox.test(ACTAC_HIPPO2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01598
res3 <- wilcox.test(ACTAC_HIPPO3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01175
res4 <- wilcox.test(ACTAC_Caudate1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.005427
res5 <- wilcox.test(ACTAC_Caudate2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.009743
res6 <- wilcox.test(ACTAC_Caudate3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.007023
res7 <- wilcox.test(tss_CP_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 4.85e-06
res8 <- wilcox.test(tss_GZ_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 2.309e-06
res9 <- wilcox.test(tss_capHIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.05499
res10 <- wilcox.test(tss_Fantom5 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.9282
res11 <- wilcox.test(tss_DLPFC_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.002539
res12 <- wilcox.test(tss_HIPPO_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.001565
res13 <- wilcox.test(HIPPO_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.04961
res14 <- wilcox.test(DLPFC_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.001369
res15 <- wilcox.test(CAUDATE_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.005996
---------------------------------------------------------------------------------------------------------------------
lncRNA$Group<- factor(lncRNA$Group, level = c("SCZ elncRNA", "eQTL lncRNA"))
res1 <- wilcox.test(ACTAC_HIPPO1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.04684
res2 <- wilcox.test(ACTAC_HIPPO2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.04979
res3 <- wilcox.test(ACTAC_HIPPO3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0498
res4 <- wilcox.test(ACTAC_Caudate1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.02668
res5 <- wilcox.test(ACTAC_Caudate2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.02271
res6 <- wilcox.test(ACTAC_Caudate3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0227
res7 <- wilcox.test(tss_CP_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01405
res8 <- wilcox.test(tss_GZ_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.006762
res9 <- wilcox.test(tss_capHIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.4267
res10 <- wilcox.test(tss_Fantom5 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.8911
res11 <- wilcox.test(tss_DLPFC_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.2011
res12 <- wilcox.test(tss_HIPPO_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.2047
res13 <- wilcox.test(HIPPO_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.5885
res14 <- wilcox.test(DLPFC_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.1337
res15 <- wilcox.test(CAUDATE_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.2204
---------------------------------------------------------------------------------------------------------------------
lncRNA$Group<- factor(lncRNA$Group, level = c("SCZ elncRNA", "Brain elncRNA"))
res1 <- wilcox.test(ACTAC_HIPPO1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.07329
res2 <- wilcox.test(ACTAC_HIPPO2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.04243
res3 <- wilcox.test(ACTAC_HIPPO3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.04433
res4 <- wilcox.test(ACTAC_Caudate1 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01115
res5 <- wilcox.test(ACTAC_Caudate2 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01803
res6 <- wilcox.test(ACTAC_Caudate3 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01389
res7 <- wilcox.test(tss_CP_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 1.844e-05
res8 <- wilcox.test(tss_GZ_HIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 9.027e-06
res9 <- wilcox.test(tss_capHIC ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.09715
res10 <- wilcox.test(tss_Fantom5 ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.9333
res11 <- wilcox.test(tss_DLPFC_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.001265
res12 <- wilcox.test(tss_HIPPO_pcHIC_po ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.0009698
res13 <- wilcox.test(HIPPO_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.07833
res14 <- wilcox.test(DLPFC_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.005671
res15 <- wilcox.test(CAUDATE_DE_P ~ Group, data = lncRNA, paired = FALSE, alternative = "greater") #p-value = 0.01724

