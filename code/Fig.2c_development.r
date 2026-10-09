module ava
module load r/4.0

library(ggplot2)
library(reshape2)
library(scales)

# 读取数据
data <- read.table("/data/home/yuanmengqin/AD_GWAS/SCZ_risk_lncRNA_prediction_2024_5_9/SREL_ANALYSIS/brain_development/average_per_sample/SREL_EL_LBL_WBL_averge_matrix.txt", header = TRUE, sep = "\t")

# 重塑数据
melted_data <- melt(data, id.vars = "Group")

# 添加新列以区分出生前和出生后
melted_data$border <- ifelse(melted_data$Group %in% c("Early prenatal", "Early mid-prenatal", "Late mid-prenatal", "Late prenatal"), "no_border", "border")

# 设置 Group 和 variable 列为因子，并按照输入顺序排列
melted_data$Group <- factor(melted_data$Group, levels = unique(data$Group))
melted_data$variable <- factor(melted_data$variable, levels = colnames(data)[-1])

# 绘制热图，使用圆圈表示；fill上限固定7，>=7全部显示为7的颜色
p <- ggplot(melted_data, aes(x = variable, y = Group, fill = value)) +
  geom_point(aes(shape = border, color = border), size = 10) +
  scale_fill_gradient(low = "white", high = "red",
                      limits = c(0,7),
                      oob = squish) +
  scale_shape_manual(values = c("no_border" = 21, "border" = 21)) +
  scale_color_manual(values = c("no_border" = NA, "border" = "black")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(x = "Median Type", y = "Development Stage", fill = "Value", shape = "Border", color = "Border")

# 保存热图
ggsave("SREL_EL_LBL_WBL_average_heatmap.pdf", plot = p, width = 6, height = 6)