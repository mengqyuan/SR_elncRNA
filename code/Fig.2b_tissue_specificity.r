data=read.table("SREL_WBL_tissue_specificity_wilcox.txt",sep="\t",header=T)
data[,3]="p< 0.05"
data[intersect(which(data[,2]<0.05),grep("Brain",data[,1])),3]="p < 0.05 in brain"
data[which(data[,2]>0.05),3]="p> 0.05"
data[,2]=-log10(data[,2])
library("ggplot2")
library(ggpubr)
p1=ggdotchart( data, x = "Tissue", y = "P_value",
           color = "V3",                                # Color by groups
		   dot.size = 4,
           palette = c("#4DAD49", "#377DB7", "#E21A1C"), # Custom color palette
           sorting = "ascending",                        # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           ggtheme = theme_pubr()                        # ggplot2 theme
)+coord_flip()

pdf("SREL_WBL_tissue_specificity_bang.pdf",width=6,height=10)
print(p1)
dev.off()