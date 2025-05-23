# 热图绘制

# 加载包
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

setwd("H:\\生信文章\\复现\\DEGanalysis\\热图和火山图")
# 加载表达数据和注释数据
DEG <- read.table("exp.txt",header=T, sep="\t", check.names=F,row.names = 1)
annotation_col <- read.table("annotation_col.txt",header=T, sep="\t", check.names=F,row.names = 1)
DEG <- t(DEG)
DEG<- DEG[match(rownames(annotation_col), rownames(DEG)), ]#让DEG按annotation_col行名排列
DEG <- t(DEG)
#转换为因子型
annotation_col[] <- lapply(annotation_col, as.factor)


# 创建注释颜色
annotation_colors <- list(
  risk=c("high"="#F05006","low"="#FCA00C"),
  Gender = c("Male"="#0AA0BF","Female"="#F36E98"),
  Fustat = c("Alive"="#78B177","Dead"="#F6114A"),
  Group = c("Metastases"="#CF0BF1","No_metastases"="#C3C4E9")
)

# 绘制热图
pdf("model_geneheatmap3.pdf",width = 8.5,height =3.9)
p <- pheatmap(DEG,
              scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(50),
         show_rownames = TRUE,
         show_colnames =FALSE,
         treeheight_row =20,  # 隐藏行的树状图连线
         treeheight_col = 0 )

dev.off()
