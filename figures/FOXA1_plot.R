#cancer in TF cancer: "BRCA", "LIHC", "LUAD", "PRAD", "THCA", "UCS"


data_file = "/NAS/luozh/CancerEnhancerDB/step12_TF_summary/total_TF/FOXA1_summary.csv"
data = read.csv(data_file)

rank(data$rank_y)
data[order(data$rank_y),]
data[order(data$rank_x),]


data_nona = na.omit(data)
nrow(data_nona)
data_nona[order(data_nona$average_rank),]
rownames(data_nona) = data_nona$cancer_type
data_sub  = data_nona[,c("rank_x", "rank_y", "average_rank")]
colnames(data_sub) = c("Primary_tissue", "Cell_line", "Average")

dat = data_sub[,c("Primary_tissue", "Cell_line")]
p = pheatmap(as.matrix(dat), scale="column",angle_col = c("45"), cluster_col = FALSE, color = colorRampPalette(c("red", "white"))(300))
p = pheatmap(as.matrix(dat),scale="column",angle_col = c("45"), cluster_col = FALSE, color = colorRampPalette(c("red", "white"))(300))

dat$Primary_tissue = rank(-dat$Primary_tissue)
dat$Cell_line = rank(-dat$Cell_line)
p = pheatmap(as.matrix(dat),angle_col = c("45"), cluster_col = FALSE, treeheight_row = 0)

write.csv(data_sub)