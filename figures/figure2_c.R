library(stringr)
library(ChIPseeker)
rm(list=ls())
atac_seq_dir = "/NAS/luozh/CancerEnhancerDB/step6_peak_figure/hg38_atac_cancer_bed/"
atac_files = list.files(atac_seq_dir)
cancer_type = sapply(str_split(atac_files, "_"), "[[", 1)
file_path = paste(atac_seq_dir, atac_files, sep="/") 
df_atac_path = data.frame(cancer_type = cancer_type, file_path)

#enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step6_peak_figure/narrow_reproduced_peak"
enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step13_final_typical_enhancer/normalized_enhancer"
enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step13_final_typical_enhancer/cell_line_normalized_enhancer"
enhancer_files = list.files(enhancer_dir)
cancer_type = sapply(str_split(enhancer_files, "_"), "[[", 1)
file_path = paste(enhancer_dir, enhancer_files, sep="/")
df_enhancer_path = data.frame(cancer_type = cancer_type, file_path)

common_cancer = intersect(df_atac_path$cancer_type, df_enhancer_path$cancer_type)

df_result  = data.frame(enhancer_cancer=character(),
                        enhancer_number=integer(), 
                        atac_cancer=character(), 
                        atac_number=integer(),
                        overlapped_enhancer = integer(),
                        stringsAsFactors=FALSE) 

calculating_overlapes = function(one_cancer, df_result){
  enhancer_bed = df_enhancer_path[df_enhancer_path$cancer_type == one_cancer, "file_path"]
  enhancer_peak = readPeakFile(enhancer_bed)
  
  #for (atac_cancer in common_cancer){
  for (atac_cancer in df_atac_path$cancer_type){
    atac_bed = df_atac_path[df_atac_path$cancer_type == atac_cancer, "file_path"]
    atac_peak = readPeakFile(atac_bed)
    #overlap
    enhancer_num = length(enhancer_peak)
    atac_num = length(atac_peak)
    res = overlapsAny(enhancer_peak, atac_peak) # value: TRUE, FALSE
    enhancers_with_atac = length(res[res])
    percent = enhancers_with_atac/enhancer_num
    #percent = enhancers_with_atac/atac_num
    df_result[nrow(df_result) + 1,] = c(one_cancer, enhancer_num, atac_cancer, atac_num, enhancers_with_atac)
    #print(paste(c(one_cancer, "enhancer number:", enhancer_num, "; overlapped in", atac_cancer, "atac_num:", atac_num, "percent:", percent), sep=""))
  }
  return(df_result)
}

#result table
for (one_cancer in common_cancer){
  df_result = calculating_overlapes(one_cancer = one_cancer, df_result)
}

write.csv(df_result, "/NAS/luozh/CancerEnhancerDB/final_figures/peak_overlap_1010_all_1_primary.csv", row.names = FALSE, quote = F)
#run
#df_atac = read.csv("/NAS/luozh/CancerEnhancerDB/final_figures/peak_overlap_1009_all_1_primary.csv")
df_atac = read.csv("/NAS/luozh/CancerEnhancerDB/final_figures/peak_overlap_1010_all_1_primary.csv")
df_atac["percent_enhancer"] = df_atac$overlapped_enhancer / df_atac$enhancer_number
df_atac["percent_atac"] = df_atac$overlapped_enhancer / df_atac$atac_number
data = matrix(1:nrow(df_atac), nrow = length(unique(df_atac$enhancer_cancer)), ncol = length(unique(df_atac$atac_cancer)))
rownames(data) = unique(df_atac$enhancer_cancer)
colnames(data) = unique(df_atac$atac_cancer)
for (one_row in rownames(data)){
  for (one_col in colnames(data)){
    data[one_row,one_col] = df_atac[df_atac$enhancer_cancer == one_row & df_atac$atac_cancer == one_col, "percent_atac"]
  }
}

data = as.data.frame(data)
data_t = t(data)
data_t = data_t[rownames(data_t) %in% common_cancer,]
data_s = scale(data_t)
p = pheatmap(data_t, scale="column",cluster_cols = F, cluster_rows = F, legend_labels = "Overlapped ATAC peak Z-score")


library(ComplexHeatmap)
col_fun = colorRamp2(seq(-2,2,1), rev(brewer.pal(n = 5, name = "RdYlBu")))
col_fun = colorRamp2(seq(-2,2,1), c("#4575B4", "#91BFDB","#FFFFBF","#FC8D59","#D73027"))
rev(brewer.pal(n = 7, name = "RdYlBu"))
p = Heatmap(data_s, col=col_fun, name = "Overlapped\nZ-score", column_title = "Typical enhancer", row_title = "ATAC peak", cluster_rows = FALSE, cluster_columns = FALSE,
            row_title_side = c("right"), column_title_side = c("bottom"), column_title_gp = gpar(fontsize = 18), 
            row_title_gp = gpar(fontsize = 18), show_heatmap_legend = T, column_names_gp = gpar(fontsize = 12), column_names_rot = 90,
            colo = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                      "RdYlBu")))(100)
            )

pdf("/NAS/luozh/CancerEnhancerDB/final_figures/figure2S_a.pdf", width = 6, height=5)
draw(p, align_heatmap_legend = "heatmap_top", legend_title_gp=gpar(fontsize = 16))
dev.off()

#cell line















































###useless
draw(p, heatmap_legend_side = "bottom")
pushViewport(viewport(width = 1, height = 1))
#grid.rect()
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "Overlapped Z-score", at = c(-4, 0, 4), 
             labels = c("low", "median", "high"), direction = "horizontal")

#draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
draw(lgd,  x = unit(18, "cm"), y = unit(0.1, "cm"),just = c("right", "bottom"))
#popViewport()

#draw(p, heatmap_legend_side="")