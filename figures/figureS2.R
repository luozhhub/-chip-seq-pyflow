rm(list=ls())
atac_seq_dir = "/NAS/luozh/CancerEnhancerDB/step6_peak_figure/hg38_atac_cancer_bed/"
atac_files = list.files(atac_seq_dir)
cancer_type = sapply(str_split(atac_files, "_"), "[[", 1)
file_path = paste(atac_seq_dir, atac_files, sep="/") 
df_atac_path = data.frame(cancer_type = cancer_type, file_path)

enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step13_final_typical_enhancer/normalized_enhancer"
enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step13_final_typical_enhancer/cell_line_normalized_enhancer"
enhancer_files = list.files(enhancer_dir)
cancer_type = sapply(str_split(enhancer_files, "_"), "[[", 1)
file_path = paste(enhancer_dir, enhancer_files, sep="/")
df_enhancer_path = data.frame(cancer_type = cancer_type, file_path)

common_cancer = intersect(df_atac_path$cancer_type, df_enhancer_path$cancer_type)
cancers = common_cancer
cancer = c()
percent_of_atac_overlap = c()
sample_of_peaks = c()
for (one_cancer in cancers){
  #one_cancer = "COAD"
  enhancer_bed = df_enhancer_path[df_enhancer_path$cancer_type == one_cancer, "file_path"]
  enhancer_peak = readPeakFile(enhancer_bed)
  
  atac_bed = df_atac_path[df_atac_path$cancer_type == one_cancer, "file_path"]
  atac_peak = readPeakFile(atac_bed)
  
  if (enhancer_peak$sample_coverage[1] >= 10){
      cancer = c(cancer, rep(one_cancer, 10))
      sample_of_peaks = c(sample_of_peaks, 1:10)
  
      for (i in 1:10){
        reduce_number = i
        enhancer_peak = enhancer_peak[enhancer_peak$sample_coverage >= reduce_number,]
        #overlap
        res = overlapsAny(enhancer_peak, atac_peak)
        
        enhancer_num = length(enhancer_peak)
        enhancers_with_atac = length(res[res])
        percent = enhancers_with_atac/enhancer_num
        percent_of_atac_overlap = c(percent_of_atac_overlap, percent)
      }
  }
}
data = data.frame(cancertype=cancer, sample_of_peaks, percent_of_atac_overlap)

library(ggplot2)
library(ggsci)
p = ggplot(data, aes(x=sample_of_peaks, y=percent_of_atac_overlap, group=cancertype, color=cancertype)) +
  geom_line() +
  scale_color_uchicago() +
  xlab("Typical enhancers covering sample number") +
  ylab("Overlapped percentage with ATAC-seq peaks") + 
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title = element_text(size=14),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
p
ggsave(filename = "/NAS/luozh/CancerEnhancerDB/final_figures/figure2s_b.pdf", width = 7, height=6)
