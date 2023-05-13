rm(list=ls())
typical_enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step13_final_typical_enhancer/normalized_enhancer"
peak_files = list.files(typical_enhancer_dir)
cancer_type = sapply(str_split(peak_files, "_"), "[[", 1)
#cancer_type = gsub(".csv", "", cancer_type)
file_path = paste(typical_enhancer_dir, peak_files, sep="/")
df_typical_path = data.frame(cancer_type = cancer_type, file_path)



super_enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step13_final_typical_enhancer/cell_line_normalized_enhancer"
super_files = list.files(super_enhancer_dir)
cancer_type = sapply(str_split(super_files, "_"), "[[", 1)
#cancer_type = gsub(".csv", "", cancer_type)
file_path = paste(super_enhancer_dir, super_files, sep="/")
df_super_path = data.frame(cancer_type = cancer_type, file_path)

cancer_type = c()
peak_number = c()
enhancer_type = c()

cancer_common = intersect(df_typical_path$cancer_type, df_super_path$cancer_type)

for (one_cancer in cancer_common){
  enhancer_bed = df_typical_path[df_typical_path$cancer_type == one_cancer, "file_path"]
  enhancer_peak = read.csv(enhancer_bed)
  peak_number = c(peak_number, nrow(enhancer_peak))
  cancer_type = c(cancer_type, one_cancer)
  enhancer_type = c(enhancer_type, "Primary tissue")
}

for (one_cancer in cancer_common){
  enhancer_bed = df_super_path[df_super_path$cancer_type == one_cancer, "file_path"]
  enhancer_peak = read.csv(enhancer_bed)
  peak_number = c(peak_number, nrow(enhancer_peak))
  cancer_type = c(cancer_type, one_cancer)
  enhancer_type = c(enhancer_type, "Cell line")
}

data = data.frame(enhancer_number=peak_number, cancer_type=cancer_type, data_type=enhancer_type)
data_sub = data[data$data_type == "Primary tissue",]
data_sub = data_sub[order(data_sub$enhancer_number,decreasing = T),]

rank = data_sub$cancer_type

library(ggsci)
data$cancer_type = factor(data$cancer_type, levels=rank)
data$data_type = factor(data$data_type, levels = c("Primary tissue", "Cell line"))
p = ggplot(data, aes(x=cancer_type, y = enhancer_number, fill=data_type))+ geom_histogram(stat='identity', position = position_dodge()) +
  scale_fill_locuszoom() +
  #xlab("Cancer types") +
  ylab("Typical enhancer number") + 
  theme_classic() +
  #theme(legend.position = "right") +
  theme(legend.position = c(0.9, 0.8)) +
  theme(axis.title = element_text(size=18),
        axis.title.x = element_blank(),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=14), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
p
ggsave("/NAS/luozh/CancerEnhancerDB/final_figures/figure2b.pdf", p, width=13.5, height=3.5)




###super
typical_enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step8_super_enhancer/super_enhancer_peaks_merged"
peak_files = list.files(typical_enhancer_dir)
cancer_type = sapply(str_split(peak_files, "_"), "[[", 1)
#cancer_type = gsub(".csv", "", cancer_type)
file_path = paste(typical_enhancer_dir, peak_files, sep="/")
df_typical_path = data.frame(cancer_type = cancer_type, file_path)



super_enhancer_dir = "/NAS/luozh/CancerEnhancerDB/step8_super_enhancer/cell_line_super_enhancer_peaks_merged"
super_files = list.files(super_enhancer_dir)
cancer_type = sapply(str_split(super_files, "_"), "[[", 1)
#cancer_type = gsub(".csv", "", cancer_type)
file_path = paste(super_enhancer_dir, super_files, sep="/")
df_super_path = data.frame(cancer_type = cancer_type, file_path)

cancer_type = c()
peak_number = c()
enhancer_type = c()

#cancer_common = intersect(df_typical_path$cancer_type, df_super_path$cancer_type)
for (one_cancer in cancer_common){
  enhancer_bed = df_typical_path[df_typical_path$cancer_type == one_cancer, "file_path"]
  enhancer_peak = read.csv(enhancer_bed)
  peak_number = c(peak_number, nrow(enhancer_peak))
  cancer_type = c(cancer_type, one_cancer)
  enhancer_type = c(enhancer_type, "Primary tissue")
}

for (one_cancer in cancer_common){
  enhancer_bed = df_super_path[df_super_path$cancer_type == one_cancer, "file_path"]
  enhancer_peak = read.csv(enhancer_bed)
  peak_number = c(peak_number, nrow(enhancer_peak))
  cancer_type = c(cancer_type, one_cancer)
  enhancer_type = c(enhancer_type, "Cell line")
}

data = data.frame(enhancer_number=peak_number, cancer_type=cancer_type, data_type=enhancer_type)
data_sub = data[data$data_type == "Primary tissue",]
data_sub = data_sub[order(data_sub$enhancer_number,decreasing = T),]

rank1 = data_sub$cancer_type

data$cancer_type = factor(data$cancer_type, levels=rank)
data$data_type = factor(data$data_type, levels = c("Primary tissue", "Cell line"))
p1 = ggplot(data, aes(x=cancer_type, y = enhancer_number, fill=data_type))+ geom_histogram(stat='identity', position = position_dodge()) +
  scale_fill_locuszoom() +
  xlab("Cancer types") +
  ylab("Super enhancer number") + 
  theme_classic() +
  #theme(legend.position = "right") +
  theme(legend.position='none') +
  theme(axis.title = element_text(size=18),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=14), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
p1



library(cowplot)


all_p = plot_grid(p, p1, ncol = 1, align = 'v')

ggsave("/NAS/luozh/CancerEnhancerDB/final_figures/figure2_b.pdf", all_p, width = 14, height = 7)
pdf("/NAS/luozh/CancerEnhancerDB/final_figures/figure2_b.pdf", all_p, width = 14, height = 7)
all_p
dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#super enhancer
file_super_enhancer_primary = "/NAS/luozh/CancerEnhancerDB/step14_final_super_enhancer/super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv"
file_super_enhancer_cell_line = "/NAS/luozh/CancerEnhancerDB/step14_final_super_enhancer/cell_line_super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv"

df_super_enhancer_primary = read.csv(file_super_enhancer_primary, sep="\t")
coln = colnames(df_super_enhancer_primary)
coln[1] = "PeakID"
colnames(df_super_enhancer_primary) = coln
cancer_type = sapply(str_split(df_super_enhancer_primary$PeakID, "_"), "[[", 2)
prim_super_df = as.data.frame(table(cancer_type))
colnames(prim_super_df) = c("Cancer_types", "number")

df_super_enhancer_cell_line = read.csv(file_super_enhancer_cell_line, sep="\t")
coln = colnames(df_super_enhancer_cell_line)
coln[1] = "PeakID"
colnames(df_super_enhancer_cell_line) = coln
cancer_type = sapply(str_split(df_super_enhancer_cell_line$PeakID, "_"), "[[", 2)
cell_super_df = as.data.frame(table(cancer_type))
colnames(cell_super_df) = c("Cancer_types", "number")

df_merge = dplyr:: inner_join(prim_super_df, cell_super_df, c("Cancer_types" = "Cancer_types"))
df_sub = df_merge[df_merge$Cancer_types %in% cancer_common,]
colnames(df_sub) = c("Cancer_types", "Primary tissue enhancer", "Cell line enhancer")

data_sup <- melt(df_sub, id=c("Cancer_types"))
colnames(data_sup) = c("cancer_type", "data_type", "enhancer_number")
data_sup$cancer_type = factor(data_sup$cancer_type, levels=rank)
data_sup$data_type = factor(data_sup$data_type, levels = c("Primary tissue enhancer", "Cell line enhancer"))


p2 = ggplot(data_sup, aes(x=cancer_type, y = enhancer_number, fill=data_type))+ geom_histogram(stat='identity', position = position_dodge()) +
  scale_fill_locuszoom() +
  xlab("Cancer types") +
  ylab("Super enhancer number") + 
  theme_classic() +
  theme(legend.position = "right") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.title = element_text(size=18),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )
p2

ggsave("/NAS/luozh/CancerEnhancerDB/final_figures/figure2b_super.pdf", p2, width=13.5, height=3.5)
