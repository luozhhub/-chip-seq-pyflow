rm(list=ls())


#primary data
super_enhancer_primary_file = "/NAS/luozh/CancerEnhancerDB/step14_final_super_enhancer/super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv"
df_primary = read.csv(super_enhancer_primary_file, sep="\t")
col_n = colnames(df_primary)
col_n[1] = "PeakID"
colnames(df_primary) = col_n 
df_primary_sub = df_primary[,c("PeakID", "Gene.Name")]
str_split_fixed(df_primary_sub$PeakID, "_", 3)
df_primary_sub = df_primary_sub %>% separate(PeakID, c("GSE", "cancer_type", "peak_num"), "_")


#all sample summary
sample_info = df_primary_sub[,c("GSE", "cancer_type")] 
sample_info = sample_info  %>% distinct(GSE, cancer_type, .keep_all = TRUE)
sample_info = as.data.frame(table(sample_info$cancer_type))
colnames(sample_info) = c("cancer_type", "sample_number")

#gene_number
gene_list = c("IER3", "LIF", "SLC7A5", "CYP2S1", "PHF19", "RNF43", "CEBPB", "TBC1D16", "TNFRSF6B", "VEGFA")
total_gene_df = purrr::map(gene_list,function(x){
gene_df = df_primary_sub[df_primary_sub$Gene.Name == x, ]
sample_info_gene = gene_df[,c("GSE", "cancer_type")] 
sample_info_gene = sample_info_gene  %>% distinct(GSE, cancer_type, .keep_all = TRUE)
sample_info_gene = as.data.frame(table(sample_info_gene$cancer_type))
colnames(sample_info_gene) = c("cancer_type", "sample_number")
gene_percent = dplyr::left_join(sample_info, sample_info_gene, c("cancer_type" = "cancer_type"))
colnames(gene_percent) = c("cancer_type", "total_sample", "gene_sample")
gene_percent$gene = x
gene_percent[is.na(gene_percent)] <- 0
return(gene_percent)
}
)%>%dplyr::bind_rows()

#compare
total_gene_df$percent = total_gene_df$gene_sample / total_gene_df$total_sample

gene_percent_df = purrr::map(gene_list,function(x){
  gene_df =  total_gene_df[total_gene_df$gene == x,]
  sum_all = sum(gene_df$total_sample) - gene_df[gene_df$cancer_type == "COAD", "total_sample"]
  sum_cancer = sum(gene_df$gene_sample) - gene_df[gene_df$cancer_type == "COAD", "gene_sample"]
  all_percent = sum_cancer/sum_all
  data.frame(gene=x, all_percent=all_percent, CRC_percent=gene_df[gene_df$cancer_type == "COAD", "percent"])
  
})%>%dplyr::bind_rows()

write.csv( gene_percent_df, "/NAS/luozh/CancerEnhancerDB/final_figures/gene_percent_CRC.csv", row.names = F, quote = F)

#figure 3a
colnames(gene_percent_df) = c("Gene", "All cancer", "COAD")
data = melt(gene_percent_df, id=c("Gene"))
colnames(data) = c("Gene", "Group", "Percent")

data$Gene = factor(data$Gene, levels=unique(data$Gene))
p = ggplot(data, aes(x=Gene, y=Percent, fill=Group)) +
  #ggtitle("Up regulated genes") +
  #geom_text(aes(x=class, y=P_value, label=P_value, vjust = -1.5), position=position_dodge(.9), size=3) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=c('#999999','#E69F00')) +
  theme(axis.text.y   = element_text(size=12),
        axis.text.x   = element_text(size=12, angle = 90, vjust = 0.5, hjust=1),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        #panel.background = element_blank(),
        panel.grid.major.y = element_line(), 
        panel.grid.minor = element_line(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p
ggsave(filename = "/NAS/luozh/CancerEnhancerDB/final_figures/figure3_a.pdf", p, width = 7, height=4)

#figure 3c
data_3c = total_gene_df[total_gene_df$total_sample >= 10,]

p2 = ggplot(data = data_3c, aes(x=cancer_type, y=percent)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(vars(gene), nrow = 10, strip.position="right") +
  xlab("Cancer type") +
  ylab("Percent") + 
  theme(
    strip.background = element_blank(),
    #strip.text.x = element_blank()
  )
p2

ggsave(filename = "/NAS/luozh/CancerEnhancerDB/final_figures/figure3s_a.pdf", p2, width = 9, height=8)
