rm(list=ls())

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

yAxis = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
n = 0
for (gene in gene_list){
#gene="RNF43"
coad_sample = total_gene_df[(total_gene_df$cancer_type == "COAD") & (total_gene_df$gene == gene),]
coad_total_sample = coad_sample$total_sample
coad_gene_sample = coad_sample$gene_sample
cancer_sample = total_gene_df[total_gene_df$gene == gene,]
cancer_total_sample = sum(cancer_sample$total_sample)
cancer_gene_sample = sum(cancer_sample$gene_sample)

genes_fisher_cluster1 <- matrix(c(coad_gene_sample, coad_total_sample - coad_gene_sample, cancer_gene_sample - coad_gene_sample, cancer_total_sample - coad_total_sample -(cancer_gene_sample - coad_gene_sample)) , nrow = 2,
                                dimnames =list(c("COAD_specific", "Non_COAD_specific"),
                                               c("cancer", "Not cancer")))
OR_1 = fisher.test(genes_fisher_cluster1, conf.level = 0.95)
n = n +1
yAxis = c(yAxis, n)
boxOdds = c(boxOdds, OR_1$estimate)
boxCILow = c(boxCILow, OR_1$conf.int[1])
boxCIHigh = c(boxCIHigh, OR_1$conf.int[2])

}




#boxLabels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")
df_primary_od <- data.frame(yAxis = yAxis, 
                 BoxOdds = boxOdds,
                 BoxCILow = boxCILow,
                 BoxCIHigh = boxCIHigh,
                 gene = gene_list
)


df_primary_od$group = "Primary tissue"













#####################################cell line############################################################################
super_enhancer_primary_file = "/NAS/luozh/CancerEnhancerDB/step14_final_super_enhancer/cell_line_super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv"
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

yAxis = c()
boxOdds = c()
boxCILow = c()
boxCIHigh = c()
n = 0
for (gene in gene_list){
  #gene="RNF43"
  coad_sample = total_gene_df[(total_gene_df$cancer_type == "COAD") & (total_gene_df$gene == gene),]
  coad_total_sample = coad_sample$total_sample
  coad_gene_sample = coad_sample$gene_sample
  cancer_sample = total_gene_df[total_gene_df$gene == gene,]
  cancer_total_sample = sum(cancer_sample$total_sample)
  cancer_gene_sample = sum(cancer_sample$gene_sample)
  
  genes_fisher_cluster1 <- matrix(c(coad_gene_sample, coad_total_sample - coad_gene_sample, cancer_gene_sample - coad_gene_sample, cancer_total_sample - coad_total_sample -(cancer_gene_sample - coad_gene_sample)) , nrow = 2,
                                  dimnames =list(c("COAD_specific", "Non_COAD_specific"),
                                                 c("cancer", "Not cancer")))
  OR_1 = fisher.test(genes_fisher_cluster1, conf.level = 0.95)
  n = n +1
  yAxis = c(yAxis, n)
  boxOdds = c(boxOdds, OR_1$estimate)
  boxCILow = c(boxCILow, OR_1$conf.int[1])
  boxCIHigh = c(boxCIHigh, OR_1$conf.int[2])
  
}




#boxLabels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")
df_cell_line_od <- data.frame(yAxis = yAxis, 
                 BoxOdds = boxOdds,
                 BoxCILow = boxCILow,
                 BoxCIHigh = boxCIHigh,
                 gene = gene_list
)

df_cell_line_od$group = "Cell line"

df = rbind(df_primary_od, df_cell_line_od)
df$group=factor(df$group)
# Plot

p <- ggplot(df, aes(x = BoxOdds, y = gene,group=group,color = group)) + 
  # geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = BoxCIHigh, xmin = BoxCILow), size = 1.5, height = 
                   .2,position=position_dodge(.9),color = "gray50") +
  geom_point(size = 3.5,position=position_dodge(.9)) +
  #scale_x_continuous(limits = c(0, 5)) +
  scale_color_manual(values=c('#999999','#E69F00'))+
  xlab("Odds ratio") +
  ylab("Genes") + 
  #coord_trans(x = "log2") +
  theme(axis.title = element_text(size=14),
        #axis.text = element_text(size=24, color="black"),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=12), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.85,0.9)
  )
ggsave(filename = "/NAS/luozh/CancerEnhancerDB/final_figures/figure3_d.pdf", p, width = 7, height=5)
