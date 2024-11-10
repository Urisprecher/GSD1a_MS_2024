#libraries
library(tidyverse)
library(psych)
library(dplyr)
library(boot)
library(GGally)
library(FactoMineR)
library(mixOmics)
library(readxl)
library("writexl")
library(stats)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(ggfortify)
library(ggplot2)
library(ggVennDiagram)
library(tidyr)
library(scales)   
library(patchwork)
###################
BiocManager::install("", force = TRUE)
#working dir
setwd("data_integration/INTEG_2024")
#load data
##met dat
methylation <- read_excel("met_integ_3.xlsx")
## rep methyl gene symbols by Lowest p val to remove duplicate gene symbols. 
methylation <- methylation %>% 
  group_by(SYMBOL) %>% 
  slice(which.min(adj.P.Val)) %>% 
  ungroup()
drop <- c("adj.P.Val")
methylation = methylation[,!(names(methylation) %in% drop)]
## atac all 
atac_seq <- read_excel("atac_integ_1.xlsx")
## atac promoters 
atac_seq <- read_excel("atac_integ_2.xlsx")
## rna dat 
rna_seq <- read_excel("rna_integ_1.xlsx")

#normalize each dataset
normalize_data <- function(df) {
  #numerical columns only following column 1 which is sample column
  df[, 2:ncol(df)] <- scale(df[, 2:ncol(df)])
  return(df)
}
#norm frames
methylation <- normalize_data(methylation)
atac_seq <- normalize_data(atac_seq)
rna_seq <- normalize_data(rna_seq)
## remove var - top genes
#filter genes based on variance top genes
filter_top_variance <- function(df, top_n = 10000) {
  feature_variances <- apply(df[, -1], 1, var)
  top_indices <- order(feature_variances, decreasing = TRUE)[1:min(top_n, length(feature_variances))]
  #filt and return
  filtered_df <- df[top_indices, ]
  
  return(filtered_df)
}

#filter each dataframe -
methylation <- filter_top_variance(methylation)
atac_seq <- filter_top_variance(atac_seq)
rna_seq <- filter_top_variance(rna_seq)
#common gene symbols
common_genes <- Reduce(intersect, list(methylation$SYMBOL, atac_seq$SYMBOL, rna_seq$SYMBOL))
#subset all datasets to common gene symbols
sub_methylation <- methylation %>% filter(SYMBOL %in% common_genes)
sub_atac_seq <- atac_seq %>% filter(SYMBOL %in% common_genes)
sub_atac_seq <- sub_atac_seq %>% 
  group_by(SYMBOL) %>% 
  summarise(across(everything(), mean),
                .groups = 'drop') %>%
  as.data.frame()

sub_rna_seq <- rna_seq %>% filter(SYMBOL %in% common_genes)
##################
common_genes <- Reduce(intersect, list(sub_methylation$SYMBOL, sub_atac_seq$SYMBOL, sub_rna_seq$SYMBOL))
write.csv(as.data.frame(common_genes), "common_genes-backround.csv", row.names = FALSE)
## transpose frames for cor
sub_methylation_t <- as.data.frame(t(sub_methylation))
names(sub_methylation_t) <- lapply(sub_methylation_t[1, ], as.character)
sub_methylation_t <- sub_methylation_t[-1,] 
sub_methylation_t <- sub_methylation_t %>% rownames_to_column("sample")
sub_atac_seq_t <- as.data.frame(t(sub_atac_seq))
names(sub_atac_seq_t) <- lapply(sub_atac_seq_t[1, ], as.character)
sub_atac_seq_t <- sub_atac_seq_t[-1,] 
sub_atac_seq_t <- cbind("sample" = row.names(sub_atac_seq_t), `row.names<-`(sub_atac_seq_t, NULL))
sub_rna_seq_t <- as.data.frame(t(sub_rna_seq))
names(sub_rna_seq_t) <- lapply(sub_rna_seq_t[1, ], as.character)
sub_rna_seq_t <- sub_rna_seq_t[-1,] 
sub_rna_seq_t <- cbind("sample" = row.names(sub_rna_seq_t), `row.names<-`(sub_rna_seq_t, NULL))
## add sufix -met
suffix <- "met"
colnames(sub_methylation_t)[-1] <- paste0(colnames(sub_methylation_t)[-1], "_", suffix)
## add sufix -atac
suffix <- "atac"
colnames(sub_atac_seq_t)[-1] <- paste0(colnames(sub_atac_seq_t)[-1], "_", suffix)
## add sufix -rna
suffix <- "rna"
colnames(sub_rna_seq_t)[-1] <- paste0(colnames(sub_rna_seq_t)[-1], "_", suffix)
#merge dataframes
## merging options - with suf
combined_df_1_3 <- merge(sub_methylation_t, sub_atac_seq_t, by = "sample")
combined_df_1_2 <- merge(sub_rna_seq_t, sub_atac_seq_t, by = "sample")
combined_df_2_3 <- merge(sub_methylation_t, sub_rna_seq_t, by = "sample")
combined_df <- merge(combined_df_1_3, sub_rna_seq_t, by = "sample")
#combined data to CSV
write.csv(combined_df, file = "combined_data_integ_reduced_cor_vs.csv", row.names = FALSE)
####
gc()
#cor segemnt
#prep 
cor_frame <- combined_df
rownames(cor_frame) <- cor_frame$sample
cor_frame = subset(cor_frame, select = -c(sample))
cor_frame <- data.frame(sapply(cor_frame, function(x) as.numeric(as.character(x))))
sapply(cor_frame, class)

#cor- boot 
#compute correlation for bootstrapping
boot_correlation <- function(data, indices) {
  resampled_data <- data[indices, ]
  if (nrow(resampled_data) < 2) {
    return(NA)  
  }
  
  cor(resampled_data[, 1], resampled_data[, 2], method = "pearson")
}

#bootstrapping for correlation
bootstrap_correlations <- function(df, n_bootstrap = 1000) {
  # base gene names
  gene_base_names <- unique(gsub("(_atac|_rna|_met)$", "", colnames(df)))
  
  results_list <- list()
  ##select correct suffix
  for (gene in gene_base_names) {
    gene_1 <- paste0(gene, "_atac")
    gene_2 <- paste0(gene, "_met")
    
    if (gene_1 %in% colnames(df) && gene_2 %in% colnames(df)) {
      #data for gene from datasets
      gene_1_data <- df[, gene_1]
      gene_2_data <- df[, gene_2]
      
      data_for_bootstrap <- data.frame(gene_1_data, gene_2_data)
      
      # bootstrapping
      boot_out <- boot(data_for_bootstrap, statistic = boot_correlation, R = n_bootstrap)
      
      # confidence intervals
      ci <- tryCatch(boot.ci(boot_out, type = "perc", conf = 0.9), error = function(e) NULL)
      
      if (!is.null(ci)) {
        results_list[[gene]] <- list(
          mean_correlation = mean(boot_out$t, na.rm = TRUE),
          lower_ci = ci$percent[4],
          upper_ci = ci$percent[5]
        )
      } else {
        results_list[[gene]] <- list(
          mean_correlation = NA,
          lower_ci = NA,
          upper_ci = NA
        )
      }
    }
  }
  
  results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  rownames(results_df) <- names(results_list)
  
  return(results_df)
}


#run per group
bootstrap_results_1_2 <- bootstrap_correlations(cor_frame, n_bootstrap = 1000)
bootstrap_results_1_3 <- bootstrap_correlations(cor_frame, n_bootstrap = 1000)
bootstrap_results_2_3 <- bootstrap_correlations(cor_frame, n_bootstrap = 1000)
###
write.csv(bootstrap_results_1_2, "bootstrap_results_1_2.csv", row.names = TRUE)
write.csv(bootstrap_results_1_3, "bootstrap_results_1_3.csv", row.names = TRUE)
write.csv(bootstrap_results_2_3, "bootstrap_results_2_3.csv", row.names = TRUE)
## prep for intersection
bootstrap_results_1_2$Gene <- rownames(bootstrap_results_1_2)
bootstrap_results_1_3$Gene <- rownames(bootstrap_results_1_3)
bootstrap_results_2_3$Gene <- rownames(bootstrap_results_2_3)
colnames(bootstrap_results_1_2)
colnames(bootstrap_results_1_3)
colnames(bootstrap_results_2_3)
## if want abs mean cor
bootstrap_results_1_2$abs_cor <- abs(bootstrap_results_1_2$mean_correlation)
bootstrap_results_1_3$abs_cor <- abs(bootstrap_results_1_3$mean_correlation)
bootstrap_results_2_3$abs_cor <- abs(bootstrap_results_2_3$mean_correlation)
##filter top cors 
#thresholds
positive_threshold <- 0.7 
negative_threshold <- -0.7  
ci_diff_threshold <- 1   

filtered_boot_df_1_2 <- bootstrap_results_1_2 %>%
  mutate(CI_Diff = upper_ci - lower_ci)  %>%   
  filter(
    (mean_correlation > positive_threshold),  
    CI_Diff < ci_diff_threshold)


filtered_boot_df_1_3 <- bootstrap_results_1_3 %>%
  mutate(CI_Diff = upper_ci - lower_ci)  %>%   
  filter(
    (mean_correlation < negative_threshold),  
    CI_Diff < ci_diff_threshold)

filtered_boot_df_2_3 <- bootstrap_results_2_3 %>%
  mutate(CI_Diff = upper_ci - lower_ci)  %>%   
  filter(
    (mean_correlation < negative_threshold),  
    CI_Diff < ci_diff_threshold)

####
write.csv(filtered_boot_df_1_2, "final-bootstrap_results_1_2.csv", row.names = TRUE)
write.csv(filtered_boot_df_1_3, "final-bootstrap_results_1_3.csv", row.names = TRUE)
write.csv(filtered_boot_df_2_3, "final-bootstrap_results_2_3.csv", row.names = TRUE)
#gene symbols
common_genes <- Reduce(intersect, list(filtered_boot_df_1_2$Gene, filtered_boot_df_2_3$Gene))
write.csv(as.data.frame(common_genes), "common_genes-for string_analysis.csv", row.names = TRUE)



### splsda
pls_frame <- combined_df
combined_df$group <- c("control", "disease", "control", "disease")
####
pls_frame <- data.frame(sapply(pls_frame, function(x) as.numeric(as.character(x))))
sapply(pls_frame, class)

#feature matrix-X and group vector-Y
X <- pls_frame[, -c(1, ncol(pls_frame))]  
Y <- combined_df$group
#run sPLS-DA
splsda_result <- splsda(X, Y, ncomp = 2, keepX = c(5000, 5000))
# sPLS-DA results
plotIndiv(splsda_result, comp = 1:2, group = Y, ind.names = TRUE, ellipse = TRUE)
plotVar(splsda_result, comp = 1:2)
# top contributing variables 
top_variables <- unique(c(
  selectVar(splsda_result, comp = 1)$name
))
### important! 
top_var_df <- as.data.frame(top_variables)
# Save 
write.csv(top_var_df, file = "top_variables-comp1-(splsada-reg-5000).csv", row.names = FALSE)

top_variables_2 <- unique(c(
  selectVar(splsda_result, comp = 2)$name
))
### important! 
top_var_df <- as.data.frame(top_variables_2)
# Save 
write.csv(top_var_df, file = "top_variables-comp2-(splsada-reg-5000).csv", row.names = FALSE)

combined_top_genes <- unique(c(
  selectVar(splsda_result, comp = 2)$name,
  selectVar(splsda_result, comp = 1)$name
  
))
### important! 
top_var_df <- as.data.frame(combined_top_genes)
# Save 
write.csv(top_var_df, file = "top_variables-comp1+2-(splsada-reg-5000).csv", row.names = FALSE)
## vis for pla-da /sparase analysis
## plot dist
#suffix 
top_var_df$Suffix <- sub(".*_(atac|rna|met)$", "\\1", top_var_df$combined_top_genes)
top_var_df$Suffix2 <- sub(".*_(atac).*", "\\1", top_var_df$Suffix)
top_var_df$Suffix3 <- sub(".*_(rna).*", "\\1", top_var_df$Suffix2)
write.csv(top_var_df, file = "combined_top_genes-arranged-(splsada-reg-5000).csv", row.names = FALSE)
#number of genes in each category 
gene_distribution <- table(top_var_df$Suffix3)
gene_distribution_df <- as.data.frame(gene_distribution)
#pie chart with percentage labels
ggplot(gene_distribution_df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  labs(title = "Distribution of Top 500 Genes by Dataset",
       fill = "Dataset",
       x = NULL, y = NULL) +
  theme_void() +
  geom_text(aes(label = paste0(round(Freq/sum(Freq)*100, 1), "%")),
            position = position_stack(vjust = 0.5))

###duplicated
# gene names without suffixes
top_var_df$GeneName <- sub("_(rna|atac|met)$", "", top_var_df$combined_top_genes)
#gene name count across datasets
gene_counts <- table(top_var_df$GeneName)
#presence in two datasets and all three datasets
top_var_df$InTwoDatasets <- ifelse(gene_counts[top_var_df$GeneName] == 2, "Yes", "No")
top_var_df$InThreeDatasets <- ifelse(gene_counts[top_var_df$GeneName] == 3, "Yes", "No")
##
write.csv(top_var_df, "combined_top_genes-analyzed_(splsada-reg-5000).csv", row.names = FALSE)


#top important features from sPLS-DA model
important_genes <- rownames(splsda_result$loadings$X)[which(abs(splsda_result$loadings$X[,1]) > 0.05)]  
# subset important genes
important_genes_data <- pls_pca_run[, important_genes]
important_genes_data_mat <- as.matrix(important_genes_data)
#heatmap plot
Heatmap(important_genes_data_mat, name = "Expression", cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE, col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
?Heatmap

###
## chr dist pie charts
gc()
## read data 
methylation_sig <- read_excel("met_sig_0.1.xlsx")
atac_seq_sig <- read_excel("atac_sig_0.1.xlsx")
rna_seq_sig <- read_excel("rna_sig_0.1.xlsx")
methylation_sig$chr <- as.character(methylation_sig$chr)
atac_seq_sig$chr <- as.character(atac_seq_sig$chr)
rna_seq_sig$chr <- as.character(rna_seq_sig$chr)
methylation_sig$chr = sort(methylation_sig$chr, decreasing = TRUE)
rna_seq_sig$chr = sort(rna_seq_sig$chr, decreasing = TRUE)
atac_seq_sig$chr = sort(atac_seq_sig$chr, decreasing = TRUE)

# chromosome distribution data preparation
prepare_chromosome_distribution <- function(df, suffix) {
  df %>%

    group_by(chr) %>%
    summarise(count = n()) %>%
    mutate(dataset = suffix) 
}
methylation_chrom <- prepare_chromosome_distribution(methylation_sig, "Methylation")
methylation_chrom$percentages <- round(methylation_chrom$count / sum(methylation_chrom$count) * 100, 1)
atac_seq_chrom <- prepare_chromosome_distribution(atac_seq_sig, "ATAC-seq")
atac_seq_chrom$percentages <- round(atac_seq_chrom$count / sum(atac_seq_chrom$count) * 100, 1)
rna_expression_chrom <- prepare_chromosome_distribution(rna_seq_sig, "RNA Expression")
rna_expression_chrom$percentages <- round(rna_expression_chrom$count / sum(rna_expression_chrom$count) * 100, 1)
nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
#pie charts for datasets
create_pie_chart <- function(data, dataset_name) {
  ggplot(data, aes(x = "", y = count, fill = chr)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    ggtitle(paste("Chromosome Distribution of Significant Genes -", dataset_name)) +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_fill_manual(values = mycolors) +
    geom_text(aes(label = paste0(percentages, "%")),
              
              
              position = position_stack(vjust = 0.5), size=2)
}

# run
pie_chart_methylation <- create_pie_chart(methylation_chrom, "Methylation")
pie_chart_atac_seq <- create_pie_chart(atac_seq_chrom, "ATAC-seq")
pie_chart_rna_expression <- create_pie_chart(rna_expression_chrom, "RNA Expression")
#pie charts vis
print(pie_chart_methylation)
print(pie_chart_atac_seq)
print(pie_chart_rna_expression)
### save frames
write.csv(methylation_chrom, "sig_o.1-methylation_chrom_dist.csv", row.names = FALSE)
write.csv(atac_seq_chrom, "sig_o.1-atac_seq_chrom_dist.csv", row.names = FALSE)
write.csv(rna_expression_chrom, "sig_o.1-rna_expression_chrom_dist.csv", row.names = FALSE)


## vens + perms
## load dat
atac = read.csv('D:/MiguelW12/Documents/data_integration/res_deseq_hcvsgsd_peaks_DF_withanno.csv', header=TRUE)
met = read.csv('D:/MiguelW12/Documents/data_integration/DMP_proc.csv', header=TRUE)
rna = read.csv('D:/MiguelW12/Documents/data_integration/male_ANO-star-40filt.csv', header=TRUE)
##atac prep
pgdata_atac<-atac %>%
  group_by(SYMBOL) 
nrow(pgdata_atac)
pgdata_atac <- na.omit(pgdata_atac)
pgdata_atac_extream<-pgdata_atac[!duplicated(pgdata_atac$SYMBOL),]
#met  prep
metil_data_melt_na<-met[met$SYMBOL!="",]
pgdata_met <- metil_data_melt_na %>% 
  group_by(SYMBOL) %>% 
  slice(which.min(adj.P.Val)) %>% 
  ungroup()
nrow(pgdata_met)
pgdata_met <- na.omit(pgdata_met)
pgdata_met_extream<-pgdata_met[!duplicated(pgdata_met$SYMBOL),]
# rna prep
pgdata_rna<-rna %>%
  group_by(SYMBOL) 
nrow(pgdata_rna)
pgdata_rna <- na.omit(pgdata_rna)
pgdata_rna_extream<-pgdata_rna[!duplicated(pgdata_rna$SYMBOL),]

### 
merge_extream_1_2_perm <- merge(pgdata_rna_extream,pgdata_atac_extream,by="SYMBOL",how="inner")
merge_extream_1_3_perm <- merge(pgdata_met_extream,pgdata_atac_extream,by="SYMBOL",how="inner")
merge_extream_2_3_perm <- merge(pgdata_rna_extream,pgdata_met_extream,by="SYMBOL",how="inner")
merge_extream_1_2_3 <- merge(merge_extream_1_2_perm,pgdata_met_extream,by="SYMBOL",how="inner")

methyl_sig_perm <-merge_extream_1_3_perm[merge_extream_1_3_perm$adj.P.Val<0.1,]
atac_sig_perm <-merge_extream_1_3_perm[merge_extream_1_3_perm$FDR<0.1,]
rna_sig_perm <- merge_extream_1_2_perm[merge_extream_1_2_perm$padj<0.1,]

methyl_sig <-pgdata_met[pgdata_met$adj.P.Val<0.1,]
atac_sig <-pgdata_atac[pgdata_atac$FDR<0.1,]
rna_sig <-pgdata_rna[pgdata_rna$padj<0.1,]

###
overlap_genes_1_2_perm<-merge(rna_sig_perm,atac_sig_perm,by="SYMBOL",how="inner")
overlap_genes_1_3_perm<-merge(methyl_sig_perm,atac_sig_perm,by="SYMBOL",how="inner")
overlap_genes_2_3_perm<-merge(rna_sig_perm,methyl_sig_perm,by="SYMBOL",how="inner")

overlap_genes_1_2_3 <- merge(overlap_genes_1_2_perm,methyl_sig_perm,by="SYMBOL",how="inner")
write.csv(overlap_genes_1_2_3, "overlap_genes_1_2_3-tes.csv", row.names = FALSE)

merge_extream_1_2_3 <- merge(merge_extream_1_2,pgdata_met_extream,by="SYMBOL",how="inner")

overlap_genes_1_2_3 <- merge(overlap_genes_1_2_perm,methyl_sig_perm,by="SYMBOL",how="inner")

overlap_genes_1_2<-merge(rna_sig_from_merge_extream_dataset,atac_sig_from_merge_extream_dataset,by="SYMBOL",how="inner")
overlap_genes_1_3<-merge(methyl_sig_from_merge_extreame_dataset,atac_sagn_from_merge_extream_dataset,by="SYMBOL",how="inner")
overlap_genes_2_3<-merge(rna_sagn_from_merge_extream_dataset,methyl_sig_from_merge_extreame_dataset,by="SYMBOL",how="inner")
overlap_genes_1_2_3 <- merge(overlap_genes_1_2,methyl_sig_from_merge_extreame_dataset,by="SYMBOL",how="inner")
write.csv(overlap_genes_1_2, "sig_o.1-overlap_genes_1_2.csv", row.names = FALSE)
write.csv(overlap_genes_1_3, "sig_o.1-overlap_genes_1_3.csv", row.names = FALSE)
write.csv(overlap_genes_2_3, "sig_o.1-overlap_genes_2_3.csv", row.names = FALSE)
write.csv(overlap_genes_1_2_3, "sig_o.1-overlap_genes_1_2_3.csv", row.names = FALSE)
gc()
##good!!!
merge_extream_1_2<-merge(pgdata_rna,pgdata_atac,by="SYMBOL",how="inner")
merge_extream_1_3<-merge(pgdata_met,pgdata_atac,by="SYMBOL",how="inner")
merge_extream_2_3<-merge(pgdata_rna,pgdata_met,by="SYMBOL",how="inner")
merge_extream_1_2_3 <- merge(merge_extream_1_2,pgdata_met,by="SYMBOL",how="inner")


## vis vens 
## rna+atac
x<-list("RNA"=rna_sig$SYMBOL,"ATAC"=atac_sig$SYMBOL)
# 2D Venn diagram
venn_diagram<-ggVennDiagram(x,label = "percent", color = "black", lwd = 0.6, lty = 1, label_alpha = 0.1)+
  scale_fill_gradient(low = "pink", high = "blue") 
#show
venn_diagram
## annotation atac+rna
overlap_genes <- intersect(rna_sig$SYMBOL, atac_sig$SYMBOL)
overlap_genes_1_2_new <- merge(atac_sig[atac_sig$SYMBOL %in% overlap_genes, ],
                               rna_sig[rna_sig$SYMBOL %in% overlap_genes, ],
                               by = "SYMBOL")

write.csv(overlap_genes_1_2_new, "sig_o.1-overlap_genes_1_2.csv", row.names = FALSE)
# Subset the ATAC-seq data for the overlapping genes
overlap_atac <- atac_sig %>%
  filter(SYMBOL %in% overlap_genes)
# Create a summary of annotations
annotation_summary <- overlap_atac %>%
  group_by(annotation) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

#######order introns annotations
atac_intron<-annotation_summary %>%
  filter(str_detect(annotation, "intron"))

atac_first_intron<-atac_intron %>%
  filter(str_detect(annotation, "1 of"))

atac_con_intron<-atac_intron %>%
  filter(!str_detect(annotation, "1 of"))

atac_first_intron[,"annotation"]<-"intron 1"

atac_con_intron[,"annotation"]<-"intron >=2"

atac_not_intron<-annotation_summary %>%
  filter(!str_detect(annotation, "intron"))

####order exon annotations
atac_not_intron_or_exon<-atac_not_intron %>%
  filter(!str_detect(annotation, "Exon"))

atac_exon<-atac_not_intron %>%
  filter(str_detect(annotation, "Exon"))
atac_exon[,"annotation"]<-'exon'
#merge all the sub datasets
atac_data_new_annotations<-rbind(atac_not_intron_or_exon,atac_exon,atac_con_intron,atac_first_intron)

#create annotation dataset
atac_annotation_chart<-atac_data_new_annotations%>%
  group_by(annotation)

names(atac_annotation_chart)<-c("annotation","count")

## plot annotations for atac+rna overlaping genes
annotation_plot <- ggplot(atac_annotation_chart, aes(x = annotation, y = count, fill = annotation)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(title = "Annotations of Overlapping Genes-rna+atac",
       x = "Annotation",
       y = "Count") +
  scale_fill_brewer(palette = "Set3")
annotation_plot

#grid for rna+atac overlap
gridExtra::grid.arrange(venn_diagram, annotation_plot, ncol = 2)
### annotation gene dists 
p <- overlap_atac
atac_intron<-p %>%
  filter(str_detect(annotation, "intron"))

atac_first_intron<-atac_intron %>%
  filter(str_detect(annotation, "1 of"))

atac_con_intron<-atac_intron %>%
  filter(!str_detect(annotation, "1 of"))

atac_first_intron[,"annotation"]<-"intron 1"

atac_con_intron[,"annotation"]<-"intron >=2"

atac_not_intron<-p %>%
  filter(!str_detect(annotation, "intron"))

####order exon annotations
atac_not_intron_or_exon<-atac_not_intron %>%
  filter(!str_detect(annotation, "Exon"))

atac_exon<-atac_not_intron %>%
  filter(str_detect(annotation, "Exon"))
atac_exon[,"annotation"]<-'exon'

atac_data_new_annotations<-rbind(atac_not_intron_or_exon,atac_exon,atac_con_intron,atac_first_intron)

#unique annotations for each gene
annotation_count <- atac_data_new_annotations %>%
  group_by(SYMBOL) %>%
  summarise(
    num_annotations = n_distinct(annotation),
    annotations_name = paste0(unique(annotation), collapse = ",")
    )
write.csv(annotation_count, "annotation_distribution_count-rna+atac.csv", row.names = FALSE)
#
#bar plot of the distribution of annotation counts
ggplot(annotation_count, aes(x = num_annotations)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Gene Annotations-rna+atac",
       x = "Number of Annotations per Gene",
       y = "Number of Genes") +
  theme_minimal()

### up/down regulation
overlap_1_2_plot <- overlap_genes_1_2_new %>%
  mutate(Regulation = case_when(
    Fold > 0 & log2FoldChange > 0 ~ "up_up",
    Fold < 0 & log2FoldChange < 0 ~ "down_down",
    Fold > 0 & log2FoldChange < 0 ~ "up_down",
    Fold < 0 & log2FoldChange > 0 ~ "down_up",
    TRUE ~ "other"
  ))

write.csv(overlap_1_2_plot, "overlap_1_2_regulation-rna+atac_overlap.csv", row.names = FALSE)

ggplot(overlap_1_2_plot, aes(x = Regulation)) +
  geom_bar(fill = "#4981BF") +
  labs(x = "Regulation Direction", y = "Count of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Adjust label size if needed


### perms   
permutation.test <- function(data, sagn, genes, n){
  distribution=c()
  result=0
  
  for(i in 1:n){
    sample_index <- sample(1:nrow(data),nrow(sagn))
    sapmle_genes<-data[sample_index,]
    distribution[i]<-nrow(merge(genes,sapmle_genes,by="SYMBOL",how="inner"))
    
  }
  result=sum(distribution >= nrow(overlap_genes))/n
  return(list(result, distribution))
}
## run and plot
test_1_2 <- permutation.test(merge_extream_1_2_perm,rna_sig_perm, atac_sig_perm, 1000)
hist(test_1_2[[2]], breaks=30, col='grey', xlim = range(25,80), main="Permutation Distribution", las=1, xlab='')
abline(v=nrow(overlap_genes_1_2_perm), lwd=3, col="red")
text(x=70, y=65, paste0("pval-", test_1_2[[1]]), col='black', cex=1.2, font=3) 
test_1_2[[1]]
###
## met+atac
x<-list("Methylation"=methyl_sig$SYMBOL,"ATAC"=atac_sig$SYMBOL)
# 2D Venn diagram
venn_diagram<-ggVennDiagram(x,label = "percent", color = "black", lwd = 0.8, lty = 1, label_alpha = 0.1)+
  scale_fill_gradient(low = "pink", high = "blue") 
#show
venn_diagram
## annotation atac+rna
overlap_genes <- intersect(methyl_sig$SYMBOL, atac_sig$SYMBOL)
overlap_genes_1_3_new <- merge(atac_sig[atac_sig$SYMBOL %in% overlap_genes, ],
                               methyl_sig[methyl_sig$SYMBOL %in% overlap_genes, ],
                               by = "SYMBOL")
write.csv(overlap_genes_1_3_new, "sig_o.1-overlap_genes_1_3.csv", row.names = FALSE)
# Subset the ATAC-seq data for the overlapping genes
overlap_atac <- atac_sig %>%
  filter(SYMBOL %in% overlap_genes)
# Create a summary of annotations

annotation_summary <- overlap_atac %>%
  group_by(annotation) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

##
#######order introns annotations
atac_intron<-annotation_summary %>%
  filter(str_detect(annotation, "intron"))

atac_first_intron<-atac_intron %>%
  filter(str_detect(annotation, "1 of"))

atac_con_intron<-atac_intron %>%
  filter(!str_detect(annotation, "1 of"))

atac_first_intron[,"annotation"]<-"intron 1"

atac_con_intron[,"annotation"]<-"intron >=2"

atac_not_intron<-annotation_summary %>%
  filter(!str_detect(annotation, "intron"))

####order exon annotations
atac_not_intron_or_exon<-atac_not_intron %>%
  filter(!str_detect(annotation, "Exon"))

atac_exon<-atac_not_intron %>%
  filter(str_detect(annotation, "Exon"))
atac_exon[,"annotation"]<-'exon'
#merge all the sub datasets
atac_data_new_annotations<-rbind(atac_not_intron_or_exon,atac_exon,atac_con_intron,atac_first_intron)

#create annotation dataset
atac_annotation_chart<-atac_data_new_annotations%>%
  group_by(annotation)

names(atac_annotation_chart)<-c("annotation","count")

## plot atac+met annotations 
annotation_plot <- ggplot(atac_annotation_chart, aes(x = annotation, y = count, fill = annotation)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  labs(title = "Annotations of Overlapping Genes-met+atac",
       x = "Annotation",
       y = "Count") +
  scale_fill_brewer(palette = "Set3")
annotation_plot
#grid for atac+met overlap 
gridExtra::grid.arrange(venn_diagram, annotation_plot, ncol = 2)
### annotation gene dists 
p <- overlap_atac
atac_intron<-p %>%
  filter(str_detect(annotation, "intron"))

atac_first_intron<-atac_intron %>%
  filter(str_detect(annotation, "1 of"))

atac_con_intron<-atac_intron %>%
  filter(!str_detect(annotation, "1 of"))

atac_first_intron[,"annotation"]<-"intron 1"

atac_con_intron[,"annotation"]<-"intron >=2"

atac_not_intron<-p %>%
  filter(!str_detect(annotation, "intron"))

####order exon annotations
atac_not_intron_or_exon<-atac_not_intron %>%
  filter(!str_detect(annotation, "Exon"))

atac_exon<-atac_not_intron %>%
  filter(str_detect(annotation, "Exon"))
atac_exon[,"annotation"]<-'exon'

atac_data_new_annotations<-rbind(atac_not_intron_or_exon,atac_exon,atac_con_intron,atac_first_intron)

#unique annotations for each gene
annotation_count <- atac_data_new_annotations %>%
  group_by(SYMBOL) %>%
  summarise(
    num_annotations = n_distinct(annotation),
    annotations_name = paste0(unique(annotation), collapse = ",")
  )

#
write.csv(annotation_count, "annotation_distribution-met+atac_overlap.csv", row.names = FALSE)

#bar plot of the distribution of annotation counts
ggplot(annotation_count, aes(x = num_annotations)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Gene Annotations-met_atac",
       x = "Number of Annotations per Gene",
       y = "Number of Genes") +
  theme_minimal()

### up/down 
overlap_1_3_plot <- overlap_genes_1_3_new %>%
  mutate(Regulation = case_when(
    Fold > 0 & logFC > 0 ~ "up_up",
    Fold < 0 & logFC < 0 ~ "down_down",
    Fold > 0 & logFC < 0 ~ "up_down",
    Fold < 0 & logFC > 0 ~ "down_up",
    TRUE ~ "other"
  ))
write.csv(overlap_1_3_plot, "overlap_1_3_regulation-met+atac_overlap.csv", row.names = FALSE)
ggplot(overlap_1_3_plot, aes(x = Regulation)) +
  geom_bar(fill = "#4981BF") +
  labs(x = "Regulation Direction", y = "Count of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Adjust label size if needed

### perms   
test_1_3 <- permutation.test(merge_extream_1_3_perm, methyl_sig_perm, atac_sig_perm, 1000)
hist(test_1_3[[2]], breaks=30, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=nrow(overlap_genes_1_3_perm), lwd=3, col="red")
text(x=68, y=60, paste0("pval-", test_1_3[[1]]), col='black', cex=1.2, font=3) 
test_1_3[[1]]

## vis vens 
## met+rna
x<-list("Methylation"=methyl_sig$SYMBOL,"RNA"=rna_sig$SYMBOL)
# 2D Venn diagram
venn_diagram<-ggVennDiagram(x,label = "percent", color = "black", lwd = 0.8, lty = 1, label_alpha = 0.1)+
  scale_fill_gradient(low = "pink", high = "blue") 
#show
venn_diagram
## annotation atac+rna
overlap_genes <- intersect(rna_sig$SYMBOL, methyl_sig$SYMBOL)
overlap_genes_2_3_new <- merge(rna_sig[rna_sig$SYMBOL %in% overlap_genes, ],
                               methyl_sig[methyl_sig$SYMBOL %in% overlap_genes, ],
                               by = "SYMBOL")
write.csv(overlap_genes_2_3_new, "sig_o.1-overlap_genes_2_3.csv", row.names = FALSE)
### up/down 
overlap_2_3_plot <- overlap_genes_2_3_new %>%
  mutate(Regulation = case_when(
    log2FoldChange > 0 & logFC > 0 ~ "up_up",
    log2FoldChange < 0 & logFC < 0 ~ "down_down",
    log2FoldChange > 0 & logFC < 0 ~ "up_down",
    log2FoldChange < 0 & logFC > 0 ~ "down_up",
    TRUE ~ "other"
  ))

write.csv(overlap_2_3_plot, "overlap_2_3_plot_regulation-met+rna_overlap.csv", row.names = FALSE)
ggplot(overlap_2_3_plot, aes(x = Regulation)) +
  geom_bar(fill = "#4981BF") +
  labs(x = "Regulation Direction", y = "Count of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  # Adjust label size if needed

### perms   
test_2_3 <- permutation.test(merge_extream_2_3_perm,methyl_sig_perm, rna_sig_perm, 1000)
hist(test_2_3[[2]], breaks=30, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=nrow(overlap_genes_2_3_perm), lwd=3, col="red")
text(x=88, y=95, paste0("pval-", test_2_3[[1]]), col='black', cex=1.2, font=3) 
test_2_3[[1]]

##all 
# List of items
x<-list("ATAC"=atac_sig$SYMBOL,"Methylation"=methyl_sig$SYMBOL, "RNA"=rna_sig$SYMBOL)
# 2D Venn diagram
venn_diagram<-ggVennDiagram(x, color = "black", lwd = 0.8, lty = 1, label_alpha = 0.1)+
  scale_fill_gradient(low = "pink", high = "blue") 
#show
venn_diagram
##
overlap_genes <- intersect(rna_sig$SYMBOL, methyl_sig$SYMBOL, atac_sig$SYMBOL)
overlap_genes_1_2_3 <- merge(overlap_genes_1_2_new,methyl_sig,by="SYMBOL",how="inner")
### up/down 
overlap_1_2_3_plot <- overlap_genes_1_2_3 %>%
  mutate(Regulation = case_when(
    Fold > 0 & log2FoldChange > 0 & logFC > 0 ~ "up_up_up",
    Fold < 0 & log2FoldChange < 0 & logFC < 0 ~ "down_down_down",
    Fold > 0 & log2FoldChange > 0 & logFC < 0 ~ "up_up_down",
    Fold > 0 & log2FoldChange < 0 & logFC > 0 ~ "up_down_up",
    Fold < 0 & log2FoldChange > 0 & logFC > 0 ~ "down_up_up",
    Fold > 0 & log2FoldChange < 0 & logFC < 0 ~ "up_down_down",
    Fold < 0 & log2FoldChange < 0 & logFC > 0 ~ "down_down_up",
    Fold < 0 & log2FoldChange > 0 & logFC < 0 ~ "down_up_down",
    TRUE ~ "other"
  ))
write.csv(overlap_1_2_3_plot, "overlap_1_2_3_plot_regulation-met+atac+RNA_overlap.csv", row.names = FALSE)

ggplot(overlap_1_2_3_plot, aes(x = Regulation)) +
  geom_bar(fill = "#4981BF") +
  labs(x = "Regulation Direction", y = "Count of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))  

## enrichr plot 
# function 1
plot_gene_pathway_heatmap <- function(data, pathways_col = "pathway", genes_col = "gene", score_col = "score") {
  
  #binary matrix for gene-pathway membership
  heatmap_data <- data %>%
    mutate(membership = 1) %>%
    select(!!sym(genes_col), !!sym(pathways_col), membership) %>%
    pivot_wider(names_from = !!sym(pathways_col), values_from = membership, values_fill = list(membership = 0)) %>%
    pivot_longer(cols = -!!sym(genes_col), names_to = pathways_col, values_to = "membership")
  
  #score normalization
  pathway_scores <- data %>%
    group_by(!!sym(pathways_col)) %>%
    summarize(avg_score = mean(!!sym(score_col), na.rm = TRUE)) %>%
    mutate(norm_score = rescale(avg_score, to = c(0, 1)))  # normalize to [0, 1]
  
  #plot
  heatmap_plot <- ggplot(heatmap_data, aes_string(x = pathways_col, y = genes_col)) +
    geom_point(aes(color = as.factor(membership)), shape = 21, size = 4) +
    scale_color_manual(values = c("0" = "grey", "1" = "darkblue"), guide = "none") +
    theme_minimal() +
    labs(x = "Pathways", y = "Genes", title = "Gene-Pathway Membership") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #normalized enrichment scores as a color gradient below the pathways
  score_plot <- ggplot(pathway_scores, aes_string(x = pathways_col, y = 1, fill = "norm_score")) +
    geom_tile() +
    scale_fill_gradient(low = "yellow", high = "purple", name = "normalized enrichment score") +
    theme_void() +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  #combine plots
  combined_plot <- heatmap_plot / score_plot + plot_layout(heights = c(18, 2))
  
  print(combined_plot)
}


# function 2
plot_gene_pathway_heatmap <- function(data, pathways_col = "pathway", genes_col = "gene", score_col = "score") {
  
  #binary matrix indicating gene-pathway membership
  heatmap_data <- data %>%
    mutate(membership = 1) %>%
    select(!!sym(genes_col), !!sym(pathways_col), membership, !!sym(score_col))
  
  #normalize scores for size scaling
  heatmap_data <- heatmap_data %>%
    mutate(norm_score = rescale(!!sym(score_col), to = c(1, 10)))  # Adjust range as needed for visibility
  
  #plot
  heatmap_plot <- ggplot(heatmap_data, aes_string(x = pathways_col, y = genes_col)) +
    geom_point(aes(size = norm_score, fill = as.factor(membership)), shape = 21, color = "black") +
    scale_size_continuous(name = "Enrichment Score") +
    scale_fill_manual(values = c("0" = "grey", "1" = "purple"), guide = "none") +
    theme_minimal() +
    labs(x = "Pathways", y = "Genes", title = "gene-pathway with enrichment scores") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(heatmap_plot)
}

##load df
df = read.csv("integ_enrich.csv")
plot_gene_pathway_heatmap(df)
##########################################################################################