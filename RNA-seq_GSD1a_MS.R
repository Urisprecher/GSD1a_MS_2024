#libs
library(data.table)  
library(ggplot2)   
library(readr)
library(airway)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(genefilter)
library(PoiClaClu)
library(AnnotationDbi)
library(dplyr)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(dplyr)
library(ggfortify)
library(pcaExplorer)
library(cli)
library(topGO)
library(markdown)
library(PCAtools)
library(DEGreport)
library(radiant.data)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(pathview)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(AnnotationHub)
library(dplyr)
library(tibble)
##install 
BiocManager::install('')
##options 
options(ggrepel.max.overlaps = Inf)
colramp = colorRampPalette(c(3,"white",2))(20)
orgainsm = 'org.Hs.eg.db'
gc()
### full data - QC check
#original 
count_dir <- "D:/MiguelW12/Documents/gsd1a_rna_2024/count_initial_all/count_matrices_text_star"

#read and merge count files frim different samples
merge_count_files <- function(count_dir) {
  #list files in dir
  count_files <- list.files(count_dir, pattern = "\\.txt$", full.names = TRUE)
  print(count_files)
  count_data <- lapply(count_files, function(file) {
    sample_name <- gsub("\\.counts\\.txt", "", basename(file))
    counts <- fread(file, header = FALSE)
    colnames(counts) <- c("GeneID", sample_name)
    return(counts)
    
  })
  
  #merge by GeneID
  print(count_data)
  merged_count_data <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), count_data)
  return(merged_count_data)
}

#run
count_data <- merge_count_files(count_dir)
write.csv(as.data.frame(count_data), file="count_data_gsd1a_2024_rna_star.csv")


###basic statistics
summary_stats <- apply(count_data[, -1], 2, summary)
write.csv(as.data.frame(summary_stats), file="star- count_data_summary_stats.csv")
# count zero counts for each sample
zero_counts <- colMeans(count_data[, -1] == 0, na.rm = TRUE) * 100
print(zero_counts)
### qc met
create_qc_metrics <- function(count_data, output_dir) {

  dir.create(output_dir, showWarnings = FALSE)
  
  #data frame for qc metrics
  qc_df <- data.frame(Sample = character(),
                      Total_Genes = numeric(),
                      Num_NAs = numeric(),
                      Num_Zeros = numeric(),
                      stringsAsFactors = FALSE)
  
  for (col_name in colnames(count_data)) {
    counts <- count_data[[col_name]]
    
    total_genes <- length(counts)
    num_nas <- sum(is.na(counts))
    num_zeros <- sum(counts == 0)
    
    qc_df <- rbind(qc_df, data.frame(Sample = col_name,
                                     Total_Genes = total_genes,
                                     Num_NAs = num_nas,
                                     Num_Zeros = num_zeros))
  }
  
  write.csv(qc_df, file.path(output_dir, "qc_metrics.csv"), row.names = FALSE)
}

create_qc_metrics(count_data, "gsd1a_rna_2024/star-count_stats")
## vis qc
# hists
create_count_histograms <- function(count_data, output_dir = "path/to/histograms") {
  dir.create(output_dir, showWarnings = FALSE)
  
  for (col_name in colnames(count_data)) {
    if (col_name == "GeneID") next
    
    counts <- as.numeric(count_data[[col_name]])
    
    if (all(!is.na(counts)) && all(counts >= 0)) {
      counts <- counts[counts > 0]
      
      log_counts <- log10(counts)
      
      png(file.path(output_dir, paste0(col_name, "_histogram.png")))
      hist(log_counts, main = paste("Histogram of Counts -", col_name),
           xlab = "Log10(Counts)", ylab = "Frequency", xlim = c(min(log_counts), max(log_counts)),
           freq = TRUE, breaks = "Sturges")
      dev.off()
      
      missing_values <- sum(is.na(counts))
      if (missing_values > 0) {
        cat(paste("Number of missing values for", col_name, ":", missing_values, "\n"))
      }
      
      cat("histogram saved-", file.path(output_dir, paste0(col_name, "_histogram.png")), "\n")
    } else {
      warning(paste("invalid count data in sample", col_name, "- skipping histogram generation."))
    }
  }
}


create_count_histograms(count_data, output_dir = "gsd1a_rna_2024/star-count_histograms")

#boxes per sample
create_count_boxplots <- function(count_data, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  
  
  for (col_name in colnames(count_data)) {
    if (col_name == "GeneID") next
    
    counts <- na.omit(count_data[[col_name]])
    counts <- counts[counts > 0]  # Filter out zero and negative counts
    
    df <- data.frame(Sample = rep(col_name, length(counts)),
                     Counts = counts)
    
    p <- ggplot(df, aes(x = Sample, y = Counts)) +
      geom_boxplot(fill = "skyblue", color = "black") +
      labs(title = paste("Boxplot of Counts -", col_name),
           x = "Sample", y = "Counts") +
      scale_y_continuous(trans = "log10")
    
    png(paste(output_dir, "/", col_name, "_boxplot.png", sep = ""), width = 800, height = 600)
    print(p)
    dev.off()
  }
}

create_count_boxplots(count_data, "gsd1a_rna_2024/star-count_BOXPLOTS")

## boxplot all 
create_count_boxplots2 <- function(count_data, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  
  
  
  all_counts <- data.frame(Sample = character(),
                           Counts = numeric(),
                           stringsAsFactors = FALSE)
  
  for (col_name in colnames(count_data)) {
    if (col_name == "GeneID") next
    
    counts <- na.omit(count_data[[col_name]])
    counts <- counts[counts > 0]  
    
    all_counts <- rbind(all_counts, data.frame(Sample = rep(col_name, length(counts)),
                                               Counts = counts,
                                               stringsAsFactors = FALSE))
  }
  
  p <- ggplot(all_counts, aes(x = Sample, y = Counts, fill = Sample)) +
    geom_boxplot(color = "black") +
    labs(title = "Boxplot of Counts",
         x = "Sample", y = "Counts") +
    scale_y_continuous(trans = "log10") +
    theme(legend.position = "top")
  
  output_path <- file.path(output_dir, "all_samples_boxplot.png")
  png(output_path, width = 800, height = 600)
  print(p)
  dev.off()
  
  cat("boxplot saved as", output_path, "\n")
}

create_count_boxplots2(count_data, "gsd1a_rna_2024/star-count_BOXPLOTS")



### differential prep
Count_Data_ALL = read.csv('D:/MiguelW12/Documents/gsd1a_rna_2024/final_res(male_ut)/star_all/data_frames/count_data_gsd1a_2024_rna_all_star.csv', header=TRUE)
Col_Data_ALL = read.csv('D:/MiguelW12/Documents/gsd1a_rna_2024/final_res(male_ut)/star_all/data_frames/sample_data_gsd1a_2024_rna_ut_star.csv', header=TRUE)
rownames(Count_Data_ALL) <- Count_Data_ALL$GeneID
Count_Data_ALL = subset(Count_Data_ALL, select = -c(GeneID))
rownames(Col_Data_ALL) <- Col_Data_ALL$X
Col_Data_ALL  = subset(Col_Data_ALL, select = -c(X))
all(colnames(Count_Data_ALL) %in% rownames(Col_Data_ALL))
######
## desqe build
dds_all <- DESeq2::DESeqDataSetFromMatrix(countData = Count_Data_ALL, colData = Col_Data_ALL, design = ~ group)
dds_all
#add filter for counts
# run deseq
keep_all <- rowMeans(counts(dds_all)) >= 40
dds_all <- dds_all[keep_all,]
dds_all
dds_all_run <- DESeq2::DESeq(dds_all)
resultsNames(dds_all_run)
### hcvsgsd
res_all <- DESeq2::results(dds_all_run, contrast = c("group","HC_UT", "GSD_UT"))
#### vis results
hist(res_all$pvalue)
hist(res_all$padj)
hist(res_all$log2FoldChange)
summary(res_all)
write.csv(as.data.frame(res_all), file="diff_analysis_res_star.csv")
###anno
orgainsm = 'org.Hs.eg.db'
head(res_all)
rownames(res_all)
anno_2 <- AnnotationDbi::select(org.Hs.eg.db, rownames(res_all), columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME", "CHR"),keytype = "ENSEMBL")
keytypes(org.Hs.eg.db)
keys(org.Hs.eg.db)
head(anno_2)
res_star_anno = cbind( ENSEMBL = rownames( res_all), res_all)
outtable_star <- left_join(as.data.frame(res_star_anno), anno_2 )
outtable_star <- outtable_star[order(outtable_star$padj),]
write.csv(as.data.frame(outtable_star), file="diff_analysis_ANO-star_res.csv")
#EV
EnhancedVolcano::EnhancedVolcano(as.data.frame(outtable_star), lab = outtable_star$SYMBOL, x = 'log2FoldChange', y = 'padj', xlim = c(-8,8), title = 'HC vs Gsd', pCutoff = 0.1, FCcutoff = 1, pointSize = 2.0, labSize = 4.0)
## EXPORT NORM DAT 
norm_data_star <- counts(dds_all_run, normalized=TRUE)
write.csv(as.data.frame(norm_data_star ), file="star-normatotaldata.csv")
## plot main genes following examining lists
## load dfs 
Imp_genes_fold = read.csv('D:/MiguelW12/Documents/gsd1a_rna_2024/final_res(male_ut)/star_all/res/filt40/list of major genes -fold.csv', header=TRUE)
Imp_genes_class = read.csv('D:/MiguelW12/Documents/gsd1a_rna_2024/final_res(male_ut)/star_all/res/filt40/list of major genes -class.csv', header=TRUE)
sorted_imp_genes_df <- merge(Imp_genes_fold, Imp_genes_class, by = "Gene", all.x = TRUE)
sorted_imp_genes_df$Gene <- factor(sorted_imp_genes_df$Gene, levels = sorted_imp_genes_df$Gene[order(sorted_imp_genes_df$Class)])
# colors for each class
class_colors <- c("Hox_genes" = "skyblue", "Metabolic_genes" = "orange", "Fibroblast_genes" = "yellow", "Wnt_pathway_genes" = "purple", "Lysosomal/Mitochondrial_genes" = "red")  
#plot
level_order <- c("Hox_genes", "Metabolic_genes", "Fibroblast_genes", "Wnt_pathway_genes", "Lysosomal/Mitochondrial_genes")
ggplot(sorted_imp_genes_df, aes(x = Gene, y = L2FC, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = class_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(title = "Genes", x = "Gene", y = "L2FC(HC/GSD1a")

### PCA 
rdat <- assay(rlog(dds_all))
p <- pca(rdat, metadata = colData(dds_all), removeVar = 0.5)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
elbow <- findElbowPoint(p$variance)
elbow
horn <- parallelPCA(rdat)
horn$n
biplot(p,
       colby = 'group',
       colLegendTitle = 'pca_gsd1a_rnaseq',
       # encircle config
       encircle = TRUE,
       encircleFill = TRUE,
       ellipseLevel = 0.95,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

#gene_set_enrichment_analysis
new_df_Ut_mean = read.csv('D:/MiguelW12/Documents/gsd1a_rna_2024/final_res(male_ut)/star_all/res/filt40/res_initial-star-40filt.csv', header=TRUE)
original_gene_list <- new_df_Ut_mean$log2FoldChange
names(original_gene_list) <- new_df_Ut_mean$X
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse_ut_n4 <- gseGO(geneList=gene_list, 
                ont ="ALL", 
                keyType = "ENSEMBL", 
                minGSSize = 20, 
                maxGSSize = 10000, 
                pvalueCutoff = 0.1,
                nPerm = 10000,
                verbose = TRUE, 
                OrgDb = 'org.Hs.eg.db',
                pAdjustMethod = "fdr")

#save file
write.csv(as.data.frame(gse_ut_n), file="gseago_reg_male_ut_all-filt40star-for sup8.csv")
##
id_list_gsea_plot <- c("GO:1902493",
                            "GO:0004402",
                            "GO:0000123",
                            "GO:0061733",
                            "GO:0010628",
                            "GO:0031248",
                            "GO:0008276",
                            "GO:0008135",
                            "GO:0090079",
                            "GO:0006418") 
gse_ut_n4@result = gse_ut_n4@result[gse_ut_n4@result$ID %in% id_list_gsea_plot, ]
dotplot(gse_ut_n4, showCategory=30, font.size = 8,label_format = 60, color = "p.adjust", title = "dotplot_GSEA", split=".sign") + facet_grid(.~.sign)

########################################################################################






















