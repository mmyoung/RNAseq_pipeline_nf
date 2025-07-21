library(dplyr)
library(tibble)
library(tidyverse)
library(pheatmap)
library(optparse)
library(stats)
library(edgeR)
library(factoextra)

option_list <- list(
    make_option(c("-i", "--count_matrix"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

## read the full count matrix
raw_count_mat <- read.table(opt$count_matrix,header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
new_colnames <- gsub("_Q20_sorted.bam","",basename(colnames(raw_count_mat)))
colnames(raw_count_mat) <- new_colnames

# Create DGEList object
gene_lengths <- raw_count_mat$Length / 1000
gene_ids <- raw_count_mat$Geneid
count_matrix <- raw_count_mat %>% 
  select(-Chr,-Start,-End,-Strand,-Length) %>%
  column_to_rownames(var="Geneid")

write.csv(count_matrix,file="all_sample_raw_count.csv",row.names=T,quote=F)

dge <- DGEList(counts = count_matrix)

# Get library sizes in millions
lib_sizes <- dge$samples$lib.size / 1e6

# Calculate FPKM
fpkm_matrix <- t(t(count_matrix) / lib_sizes) / gene_lengths
rownames(fpkm_matrix) <- gene_ids
colnames(fpkm_matrix) <- colnames(count_matrix)

write.csv(fpkm_matrix,file="all_sample_fpkm_matrix.csv",row.names=T,quote=F)

# Log-transform for visualization
log_fpkm <- log2(fpkm_matrix + 1)

# Compute correlation matrix
cor_matrix <- cor(log_fpkm, method = "pearson")

# Plot heatmap
png("sample_corr_heatmap.png")
pheatmap(cor_matrix, display_numbers = TRUE, main = "Sample Correlation (FPKM)")
dev.off()

pca <- prcomp(t(log_fpkm))

png("sample_pca_dotplot.png")
fviz_pca_ind(pca,
             #palette = c("#00AFBB",  "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
dev.off()