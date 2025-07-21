#!/usr/bin/env Rscript
## the script performs 

library(optparse)
library(DESeq2)
library(ggplot2)
library(stringr)

## PARSE COMMAND-LINE PARAMETERS ##

option_list <- list(
    make_option(c("-i", "--count_file"), type="character", default=NULL, metavar="path", help="Count file matrix where rows are genes and columns are samples."),
    make_option(c("-f", "--col_data"), type="character"  , default=3, metavar="integer", help="First column containing sample count data."),
    #make_option(c("-o", "--outdir"), type="character", default='./', metavar="path", help="Output directory."),
    make_option(c("-c", "--cores"), type="integer"  , default=1, metavar="integer", help="Number of cores.")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)


## read full count matrix
full_count_mat <- read.csv(file=opt$count_file,header=T,row.names=1)

## read meta information
col_meta <- read.csv(file=opt$col_data,header=T,row.names=NULL)
## col_meta has five columns: comma separated control sample IDs, comma separated control sample group, same for treated samples, name of comparison (used for prefix of output)
## if multiple comparisons included, there'll be multiple rows

for(n_row in seq_len(nrow(col_meta))){

    ## build colData for each comparison
    coldata=data.frame(sample=c(str_split(col_meta[n_row,"control_id"],";",simplify=T),str_split(col_meta[n_row,"treated_id"],";",simplify=T)),
            condition=c(str_split(col_meta[n_row,"control_group"],";",simplify=T),str_split(col_meta[n_row,"treated_group"],";",simplify=T)))
    row.names(coldata) <- coldata$sample

    sample_count_matrix <- full_count_mat[,coldata[,"sample"]]

    dds <- DESeqDataSetFromMatrix(countData=sample_count_matrix, colData=coldata, design= ~condition)
    #dds <- estimateSizeFactors(dds)

    dds$condition <- factor(dds$condition, levels = c("ctrl","treated"))

    # Differential expression analysis
    dds <- DESeq(dds)
    res <- results(dds)

    write.table(as.data.frame(res),
                file=paste0(col_meta[n_row,"comparison_name"],".tsv"),
                quote=F,row.names=T,sep="\t")
    
}
