# DAPseq_pipeline_nf
A workflow for RNA-seq analysis.

## The workflow includes the following steps:
1. Reads trimming (trim_galore)
2. Clean reads mapping (bowtie2)
3. reads count (featurecounts)
4. Quality control (sample clustering and PCA analysis)
5. differential expression analysis (DESeq2)

## Prerequisite
```
nextflow program
conda environment: DAPseq_env (trim_galore and bedtools and bowtie2 and MACS3 and HOMER and MEME suite installed)
R packages: pheatmap, optparse, factoextra, dplyr, DT, kableExtra, base64enc, edgeR, DESeq2
Indexed genome

```

## Test the workflow
```
git clone git@github.com:mmyoung/RNAseq_pipeline_nf.git
conda activate DAPseq_env
nextflow run /project/gzy8899/lyang/RNAseq_pipeline_nf -params-file params.yml  ## yml file saving all parameters, refer to the file in ./test folder for format

```

## Parameters

Parameters can be passed to the pipeline in the command line or be put in a .yml file and passed to the pipeline in a whole (-params-file params.yml)
```
--fq_sheet A tab-delimited file storing the samples information, with five columns: sample,fq1,fq2,single_end,control
--fasta    Genome fasta file for the analyzing species.
--gtf    Genome gtf file for the analyzing species.
--output_dir    Name for directory for saving the results. Default: ./results
--fq_dir  The folder where the raw .fastq files are.
--gsize The size of analyzing genome.
--bowtie_idx The bowtei2 index directory. ## built with bowtie2-build command
--DESeq2_coldata Meta information of the treated and control samples used for differential expression analysis.
--strandness Library strandness for featureCounts: 0 (unstranded), 1 (stranded), 2 (reverse stranded). Default: 1
--pairend Whether the library is pair-end. 1 (pair-end), 0 (single-end). Default: 1
```

**Caveats**:
* check the required packages are all installed.
* need to go through the scripts to make sure the path to softwares are executable for current user.
* pay attention to the path to cutadapt tool which is required to execute trim_galore and installation of it isn't compatibale with the DAPseq_env environment

## example of sample sheet passed to --fq_sheet
fq_sheet.csv (comma separated meta data)
```
sample,fq1,fq2,single_end,control
IP,SRR27496336_1.fastq,SRR27496336_2.fastq,0,Input
Input,SRR27496337_1.fastq,SRR27496337_2.fastq,0,
```
sample: sample name or id, need to be unique
fq1: read1 fastq name for pair-end sequencing (the only fastq file name for single-end sequencing)
fq2: read2 fastq name for pair-end sequencing (empty for single-end sequencing)
single_end: indicate whether the data is from single-end sequencing (1: single-end; 0: pair-end)
control; the corresponding control sample name/id for the DAP library (empty if there is no control)

## example of DESeq2_coldata passed to --DESeq2_coldata
```
control_id,control_group,treated_id,treated_group,comparison_name
ctr1;ctr2;ctr3,ctrl;ctrl;ctrl,treated1;treated2;treated3,treated;treated;treated,treated_vs_ctrl
ctr1;ctr2;ctr3,ctrl;ctrl;ctrl,group1;group2;group3,treated;treated;treated,group_vs_ctrl
```
control_id: samples IDs that are in the control (untreated) group, need to be unique and match the IDs in sample sheet.
control_group: the group labels for control samples, should be ctrl, the count should be consistent with the control_id column.
treated_id: the same operation as control samples.
treated_group: should be treated, count matches the treated_id column.
comparison_name: the names for corresponding comparison, mainly used for the file names saved from DESeq2.

## Results 
The intermediate output from all procedures are saved in the output directory the user defined, the structure looks like the following:
test_out/
├── alignment
│   ├── ctr1.bowtie2.log
│   ├── ctr1_Q20_sorted.bam
│   ├── ctr1_Q20_sorted.bam.bai
│   ├── ctr2.bowtie2.log
│   ├── ctr2_Q20_sorted.bam
│   ├── ctr2_Q20_sorted.bam.bai
│   ├── ctr3.bowtie2.log
│   ├── ctr3_Q20_sorted.bam
│   ├── ctr3_Q20_sorted.bam.bai
│   ├── group1.bowtie2.log
│   ├── group1_Q20_sorted.bam
│   ├── group1_Q20_sorted.bam.bai
│   ├── group2.bowtie2.log
│   ├── group2_Q20_sorted.bam
│   ├── group2_Q20_sorted.bam.bai
│   ├── group3.bowtie2.log
│   ├── group3_Q20_sorted.bam
│   ├── group3_Q20_sorted.bam.bai
│   ├── treated1.bowtie2.log
│   ├── treated1_Q20_sorted.bam
│   ├── treated1_Q20_sorted.bam.bai
│   ├── treated2.bowtie2.log
│   ├── treated2_Q20_sorted.bam
│   ├── treated2_Q20_sorted.bam.bai
│   ├── treated3.bowtie2.log
│   ├── treated3_Q20_sorted.bam
│   └── treated3_Q20_sorted.bam.bai
├── bw_output
│   ├── ctr1_sorted_bam.bw
│   ├── ctr2_sorted_bam.bw
│   ├── ctr3_sorted_bam.bw
│   ├── group1_sorted_bam.bw
│   ├── group2_sorted_bam.bw
│   ├── group3_sorted_bam.bw
│   ├── treated1_sorted_bam.bw
│   ├── treated2_sorted_bam.bw
│   └── treated3_sorted_bam.bw
├── clean_reads_fastQC
│   ├── ctr1
│   │   ├── ctr1_1_fastqc.html
│   │   ├── ctr1_1_fastqc.zip
│   │   ├── ctr1_2_fastqc.html
│   │   └── ctr1_2_fastqc.zip
│   ├── ctr2
│   │   ├── ctr2_1_fastqc.html
│   │   ├── ctr2_1_fastqc.zip
│   │   ├── ctr2_2_fastqc.html
│   │   └── ctr2_2_fastqc.zip
│   ├── ctr3
│   │   ├── ctr3_1_fastqc.html
│   │   ├── ctr3_1_fastqc.zip
│   │   ├── ctr3_2_fastqc.html
│   │   └── ctr3_2_fastqc.zip
│   ├── group1
│   │   ├── group1_1_fastqc.html
│   │   ├── group1_1_fastqc.zip
│   │   ├── group1_2_fastqc.html
│   │   └── group1_2_fastqc.zip
│   ├── group2
│   │   ├── group2_1_fastqc.html
│   │   ├── group2_1_fastqc.zip
│   │   ├── group2_2_fastqc.html
│   │   └── group2_2_fastqc.zip
│   ├── group3
│   │   ├── group3_1_fastqc.html
│   │   ├── group3_1_fastqc.zip
│   │   ├── group3_2_fastqc.html
│   │   └── group3_2_fastqc.zip
│   ├── treated1
│   │   ├── treated1_1_fastqc.html
│   │   ├── treated1_1_fastqc.zip
│   │   ├── treated1_2_fastqc.html
│   │   └── treated1_2_fastqc.zip
│   ├── treated2
│   │   ├── treated2_1_fastqc.html
│   │   ├── treated2_1_fastqc.zip
│   │   ├── treated2_2_fastqc.html
│   │   └── treated2_2_fastqc.zip
│   └── treated3
│       ├── treated3_1_fastqc.html
│       ├── treated3_1_fastqc.zip
│       ├── treated3_2_fastqc.html
│       └── treated3_2_fastqc.zip
├── differential_expression
│   ├── group_vs_ctrl.tsv
│   ├── treated_vs_ctrl.tsv
│   └── versions.yml
├── featurecounts
│   ├── all.featureCounts.tsv
│   ├── all.featureCounts.tsv.summary
│   └── versions.yml
├── report
│   ├── all_sample_fpkm_matrix.csv
│   ├── all_sample_raw_count.csv
│   ├── read_peak.num.summary
│   ├── report.html
│   ├── sample_corr_heatmap.png
│   └── sample_pca_dotplot.png
└── trimm
    ├── ctr1_1.fastq_trimming_report.txt
    ├── ctr1_2.fastq_trimming_report.txt
    ├── ctr1_val_1.fq.gz
    ├── ctr1_val_2.fq.gz
    ├── ctr2_1.fastq_trimming_report.txt
    ├── ctr2_2.fastq_trimming_report.txt
    ├── ctr2_val_1.fq.gz
    ├── ctr2_val_2.fq.gz
    ├── ctr3_1.fastq_trimming_report.txt
    ├── ctr3_2.fastq_trimming_report.txt
    ├── ctr3_val_1.fq.gz
    ├── ctr3_val_2.fq.gz
    ├── group1_1.fastq_trimming_report.txt
    ├── group1_2.fastq_trimming_report.txt
    ├── group1_val_1.fq.gz
    ├── group1_val_2.fq.gz
    ├── group2_1.fastq_trimming_report.txt
    ├── group2_2.fastq_trimming_report.txt
    ├── group2_val_1.fq.gz
    ├── group2_val_2.fq.gz
    ├── group3_1.fastq_trimming_report.txt
    ├── group3_2.fastq_trimming_report.txt
    ├── group3_val_1.fq.gz
    ├── group3_val_2.fq.gz
    ├── treated1_1.fastq_trimming_report.txt
    ├── treated1_2.fastq_trimming_report.txt
    ├── treated1_val_1.fq.gz
    ├── treated1_val_2.fq.gz
    ├── treated2_1.fastq_trimming_report.txt
    ├── treated2_2.fastq_trimming_report.txt
    ├── treated2_val_1.fq.gz
    ├── treated2_val_2.fq.gz
    ├── treated3_1.fastq_trimming_report.txt
    ├── treated3_2.fastq_trimming_report.txt
    ├── treated3_val_1.fq.gz
    └── treated3_val_2.fq.gz
    
```
The report.html output is the one to check with some important statistics report about the data.

## Report explanation
```
sample: sample names
raw reads pairs#: total raw read (pairs for pair-end sequencing)
clean read pairs#: total read pairs passing trim_galore filters (-q 20, --length 20)
mapping ratio: % of clean reads that are mapped to genome
unique mapped pairs#: total reads number that are uniquely mapped to genome (filtered with samtools view -q 20)

heatmap: sample clustering results based on FPKM normalized expression values.
PCA plot: the first two PCs calculated based on FPKM values.
```

