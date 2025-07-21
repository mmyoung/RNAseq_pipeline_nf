process REPORT {
    label 'report'
    publishDir "${params.output_dir}/report", mode: 'copy'

    input:
    path(sample_sheet)
    val(output_dir)
    val(ready_signal)
    path(count_matrix)

    output:
    path("sample_corr_heatmap.png"), emit: corrplot
    path("sample_pca_dotplot.png"), emit: pcaplot
    path("all_sample_raw_count.csv"), emit: counts
    path("all_sample_fpkm_matrix.csv"), emit: fpkm

    when:
    task.ext.when == null || task.ext.when

    script:

    """
        module load samtools
        echo -e "Sample_ID\\tRaw_reads\\tClean_reads\\tMapping_ratio\\tMapped_reads" > read_peak.num.summary
            
        sed "1d" ${sample_sheet} | while IFS=',' read ID fq1 _ _ _ || [ -n "\$ID" ]; do
            [ -z "\$ID" ] && continue
            
            raw_num=\$(cat ${output_dir}/trimm/\$fq1"_trimming_report.txt" | grep "Total reads processed:" | sed 's/ //g' | cut -d ":" -f 2 || echo "NA")
            clean_num=\$(cat ${output_dir}/trimm/\$fq1"_trimming_report.txt" | grep "Reads written (passing filters):" | sed 's/(/:/g' |sed 's/ //g' | cut -d ":" -f 3 || echo "NA")

            mapping_ratio=\$(cat ${output_dir}/alignment/\$ID".bowtie2.log" | grep "overall alignment rate" | sed 's/overall alignment rate//g' || echo "NA")
            
            unique_mapping_reads=\$(samtools view ${output_dir}/alignment/\$ID"_Q20_sorted.bam" |cut -f 1|sort| uniq|wc -l)

            printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\$ID" "\$raw_num"  "\$clean_num" "\$mapping_ratio" "\$unique_mapping_reads"
        done >> read_peak.num.summary

        mkdir -p ${params.output_dir}/report/
        cp read_peak.num.summary ${params.output_dir}/report/
        
        Rscript ${projectDir}/bin/integrate_N_clustering.r --count_matrix ${count_matrix} --outdir ${params.output_dir}

        WORK_DIR=\$(pwd)

        Rscript -e "rmarkdown::render(input = '${projectDir}/bin/report.Rmd', 
                                        output_file='report.html', 
                                        params=list(summary_table='${params.output_dir}/report/read_peak.num.summary',
                                                    corr_plot = file.path('\$WORK_DIR', 'sample_corr_heatmap.png'),
                                                    pca_plot = file.path('\$WORK_DIR', 'sample_pca_dotplot.png')))"

        mv ${projectDir}/bin/report.html ${params.output_dir}/report/

    
    """

}