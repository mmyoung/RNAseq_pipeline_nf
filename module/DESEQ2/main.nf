process DESEQ {
    label "DESeq2"
    publishDir "${params.output_dir}/differential_expression", mode: 'copy'

    input:
    path counts

    output:
    path "*.tsv"
    path "versions.yml"

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    Rscript ${projectDir}/bin/DESeq2.r \\
        --count_file $counts \\
        --col_data ${params.DESeq2_coldata} \\
        --cores $task.cpus 
    cat <<-END_VERSIONS > versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | cut -d' ' -f3)
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}