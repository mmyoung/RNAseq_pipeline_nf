process FEATURECOUNTS {
    label 'featurecounts'
    publishDir "${params.output_dir}/featurecounts", mode: 'copy'

    input:
    path(bam_files)
    path(annotation)
    val(strandness)
    val(pairend)

    output:
    path("*featureCounts.tsv"), emit: counts
    path("*featureCounts.tsv.summary"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def strand_param = strandness ? "-s ${strandness}" : ""
    def bam_list = bam_files.join(' ')
    def pair_end = pairend==0 ? '' : '-p'

    """
    /project/gzy8899/softwares/subread-2.1.1-Linux-x86_64/bin/featureCounts \\
        ${pair_end} \\
        -T ${task.cpus} \\
        -a ${annotation} \\
        ${strand_param} \\
        -o all.featureCounts.tsv \\
        ${bam_list}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """

}