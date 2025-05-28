
process JOIN_TABLES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::humann=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.9--py312hdfd78af_0':
        'biocontainers/humann:3.9--py312hdfd78af_0' }"

    input:
    path gene_files
    path abund_files
    path coverage_files

    output:
    path ("all_genefamilies_joined.tsv")
    path ("all_pathabundance_joined.tsv")
    path ("all_pathcoverage_joined.tsv")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p gene_join_dir
    cp ${gene_files.join(' ')} gene_join_dir/
    humann_join_tables \\
        --input gene_join_dir \\
        --output all_genefamilies_joined.tsv

    mkdir -p abund_join_dir
    cp ${abund_files.join(' ')} abund_join_dir/
    humann_join_tables \\
        --input abund_join_dir \\
        --output all_pathabundance_joined.tsv

    mkdir -p coverage_join_dir
    cp ${coverage_files.join(' ')} coverage_join_dir/
    humann_join_tables \\
        --input coverage_join_dir \\
        --output all_pathcoverage_joined.tsv

    """

}
