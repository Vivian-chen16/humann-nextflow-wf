
process HUMANN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::humann=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.9--py312hdfd78af_0':
        'biocontainers/humann:3.9--py312hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path uniref_db
    path chocophlan_db
    path mapping_db
    path metaphlan_db
    val metaphlan_db_index

    output:
    tuple val(meta), path("${meta.id}_genefamilies.tsv")    , emit: gene_family
    tuple val(meta), path("${meta.id}_pathabundance.tsv")   , emit: path_abundance
    tuple val(meta), path("${meta.id}_pathcoverage.tsv")    , emit: path_coverage
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ( meta.type == "PE") {
        read1 = "${reads[0]}"
        read2 = "${reads[1]}"

        """
        cat $read1 $read2 > merged.fastq.gz

        export PYTHONWARNINGS="ignore::SyntaxWarning"
        humann_config --update database_folders nucleotide $chocophlan_db
        humann_config --update database_folders protein $uniref_db
        humann_config --update database_folders utility_mapping $mapping_db

        humann \\
            --input merged.fastq.gz \\
            --input-format fastq.gz \\
            --output ./ \\
            --output-basename ${prefix} \\
            --search-mode uniref90 \\
            --pathways metacyc \\
            --metaphlan-options "--bowtie2db $metaphlan_db --index $metaphlan_db_index" \\
            --threads $task.cpus


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            : \$(echo \$(humann --version 2>&1) | sed 's/^.*humann //; s/Using.*\$//' ))
        END_VERSIONS
        """

    } else {

        """
        export PYTHONWARNINGS="ignore::SyntaxWarning"
        humann_config --update database_folders nucleotide $chocophlan_db
        humann_config --update database_folders protein $uniref_db
        humann_config --update database_folders utility_mapping $mapping_db

        humann \\
            --input $reads \\
            --input-format fastq.gz \\
            --output ./ \\
            --output-basename ${prefix} \\
            --search-mode uniref90 \\
            --pathways metacyc \\
            --metaphlan-options "--bowtie2db $metaphlan_db --index $metaphlan_db_index" \\
            --threads $task.cpus


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            : \$(echo \$(humann --version 2>&1) | sed 's/^.*humann //; s/Using.*\$//' ))
        END_VERSIONS
        """
    }

}
