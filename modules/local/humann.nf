
process HUMANN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::humann=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.9--py312hdfd78af_0':
        'biocontainers/humann:3.9--py312hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path uniref_db

    output:
    tuple val(meta), path("${prefix}_genefamilies.tsv")    , emit: gene_family
    tuple val(meta), path("${prefix}_pathabundance.tsv")   , emit: path_abundance
    tuple val(meta), path("${prefix}_pathcoverage.tsv")    , emit: path_coverage
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ( meta.type == "PE") {
        read1 = "${reads[0]}"
        read2 = "${reads[1]}"

        """
        humann \\
            --input $read1 \\
            --input $read2 \\
            --output ${prefix} \\
            --bypass-nucleotide-search \\
            --protein-database $uniref_db \\
            --threads $task.cpus


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            : \$(echo \$(humann --version 2>&1) | sed 's/^.*humann //; s/Using.*\$//' ))
        END_VERSIONS
        """

    } else {
        """
        humann \\
            --input $reads \\
            --output ${prefix} \\
            --bypass-nucleotide-search \\
            --protein-database $uniref_db \\
            --threads $task.cpus


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            : \$(echo \$(humann --version 2>&1) | sed 's/^.*humann //; s/Using.*\$//' ))
        END_VERSIONS
        """
    }

}
