
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
    path chocophlan_db
    path mapping_db
    path metaphlan_db

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
        index=$(basename "$metaphlan_db")
        humann_config --update database_folders nucleotide $chocophlan_db
        humann_config --update database_folders protein $uniref_db
        humann_config --update database_folders utility_mapping $mapping_db

        humann \\
            --input $read1 \\
            --input $read2 \\
            --input-format fastq.gz \\
            --output ${prefix} \\
            --output-basename ${prefix} \\
            --search-mode uniref90 \\
            --pathways metacyc \\
            --metaphlan-options "--bowtie2db $metaphlan_db --index $index" \\
            --threads $task.cpus


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            : \$(echo \$(humann --version 2>&1) | sed 's/^.*humann //; s/Using.*\$//' ))
        END_VERSIONS
        """

    } else {

        """
        index=$(basename "$metaphlan_db")
        humann_config --update database_folders nucleotide $chocophlan_db
        humann_config --update database_folders protein $uniref_db
        humann_config --update database_folders utility_mapping $mapping_db

        humann \\
            --input $reads \\
            --input-format fastq.gz \\
            --output ${prefix} \\
            --output-basename ${prefix} \\
            --search-mode uniref90 \\
            --pathways metacyc \\
            --metaphlan-options "--bowtie2db $metaphlan_db --index $index" \\
            --threads $task.cpus


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            : \$(echo \$(humann --version 2>&1) | sed 's/^.*humann //; s/Using.*\$//' ))
        END_VERSIONS
        """
    }

}
