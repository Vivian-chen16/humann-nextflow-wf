//
// Process reads.
//

include { FASTQC                           } from '../../modules/nf-core/fastqc/main'

workflow PROCESS_READS {
    take:
    samples // [[meta], reads]


    main:
    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()

    ch_PE = samples
        .filter { meta, reads -> meta.type == 'PE' }
        .map { meta, reads -> [
                [ id:meta.id ,
                type:meta.type,
                single_end:false],
            reads ]
        }

    ch_not_PE = samples
        .filter { meta, reads -> meta.type != 'PE' }
        .map { meta, reads -> [
                [ id:meta.id ,
                type:meta.type,
                single_end:true],
            reads ]
        }

    ch_with_nfcore_SE_meta = ch_PE.concat(ch_not_PE)

    FASTQC ( ch_with_nfcore_SE_meta )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    ch_processed_reads = FASTP.out.reads
        .map { meta, reads -> [
                [ id:meta.id ,
                type:meta.type,
                background:meta.background],
            reads ]
        }

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PROCESSED.out.zip.collect{it[1]}.ifEmpty([]))

    emit:
    reads    = ch_processed_reads // channel: [ val(meta), [ coverage_depth ] ]
    versions = ch_versions        // channel: [ versions.yml ]
    mqc      = ch_multiqc_files

}