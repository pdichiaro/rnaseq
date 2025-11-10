//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT          } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT        } from '../../../modules/nf-core/kallisto/quant'
include { KALLISTO_STATS        } from '../../../modules/local/kallisto_stats'
include { KALLISTO_STATS_MERGE  } from '../../../modules/local/kallisto_stats_merge'
include { CUSTOM_TX2GENE        } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT      } from '../../../modules/nf-core/tximeta/tximport'
include { TXIMPORT_MERGE_COUNTS } from '../../../modules/local/tximport_merge_counts'



workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet               // channel: [ val(meta), /path/to/samplsheet ]
    reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: /path/to//index/
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf                       // channel: /path/to/genome.gtf
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner            //     val: kallisto or salmon
    alignment_mode            //    bool: Run Salmon in alignment mode
    lib_type                  //     val: String to override Salmon library type
    kallisto_quant_fraglen    //     val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd //     val: Estimated standard error for fragment length required by Kallisto in single-end mode
    annotation_matrix         // channel: /path/to/annotation_matrix.txt
    chrom_sizes               // channel: /path/to/genome.sizes

    main:
    ch_versions = Channel.empty()

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        SALMON_QUANT (
            reads,
            index,
            gtf,
            transcript_fasta,
            alignment_mode,
            lib_type
        )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results
        ch_pseudo_bam = Channel.empty()  // Salmon doesn't generate BAM files in pseudoalignment mode
        ch_pseudo_bai = Channel.empty()  // Salmon doesn't generate BAI files in pseudoalignment mode
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
    } else {
        KALLISTO_QUANT (
            reads,
            index,
            gtf,
            chrom_sizes,
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_pseudo_results = KALLISTO_QUANT.out.results
        ch_pseudo_multiqc = KALLISTO_QUANT.out.log
        ch_pseudo_bam = KALLISTO_QUANT.out.bam
        ch_pseudo_bai = KALLISTO_QUANT.out.bai
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
        
        // Parse Kallisto run_info.json for alignment statistics bar plots
        KALLISTO_STATS (
            KALLISTO_QUANT.out.json_info
        )
        
        // Merge per-sample stats into single .txt files (no headers, no MultiQC)
        KALLISTO_STATS_MERGE (
            KALLISTO_STATS.out.multiqc_files.collect(),
            KALLISTO_STATS.out.multiqc_txt_files.collect()
        )
        ch_versions = ch_versions.mix(KALLISTO_STATS.out.versions.first())
        ch_versions = ch_versions.mix(KALLISTO_STATS_MERGE.out.versions)
    }

    CUSTOM_TX2GENE (
        gtf.map { [ [:], it ] },
        ch_pseudo_results.collect{ it[1] }.map { [ [:], it ] },
        pseudo_aligner,
        gtf_id_attribute,
        gtf_extra_attribute
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    TXIMETA_TXIMPORT (
        ch_pseudo_results.collect{ it[1] }.map { [ ['id': 'all_samples'], it ] },
        CUSTOM_TX2GENE.out.tx2gene,
        pseudo_aligner
    )
    ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)

    //
    // Merge tximport counts across samples with annotation
    //
    TXIMPORT_MERGE_COUNTS (
        TXIMETA_TXIMPORT.out.counts_gene,
        TXIMETA_TXIMPORT.out.tpm_gene,
        TXIMETA_TXIMPORT.out.lengths_gene,
        TXIMETA_TXIMPORT.out.counts_transcript,
        TXIMETA_TXIMPORT.out.tpm_transcript,
        TXIMETA_TXIMPORT.out.lengths_transcript,
        annotation_matrix
    )
    ch_versions = ch_versions.mix(TXIMPORT_MERGE_COUNTS.out.versions)



    emit:
    results                       = ch_pseudo_results                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                              // channel: [ val(meta), files_for_multiqc ]
    bam                           = ch_pseudo_bam                                  // channel: [ val(meta), path(*.bam) ] - Kallisto pseudoBAM files
    bai                           = ch_pseudo_bai                                  // channel: [ val(meta), path(*.bam.bai) ] - Kallisto pseudoBAM indices

    // Merged outputs with annotation (for gene-level data)
    tpm_gene                      = TXIMPORT_MERGE_COUNTS.out.tpm_gene             //    path: tximport.merged.gene_tpm.tsv
    counts_gene                   = TXIMPORT_MERGE_COUNTS.out.counts_gene          //    path: tximport.merged.gene_counts.tsv
    lengths_gene                  = TXIMPORT_MERGE_COUNTS.out.lengths_gene         //    path: tximport.merged.gene_lengths.tsv
    
    // Merged outputs without annotation (for transcript-level data)
    tpm_transcript                = TXIMPORT_MERGE_COUNTS.out.tpm_transcript       //    path: tximport.merged.transcript_tpm.tsv
    counts_transcript             = TXIMPORT_MERGE_COUNTS.out.counts_transcript    //    path: tximport.merged.transcript_counts.tsv
    lengths_transcript            = TXIMPORT_MERGE_COUNTS.out.lengths_transcript   //    path: tximport.merged.transcript_lengths.tsv



    versions                      = ch_versions                                    // channel: [ versions.yml ]
}
