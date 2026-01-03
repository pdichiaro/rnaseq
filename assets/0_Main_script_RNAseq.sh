#!/bin/bash
CONDA_BASE_PATH=$(conda info --base)
source $CONDA_BASE_PATH/etc/profile.d/conda.sh

#####################################
# The main script will run on the frontend (it is not consuming many resources).
# All processes will run in parellel.
####################################


# paths
NF_FOLDER=/mnt/ngs_ricerca/NEXTFLOW/
UsefulData=/mnt/ngs_ricerca/Software/
work_dir=/mnt/ngs_ricerca/NEXTFLOW/nextflow_temp/
conf_file=/mnt/ngs_ricerca/NEXTFLOW/local.config

######## PARAMS #########
sample_file=/mnt/ngs_ricerca/dichiarop/Master_scripts/202505_Mesothelioma_OMIM_study_2/samplesheet.csv
outdir=/mnt/ngs_ricerca/dichiarop/Data/202505_Mesothelioma_OMIM_study_2/
star_index=$UsefulData/reference_genome/STAR_index/star.v2.7.11b_gencode.v47_hg38/
index=$UsefulData/reference_transcriptome/kallisto_index/kallisto.v0.48_gencode.v47_hg38_k31/kallisto  
GTF=$UsefulData/reference_genome/gencode.v47.basic.annotation.gtf    #gencode.v47.primary_assembly.annotation.gtf   gencode.v47.basic.annotation.gtf
bed_file=$UsefulData/reference_genome/gencode.v47.bed
genome_fasta=$UsefulData/reference_genome/GRCh38.primary_assembly.genome.fa
transcriptome_fasta=$UsefulData/reference_transcriptome/Gencode_v47_GRCh38_p14/gencode.v47.transcripts.fa
########################

# Where will be stored work files, it can be deleted once the script is finished
mkdir -p $work_dir

project_name=BulkRNAseq_test
work_dir_project=$work_dir/$project_name/

mkdir -p $work_dir_project
cd $work_dir_project


###--- Run the script ---###
### NOTE ###
#--email pierluigi.dichiaro@ausl.re.it \   #if MTA (Mail Transfer Agent) is installed
#--genome GRCh38  \  #GRCh38 --> NCBI  #hg38 --> UCSC   ## if do not you have genome .fasta (not raccomended)
#--extra_star_align_args '--outFilterMultimapNmax 20' \
conda activate nextflow

export NXF_ASSETS=/mnt/ngs_ricerca/NEXTFLOW/Nextflow_pipeline

NXF_VER=25.04.7 nextflow run pdichiaro/rnaseq -r main \
            --publish_dir $outdir \
            --input $sample_file \
            --outdir $outdir \
            --fasta $genome_fasta \
            --transcript_fasta $transcriptome_fasta \
            --star_index $star_index \
            --kallisto_index $index \
            --gene_bed $bed_file \
            --gtf $GTF \
            --gencode True \
            --gtf_extra_attributes 'gene_name' \
            --trimmer trimgalore \
            --extra_trimgalore_args '--quality 20 --stringency 3 --length 20' \
            --min_trimmed_reads 0 \
            --aligner star \
            --pseudo_aligner kallisto \
            --min_mapped_reads 0 \
            --extra_kallisto_quant_args '-b 50 --single-overhang' \
            --pseudo_aligner_kmer_size 31 \
            --with_umi False \
            --remove_ribo_rna False \
            --quantification genome \
            --normalization_method 'all_genes, invariant_genes' \
            --deseq2_vst True \
            --skip_linting True \
            --skip_gtf_filter True \
            --skip_fastqc False \
            --skip_rseqc False \
            --skip_trimming False \
            --skip_alignment False \
            --skip_pseudo_alignment True \
            --skip_markduplicates True \
            --skip_quantification_method False \
            --skip_bigwig False \
            --skip_deseq2_qc False \
            --skip_deeptools_norm False \
            --save_merged_fastq False \
            --save_trimmed False \
            --save_reference False \
            --save_align_intermeds False \
            --multiqc_title $project_name \
            -profile singularity \
            -c $conf_file \
            -process.echo \
            -resume
        
if [ $? -eq 0 ]; then
    echo "Nextflow finished successfully"
    scp -r {./report*,./timeline*,./trace*,./flowchart*} $outdir
    rm -rf $work_dir_project
else
    echo "Nextflow encountered an error"
fi

conda deactivate
