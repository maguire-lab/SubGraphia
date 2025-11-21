// Nextflow processes for removing path redundancy and annotating extracted paths

process MINIMAP_REDUNDANCY_REMOVER {
    conda "${projectDir}/path_filtering/path_filtering.yml"

    tag "${subgraphID}"

    // publishDir "$params.outdir/path_filtering/${subgraphID}/", mode: 'symlink', pattern : '*.fasta'

    input:
    tuple val(subgraphID), path(fasta)

    output:
    path('*.fasta'), emit: representative_fasta

    script:
    """
    # if the input $fasta has only one sequence
    if [ \$(grep -c '^>' $fasta) -eq 1 ]; then
        #rename fasta to have same name as its > header
        basename=\$(grep '^>' $fasta | head -n 1 | sed 's/>//')
        mv $fasta "\${basename}.fasta"
    else
        # run minimap2 to find redundant paths
        minimap2 -x asm5 -t 8 $fasta $fasta > tmp.tsv
        cut -f1-11 tmp.tsv > ava_mm_out.tsv
        rm tmp.tsv

        # run redundancy remover script on output 
        python3 ${projectDir}/path_filtering/redundancy_remover.py ./ava_mm_out.tsv 95 75

        #using samtools extract the representative sequences listed in representative_sequences.tsv from $fasta
        while read -r line; do
            samtools faidx $fasta "\$line" > "\$line.fasta"
        done < representative_sequences.tsv
    fi
    """
    stub:
    """
    touch representative_path1.fasta
    touch representative_path2.fasta
    """
}

process BWAMEM2 {
    conda "${projectDir}/path_filtering/path_filtering.yml"

    tag "${readID}"

    publishDir "$params.outdir/", mode: 'copy', pattern : '*.bam'

    input:
    //input only the fasta files from the channel
    path(representative_fasta)
    tuple val(readID), path(reads)

    output:
    tuple val(readID), path('*.bam'), emit: bam

    script:
    """
    # concatenate all fasta files into single fasta file
    cat $representative_fasta > all_paths.fasta

    # index the fasta file
    bwa-mem2 index all_paths.fasta

    # align reads to all paths
    bwa-mem2 mem -a -K 100000000 -t 12 all_paths.fasta ${reads[0]} ${reads[1]} > out.sam

    # convert sam to bam, sort and index
    samtools sort --threads 12 -o all_paths_readaln.bam out.sam
    rm out.sam
    """
    stub:
    """
    touch all_paths_readaln.bam
    """
}

process BAM_FILTER {

    conda "${projectDir}/path_filtering/path_filtering.yml"

    publishDir "$params.outdir/", mode: 'copy', pattern : 'bam_filtering_summary.csv'

    tag "${readID}"

    input:
    tuple val(readID), path(bam)

    output:
    tuple val(readID), path('*_summary.csv'), emit: bam_summary

    script:
    """
    python3 ${projectDir}/path_filtering/bam_filtering.py $bam
    """
    stub:
    """
    touch ${readID}_summary.csv
    """
}

process PATH_READ_ALN_FILTER {

    conda "${projectDir}/path_filtering/path_filtering.yml"

    publishDir "$params.outdir/${readID}/", mode: 'symlink', pattern : '*.fa'
    publishDir "$params.outdir/", mode: 'copy', pattern : '*final_metadata.tsv'
    publishDir "$params.outdir/", mode: 'copy', pattern : 'AMR_genes_summary.tsv'

    tag "${readID}"

    input:
    tuple val(readID), path(bam_summary)
    path(metadata)
    path(representative_fasta)
    val(strictness)

    output:
    tuple val(readID), path('*.fa'), emit: filtered_fasta
    tuple val(readID), path('*final_metadata.tsv'), emit: filtered_metadata
    tuple val(readID), path('AMR_genes_summary.tsv'), emit: amr_summary

    script:
    """
    #concatenate all metadata files into single metadata file
    cat $metadata | grep -v path_id > all_metadata.tsv
    #run the filtering script
    python3 ${projectDir}/path_filtering/final_reconciliation.py all_metadata.tsv $bam_summary $strictness

    #concatenate all fasta files into single fasta file
    cat $representative_fasta > all_rep.fasta

    # extract the IDs of the filtered paths from the final_metadata.tsv, exclude header
    cut -f1 final_metadata.tsv | tail -n +2 > filtered_path_ids.txt

    # extract the filtered paths from all_rep.fasta using samtools
    while read -r line; do
        samtools faidx all_rep.fasta "\$line" > "\$line.fa"
    done < filtered_path_ids.txt
    """
    stub:
    """
    touch filtered_path1.fa
    touch filtered_path2.fa
    touch final_metadata.tsv
    touch AMR_genes_summary.tsv
    """
}