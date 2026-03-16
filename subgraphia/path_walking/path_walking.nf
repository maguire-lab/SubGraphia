// Nextflow script for processing input graphs, aligning them, extracting subgraphs, and walking paths.
process GRAPHALIGNER {
    conda "${projectDir}/path_walking/graphaligner.yml"
    cpus params.max_cpus/3
    tag "${graphID}"
    
    input:
    tuple val(graphID), path(graph)
    path genes

    output:
    tuple val(graphID), path('*.gaf'), emit: align

    script:
    """
    GraphAligner -g $graph -f $genes -a ${graphID}_GA_.gaf -x dbg -t ${task.cpus}
    """
    stub:
    """
    touch ${graphID}_GA_.gaf
    """

}

process GRAPHALIGNER_OUTPUT_PROCESSING {
    conda "${projectDir}/path_walking/path_walking.yml"
    tag "${graphID}"
    publishDir "$params.outdir/", mode: 'copy', pattern : '*GA_locus.tsv'

    input:
    tuple val(graphID), path(align)
    tuple val(graphID), path(graph)

    output:
    tuple val(graphID), path('*.tsv'), emit: align_processed
    tuple val(graphID), path('*_gene_sequences.fasta'), emit: gene_sequences

    script:
    """
    python3 ${projectDir}/path_walking/graphaligner_hit_locus.py ${graphID}_GA_.gaf ${graphID}_GA_locus.tsv ${graph}
    """
    stub:
    """
    touch ${graphID}_GA_locus.tsv
    touch ${graphID}_gene_sequences.fasta
    """

    }

process GFAKRAKEN2 {
    conda "${projectDir}/path_walking/path_walking.yml"
    cpus params.max_cpus
    tag "${graphID}"

    publishDir "$params.outdir/", mode: 'copy', pattern : '*kraken_out.txt'

    input:
    tuple val(graphID), path(graph)
    path kraken_db

    output:
    tuple val(graphID), path('*kraken_out.txt'), emit: kraken_out

    script:
    """
    python3 ${projectDir}/path_walking/gfaKraken.py $graph $kraken_db ${task.cpus}
    """

    stub:
    """
    touch ${graphID}_kraken_out.txt
    """
}


process SUBGRAPH_EXTRACT {
    conda "${projectDir}/path_walking/path_walking.yml"
    tag "${graphID}"

    input:
    // input joined tuple of graph and graphaligner output
    tuple val(graphID), path(graph), path(align)
    // input the extraction radius
    val(radius)

    output:
    // define an output tuple with the output file basename and the output file
    tuple val(graphID), path('*.gfa'), emit: subgraphs

    script:
    """
    python3 ${projectDir}/path_walking/subgraph_extract.py $graph $align $radius
    """
    stub:
    """
    touch ARO_subgraph.gfa
    """

}

process PATH_WALK {
    conda "${projectDir}/path_walking/path_walking.yml"
    maxForks params.max_cpus
    tag "${graphID}"

    input:
    tuple val(graphID), path(subgraphs), path(align), path(kraken_out)

    output:
    tuple val(graphID), path('*congruent_paths.csv'), emit: paths
    path('*metadata.tsv'), emit: metadata
    tuple val(subgraphs.baseName), path('*.fasta'), emit: fasta

    script:
    """
    #! /bin/bash
    # extract the specific target path from the graphaligner output corresponding to the subgraph
    target_path=\$(grep -w ${subgraphs.baseName} $align | cut -f 6)
    # run the path walk script
    python3 ${projectDir}/path_walking/path_walk.py $subgraphs \$target_path $kraken_out $align
    """
    stub:
    """
    touch ${subgraphs.baseName}_congruent_paths.csv
    touch ${subgraphs.baseName}_metadata.tsv
    touch ${subgraphs.baseName}.fasta
    """
}