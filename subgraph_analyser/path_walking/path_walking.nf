// Nextflow script for processing input graphs, aligning them, extracting subgraphs, and walking paths.

process DATABASE_CHECK {
    // Check if files exist at the specified paths for kraken_db, genes, and bakta_db, else throw an error.
    input:
    path kraken_db
    path genes
    path bakta_db

    script:
    """
    if [ ! -d "$kraken_db" ]; then
        echo "Error: kraken2 database path $kraken_db does not exist."
        exit 1
    fi
    if [ ! -f "$genes" ]; then
        echo "Error: genes fasta file path $genes does not exist."
        exit 1
    fi
    if [ ! -d "$bakta_db" ]; then
        echo "Error: bakta database path $bakta_db does not exist."
        exit 1
    fi
    """

}

process GRAPHALIGNER {
    conda "${projectDir}/path_walking/graphaligner.yml"

    tag "${graphID}"

    publishDir "$params.outdir/", mode: 'copy', pattern : '*GA_locus.tsv'

    input:
    tuple val(graphID), path(graph)
    path genes

    output:
    tuple val(graphID), path('*.tsv'), emit: align

    script:
    """
    GraphAligner -g $graph -f $genes -a ${graphID}_GA_.gaf -x dbg
    python3 ${projectDir}/path_walking/graphaligner_hit_locus.py ${graphID}_GA_.gaf ${graphID}_GA_locus.tsv
    """
    stub:
    """
    touch ${graphID}_GA_locus.tsv
    """

}

process GFAKRAKEN2 {
    conda "${projectDir}/path_walking/path_walking.yml"

    tag "${graphID}"

    publishDir "$params.outdir/", mode: 'copy', pattern : '*kraken_out.txt'

    input:
    tuple val(graphID), path(graph)
    path kraken_db

    output:
    tuple val(graphID), path('*kraken_out.txt'), emit: kraken_out

    script:
    """
    python3 ${projectDir}/path_walking/gfaKraken.py $graph $kraken_db
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

    // publishDir "$params.outdir/${graphID}/${subgraphs.baseName}/", mode: 'copy', pattern : '*congruent_paths.csv'
    // publishDir "$params.outdir/${graphID}/${subgraphs.baseName}/", mode: 'copy', pattern : '*metadata.tsv'
    // publishDir "$params.outdir/${graphID}/${subgraphs.baseName}/", mode: 'symlink', pattern : '*.fasta'

    tag "${graphID}"

    input:
    //input joined tuple of subgraph, graphaligner output, and kraken output
    tuple val(graphID), path(subgraphs), path(align), path(kraken_out)

    val(overlap)

    output:
    tuple val(graphID), path('*congruent_paths.csv'), emit: paths
    path('*metadata.tsv'), emit: metadata
    tuple val(subgraphs.baseName), path('*.fasta'), emit: fasta

    script:
    """
    #! /bin/bash
    # extract the specific target path from the graphaligner output corresponding to the subgraph
    target_path=\$(grep ${subgraphs.baseName} $align | cut -f 6)
    # run the path walk script
    python3 ${projectDir}/path_walking/path_walk.py $subgraphs \$target_path $kraken_out $overlap $align
    """
    stub:
    """
    touch ${subgraphs.baseName}_congruent_paths.csv
    touch ${subgraphs.baseName}_metadata.tsv
    touch ${subgraphs.baseName}.fasta
    """
}