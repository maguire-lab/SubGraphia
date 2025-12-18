// Nextflow script for processing input graphs, aligning them, extracting subgraphs, and walking paths.

process DATABASE_CHECK {
    conda "${projectDir}/path_walking/path_walking.yml"
    // Check if files exist at the specified paths for kraken_db, genes, and bakta_db, else download the databases.
    input:
    val kraken_db
    val genes
    val bakta_db

    output:
    path 'DATABASE_CHECK.done', emit: db_check

    script:
    """
    if [ ! -d "$kraken_db" ]; then
        kraken2-build --standard --threads 24 --db $kraken_db
    fi
    if [ ! -f "$genes" ]; then
        wget https://card.mcmaster.ca/latest/data
        tar -xf data ./nucleotide_fasta_protein_homolog_model.fasta
        mv nucleotide_fasta_protein_homolog_model.fasta $genes
        rm data
    fi
    if [ ! -d "$bakta_db" ]; then
        bakta_db download --output $bakta_db --type light
    fi
    touch DATABASE_CHECK.done
    """
    stub:
    """
    touch DATABASE_CHECK.done
    """
}

process GRAPHALIGNER {
    conda "${projectDir}/path_walking/graphaligner.yml"

    tag "${graphID}"

    publishDir "$params.outdir/", mode: 'copy', pattern : '*GA_locus.tsv'

    input:
    tuple val(graphID), path(graph)
    path genes
    path db_check

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
    path db_check

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

    // publishDir "$params.outdir/${graphID}/subgraphs/", mode: 'copy', pattern : '*.gfa'

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

    maxForks 8

    publishDir "$params.outdir/${graphID}/${subgraphs.baseName}/", mode: 'copy', pattern : '*congruent_paths.csv'
    publishDir "$params.outdir/${graphID}/${subgraphs.baseName}/", mode: 'copy', pattern : '*metadata.tsv'
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
    target_path=\$(grep -w ${subgraphs.baseName} $align | cut -f 6)
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