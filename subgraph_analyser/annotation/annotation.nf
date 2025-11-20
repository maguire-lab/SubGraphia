process BAKTA {
    conda "${projectDir}/annotation/annotation.yml"

    tag "${readID}"

    publishDir "$params.outdir/annotations/", mode: 'symlink', pattern : '*.gbff'

    input:
    tuple val(readID), path(filtered_fasta)
    path bakta_db

    output:
    tuple val(readID), path('*.gbff'), emit: bakta_output

    script:
    """
    bakta --db $bakta_db --skip-plot $filtered_fasta
    """
}