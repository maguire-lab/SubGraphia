#!/usr/bin/env nextflow
// main script to exceute the subgraph analysis pipeline using the modules created for each tool. 

// establish input channel for graphs 
Channel
    .fromPath(params.graph)
    .map { tuple(it.baseName, it) }
    .set { graph }
// establish input channel for reads
Channel
    .fromFilePairs(params.reads)
    .map { tuple(it[0].replaceAll(/_R1$/,''), it[1]) }
    .set { reads }

// Import processes
include {DATABASE_CHECK; GRAPHALIGNER; GFAKRAKEN2; SUBGRAPH_EXTRACT; PATH_WALK} from "${projectDir}/path_walking/path_walking.nf"
include {MINIMAP_REDUNDANCY_REMOVER; BWAMEM2; BAM_FILTER; PATH_READ_ALN_FILTER; LR_MINIMAP2} from "${projectDir}/path_filtering/path_filtering.nf"
include {BAKTA} from "${projectDir}/annotation/annotation.nf"

workflow {
    // check that databases have been provided
    if (params.kraken_db == "specify_path") {
        error "ERROR: kraken2 database path must be specified in nextflow.config, expect large disk space requirements for standard refseq bacteria db"
    }
    if (params.genes == "specify_path") {
        error "ERROR: gene database path must be specified in nextflow.config, default is CARD protein homolog fasta"
    }
    if (params.bakta_db == "specify_path") {
        error "ERROR: bakta database path must be specified in nextflow.config"
    }
    // check that reads have been provided, R1 and R2 fastq files
    if (params.reads == "specify_path") {
        error "ERROR: reads must be specified in nextflow.config as a glob pattern e.g. '/path/to/reads/*_R{1,2}.fastq'"
    }
    // run the subgraph extraction process if a graph is supplied
    if (params.graph != "specify_path") {
        DATABASE_CHECK(params.kraken_db, params.genes, params.bakta_db)
        GRAPHALIGNER(graph, params.genes, DATABASE_CHECK.out.db_check)
        GFAKRAKEN2(graph, params.kraken_db, DATABASE_CHECK.out.db_check)
        SUBGRAPH_EXTRACT(graph.join(GRAPHALIGNER.out.align), params.radius)
        subgraph_output = SUBGRAPH_EXTRACT.out.subgraphs.join(GRAPHALIGNER.out.align).join(GFAKRAKEN2.out.kraken_out).transpose()
        PATH_WALK(subgraph_output)
        MINIMAP_REDUNDANCY_REMOVER(PATH_WALK.out.fasta)
        if (files(params.reads).size() > 1) {
            BWAMEM2(MINIMAP_REDUNDANCY_REMOVER.out.representative_fasta.collect(), reads)
            BAM_FILTER(BWAMEM2.out.bam)
            PATH_READ_ALN_FILTER(BAM_FILTER.out.bam_summary, PATH_WALK.out.metadata.collect(), MINIMAP_REDUNDANCY_REMOVER.out.representative_fasta.collect(), params.strictness)
            // BAKTA(PATH_READ_ALN_FILTER.out.filtered_fasta.transpose(), params.bakta_db)
        }
        else if (files(params.reads).size() == 1) {
            LR_MINIMAP2(MINIMAP_REDUNDANCY_REMOVER.out.representative_fasta.collect(), reads)
        }
                
    }
    else {
        error "Usage: nextflow run main.nf --graph <path_to_graph> --reads <path_to_reads>"
    }
}