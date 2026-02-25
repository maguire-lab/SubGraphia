# Datasets used for SubGraphia evaluation 
This document describes the datasets used to evaluate SubGraphia, simulated and mock metagenomic datasets are needed to have a ground truth for comparison.
The datasets aim to represent real world metagenomic samples containing antimicrobial resistance (AMR) genes in a variety of complexities.
Real examples or artifical LGTs between reference genomes used to create the datasets. These will include integrated MGEs and plasmids for each metagenome.

## Simulated LGTs:
    1) Klebsiella pneumoniae to Escherichia coli - in silico 
        - "Simplest" data, VIM-4 bearing integron artificially put into both genomes. 
        - integrated
    2) Streptococcus Anginosus to Enterococcus faecalis - experimental PB
        - integrated
        - Integrative and conjugative element (ICE) with ErmB gene
        - MGE: PP359593.1
        - Donor: GCF_036409205.1
        - Recipient: GCF_036409255.1
    3) Acinetobacter johnsonii to E. coli 25DN
        - Plasmid transfer of a BlaOXA on a transposon 
        - MGE: 	NZ_CP037425.1
        - Donor: GCF_004337595.1
        - Recipient: GCF_015290725.1 + NZ_CP037425.1
    4) Streptococcus dysgalactiae to streptococcus suis bacteriophage 
        - integrated 
        - ermB, optrA, aph(3‚Äô)-IIIa, aac(6‚Äô)-Ie-aph(2‚Äù)-Ia
        - MGE: GCF_019856435.1 from  623692 - 679086, insertion site = SAG0725
        - Donor: GCF_019856435.1
        - recipient: CP082200.1
    5) Salmonella enterica to Citrobacter freundii transposon
        - Tn1721, both in chromosome
        - MGE: X61367.1
        - Donor: GCF_030517175.1
        - Recipient: GCF_036289245.1
    6) Pseudomonas aeruginosa to Pseusomonas aeruginosa ICE
        - integrated
        - MGE: GCF_016745135.1 position 5210218-5310913
        - Donor: GCF_016745135.1
        - Recipient: GCA_028608505.1

    7) Raoultella planticola to Enterobacter hormaechei plasmid
        - plasmid
        - MGE: NZ_CP162951.1
        - Donor: GCF_041002795.1
        - Recipient: GCF_041004125.1

    8) L. monocytogenes to L. welshimeri transformation of fosX
        - https://www.nature.com/articles/s41467-024-54459-9#MOESM3
        - integrated 
        - MGE: fosX gene
        - Donor: GCA_014223565.1 
        - Recipient: GCF_008270765.1

    9) Enterococcus faecalis to Staphylococcus aureus plasmid
        - plasmid: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.894241/full#supplementary-material
        - MGE: NZ_CP084887.1
        - Donor: GCF_000411655.2
        - Recipient: GCF_000013425.1 + NZ_CP084887.1

    10) Actinobacillus pleuropneumoniae to Vibrio cholerae ICE - sxt element
        - integrated
        - MGE: 	KT151660.1
        - Donor: GCF_023100885.1
        - Recipient: GCA_050154015.1

    11) Clostridium innocuum to Clodtridioides difficile (Tn6189)
        - ErmB
        - MGE: 	MK895712.1
        - Donor: GCA_022845735.1
        - Recipient: GCA_047551115.1

    12) Bacteroides hominis to Bacteroides fragilis transposon
        - integrated, Tn7528, CfxA2
        - MGE: OP204844.1
        - Donor: GCF_033042715.2
        - Recipient: GCA_001816225.4


## Human gut microbiome
- Well studied microbiome
- research focused sample
- Deeply sequenced but with lots of human 
- Will simulate 30M short reads which should cover even lowest abundance species (albeit at low depth)
- Will simulate 4M long reads, should cover everything, even at low depth
- Based off this paper: https://www.science.org/doi/10.1126/science.abm1483#sec-1
    - Single cell sequencing of human gut microbiome samples
    - 76 species with abundance inference and accessions (Table S3)
    - A few unidentified genomes I can plug in for LGTs
- Simulation command SR: ```iss generate --draft ../assemblies/gut_ss/*.fna --abundance_file gut_ss_iss_abundance_infile.tsv --model miseq --output miseq_reads --cpus 8 --n_reads 30M```
    - iss version 1.6.0
- fastp v1.1.0 command ```fastp -i miseq_reads_R1.fastq -I miseq_reads_R2.fastq -o gut_ss_fastpOut_reads_R1.fastq -O gut_ss_fastpOut_reads_R2.fastq```
- metaspades v3.15.5 command ```spades.py -1 gut_ss_fastpOut_reads_R1.fastq -2 gut_ss_fastpOut_reads_R2.fastq --meta -m 300 -o gut_ss_metaspades_out```
- Simulation command LR: ```simulator.py metagenome -gl gut_ss_nanosim_assembly_paths_infile.tsv -a gut_ss_nanosim_abundance_infile.tsv -dl gut_ss_nanosim_dna_type_infile.tsv -c ../nano_sim_model/metagenome_ERR3152364_Even_v3.2.2/training --chimeric -t 8 --fastq```
    - nanosim version 3.2.2
- filtlong v0.3.1 command ```filtlong --min_length 1000 --keep_percent 90 simulated_sample0_aligned_reads.fastq```
- rename reads ```seqkit replace -p ".+" -r "read_{nr}" input.fastq.gz -o output.fastq.gz ```
- save old read names ```sed -n '1~4p' input.fastq | cut -c 2- > read_names.txt```
- flye v2.9.6 command ``` flye --nano-corr gut_ss_filtlong_out.fastq --meta --out-dir gut_ss_flye_out --threads 16```


### Included LGTs:
    3) Acinetobacter johnsonii to E. coli 25DN
        - Plasmid transfer of a BlaOXA on a transposon 
        - MGE: 	NZ_CP037425.1
        - Donor: GCF_004337595.1
        - Recipient: GCF_015290725.1 + NZ_CP037425.1

    8) L. welshimeri to L. monocytogenes transformation of fosX
        - https://www.nature.com/articles/s41467-024-54459-9#MOESM3
        - integrated 
        - MGE: fosX gene
        - Donor: GCA_014223565.1 
        - Recipient: GCF_008270765.1

    9) Enterococcus faecalis to Staphylococcus aureus plasmid
        - plasmid: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.894241/full#supplementary-material
        - MGE: NZ_CP084887.1
        - Donor: GCF_000411655.2
        - Recipient: GCF_000013425.1 + NZ_CP084887.1

    11) Clostridium innocuum to Clodtridioides difficile (Tn6189)
        - ErmB
        - MGE: 	MK895712.1
        - Donor: GCA_022845735.1
        - Recipient: GCA_047551115.1

    12) Bacteroides hominis to Bacteroides fragilis transposon
        - integrated, Tn7528, CfxA2
        - MGE: OP204844.1
        - Donor: GCF_033042715.2
        - Recipient: GCA_001816225.4

## Simulated wastewater metagenome 
- Classic One Health sample, high complexity. 
- simulating depth on the lower side, to represent real world sequencing effort keeping a recurrent sequence cost in mind.
- 26M short reads, ~$350 at the IMR
- 3.5M long reads, not clear on cost. $2-3K
- [ ] Need an example metagenome to base off of
    - Abundances
    - reference genomes
- Global wastewater metagenome paper: https://www.nature.com/articles/s41467-025-66070-7#Sec7
    - Halifax sample: ERS22875544
    - run: ERR14173773
    - Number of reads: 85,305,236
    - Number of human reads: 283,911
- Simulation command: ```iss generate --draft ../assemblies/wastewater/*.fna --abundance_file ww_iss_abundance_infile.tsv --model miseq --output miseq_reads --cpus 8 --n_reads 26M```
    - iss version 1.6.0
- fastp v1.1.0 command ```fastp -i miseq_reads_R1.fastq -I miseq_reads_R2.fastq -o ww_fastpOut_reads_R1.fastq -O ww_fastpOut_reads_R2.fastq```
- metaspades v3.15.5 command ```spades.py -1 ww_fastpOut_reads_R1.fastq -2 ww_fastpOut_reads_R2.fastq --meta -m 300 -o ww_metaspades_out```
- Simulation command LR: ```simulator.py metagenome -gl ww_nanosim_assembly_paths_infile.tsv -a ww_nanosim_abundance_infile.tsv -dl ww_nanosim_dna_type_infile.tsv -c ../nano_sim_model/metagenome_ERR3152364_Even_v3.2.2/training --chimeric -t 12 --fastq```
    - nanosim version 3.2.2
- filtlong v0.3.1 command ```filtlong --min_length 1000 --keep_percent 90 simulated_sample0_aligned_reads.fastq```
- rename reads ```seqkit replace -p ".+" -r "read_{nr}" input.fastq.gz -o output.fastq.gz ```
- save old read names ```sed -n '1~4p' input.fastq | cut -c 2- > read_names.txt```
- flye v2.9.6 command ``` flye --nano-corr ww_filtlong_out.fastq --meta --out-dir ww_flye_out --threads 16```

### Included LGTs:
    1) Klebsiella pneumoniae to Escherichia coli - in silico 
        - "Simplest" data, VIM-4 bearing integron artificially put into both genomes. 
        - integrated
    2) Streptococcus Anginosus to Enterococcus faecalis - experimental PB
        - integrated
        - Integrative and conjugative element (ICE) with ErmB gene
        - MGE: PP359593.1
        - Donor: GCF_036409205.1
        - Recipient: GCF_036409255.1
    3) Acinetobacter johnsonii to E. coli 25DN
        - Plasmid transfer of a BlaOXA on a transposon 
        - MGE: 	NZ_CP037425.1
        - Donor: GCF_004337595.1
        - Recipient: GCF_015290725.1 + NZ_CP037425.1
    4) Streptococcus dysgalactiae to streptococcus suis bacteriophage 
        - integrated 
        - ermB, optrA, aph(3‚Äô)-IIIa, aac(6‚Äô)-Ie-aph(2‚Äù)-Ia
        - MGE: GCF_019856435.1 from  623692 - 679086, insertion site = SAG0725
        - Donor: GCF_019856435.1
        - recipient: GCF_019856415.1
    5) Salmonella enterica to Citrobacter freundii transposon
        - Tn1721, both in chromosome
        - MGE: X61367.1
        - Donor: GCF_030517175.1
        - Recipient: GCF_036289245.1
    6) Pseudomonas aeruginosa to Pseusomonas aeruginosa ICE
        - integrated
        - MGE: GCF_016745135.1 position 5210218-5310913
        - Donor: GCF_016745135.1
        - Recipient: GCA_028608505.1

    7) Raoultella planticola to Enterobacter hormaechei plasmid
        - plasmid
        - MGE: NZ_CP162951.1 in donor, NZ_CP163027.1 in recipient
        - Donor: GCF_041002795.1
        - Recipient: GCF_041004125.1

    8) L. welshimeri to L. monocytogenes transformation of fosX
        - https://www.nature.com/articles/s41467-024-54459-9#MOESM3
        - integrated 
        - MGE: fosX gene
        - Donor: GCA_014223565.1 
        - Recipient: GCF_008270765.1

    9) Enterococcus faecalis to Staphylococcus aureus plasmid
        - plasmid: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.894241/full#supplementary-material
        - MGE: NZ_CP084887.1
        - Donor: GCF_000411655.2
        - Recipient: GCF_000013425.1 + NZ_CP084887.1

    10) Actinobacillus pleuropneumoniae to Vibrio cholerae ICE - sxt element
        - integrated
        - MGE: 	KT151660.1
        - Donor: GCF_023100885.1
        - Recipient: GCA_050154015.1

## Simulated wastewater complete references only
- 26M short reads
- exponential distribution
- Simulation command: ```iss generate --genomes ../assemblies/complete_ww/*.fna --abundance exponential --model miseq --output miseq_reads --cpus 4 --n_reads 26M```
    - iss version 1.6.0



## Simulated manure pit metagenome
- High complexity, animal health associated sample.
- Simulating higher depth than wastewater to represent a more targeted sequencing effort.
- 90M reads, similar to WW metagenome depth above, ~$1000 at the IMR
- 12M long reads, three gridion flow cells
- expecting some strain level diversity and genus level diversity.
- [ ] Need an example metagenome to base off of
    - Abundances
    - reference genomes
- Manure pit AMR metagenome paper with many samples at various depths: https://doi.org/10.1128/mbio.00798-21
    - Would need to run kraken2/bracken then find reference genomes https://nf-co.re/taxprofiler/1.2.5/
    - 6" depth manure pit: SAMN16546522
        - chosen because lots of AMR genes in the paper in this sample
    - SRA: SRS7598713
    - Run: SRR12915125
    - Number of reads: 15,993,250
    - Number of human reads: 22,338 * ideally should have bos taurus in the database
- Simulation command: ```iss generate --draft ../assemblies/manure_pit/*.fna --abundance_file mp_iss_abundance_infile.tsv --model miseq --output miseq_reads --cpus 8 --n_reads 90M```
    - iss version 1.6.0
- fastp v1.1.0 command ```fastp -i miseq_reads_R1.fastq -I miseq_reads_R2.fastq -o mp_fastpOut_reads_R1.fastq -O mp_fastpOut_reads_R2.fastq```
- metaspades v3.15.5 command ```spades.py -1 mp_fastpOut_reads_R1.fastq -2 mp_fastpOut_reads_R2.fastq --meta -m 300 -o mp_metaspades_out```
    - Failed, tried to allocate 419Gb
- Simulation command LR: ```simulator.py metagenome -gl mp_nanosim_assembly_paths_infile.tsv -a mp_nanosim_abundance_infile.tsv -dl mp_nanosim_dna_type_infile.tsv -c ../nano_sim_model/metagenome_ERR3152364_Even_v3.2.2/training --chimeric -t 12 --fastq```
    - nanosim version 3.2.2
- filtlong v0.3.1 command ```filtlong --min_length 1000 --keep_percent 90 simulated_sample0_aligned_reads.fastq```
- rename reads ```seqkit replace -p ".+" -r "read_{nr}" input.fastq.gz -o output.fastq.gz ```
- save old read names ```sed -n '1~4p' input.fastq | cut -c 2- > read_names.txt```
- flye v2.9.6 command ``` flye --nano-corr mp_filtlong_out.fastq --meta --out-dir mp_flye_out --threads 16```

### Included LGTs:
    1) Klebsiella pneumoniae to Escherichia coli - in silico 
        - "Simplest" data, VIM-4 bearing integron artificially put into both genomes. 
        - integrated
    2) Streptococcus Anginosus to Enterococcus faecalis - experimental PB
        - integrated
        - Integrative and conjugative element (ICE) with ErmB gene
        - MGE: PP359593.1
        - Donor: GCF_036409205.1
        - Recipient: GCF_036409255.1
    3) Acinetobacter johnsonii to E. coli 25DN
        - Plasmid transfer of a BlaOXA on a transposon 
        - MGE: 	NZ_CP037425.1
        - Donor: GCF_004337595.1
        - Recipient: GCF_015290725.1 + NZ_CP037425.1
    4) Streptococcus dysgalactiae to streptococcus suis bacteriophage 
        - integrated 
        - ermB, optrA, aph(3‚Äô)-IIIa, aac(6‚Äô)-Ie-aph(2‚Äù)-Ia
        - MGE: GCF_019856435.1 from  623692 - 679086, insertion site = SAG0725
        - Donor: GCF_019856435.1
        - recipient: GCF_019856415.1
    5) Salmonella enterica to Citrobacter freundii transposon
        - Tn1721, both in chromosome
        - MGE: X61367.1
        - Donor: GCF_030517175.1
        - Recipient: GCF_036289245.1
    6) Pseudomonas aeruginosa to Pseusomonas aeruginosa ICE
        - integrated
        - MGE: GCF_016745135.1 position 5210218-5310913
        - Donor: GCF_016745135.1
        - Recipient: GCA_028608505.1

    7) Raoultella planticola to Enterobacter hormaechei plasmid
        - plasmid
        - MGE: NZ_CP162951.1
        - Donor: GCF_041002795.1
        - Recipient: GCF_041004125.1

    8) L. welshimeri to L. monocytogenes transformation of fosX
        - https://www.nature.com/articles/s41467-024-54459-9#MOESM3
        - integrated 
        - MGE: fosX gene
        - Donor: GCA_014223565.1 
        - Recipient: GCF_008270765.1

    9) Enterococcus faecalis to Staphylococcus aureus plasmid
        - plasmid: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.894241/full#supplementary-material
        - MGE: NZ_CP084887.1
        - Donor: GCF_000411655.2
        - Recipient: GCF_000013425.1 + NZ_CP084887.1

    10) Actinobacillus pleuropneumoniae to Vibrio cholerae ICE - sxt element
        - integrated
        - MGE: 	KT151660.1
        - Donor: GCF_023100885.1
        - Recipient: GCA_050154015.1

    11) Clostridium innocuum to Clodtridioides difficile (Tn6189)
        - ErmB
        - MGE: 	MK895712.1
        - Donor: GCA_022845735.1
        - Recipient: GCA_047551115.1

    12) Bacteroides hominis to Bacteroides fragilis transposon
        - integrated, Tn7528, CfxA2
        - MGE: OP204844.1
        - Donor: GCF_033042715.2
        - Recipient: GCA_001816225.4

## [Zymobiomics complex mock community](https://www.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.001469#R28) 
- Simulation independent mock community
- kind of more environmental?
- high depth sequencing
- Created by pooling even amounts (by mass) of genomic DNA from 227 bacterial strains belonging to eight phyla (Actinobacteria, Bacteroidetes, Deinococcus-Thermus, Firmicutes, Fusobacteria, Planctomycetes, Proteobacteria and Verrucomicrobia), 19 classes, 47 orders, 85 families, 175 genera and 197 species
- Covers a wide diversity of taxa not necessarily found together in real world samples
- Available in illumina, promethion, and pacbio
- Reads:
    - Illumina: ERR5321934
    - Promethion: ERR6109821
    - Pacbio: PRJEB45869
- Reference genomes:
    - PRJEB32402
    - read accessions: /data/raid1_HDD/David/lgtfinder/simulations/assemblies/hc227/HCmockAccns.csv
    - Assembly accessions for ENA from ATB: /home/david/dm-lab/listeria_evol/atb/atb.metadata.202408.sqlite.assembly.tsv
    - One assembly rejected from ATB (SAMEA5613686) assembled with spades 3.15.5
        - `spades.py -1 ERR3330973_1.fastq -2 ERR3330973_2.fastq -o SAMEA5613686_out --isolate`
            -  `* 0:00:09.291    55M / 940M  WARN    General                 (kmer_coverage_model.cpp   : 218)   Too many erroneous kmers, the estimates might be unreliable`
    - fastp v1.1.0 command ```fastp -i ERR5321934_1.fastq -I ERR5321934_2.fastq -o hc227_fastpOut_reads_R1.fastq -O hc227_fastpOut_reads_R2.fastq```
    - metaspades v3.15.5 command ```spades.py -1 hc227_fastpOut_reads_R1.fastq -2 hc227_fastpOut_reads_R2.fastq --meta -m 300 -o hc227_metaspades_out```
    - filtlong v0.3.1 command ```filtlong --min_length 1000 --keep_percent 90 ERR6109821.fastq```
    - flye v2.9.6 command ``` flye --nano-corr hc227_filtlong_out.fastq --meta --out-dir hc227_flye_out --threads 16```


# Simulation process
## Kraken and bracken 
- taxprofiler nf-core pipeline see ```/data/raid1_HDD/David/lgtfinder/simulations/base_metagenomes```

## Genus taxID to reference genome
- NCBI datasets to het first refseq hit ```datasets summary genome taxon 286 --reference --assembly-source RefSeq --assembly-level chromosome --report ids_only --limit 1 --as-json-lines```
- get curate accession list to include LGTs
- Download genomes with ncbi datasets

## Limiting diversity to taxa making up 98% of the abundance, the remaining 2% being used for spike in LGTs
- This will limit some complexity from very sparsely sqeuenced genomes