# pre-mRNA_SecStruct
Modeling sequence and structural features for yeast pre-mRNA

# Overview
This repository contains:
* **Database** with gene annotations, RNA-seq RPKMs, intron sets in S. cerevisiae
    * Intron sets: standard introns, controls, decoys, proto introns
* **Feature calculation**: 
    * transcription levels
    * stop-codon features
    * sequence-based features: splicing motifs PWMs, length distributions
    * secondary structure features: protection, stem-formation
    * 3D motifs
* Basic **utilities**: 
    * Conversions between Refseq / Ensembl / symbolic gene names
    * Secondary structure and ensemble generation
* **Cache** expensive-to-compute secondary structure ensembles and features
* **Analyses**: 
    * Statistical comparisons between intron and control / decoy sets
    * Identifying decoy introns from whole transcriptome
    * Classification of introns vs decoys


# Requirements
* bedtools
* mygene
* RNAstructure, Vienna, Contrafold
