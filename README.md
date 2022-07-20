# pre-mRNA_SecStruct
Computational modeling for sequence and structural features of yeast pre-mRNA

# Overview
This repository contains:
* **Database** with gene annotations, RNA-seq RPKMs, intron sets in S. cerevisiae, alignments across Saccharomyces genus
    * Each intron set (e.g. standard introns, controls, decoys, proto introns) includes:  
         * sequences
         * splice site positions
         * cached secondary strucure features
         * minimum free energy secondary structures
         * stochastically sampled secondary structure ensembles (will be stored here but not uploaded to repo)
* **Feature calculation**: 
    * transcription levels
    * stop-codon features
    * sequence-based features: splicing motifs PWMs, length distributions
    * secondary structure features: stem dGs, longest stem, maximum path length, sequence protection
* Basic **utilities**: 
    * Conversions between Refseq / Ensembl / symbolic gene names
    * Secondary structure and ensemble generation
    * Processing sequence alignments
* **Cache** expensive-to-compute secondary structure ensembles and features
* **Analyses**: 
    * Statistical comparisons between intron and control / decoy sets
    * Analyzing sequence and secondary structure properties across yeast species
    * Identifying decoy introns from whole transcriptome
    * Classification of introns vs decoys

# Generating figures
A description for how to generate figures is included in `src/README.md`.

# Requirements
* bedtools
* mygene
* arnie, RNAstructure, Vienna, Contrafold

# Running
* Copy `src/sample_config.py` to `src/config.py` and update paths.
