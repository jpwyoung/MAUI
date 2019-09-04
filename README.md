Using PCR amplicons to describe microbial diversity

MAUI-seq: Metabarcoding using amplicons with unique molecular identifiers to improve error correction


This repository includes scripts for the analysis of amplicon sequences that incorporate a random sequence tag (seqid or Unique Molecular Identifier). An early version of the method was described in the following preprint: https://www.biorxiv.org/content/10.1101/538587v2.abstract, but the current version is substantially different.

MAUIsortgenes is a simple script that sorts multiplexed fastq sequences into gene-specific files based on a short motif in the forward primer. It was written for our own data but is easily modified.

MAUIcount analyses amplicon data for one gene across multiple samples, producing a table of allele frequencies and a fasta file of allele sequences.

MC_parameters is a user-edited file that specifies the data source and parameters for MAUIcount.