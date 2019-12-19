
MAUI-seq: Metabarcoding using amplicons with unique molecular identifiers to improve error correction
===

This is a method for using PCR amplicons to describe microbial diversity. A manuscript has been submitted for publication, and a reference will be provided here once it has been published. An early version of the method was described in the following preprint: https://www.biorxiv.org/content/10.1101/538587v2.abstract, but the current version is substantially different.

This repository includes scripts for the analysis of amplicon sequences that incorporate a random sequence tag (seqid or Unique Molecular Identifier).

===

MAUIsortgenes is a simple script that sorts multiplexed fastq sequences into gene-specific files based on a short motif in the forward primer. It was written for our own data but is easily modified.

MAUIcount analyses amplicon data for one gene across multiple samples, producing a table of allele frequencies and a fasta file of allele sequences.

MC_parameters is a user-edited file that specifies the data source and parameters for MAUIcount.

Sample data files for testing are provided in the 'Test' folder.

===

MAUIsortgenes.py 

This script uses BioPython, but is otherwise written in standard Python 3. The input files contain amplicon sequences in fastq format that have been assembled from paired-end reads using PEAR [Zhang, J, Kobert, K, Flouri, T, and Stamatakis, A. (2014). PEAR: a fast and accurate Illumina paired-end read merger. Bioinformatics. 30(5); 614.]. 

It sorts sequences into separate files for each gene, based on a short unambiguous motif (tag) in the forward primer. Sequences that do not match any tag are filed separately as "short" or "unknown". It is very simple and does not allow mismatches. It produces files that are suitable as input for MAUIcount, but users may prefer to use other methods to create these single-gene fastq files.

===

MAUIcount.py

This is the core script for the MAUI method. It is written in standard Python 3 and uses only modules in the Python standard library. It reads amplicon sequences from a set of fastq files.

The first seqid_len bases of each read are a random tag, called a seqid or Unique Molecular Identifier (UMI). The script keeps track of how many times each UMI is used with each unique sequence. For the set of samples, outputs are files with a list of fasta sequences in descending rank order of abundance, and corresponding tables with the counts of each sequence in each sample.

Three output sets are produced by default:

"accepted_sequences": the main MAUI-seq output, using secondary sequences to filter out errors

"all_primary_sequences": the same UMI counts before error filtering

"read_sequences": 'conventional' analysis of the sequences, ignoring UMIs

An additional output set can be produced by setting output_types = 4. This has the same filtering as accepted_sequences, but applied on a per-sample basis, rather than on totals across all samples. If allele frequencies vary greatly across samples, this would be preferable in principle, but it is not reliable unless read counts are very high. Otherwise, it can lead to sequences being stochastically deleted from some samples but not others.

In all cases, the outputs are truncated to discard rare sequences that would have frequencies below add_limit in the overall set of samples. 

If there is a file fastq_file_list.txt in the same folder as the data, only the files listed in this file will be included in the analysis. If this file is not present, it will be created with a list of all files that have the extension .fastq.

===

MC_parameters.py

This file is imported by MAUIcount.py and specifies the files to be analysed and gene-specific parameters. It should be in the same folder as the MAUIcount.py script (or another path that will be found). It needs to be edited appropriately for your project, and avoids the need for a long list of command line arguments.

Sets of parameters for different genes and data sets can be 'stored' as commented-out lines in MC_parameters.

===

Test

The folder Test provides some sample input data for MAUIcount.py. The supplied version of MC_parameters is set up to run an analysis of these data, provided that MAUIcount is run from the folder immediately above the Test folder. The command "python MAUIcount.py" should be sufficient, since all parameters are already set up in MC_parameters.py. If it is working correctly, MAUIcount should generate a new folder called MAUIcount_output within Test that has exactly the same files as the existing folder MAUIcount_expected_output.

===


