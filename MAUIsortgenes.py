#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
This script sorts sequences (previously assembled by PEAR) into separate files for 
each gene, based on a short unambiguous motif (tag) in the forward primer. 
Sequences that do not match any tag are filed separately as "short" or "unknown".
"""

#===Project-specific details that will need to be modified for another project:
working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/"
tag_to_gene = {"TCGCAG":"rpob","AGAATG":"reca", "GGATCT":"noda", "GCGTTT":"nodd"}
tag_start,tag_end = 15,21
min_length = 290
#===

def paf(print_string):
    """
    Print and file: everything printed to stdout is also written in logfile.
    """
    print(print_string)
    with open(logfile, "a") as outfile:
        outfile.write(str(print_string) + "\n")

def sort_by_tag(seq):
    """
    Identify the gene for each read using a short motif in the forward primer and append
    the sequence to the appropriate list in genes_dict
    """
    tag = str(seq[tag_start:tag_end].seq)
    if tag in tag_to_gene:
        genes_dict[tag_to_gene[tag]].append(seq)
    elif len(seq) < min_length:  genes_dict["short"].append(seq)
    else: genes_dict["unknown"].append(seq)
    return()

padding_seq = Seq("N") 
padding = SeqRecord(padding_seq)

#Process the input files one at a time
with open(working_folder + "list.txt") as list_file:
    for line in list_file:
    
        #sort out file names and get ready to read input
        file_name = line.rstrip("\n")
        if "/" in file_name:
            file_name = file_name[file_name.rindex("/")+1:] 
        infile = working_folder + file_name
        sample = file_name[:file_name.find("_")]
        logfile = working_folder + sample + "_logfile.txt"
        seq_list = SeqIO.parse(infile, "fastq")  

        #initialise genes_dict
        genes_dict = {"unknown":["null"], "short":["null"]}
        for gene in tag_to_gene.values():
            genes_dict[gene] = ["null"]
        
        #Add each sequence to the appropriate list in genes_dict
        n=0
        for seq in seq_list: 
            sort_by_tag(seq)
            n += 1
            
        #print information on sequence counts
        paf("\nSample: " + file_name[:-6])
        paf("total input sequences: " + str(n))
        paf("counts on first pass:")
        for each_gene_list in genes_dict:
            paf(each_gene_list + " " + str(len(genes_dict[each_gene_list])-1))

        #The most common error is a missing base in the seqid, so try the unidentified sequences
        #again with a 1-base offset    
        if len(genes_dict["unknown"]) > 1:
            unknown_list = genes_dict["unknown"][1:]    
            genes_dict["unknown"] = ["null"]
            for seq in unknown_list[1:]:
                annotations = seq.letter_annotations["phred_quality"]
                seq.letter_annotations = {}
                new_seq = seq
                new_seq.seq = padding.seq + seq.seq
                new_seq.letter_annotations["phred_quality"] = [0] + annotations
                sort_by_tag(new_seq)

        #print the new counts
        paf("\ncounts allowing up to 1 missing base:")
        for each_gene_list in genes_dict:
            paf(each_gene_list + " " + str(len(genes_dict[each_gene_list])-1))

        #write the sorted sequences to the appropriate files    
        for gene_type in genes_dict:
            genes_dict[gene_type].pop(0)
            gene_file = working_folder + sample + "_" + gene_type + ".fastq"
            with open(gene_file, "w") as outfile:
                SeqIO.write(genes_dict[gene_type], outfile, "fastq")
