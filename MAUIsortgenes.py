#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

#This script sorts sequences (previously assembled by PEAR) into separate files for 
#rpoB, recA, nodA, nodD genes, based on a short unambiguous motif in the forward primer.

working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/"


def paf(print_string):
    """
    Print and file: everything printed to stdout is also written in logfile.
    """
    print(print_string)
    with open(logfile, "a") as outfile:
        outfile.write(str(print_string) + "\n")

def sort_by_tag(seq):
    """
    Identify the gene for each read using a short motif in the forward primer 
    """
    genetag = str(seq[15:21].seq)
    if genetag == "TCGCAG": genes_dict["rpob"].append(seq)
    elif genetag == "AGAATG": genes_dict["reca"].append(seq)
    elif genetag == "GGATCT": genes_dict["noda"].append(seq)
    elif genetag == "GCGTTT": genes_dict["nodd"].append(seq)
    elif len(seq) < 290:  genes_dict["short"].append(seq)
    else: genes_dict["unknown"].append(seq)
    return()



with open(working_folder + "list.txt") as list_file:
    for line in list_file:
        file_name = line.rstrip("\n")
        if "/" in file_name:
            file_name = file_name[file_name.rindex("/")+1:]

 
        infile = working_folder + file_name
        sample = file_name[:file_name.find("_")]
        logfile = working_folder + sample + "_logfile.txt"



        seq_list = SeqIO.parse(infile, "fastq")  

        genes_dict = {"rpob":["null"],"reca":["null"], "noda":["null"],"nodd":["null"],"unknown":["null"], "short":["null"]}
        padding_seq = Seq("N", IUPAC.ambiguous_dna) 
        padding = SeqRecord(padding_seq)
        n=0
        for seq in seq_list: 
            sort_by_tag(seq)
            n += 1
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

        paf("\ncounts allowing up to 1 missing base:")

        for each_gene_list in genes_dict:
            paf(each_gene_list + " " + str(len(genes_dict[each_gene_list])-1))
    
        for gene_type in genes_dict:
            genes_dict[gene_type].pop(0)
            gene_file = working_folder + sample + "_" + gene_type + ".fastq"
            with open(gene_file, "w") as outfile:
                SeqIO.write(genes_dict[gene_type], outfile, "fastq")
