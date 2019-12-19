'''
This file MC_parameters.py is required by MAUIcount.py and should be in the same folder 
as the script (or another path that will be found).
It specifies gene-specific parameters and the sequence files to be analysed and needs to
be edited appropriately for your project.
'''
 
 
#Parameters that control the stringency of the analysis
#These values will replace the defaults - comment them out to use the defaults

read_diff = 2 #default = 2
#Count seqid only if most abundant sequence has at least read_diff more reads than the next

reject_threshold = 0.7 #default = 0.7
#Reject sequences that occur as second sequences with seqids at least reject_threshold
#times as often as they occur as the primary sequence

add_limit = 0.001 #default = 0.001
#sequences are included in rank order until the next would add a fraction less than add_limit
#i.e. this discards sequences with an overall relative abundance less than add_limit

output_types = 3 #default = 3
#add output filtered on a per-sample basis if set to 4


#Parameters that are specific for each amplicon, and the location of the data files.
#Comment out all except one set.

#For recA
seqid_len = 13
f_primer_len = 26
r_primer_len = 23
total_len = 313
working_folder = "test/" #Use the example data provided in the 'test' folder
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2018-05-31_NCHAIN/assembled/recA/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/PEAR_assembled_synthseq_mix/assembled/recA/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/synth_seq_mix/recA/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/clover_geno/recA/"
'''
#For rpoB
seqid_len = 13
f_primer_len = 19
r_primer_len = 20
total_len = 306
working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2018-05-31_NCHAIN/assembled/rpoB/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/PEAR_assembled_synthseq_mix/assembled/rpoB/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/synth_seq_mix/rpoB/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/clover_geno/rpoB/"

#For nodA
seqid_len = 13
f_primer_len = 17
r_primer_len = 24
total_len = 314
working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2018-05-31_NCHAIN/assembled/nodA/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/PEAR_assembled_synthseq_mix/assembled/nodA/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/synth_seq_mix/nodA/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/clover_geno/nodA/"

#For nodD
seqid_len = 13
f_primer_len = 23
r_primer_len = 18
total_len = 308
working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2018-05-31_NCHAIN/assembled/nodD/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/PEAR_assembled_synthseq_mix/assembled/nodD/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/synth_seq_mix/nodD/"
#working_folder = "/Users/peter/Documents/Sequence_data/MiSeq_2019-01_NCHAIN/clovergeno_andphusion_samples/clover_geno/nodD/"

#For viciae nodD
seqid_len = 13
f_primer_len = 23
r_primer_len = 21
total_len = 193
working_folder = ""
'''


