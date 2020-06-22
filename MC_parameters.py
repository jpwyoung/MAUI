'''
This file MC_parameters.py is required by MAUIcount.py and should be in the same folder 
as the script (or another path that will be found).
It specifies gene-specific parameters and the sequence files to be analysed and needs to
be edited appropriately for your project.
'''
 
 
#Parameters that control the stringency of the analysis and the output types
#These values will replace the defaults - comment them out to use the defaults
read_diff = 2 #default = 2
#Count UMI only if most abundant sequence has at least read_diff more reads than the next
reject_threshold = 0.7 #default = 0.7
#Reject sequences that occur as second sequences with UMIs at least reject_threshold
#times as often as they occur as the primary sequence
add_limit = 0.001 #default = 0.001
#sequences are included in rank order until the next would add a fraction less than add_limit
#i.e. this discards sequences with an overall relative abundance less than add_limit
output_types = 3 #default = 3
#add output filtered on a per-sample basis if set to 4
output_read_counts = False #default False
#If set to True, output the total number of reads contributing to the UMI primary and secondary counts


#Parameters that are specific for each amplicon, and the location of the data files.
#Comment out all except one set.

#For recA
UMI_len = 13
f_primer_len = 26
r_primer_len = 23
total_len = 313
working_folder = "./Test/"
'''
#For rpoB
UMI_len = 13
f_primer_len = 19
r_primer_len = 20
total_len = 306
working_folder = ""

#For nodA
UMI_len = 13
f_primer_len = 17
r_primer_len = 24
total_len = 314
working_folder = ""

#For nodD
UMI_len = 13
f_primer_len = 23
r_primer_len = 18
total_len = 308
working_folder = ""

#For viciae nodD (see doi.org/10.1111/nph.16392)
UMI_len = 13
f_primer_len = 23
r_primer_len = 21
total_len = 193
working_folder = ""
'''


