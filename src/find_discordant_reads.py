# Author: Lina Sieverling

#!/usr/bin/python

# Usage: source activate telomereEnv
#        python /home/sieverli/Code/telomere_insertion_analysis/snakemake_telomere_insertions/src/find_discordant_reads.py \
#                -i <*_filtered_intratelomeric.bam resulting from TelomereHunter> \
#                -o <*_tumor_telomere_insertions.tsv output file>
# Description: extracts read names and positions of mates for reads that fullfill the following criteria:
#               - 1 mate is intratelomeric, the other mate is not
#               - mate that is not intratelomeric needs to be mapped (attention: no information about mapping quality)



import os
import sys, getopt
import re
import pysam

# ----------------------------------------------------------------
# read command line args 
# ----------------------------------------------------------------
myopts, args = getopt.getopt(sys.argv[1:],"i:o:")
 
for opt, arg in myopts:
  if opt == '-i':
    intratel_bam = arg
  elif opt == '-o':
    outfile_path = arg
  else:
    print("Usage: %s -i input_path (TelomereHunter intratelomeric bam file) -o outfile_path" % sys.argv[0])

#####################################################################################################################################


# ----------------------------------------------------------------
# make a dictionary with the number of reads with the same name
# ----------------------------------------------------------------
read_name_dict = {}

bamfile = pysam.Samfile( intratel_bam, "rb" )

for read in bamfile.fetch(until_eof=True):

  read_name = read.qname

  try:
    read_name_dict[read_name] += 1
  except:
    read_name_dict[read_name] = 1

bamfile.close()



# ---------------------------------------------------------------------------
# go through bam file again and extract mate mapping positions of reads if
# the mate is mapped and not intratelomeric
# ---------------------------------------------------------------------------   
bamfile2 = pysam.Samfile( intratel_bam, "rb" )

output = "read_name\tmate_chr\tmate_position\n"

for read in bamfile2.fetch(until_eof=True):

  #skip reads where the mate is unmapped
  if read.mate_is_unmapped:
    continue

  #skip reads where the reference ID of the mate is not known ('*', this can happen when the mapq of the mate is 255='not known')
  if read.next_reference_id==-1:
    continue

  #skip reads where the mate is also intratelomeric
  read_name = read.qname
  if read_name_dict[read_name] == 2:
    continue

  #get chromosome of mate
  mate_chr = read.next_reference_name

  # get 0-based mapping position of mate, adding 1 to get it 1-based like in SAM file
  mate_position = read.next_reference_start + 1 

  output += read_name + "\t" + str(mate_chr) + "\t" + str(mate_position) + "\n"


# ----------------------------------------------------------------
# write read name and mapping position of mate to output table
# ----------------------------------------------------------------
outfile = open( outfile_path, "w")
outfile.write(output)
outfile.close





