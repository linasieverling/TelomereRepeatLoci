# Author: Lina Sieverling

# Usage: R --no-save --slave --args <candidate_region_file> <bamfile> <outfile> <function_file> < neotelomeres_extract_fusion_reads.R
# Description: - goes through each candidate region and searches for soft- and hardclipped reads
#              - extracts clipped sequence and checks if it contains telomere repeats
#              - saves all reads in a table (read_name, clipped_sequence, position_of_clip, ...)



# get commandline arguments
commandArgs = commandArgs()
candidate_region_file = commandArgs[5]
bamfile = commandArgs[6]
outfile = commandArgs[7]
function_file = commandArgs[8]

library(GenomicAlignments)
library(stringr)

source(function_file)

#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

#----------------------------------------------------------------------
# get candidate region
#----------------------------------------------------------------------

candidate_regions = read.table(candidate_region_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE)
row.names(candidate_regions) = candidate_regions$window

clipped_reads_list = list()


for(window in candidate_regions$window){
  
  print(window) # for debugging
  
  chrom = candidate_regions[window, "chrom"]
  chromStart = candidate_regions[window, "chromStart"]
  chromEnd = candidate_regions[window, "chromEnd"]
  
  
  #----------------------------------------------------------------------
  # get soft-clipped reads in candidate region
  #----------------------------------------------------------------------
  view_window_command = paste0("samtools view ", bamfile, " ", chrom, ":", chromStart-300, "-", chromEnd+300)
  
  # extract read names
  read_names = system(paste0(view_window_command, " | cut -f 1"), intern=TRUE)

  # extract sequence
  read_seq = system(paste0(view_window_command, " | cut -f 10"), intern=TRUE)
  
  # extract position
  start = as.numeric(system(paste0(view_window_command, " | cut -f 4"), intern=TRUE))
  
  # extract cigar
  cigar = system(paste0(view_window_command, " | cut -f 6"), intern=TRUE)
  
  #extract flag
  flag = system(paste0(view_window_command, " | cut -f 2"), intern=TRUE)

  # make table with soft clipped reads
  reads_in_window_df = data.frame(read_name = read_names, flag = flag, start = start, cigar = cigar, sequence = read_seq, stringsAsFactors=FALSE)    
  reads_soft_clipped = reads_in_window_df[grep("S", reads_in_window_df$cigar), ]   
  
  # get end position
  reads_soft_clipped$end = reads_soft_clipped$start + 
    unlist(lapply(regmatches(reads_soft_clipped$cigar,gregexpr("\\d+[MD]",reads_soft_clipped$cigar)), function(x) 
      sum(as.numeric(gsub("[A-z]", "", x)))))
  
  # extract if read is first or second in pair
  reads_soft_clipped$read_1_2 = sapply(as.numeric(reads_soft_clipped$flag), function(x) 
    if(bitwAnd(x,0x40)){"READ1"}else if(bitwAnd(x,0x80)){"READ2"}else{NA})
  
  if(dim(reads_soft_clipped)[1]!=0){
    #add columns and sort
    reads_soft_clipped$window = window
    reads_soft_clipped$chr_primary_align = NA
    reads_soft_clipped$coord_primary_align = NA
    reads_soft_clipped$strand_primary_align = NA
    
    reads_soft_clipped = reads_soft_clipped[ , c("window", "read_name", "read_1_2", "start", "end", "cigar", "chr_primary_align", "coord_primary_align", "strand_primary_align", "sequence")]
    
  } else{
    reads_soft_clipped = data.frame(window=NA, read_name=NA, read_1_2=NA, start=NA, end=NA, cigar=NA, chr_primary_align=NA, coord_primary_align=NA, 
               strand_primary_align=NA, sequence=NA)[numeric(0), ]
  }
  

  
  #----------------------------------------------------------------------
  # get supplementary alignments in candidate region (hard-clipped reads)
  #----------------------------------------------------------------------
  
  SA_table_list = list()
  
  for (strandSupp in c("+", "-")){
    
    if (strandSupp == "+"){
      flag_expr="-f 2048 -F 16"
    }else if(strandSupp == "-"){
      flag_expr="-f 2064"
    }
    
    view_window_command = paste0("samtools view ", flag_expr, " ", bamfile, " ", chrom, ":", chromStart-300, "-", chromEnd+300)
    
    # extract read names
    read_names = system(paste0(view_window_command, " | cut -f 1"), intern=TRUE)
    
    # extract position
    pos = as.numeric(system(paste0(view_window_command, " | cut -f 4"), intern=TRUE))
    
    #extract cigar
    cigar = system(paste0(view_window_command, " | cut -f 6"), intern=TRUE)
    
    #extract if read is first or second in pair
    flag = system(paste0(view_window_command, " | cut -f 2"), intern=TRUE)
    read_1_2  = sapply(as.numeric(flag), function(x) if(bitwAnd(x,0x40)){"READ1"}else if(bitwAnd(x,0x80)){"READ2"}else{NA})
    
    
    # end of supp alignments (=sum of matches and deletions)
    end = pos + unlist(lapply(regmatches(cigar,gregexpr("\\d+[MD]",cigar)), function(x) sum(as.numeric(gsub("[A-z]", "", x)))))
    
    #extract SA Tag with position and strand of primary alignment
    tags = system(paste0(view_window_command, "| cut -f 12-"), intern=TRUE)    
    SA_tag = unlist(regmatches(tags,gregexpr("SA:Z:[^,]+,\\d+,[-\\+]",tags)))
    
    chr_primary_align = gsub(",\\d+,[-\\+]", "", gsub("SA:Z:", "", SA_tag))
    
    coord_primary_align = gsub(",[-\\+]", "", gsub("SA:Z:[^,]+,", "", SA_tag))
    
    strand_primary_align = gsub("SA:Z:[^,]+,\\d+,", "", SA_tag)

    
    #----------------------------------------------------------------------
    # make table with info on supplementary alignments in window
    #----------------------------------------------------------------------
    if(length(read_names)==0){
      SA_table_list[[strandSupp]] = data.frame(window=NA, 
                                               read_name = NA,
                                               read_1_2 = NA,
                                               start = NA,
                                               end = NA, 
                                               cigar = NA, 
                                               chr_primary_align = NA, 
                                               coord_primary_align = NA,
                                               strand_primary_align = NA, 
                                               stringsAsFactors=NA,
                                               sequence=NA)[numeric(0), ]
      next
    }
    SA_table = data.frame(window=window, 
                          read_name = read_names,
                          read_1_2 = read_1_2,
                          start = pos,
                          end = end, 
                          cigar = cigar, 
                          chr_primary_align = chr_primary_align, 
                          coord_primary_align = coord_primary_align,
                          strand_primary_align = strand_primary_align, 
                          stringsAsFactors=FALSE)
    
    #----------------------------------------------------------------------
    # get complete read sequences of hard-clipped reads
    #----------------------------------------------------------------------
    
    for (i in 1:dim(SA_table)[1]){
      
      coord_primary = as.numeric(SA_table$coord_primary_align[i])
      
      coord_suppl = SA_table$start[i]
      
      strand_primary = SA_table$strand_primary_align[i]
      
      read_1_2 = SA_table$read_1_2[i]
      
      if (strand_primary=="+"){
        strand_primary_flag = "-F 2064"
      }else{
        strand_primary_flag = "-F 2048 -f 16"
      }
      
      if (read_1_2=="READ1"){
        read_1_2_flag = "-f 64"
      }else if(read_1_2=="READ2"){
        read_1_2_flag = "-f 128"
      }
      
      seq = system(paste0("samtools view ", strand_primary_flag, " ", read_1_2_flag, " ", bamfile, " ", SA_table$chr_primary_align[i], ":", coord_primary, "-", coord_primary + 1, " | grep ", SA_table$read_name[i], " | awk '$4 == ", coord_primary, " {print}' | grep SA:Z: | grep ", coord_suppl, " | cut -f 10"), intern=TRUE)
     
      # if read is mapped to different strand than supplementary alignment: get reverse complement 
      if (strand_primary==strandSupp){
        SA_table[i, "sequence"] = seq
      }else{
        SA_table[i, "sequence"] = toString(reverseComplement(DNAString(seq)))
      }        
    }   
    SA_table_list[[strandSupp]] = SA_table
  }
  
  SA_table_all = do.call("rbind", SA_table_list)
  SA_table_all = SA_table_all[! is.na(SA_table_all$window), ]
  
  #--------------------------------------------------------------------------------------------------
  # merge soft- and hardclipped reads
  #--------------------------------------------------------------------------------------------------
  clipped_reads = rbind(reads_soft_clipped, SA_table_all)
  
  
  #--------------------------------------------------------------------------------------------------
  # check if clipped sequences contain telomere repeats
  #--------------------------------------------------------------------------------------------------
  if(dim(clipped_reads)[1]!=0){
    #get clipped sequence
    clipped_ranges = cigarRangesAlongQuerySpace(clipped_reads$cigar, ops=c("S", "H"), before.hard.clipping=TRUE)
    total_sequence = DNAStringSet(x=as.character(clipped_reads$sequence), start=NA, end=NA, width=NA, use.names=TRUE)
    
    clipped_sequence_DNA_string_set = extractAt(total_sequence, clipped_ranges)
    
    clipped_reads$clipped_sequence = unlist(lapply(clipped_sequence_DNA_string_set, toString))
    
    # check if sequence contains telomere repeats 
    clipped_reads$part_telomere = grepl("TTAGGG|CCCTAA", clipped_reads$clipped_sequence)   
    
    clipped_reads$TTAGGG_count = str_count(clipped_reads$clipped_sequence, "TTAGGG")
    clipped_reads$CCCTAA_count = str_count(clipped_reads$clipped_sequence, "CCCTAA")
    

    # determine whether clipped sequence (=telomere part) is upstream or downstream from read
    clipped_reads[grepl("^\\d+M.*\\d+[HS]$", clipped_reads$cigar), "expected_pos_fusion"] = "downstream"
    clipped_reads[grepl("^\\d+[HS].*\\d+M$", clipped_reads$cigar), "expected_pos_fusion"] = "upstream"
  }else{
    clipped_reads = data.frame(window=NA, read_name=NA, read_1_2=NA, start=NA, end=NA, cigar=NA, chr_primary_align=NA, coord_primary_align=NA, 
                               strand_primary_align=NA, sequence=NA, clipped_sequence=NA, part_telomere=NA, 
                               TTAGGG_count=NA, CCCTAA_count=NA, expected_pos_fusion=NA)[numeric(0), ]
  }
 
  clipped_reads_list[[window]] = clipped_reads
} 


#--------------------------------------------------------------------------------------------------
# save clipped reads in table
#--------------------------------------------------------------------------------------------------

clipped_reads_all = do.call("rbind", clipped_reads_list)

if(is.null(clipped_reads_all)){
  clipped_reads_all = data.frame(window=NA, read_name=NA, read_1_2=NA, start=NA, end=NA, cigar=NA, chr_primary_align=NA, coord_primary_align=NA, 
                                 strand_primary_align=NA, sequence=NA, clipped_sequence=NA, part_telomere=NA, 
                                 TTAGGG_count=NA, CCCTAA_count=NA, expected_pos_fusion=NA)[numeric(0), ]
}

write.table(clipped_reads_all, file=outfile, quote=FALSE, row.names = FALSE, sep="\t") 


#############################################################################################################
# > sessionInfo()
# R version 3.2.2 (2015-08-14)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: openSUSE 13.1 (Bottle) (x86_64)
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats4    parallel  stats     graphics  grDevices utils    
# [8] datasets  methods   base     
# 
# other attached packages:
#   [1] scales_0.4.1               RColorBrewer_1.1-2        
# [3] gridExtra_2.2.1            reshape2_1.4.1            
# [5] ggplot2_2.2.1              stringr_1.0.0             
# [7] GenomicAlignments_1.6.3    Rsamtools_1.22.0          
# [9] Biostrings_2.38.4          XVector_0.10.0            
# [11] SummarizedExperiment_1.0.2 Biobase_2.30.0            
# [13] GenomicRanges_1.22.4       GenomeInfoDb_1.6.3        
# [15] IRanges_2.4.8              S4Vectors_0.8.11          
# [17] BiocGenerics_0.16.1       
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.5          magrittr_1.5         zlibbioc_1.16.0     
# [4] munsell_0.4.3        BiocParallel_1.4.3   colorspace_1.2-6    
# [7] plyr_1.8.4           tools_3.2.2          gtable_0.2.0        
# [10] lambda.r_1.1.7       futile.logger_1.4.1  assertthat_0.1      
# [13] lazyeval_0.2.0       tibble_1.2           futile.options_1.0.0
# [16] bitops_1.0-6         stringi_1.0-1   

