# Author: Lina Sieverling

# Usage: R --no-save --slave --args <candidate_region_file> <clipped_reads_file> <discordant_read_file> <outfile> <function_file> < predict_insertion_sites.R
# Description: trys to predict a telomere insertion site for each candidate region from clipped reads of tumor sample
#              - takes the position where most clipped sequences start/end (if this is not unique it returns NA)
#              - add the result to the extended candidate region table


# get commandline arguments
commandArgs = commandArgs()
candidate_region_file = commandArgs[5]
clipped_reads_file = commandArgs[6]
discordant_read_file = commandArgs[7]
outfile = commandArgs[8]
function_file = commandArgs[9]

source(function_file)

#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

##########################################################################################################################


candidate_regions = read.table(candidate_region_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE)
row.names(candidate_regions) = candidate_regions$window

clipped_reads_all = read.table(clipped_reads_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE, comment.char='')

discordant_read_table = read.table(discordant_read_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE, comment.char='')


for(window in unique(clipped_reads_all$window)){

  #-------------------------------------------------------------------------------------------------------------------
  # get most likely insertion site (where do the most clipped sequences start or end?)
  #-------------------------------------------------------------------------------------------------------------------

  clipped_reads = clipped_reads_all[clipped_reads_all$window==window, ]
  
  if(sum(clipped_reads$part_telomere)==0){
    
    candidate_regions[window, "sum_TTAGGG_count"] = NA
    candidate_regions[window, "sum_CCCTAA_count"] = NA
    candidate_regions[window, "repeat_forward"] = NA
    candidate_regions[window, "insertion_site"] = NA
    candidate_regions[window, "pos_telomeres_from_insertion"] = NA
    candidate_regions[window, "reads_supporting_insertion_pos"] = NA
  }else{
    
    clipped_reads_filtered = clipped_reads[clipped_reads$part_telomere, ]
    
    
    #-----------------------------------------------------------------------------------
    # discordant reads on plus strand => clipped reads should end at same position
    #
    # 5' ------ chromosome ------ telomere 3'
    #
    #-----------------------------------------------------------------------------------
    
    if (compareNA(candidate_regions[window, "strand"], "+")){
      clipped_expected_pos_fusion = "downstream"
      clipped_start_end = "end"
    }
    

    
    #-----------------------------------------------------------------------------------
    # discordant reads on minus strand => clipped reads should start at same position
    #
    # 5' telomere ------ chromosome ------ 3'
    #
    #-----------------------------------------------------------------------------------
    
    if (compareNA(candidate_regions[window, "strand"], "-")){
      clipped_expected_pos_fusion = "upstream"
      clipped_start_end = "start"
    }
    
    
    #-------------------------------------------------------------------------
    # only keep clipped reads that match orientation of discordant reads
    #-------------------------------------------------------------------------

    #attention: does not consider reads that are clipped on both ends (these are unlikely to be indicative of telomere insertion)
    clipped_reads_filtered_matching_discordant = clipped_reads_filtered[compareNA(clipped_reads_filtered$expected_pos_fusion, clipped_expected_pos_fusion), ]

    #-------------------------------------------------------------------------
    # only keep clipped reads that match position of discordant reads
    #-------------------------------------------------------------------------
    
    # currently taking median to prevent missing insertions where there is another discordant read nearby
    chrom = candidate_regions[window, "chrom"]
    window_start = candidate_regions[window, "chromStart"]
    window_end = candidate_regions[window, "chromEnd"]
    
    discordant_reads = discordant_read_table[discordant_read_table$mate_chr==chrom & 
                                             discordant_read_table$mate_position>= window_start &
                                             discordant_read_table$mate_position<= window_end &
                                             discordant_read_table$mate_strand==candidate_regions[window, "strand"], ]
    
    discordant_read_pos_median = median(discordant_reads$mate_position) + 50   #plus 50 to get middle of read => also not ideal
    
    if (compareNA(candidate_regions[window, "strand"], "+")){
      clipped_reads_filtered_matching_discordant = clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant$end>discordant_read_pos_median, ]
      
    }else if (compareNA(candidate_regions[window, "strand"], "-")){
      clipped_reads_filtered_matching_discordant = clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant$start<discordant_read_pos_median, ]
    }
    
 
    
    #------------------------------------------------------
    # where do most reads start/end?
    #------------------------------------------------------
    
    table_pos_insertion = as.data.frame(table(clipped_reads_filtered_matching_discordant[, clipped_start_end]),
                                          stringsAsFactors=FALSE)
    
    if (dim(table_pos_insertion)[1] != 0){
      colnames(table_pos_insertion) = c("pos", "Freq")
      
      # go through insertion positions again and count number of unique cigars (this prevents that we only count reads that map at exactly the same position)
      for(pos in table_pos_insertion$pos){
        unique_cigars = length(unique(clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant[, clipped_start_end]==pos, "cigar"]))
        table_pos_insertion[table_pos_insertion$pos==pos, "unique_cigars"] = unique_cigars
      }
      
    } else{
      table_pos_insertion = data.frame(pos=NA, Freq=NA, unique_cigars=NA)
    }


    insertion_pos = table_pos_insertion[table_pos_insertion$unique_cigars==max(table_pos_insertion$unique_cigars), "pos"]

    
    if(length(insertion_pos)==1){
      candidate_regions[window, "insertion_site"] = insertion_pos
      candidate_regions[window, "pos_telomeres_from_insertion"] = clipped_expected_pos_fusion 
      candidate_regions[window, "reads_supporting_insertion_pos"] = max(table_pos_insertion$unique_cigars)  
      
      #-----------------------------------------------------------------------------------
      # get total TTAGGG and CCCTAA counts in telomere fusion reads at insertion site
      # and determine most likely repeat on forward strand
      #-----------------------------------------------------------------------------------
      
      clipped_reads_filtered_at_insertion = clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant[, clipped_start_end]==insertion_pos,]
      
      sum_TTAGGG_count = sum(clipped_reads_filtered_at_insertion$TTAGGG_count)
      sum_CCCTAA_count = sum(clipped_reads_filtered_at_insertion$CCCTAA_count)
      
      candidate_regions[window, "sum_TTAGGG_count"] = sum_TTAGGG_count
      candidate_regions[window, "sum_CCCTAA_count"] = sum_CCCTAA_count
      
      if(sum_TTAGGG_count > sum_CCCTAA_count){
        candidate_regions[window, "repeat_forward"] = "TTAGGG"
      }else if (sum_CCCTAA_count > sum_TTAGGG_count){
        candidate_regions[window, "repeat_forward"] = "CCCTAA"
      }else{
        candidate_regions[window, "repeat_forward"] = NA
      }
    }else{
      candidate_regions[window, "insertion_site"] = NA
      candidate_regions[window, "pos_telomeres_from_insertion"] = NA
      candidate_regions[window, "reads_supporting_insertion_pos"] = NA
      candidate_regions[window, "sum_TTAGGG_count"] = NA
      candidate_regions[window, "sum_CCCTAA_count"] = NA
      candidate_regions[window, "repeat_forward"] = NA
    }

  }
}

if (dim(candidate_regions)[1]==0){
  
  candidate_regions = data.frame(PID=NA, window=NA, chrom=NA, chromStart=NA, chromEnd=NA, strand=NA, tumor_discordant_read_count=NA, control_discordant_read_count=NA, blacklisted=NA,
                                 insertion_site=NA, pos_telomeres_from_insertion=NA, reads_supporting_insertion_pos=NA,
                                 sum_TTAGGG_count=NA, sum_CCCTAA_count=NA, repeat_forward=NA)[numeric(0), ]
}


write.table(candidate_regions, file=outfile,
            quote=FALSE, row.names = FALSE, sep="\t") 


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
#   [1] stats     graphics  grDevices utils     datasets  methods   base  