# Author: Lina Sieverling

# usage: R --no-save --slave --args <window_file> <candidate_region_file> <tumor_discordant_read_lower_limit> <control_discordant_read_upper_limit> <consider_blacklist> <function_file> < get_candiate_regions.R
# description: filters windows to get telomere insertion candidate regions:
#              - window must contain at least <tumor_discordant_read_lower_limit> discordant reads in tumor sample
#                and at most <control_discordant_read_upper_limit> in control sample
#              - if specified, windows in blacklist are filtered out
#              output: table with filtered results (= telomere insertion candidate regions)


# get commandline arguments
commandArgs = commandArgs()
window_file = commandArgs[5]
candidate_region_file = commandArgs[6]
tumor_discordant_read_lower_limit = as.numeric(commandArgs[7])
control_discordant_read_upper_limit = as.numeric(commandArgs[8])
consider_blacklist = commandArgs[9]
function_file = commandArgs[10]

source(function_file)

#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

#########################################################################################################################


windowTable = read.table(window_file, header=TRUE, sep="\t", comment.char='')

#------------------------------------------------------------------------------------
# filter by tumor and control discordant read thresholds
#------------------------------------------------------------------------------------
candidate_regions = windowTable[windowTable$tumor_discordant_read_count >= tumor_discordant_read_lower_limit & 
                                                windowTable$control_discordant_read_count <= control_discordant_read_upper_limit, ]


#------------------------------------------------------------------------------------
# filter by blacklist
#------------------------------------------------------------------------------------

if(consider_blacklist=="True"){
  candidate_regions = candidate_regions[!compareNA(candidate_regions$blacklisted, "yes"), ]
}


#------------------------------------------------------------------------------------
# save raw and filtered results
#------------------------------------------------------------------------------------

write.table(candidate_regions, candidate_region_file, quote=FALSE, row.names = FALSE, sep="\t")


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

