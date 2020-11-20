# Author: Lina Sieverling

# usage: R --no-save --slave --args -t <discordantReadFileTumor> -c <discordantReadFileControl> -b <blacklist_file> -o <outFile> -f <function_file> < count_discordant_reads.R
# description: - uses the output tables of add_mate_mapq.py 
#              - gets the number of discordant reads in 1 kb windows
#              - merges adjacent windows with discordant reads in the tumor sample
#              - output: table for each pid containing raw discordant read counts for tumor and control


library("optparse")
library("data.table")

option_list = list(
  make_option(c("-t", "--discordantReadFileTumor"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("-c", "--discordantReadFileControl"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("-b", "--blacklist_file"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("-o", "--outFile"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("-f", "--function_file"), type="character", default="/abi/data/sieverling/global_scripts/functions.R", 
              help="file with R functions [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

discordantReadFileTumor = opt$discordantReadFileTumor
discordantReadFileControl = opt$discordantReadFileControl
blacklist_file = opt$blacklist
outFile = opt$outFile
function_file = opt$function_file

source(function_file)


#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

#######################################################################################################################################################

#------------------------------------------------------------------------------------
# get 1 kb windows with discordant read count
#------------------------------------------------------------------------------------


windowCountList = list()

discordantReadFileList = list(tumor=discordantReadFileTumor,
                            control=discordantReadFileControl)

for (sample in c("tumor", "control")){

  # get discordant reads
  discordantReadFile = discordantReadFileList[[sample]]
  
  # check if file exists, if not: skip sample
  if(!file.exists(discordantReadFile)){
    windowCount = data.frame(window=NA)
    windowCount[, paste0(sample, '_discordant_read_count')] = NA
    windowCountList[[sample]] = windowCount[numeric(0), ]
    next
    }

  discordantReads_all = read.table(discordantReadFile, header=TRUE, sep="\t", comment.char='')

  # only keep those with mate mapping quality > 30
  discordantReads = subset(discordantReads_all, mate_mapq>30)

  # get 1kb window
  discordantReads$mate_position_1kb = floor(discordantReads$mate_position/1000) * 1000
  #discordantReads$window = paste0(discordantReads$mate_chr, '_', as.character(discordantReads$mate_position_1kb))
  discordantReads$window = paste0(discordantReads$mate_chr, '_', as.character(discordantReads$mate_position_1kb), '_', discordantReads$mate_strand)
  

  # get read count per window
  windowCount = as.data.frame(table(discordantReads$window), stringsAsFactors=FALSE)
  colnames(windowCount) = c('window', paste0(sample, '_discordant_read_count'))

  #save in list
  windowCountList[[sample]] = windowCount
}


# merge window counts into 1 table and set missing values to 0
windowCountMerged = merge(windowCountList$tumor, windowCountList$control, by="window", all=TRUE)
windowCountMerged[is.na(windowCountMerged)] = 0


# add chromosome, start and end coordinates
windowCountMerged$chrom = gsub("_.*", "", windowCountMerged$window)
# windowCountMerged$chromStart = as.numeric(gsub(".*_", "", windowCountMerged$window))
windowCountMerged$chromStart = as.numeric(gsub("_", "", regmatches(windowCountMerged$window,regexpr("_.*_",windowCountMerged$window))))
windowCountMerged$chromEnd = windowCountMerged$chromStart + 1000
windowCountMerged$strand = gsub(".*_", "", windowCountMerged$window)

#sort colnames and set rownames
windowCountMerged$PID = gsub("_discordant_reads_1_kb_windows.tsv", "", basename(outFile))
windowCountMerged = windowCountMerged[ , c("PID", "window", "chrom", "chromStart", "chromEnd", "strand", "tumor_discordant_read_count", "control_discordant_read_count")]
row.names(windowCountMerged) = windowCountMerged$window



#------------------------------------------------------------------------------------
# add blacklist
#------------------------------------------------------------------------------------

if(file.exists(blacklist_file)){
  blacklist = read.table(blacklist_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  windowCountMerged[windowCountMerged$window %in% blacklist$window, "blacklisted"] = "yes"
  windowCountMerged[! windowCountMerged$window %in% blacklist$window, "blacklisted"] = "no"
}else{
  print("No blacklist file provided, or specified file does not exist. Continuing without blacklist")
  windowCountMerged$blacklisted = NA
}



#------------------------------------------------------------------------------------
# merge adjacent 1 kb window with discordant reads in the tumor sample
#------------------------------------------------------------------------------------

for (window1 in windowCountMerged[windowCountMerged$tumor_discordant_read_count!=0, "window"]){

  #skip window if row has already been removed
  if (! window1 %in% windowCountMerged$window){next}

  while(1){

    chromEnd = windowCountMerged[window1, "chromEnd"]
    chrom = windowCountMerged[window1, "chrom"]
    strand = windowCountMerged[window1, "strand"]

    # get window where chromosome is the same and whose starting point is the same as the end point of the first window
    window2 = windowCountMerged[windowCountMerged$chrom==chrom & 
                                  windowCountMerged$chromStart == chromEnd & 
                                  windowCountMerged$strand == strand, "window"]

    if(identical(window2, character(0))){break}

    #skip if read counts in window 1 or 2 are zero in tumor sample
    if(windowCountMerged[window1, "tumor_discordant_read_count"]==0 | windowCountMerged[window2, "tumor_discordant_read_count"]==0){break}

    # add counts
    windowCountMerged[window1, "chromEnd"] = windowCountMerged[window2, "chromEnd"]
    windowCountMerged[window1, "tumor_discordant_read_count"] = windowCountMerged[window1, "tumor_discordant_read_count"] + windowCountMerged[window2, "tumor_discordant_read_count"]
    windowCountMerged[window1, "control_discordant_read_count"] = windowCountMerged[window1, "control_discordant_read_count"] + windowCountMerged[window2, "control_discordant_read_count"]

    # if any of the windows are blacklisted, merged window is also blacklisted

    if (is.na(windowCountMerged[window1, "blacklisted"]) && is.na(windowCountMerged[window2, "blacklisted"])){
      windowCountMerged[window1, "blacklisted"] = NA
    }else if(compareNA(windowCountMerged[window1, "blacklisted"],"yes") || compareNA(windowCountMerged[window2, "blacklisted"],"yes")){
      windowCountMerged[window1, "blacklisted"] = "yes"
    }else{
      windowCountMerged[window1, "blacklisted"] = "no"
    }

    
    # remove row with window 2
    windowCountMerged = windowCountMerged[windowCountMerged$window!=window2, ]

  }
}

#------------------------------------------------------------------------------------
# save results
#------------------------------------------------------------------------------------
write.table(windowCountMerged, file=outFile, quote=FALSE, row.names = FALSE, sep="\t")


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
# 
# other attached packages:
#   [1] optparse_1.3.2   data.table_1.9.6
# 
# loaded via a namespace (and not attached):
#   [1] getopt_1.20.0 chron_2.3-47 
