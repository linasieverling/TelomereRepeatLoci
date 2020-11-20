# Author: Lina Sieverling

# Usage: R --no-save --slave --args <candidate_regions_file> <clipped_read_file> <outfile> < get_consensus.R
# Description: - gets consensus sequences (consensus base must make up 65% of all clipped bases, otherwise it is "N")
#              - get flanking reference sequence
#              - determines how many bases of sequence microhomology there is



# get commandline arguments
commandArgs = commandArgs()
candidate_regions_file = commandArgs[5]
clipped_read_file = commandArgs[6]
outfile = commandArgs[7]

library('BSgenome.Hsapiens.UCSC.hg19')

########################################################################################################################################

candidate_regions = read.table(candidate_regions_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(candidate_regions) = candidate_regions$window

clipped_reads = read.table(clipped_read_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="")


for (window in candidate_regions$window){
  
  insertion_site = candidate_regions[window, "insertion_site"]
  
  if (is.na(insertion_site)){
    candidate_regions[window, "consensus"] = NA
    candidate_regions[window, "flanking_seq"] = NA
    candidate_regions[window, "bp_microhomology"] = NA
    next
  }
  
  if (candidate_regions[window, "strand"]=="+"){
    start_end = "end"
  }else if (candidate_regions[window, "strand"]=="-"){
    start_end = "start"
  }
  
  #------------------------------------------------------------
  # read in clipped sequences
  #------------------------------------------------------------
  
  clipped_reads_window = clipped_reads[clipped_reads$window == window, ]

  clipped_reads_window = clipped_reads_window[clipped_reads_window[ , start_end]==insertion_site, ]
  
  if(start_end=="end"){
    sequences = gsub(".*, ", "", clipped_reads_window$clipped_sequence)
  }else if (start_end=="start"){
    sequences = gsub(",.*", "", clipped_reads_window$clipped_sequence)
    sequences = sapply(lapply(strsplit(sequences, NULL), rev), paste, collapse="")
  }
  
  #------------------------------------------------------------
  # get consensus
  #------------------------------------------------------------
  
  consensus = c()
  
  for (i in 1:max(nchar(sequences))){
    bases = substr(sequences, i, i)
    
    read_count = sum(bases != "")
    
    freqs = c(A = sum(bases=="A")/read_count,
              C = sum(bases=="C")/read_count,
              G = sum(bases=="G")/read_count,
              T = sum(bases=="T")/read_count,
              N = sum(bases=="N")/read_count)
    
    max_freq = max(freqs)
    
    if(max_freq>=0.65){
      base = names(which.max(freqs))
    }else{
      base = "N"
    }
    
    consensus = c(consensus, base)
  }
  
  consensus = paste0(consensus, collapse = '')
  
  if (candidate_regions[window, "strand"]=="-"){
    consensus = sapply(lapply(strsplit(consensus, NULL), rev), paste, collapse="")
  }
  
  candidate_regions[window, "consensus"]  = consensus
  
  #------------------------------------------------------------
  # get flanking reference sequence         
  #------------------------------------------------------------
  
  chrom_chr = paste0("chr", candidate_regions[window, "chrom"])
  
  if (candidate_regions[window, "strand"]=="+"){
    seq = as.character(getSeq(Hsapiens, chrom_chr, insertion_site-20, insertion_site-1))
  }else if (candidate_regions[window, "strand"]=="-"){
    seq = as.character(getSeq(Hsapiens, chrom_chr, insertion_site, insertion_site+19)) 
  }
  
  candidate_regions[window, "flanking_seq"]  = seq
  
  
  #------------------------------------------------------------
  # check microhomology         
  #------------------------------------------------------------
  
  #set flag for wrong junction repeat on telomere side
  repeat_junction_wrong = FALSE
  
  #get indices of telomere repeats in telomere insertion consensus sequence
  pos_match = gregexpr("TTAGGG|CCCTAA", consensus)[[1]]

  if (pos_match[1]==-1){
    candidate_regions[window, "bp_microhomology"] = "?"
    next
  }
  
  if (is.na(candidate_regions[window, "repeat_forward"])){
    candidate_regions[window, "bp_microhomology"] = "?"
    next
  }
  
  
  if (candidate_regions[window, "strand"]=="+"){
    
    #first telomere repeat match
    indice_match = pos_match[1]
    
    #get split telomere repeat bases on reference sequence and telomere insertion
    junction_repeat_ref = substr(candidate_regions[window, "repeat_forward"], 1, 6-indice_match+1)
    junction_repeat_tel = substr(candidate_regions[window, "repeat_forward"], 6-indice_match+2, 6)
    
    #check if telomere insertion consensus sequence and expected junction repeat match
    if(junction_repeat_tel!=substr(consensus, 1, indice_match-1)){
      repeat_junction_wrong = TRUE
    }
    
    #add additional telomere repeats to sequence expected on reference side (for perfect homology)
    repeats_ref = paste0(paste0(rep(candidate_regions[window, "repeat_forward"], 3), collapse=''), junction_repeat_ref)
    
    #reverse sequences (for counting from 1)
    repeats_ref = sapply(lapply(strsplit(repeats_ref, NULL), rev), paste, collapse="")
    seq_homology = sapply(lapply(strsplit(seq, NULL), rev), paste, collapse="")

  }else if (candidate_regions[window, "strand"]=="-"){
    #last telomere repeat match
    indice_match = nchar(consensus)+1 - (pos_match[length(pos_match)] + 6) 
    
    #get split telomere repeat bases on reference sequence and telomere insertion
    junction_repeat_tel = substr(candidate_regions[window, "repeat_forward"], 1, indice_match)
    junction_repeat_ref = substr(candidate_regions[window, "repeat_forward"], indice_match+1, 6)
    
    #check if telomere insertion consensus sequence and expected junction repeat match
    if(junction_repeat_tel!=substr(consensus, (pos_match[length(pos_match)] + 6), nchar(consensus)+1)){
      repeat_junction_wrong = TRUE
    }
    
    #add additional telomere repeats to sequence expected on reference side (for perfect homology)
    repeats_ref = paste0(junction_repeat_ref, paste0(rep(candidate_regions[window, "repeat_forward"], 3), collapse=''))
    
    seq_homology = seq
  }

  #check how many bases in reference genome match with the expected telomere repeats for perfect homology
  if(indice_match>6){
    candidate_regions[window, "bp_microhomology"] = "?"
  }else if (repeat_junction_wrong==TRUE){
    candidate_regions[window, "bp_microhomology"] = "?"
  }else{
    cntr = 0  
    for (i in 1:20){
      if(substr(seq_homology,i, i) == substr(repeats_ref,i, i)){
        cntr = cntr + 1
      }else{
        break
      }
    }
    
    candidate_regions[window, "bp_microhomology"] = cntr
  }
}

if (dim(candidate_regions)[1]==0){
  candidate_regions = data.frame(PID=NA, window=NA, chrom=NA, chromStart=NA, chromEnd=NA, strand=NA, tumor_discordant_read_count=NA, control_discordant_read_count=NA,
                                 blacklisted=NA, insertion_site=NA, pos_telomeres_from_insertion=NA, reads_supporting_insertion_pos=NA,
                                 sum_TTAGGG_count=NA, sum_CCCTAA_count=NA, repeat_forward=NA,
                                 consensus=NA, flanking_seq=NA, bp_microhomology=NA)[numeric(0), ]
}


#------------------------------------------------------------------------------------
# save results
#------------------------------------------------------------------------------------
write.table(candidate_regions, file=outfile, quote=FALSE, row.names = FALSE, sep="\t")



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
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.38.0                  
# [3] rtracklayer_1.30.4                Biostrings_2.38.4                
# [5] XVector_0.10.0                    GenomicRanges_1.22.4             
# [7] GenomeInfoDb_1.6.3                IRanges_2.4.8                    
# [9] S4Vectors_0.8.11                  BiocGenerics_0.16.1              
# 
# loaded via a namespace (and not attached):
#   [1] XML_3.98-1.4               Rsamtools_1.22.0          
# [3] GenomicAlignments_1.6.3    bitops_1.0-6              
# [5] futile.options_1.0.0       zlibbioc_1.16.0           
# [7] futile.logger_1.4.1        lambda.r_1.1.7            
# [9] BiocParallel_1.4.3         tools_3.2.2               
# [11] Biobase_2.30.0             RCurl_1.95-4.8            
# [13] SummarizedExperiment_1.0.2
