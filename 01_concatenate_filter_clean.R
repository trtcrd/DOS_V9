##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Concatenate script
##################################################################################

# The script successively go through each dataset, make an ASV table, perform chimera filtering and export the table
# It uses local paths, so you will need to adjust to your system. 

# library
library(dada2) 
library(seqinr) 

# export path  **to be adjusted**
path <- "/shared5/tristan/DOS/"
path.export <- paste0(path, "/tables_and_fastas/")
if(!dir.exists(path.export)) dir.create(path.export)

# list of dataset to process (and make paths)
dataset <- c("DS_first","DS_second","eDNAbyss","malaspina","tara","tara_polar","plankton_planet","lie","paw454","xu")

for (i in dataset) {
  ## path **to be adjusted**
  ## the merger files must be located in folders of the names of datasets
  path.rds <- paste0(path, "/RDS_merger/",i)
  
  # reimport the rds files (one file / sample)
  all_samples <- sort(list.files(path.rds, pattern="merger.rds", full.names = TRUE))
  sample.names <- sapply(strsplit(basename(all_samples), "_merger"), `[`,1)
  
  ### Making the table and remove chimeras
  # Read the merged data back in
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(sam in sample.names) {
    mergers[[sam]] <- readRDS(file.path(path.rds, paste0(sam, "_merger.rds")))
  }
  
  ## table and export to parent folder
  sta <- makeSequenceTable(mergers)
  saveRDS(sta, paste0(path.export, i, "_sta.rds"))
  
  ## extract the ASV sequences to compare with chimera free, and export the chimeric seqs 
  seqs_all <- getSequences(sta)
  # export reads count by sample (with chimeras...)
  withChim_counts <- cbind(samplesID = rownames(sta), Counts = rowSums(sta))
  write.table(withChim_counts, paste0(path.export, i, "_sta_counts.tsv"), quote = F, row.names = F, sep="\t")
  
  # chimera filtering
  sta_chim <- removeBimeraDenovo(sta, method = "consensus", multithread=TRUE)
  
  ## for eDNAbyss: filter out the ASVs from negative controls
  if (i == "eDNAbyss") {
    # filter from control EssNaut
    ctrl.samples <- c("EssNaut_DNA_extraction_blank") 
    real.samples <- c(rownames(sta_chim)[grep("EssNaut_P", rownames(sta_chim))])
    #found.in.ctrls <- colSums(sta_chim[ctrl.samples,])>0
    found.in.ctrls <- sta_chim[ctrl.samples,]>0
    sum(sta_chim[real.samples,found.in.ctrls]) / sum(sta_chim[real.samples,])
    # remove controls reads found in samples (but leaving the ASV there..)
    sta_chim[real.samples,found.in.ctrls] <- 0
    # finally remove the control sample from the table
    sta_chim <- sta_chim[!ctrl.samples == rownames(sta_chim),]
    
    # filter from control MARMINE
    ctrl.samples <- c("MARMINE_DNA_BLANK") 
    real.samples <- c(rownames(sta_chim)[grep("MARMINE_S", rownames(sta_chim))])
    #found.in.ctrls <- colSums(sta_chim[ctrl.samples,])>0
    found.in.ctrls <- sta_chim[ctrl.samples,]>0
    sum(sta_chim[real.samples,found.in.ctrls]) / sum(sta_chim[real.samples,])
    # remove controls reads found in samples (but leaving the ASV there..)
    sta_chim[real.samples,found.in.ctrls] <- 0
    # finally remove the control sample from the table
    sta_chim <- sta_chim[!ctrl.samples == rownames(sta_chim),]
    
    # filter from control MDW
    ctrl.samples <- c("Temoin_extraction_MedWaves_2", "Temoin_extraction_MedWaves_3") 
    real.samples <- c(rownames(sta_chim)[grep("MDW_S", rownames(sta_chim))])
    found.in.ctrls <- colSums(sta_chim[ctrl.samples,])>0
    #found.in.ctrls <- sta_chim[ctrl.samples,]>0
    sum(sta_chim[real.samples,found.in.ctrls]) / sum(sta_chim[real.samples,])
    # remove controls reads found in samples (but leaving the ASV there..)
    sta_chim[real.samples,found.in.ctrls] <- 0
    # finally remove the control sample from the table
    sta_chim <- sta_chim[!ctrl.samples[1] == rownames(sta_chim),]
    sta_chim <- sta_chim[!ctrl.samples[2] == rownames(sta_chim),]
  
    # remove other blank
    sta_chim <- sta_chim[-c(grep("Temoin_", rownames(sta_chim))),]
    # remove empty ASVs
    sta_chim <- sta_chim[,colSums(sta_chim)>0]
  }
  
  # extract the ASV sequences and export a fasta
  ASV_headers <- paste0("ASV", c(1:ncol(sta_chim)))
  seqs <- getSequences(sta_chim)
  write.fasta(as.list(seqs), ASV_headers, paste0(path.export, i, "_sta_noChim.fasta"))
  # get the chimeric seqs and export them
  chim <- seqs_all[!seqs_all %in% seqs]
  chim_headers <- paste0("ChimericASV_", c(1:length(chim)))
  write.fasta(as.list(chim), chim_headers, paste0(path.export, i, "_sta_chimSeqs.fasta"))
  
  # export the chimera-free table (to be merged with others)
  saveRDS(sta_chim, paste0(path.export, i, "_sta_noChim.rds"))
  # export reads count by sample
  noChim_counts <- cbind(samplesID = rownames(sta_chim), Counts = rowSums(sta_chim))
  write.table(noChim_counts, paste0(path.export, i, "_sta_noChim.tsv"), quote = F, row.names = F, sep="\t")
}





