##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Aggregate sequencing replicates per BioSamples - and export FULL table
##################################################################################

# Some samples were sequenced multiple time across differents runs (deep_sea, tara and tara_polar).
# This script aggregates the reads per BioSamples, clean some samples names and export the full table.


## Path **to be adjusted**
path <- "/shared5/tristan/DOS/tables_and_fastas/"
setwd(path)

tara <- readRDS("tara_sta_noChim.rds")
polar <- readRDS("tara_polar_sta_noChim.rds")

library(dada2)
library(seqinr)

########### aggregate by BioSample
#### fetch the metadata 
comp <- read.table("../all_libs_biosample_mapping.txt", header=TRUE, sep="\t", dec = ".")
s <- "sampleID"
rownames(comp) <- comp[,s]

### TARA 
TARA_comp <- subset(comp, comp[,s] %in% rownames(tara))
TARA_samples_comp <- as.character(TARA_comp[,s])
# then check order
if(table(TARA_samples_comp == rownames(tara))["TRUE"] != nrow(TARA_comp) || table(TARA_samples_comp == rownames(tara))["FALSE"] == nrow(TARA_comp))
{
  message("The order of samples does not match with the metadata file, trying to sort it...")
  tara <- tara[TARA_samples_comp,]
  if(table(TARA_samples_comp == rownames(tara))["TRUE"] == nrow(TARA_comp))
  {
    message("Ok, done...")
  } else { stop("Could not sort") }
} else { message("Order of samples is fine") }
## check exact match
table(rownames(tara) == TARA_comp[,s])

## now agg by BioSample column
tara_AGG <- aggregate(as.data.frame(tara), list(TARA_comp$bioSampleID), FUN = sum)
rownames(tara_AGG) <- tara_AGG[,1]
tara_AGG <- tara_AGG[,2:ncol(tara_AGG)]
saveRDS(tara_AGG, "tara_sta_noChim_agg.rds")

##################
### POLAR 

POLAR_comp <- subset(comp, comp[,s] %in% rownames(polar))
POLAR_samples_comp <- as.character(POLAR_comp[,s])
# then check order
if(table(POLAR_samples_comp == rownames(polar))["TRUE"] != nrow(POLAR_comp))
{
  message("The order of samples does not match with the metadata file, trying to sort it...")
  polar <- polar[POLAR_samples_comp,]
  if(table(POLAR_samples_comp == rownames(polar))["TRUE"] == nrow(POLAR_comp))
  {
    message("Ok, done...")
  } else { stop("Could not sort") }
} else { message("Order of samples is fine") }
## check exact match
table(rownames(polar) == POLAR_comp[,s])

## now agg by BioSample column
POLAR_AGG <- aggregate(as.data.frame(polar), list(POLAR_comp$bioSampleID), FUN = sum)
rownames(POLAR_AGG) <- POLAR_AGG[,1]
POLAR_AGG <- POLAR_AGG[,2:ncol(POLAR_AGG)]
saveRDS(POLAR_AGG, "tara_polar_sta_noChim_agg.rds")

######### and concatenate all into full table

ds1 <- readRDS("DS_first_sta_noChim.rds")
ds2 <- readRDS("DS_second_sta_noChim.rds")
mala <- readRDS("malaspina_sta_noChim.rds")
# rename malaspina samples
rownames(mala) <- gsub("5477-", "", rownames(mala))
rownames(mala) <- gsub("-Euk1389F-1510R", "", rownames(mala))
pplanet <- readRDS("plankton_planet_sta_noChim.rds")
abyss <- readRDS("eDNAbyss_sta_noChim.rds")
paw <- readRDS("paw454_sta_noChim.rds")
xu <- readRDS("xu_sta_noChim.rds")
lie <- readRDS("lie_sta_noChim.rds")
## AGG tara and polar
polar <- readRDS("tara_polar_sta_noChim_agg.rds")
tara <- readRDS("tara_sta_noChim_agg.rds")
##
polar <- as.matrix(polar)
tara <- as.matrix(tara)

# merge the two runs deep sea and sum the reads of identical samples
combined_ds <- mergeSequenceTables(ds1, ds2, repeats = "sum")

# merge the deep sea, malaspina and the tara (if encountering samples with identical names: throw an error)
combined_all <- mergeSequenceTables(combined_ds, mala, pplanet, abyss, paw, xu, lie, tara, polar, repeats = "error")

## And export the full table (hooray :)
saveRDS(combined_all, "Full_ASVs_table.rds")

# export a fasta with all ASVs (identical order as in the table; named 1:n)
# Taxonomic annotations done with SLIM (Dufresne et al., 2019, BMC bioinformatics) module 'asssignment-fasta-vsearch'
# see material and methods for details on parameters and reference databases used.
ASV_headers <- paste0("ASV", c(1:ncol(combined_all)))
seqs <- getSequences(combined_all)
write.fasta(as.list(seqs), ASV_headers, "Full_ASVs_table.fasta")









