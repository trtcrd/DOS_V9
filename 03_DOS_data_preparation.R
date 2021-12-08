##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Data preparation for analysis 
##################################################################################

# This script prepares and export the objects needed for data analysis
# Scripts for figures are available separately 


## setting path **to be adjusted**
path <- "/to/be/adjusted/DOS/"
setwd(path)

## output
path.out <- paste0(path, "/output")
if(!dir.exists(path.out)) dir.create(path.out)

## custom functions **to be adjusted to the companion_functions.R path**
source("00_companionFunctions.R")


# loading libraries
library('doMC')
registerDoMC(cores = 16)
library(ShortRead)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(treemapify)

####################################################################################################################################
###### Reimport the concatenated ASV table in RDS format

asv <- readRDS("Full_ASVs_table.rds")
dna <- readFasta("Full_ASVs_table.fasta")
dna_seq <- as.character(sread(dna))

#### metadata 
comp <- read.table("full_metadata_aggregatedReplicates.txt", header=TRUE, sep="\t", dec = ".")
s <- "sampleID"
rownames(comp) <- comp[,s]

## make a vector of biome: euphotic, aphotic and sediment 
biome <- rep(NA, nrow(comp))
for (i in 1:nrow(comp))
{
  if (comp$Biome[i] == "Surface" | comp$Biome[i] == "Epipelagic" | comp$Biome[i] == "DCM") biome[i] <- "Euphotic"
  if (comp$Biome[i] == "Mesopelagic" | comp$Biome[i] == "Bathypelagic" | comp$Biome[i] == "Bathypelagic" | comp$Biome[i] == "Abyssopelagic") biome[i] <- "Aphotic"
  if (comp$Biome[i] == "Sediment") biome[i] <- "Sediment"
}
comp$biome <- biome

## make a vector for the Watling et al, 2013, prog. in oceanography abyssal provinces 
watling <- rep("", nrow(comp))
for (i in 1:nrow(comp))
{
  if (comp$Biome[i] == "Sediment" & comp$area[i] == "North Pacific Ocean" & comp$expedition[i] != "KuramBio") watling[i] <- "AB11"
  if (comp$Biome[i] == "Sediment" & comp$area[i] == "North Pacific Ocean" & comp$expedition[i] == "KuramBio") watling[i] <- "AB13"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "VEMA") watling[i] <- "AB2"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "SYSTCOII") watling[i] <- "AB6"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "MSM39") watling[i] <- "AB2"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "DIVA3" & (comp$station[i] == "M79_1/1/519" | comp$station[i] == "M79_1/1/518" |
                                                                     comp$station[i] == "M79_1/2/553" | comp$station[i] == "M79_1/2/550" | comp$station[i] == "M79_1/2/552" )) watling[i] <- "AB5"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "DIVA3" & 
      (comp$station[i] == "M79_1/4/598" | comp$station[i] == "M79_1/4/599" | comp$station[i] == "M79_1/4/600" | comp$station[i] == "M79_1/4/603" )) watling[i] <- "AB3"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "eDNAbyss" & comp$area[i] == "Artic Ocean") watling[i] <- "AB1"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "eDNAbyss" & comp$area[i] == "North Atlantic Ocean") watling[i] <- "AB2"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "eDNAbyss" & comp$area[i] == "Mediterranean Sea") watling[i] <- "Mediterranean"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "ANDEEP II/ARK XXIV" & comp$area[i] == "Southern Ocean") watling[i] <- "AB6"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "ANDEEP II/ARK XXIV" & comp$area[i] == "Artic Ocean") watling[i] <- "AB1"
  if (comp$Biome[i] == "Sediment" & comp$expedition[i] == "Global Protistan Survey") watling[i] <- "AB11"
}
comp$watling <- watling

###### taxonomic annotations using both SILVA v138 and PR2-V9 (https://doi.org/10.5281/zenodo.3768951)
### We used the module "assignment-fasta-vsearch" of SLIM with three cutoffs:
# 95: >=95% & <99% similarity: LCA with 3 candidates, >=99% direct annotation
# 85: >=85% & <99% similarity: LCA with 3 candidates, >=99% direct annotation --- (for SILVA we used 80%)
# bestMatch: Best match above 0.02 % similarity (1 candidate)
# PR2 custom V9
taxo <- read.table("Full_ASVs_table_assigned-vsearch_customPR2_95.tsv", header=TRUE, sep="\t", dec = ".", row.names = 1)
taxo85 <- read.table("Full_ASVs_table_assigned-vsearch_customPR2_85.tsv", header=TRUE, sep="\t", dec = ".", row.names = 1)
taxo02 <- read.table("Full_ASVs_table_assigned-vsearch_customPR2_bestMatch.tsv", header=TRUE, sep="\t", dec = ".", row.names = 1)
# SILVA v138
taxo_silva <- read.table("Full_ASVs_table_assigned-vsearch_SILVAv138_95.tsv", header=TRUE, sep="\t", dec = ".", row.names = 1)
taxo_silva80 <- read.table("Full_ASVs_table_assigned-vsearch_SILVAv138_80.tsv", header=TRUE, sep="\t", dec = ".", row.names = 1)

### reorder each taxonomic annotation files, and paste the ASVs seq as rownames
taxo_tmp <- cbind(taxo, ASV_ID = as.numeric(as.character(sub("ASV", "", rownames(taxo)))))
taxo <- taxo_tmp[c(paste0("ASV",sort(taxo_tmp$ASV_ID, decreasing = F))),]
## re-order the other taxo files (based on ASV_ID which is identical between files...)
taxo85 <- taxo85[c(paste0("ASV",sort(taxo_tmp$ASV_ID, decreasing = F))),]
taxo02 <- taxo02[c(paste0("ASV",sort(taxo_tmp$ASV_ID, decreasing = F))),]
taxo_silva <- taxo_silva[c(paste0("ASV",sort(taxo_tmp$ASV_ID, decreasing = F))),]
taxo_silva80 <- taxo_silva80[c(paste0("ASV",sort(taxo_tmp$ASV_ID, decreasing = F))),]
## from factors to character.... 
taxo$ASV_ID <- taxo85$ASV_ID <- taxo02$ASV_ID <- taxo_silva$ASV_ID <- taxo_silva80$ASV_ID <- rownames(taxo)
## for each for the taxon vector
taxo$taxon <- as.character(taxo$taxon)
taxo85$taxon <- as.character(taxo85$taxon)
taxo02$taxon <- as.character(taxo02$taxon)
taxo_silva$taxon <- as.character(taxo_silva$taxon)
taxo_silva80$taxon <- as.character(taxo_silva80$taxon)
## pasting the seq as rownames (sequences are in order ASV1 to ASVn, and otu table also)
identical(dna_seq, colnames(asv)) # TRUE
rownames(taxo) <- rownames(taxo85) <- rownames(taxo02) <- rownames(taxo_silva) <- rownames(taxo_silva80) <- colnames(asv)
# export taxo02 and taxo (95%..)
saveRDS(taxo02, "output/taxo02_annotations.rds")
saveRDS(taxo,   "output/taxo95_annotations.rds")

## adjust the order of samples in the ASVs table
## keep only samples from compo file that are in the OTU table (regardless of the order)
table(rownames(asv) == comp[,s])
table(rownames(asv) %in% comp[,s])

### if some not print them
rownames(asv)[!rownames(asv) %in% comp[,s]]

# keep compo for samples of the table
comp <- subset(comp, comp[,s] %in% rownames(asv))
# do a vector of character to fetch them
samples_comp <- as.character(comp[,s])

# then check order
if(table(samples_comp == rownames(asv))["TRUE"] != nrow(comp) || table(samples_comp == rownames(asv))["FALSE"] == nrow(comp))
{
  message("The order of samples does not match with the metadata file, trying to sort it...")
  asv <- asv[samples_comp,]
  if(table(samples_comp == rownames(asv))["TRUE"] == nrow(comp))
  {
    message("Ok, done...")
  } else { stop("Could not sort") }
} else { message("Order of samples is fine") }
## check exact match
table(rownames(asv) == comp[,s]) # TRUE

########################################################################################################
## Now we have matching table and metadata
########################################################################################################

## We remove the sediment samples shallower than 1900 (few coastal stations of eDNAbyss, but one in the north atlantic at 1920, so we keep this one only...)
# get the list of sediment samples shallower than 1900
tmp <- rownames(subset(comp, comp$biome == "Sediment" & comp$depth < 1900))
## show them 
comp[tmp,c("dataset", "basin", "depth")]

## remove them from asv and comp 
asv <- asv[!rownames(asv) %in% tmp,]
asv <- asv[,colSums(asv)>0]  
comp <- comp[rownames(asv),]

# check 
identical(rownames(asv), rownames(comp)) # TRUE
table(comp$dataset)


##########################################################################################
##########################################################################################
############### RAW seq depth

pdf("output/seq_depth.pdf", useDingbats = F)
row_S <- rowSums(asv)
coll <- paste0(comp$dataset, comp$molecule)
under <- length(subset(row_S, row_S < 10000))
plot(log(row_S+1), ylab = "Sequencing depth - log (Abundance)", yaxt="n", xlab = "Samples ID", pch=1,
     cex=0.5, col= as.numeric(as.factor(coll)), main=paste0(under, " out of ", length(row_S), " samples under 10000 reads"))
axis(2, at=log(c(1,10,100,1000,10000,100000,1000000)),labels=c(0,10,100,1000,10000,100000,1000000),cex.axis=0.7, las=2)
legend("bottomright", legend = unique(coll), col=unique(as.numeric(as.factor(coll))), pch=1, box.lty=0)
abline(h=log(10000+1), col="blue")
dev.off()

##########################################################################################
##########################################################################################
### Distribution bacteria / archaea / eukaryotes / unassigned
### plot seq length
identical(as.character(sread(dna)), rownames(taxo)) # TRUE

### get kingdoms infos from SILVA 80 and by length
bac   <- rownames(taxo_silva80[grep("Bacteria", taxo_silva80$taxon),]) ## SILVA is more extensive for prokaryotes
arc   <- rownames(taxo_silva80[grep("Archaea", taxo_silva80$taxon),])  ## SILVA is more extensive for prokaryotes
no18S <- rownames(subset(taxo02, taxo02$mean.similarity < 0.20))       ## anything that loosely match anything in PR2 is removed (because primers also match V9 of mitochondriale 16S and GC rich genomes)
euk <- rownames(taxo_silva80[grep("Eukaryota", taxo_silva80$taxon),])                               ## just for the plot to set the length cutoff here
una <- rownames(taxo_silva80[grep("unassigned", taxo_silva80$taxon, ignore.case = F, fixed = T),])  ## just for the plot to set the length cutoff here

# df for ggplot
hist_all <- data.frame(Kingdom = c(rep("Bacteria", length(bac)), rep("Archaea", length(arc)),
                                   rep("Eukaryotes", length(euk)), rep("Unassigned", length(una)), rep("no18S", length(no18S))),
                       Length = c(nchar(bac), nchar(arc), nchar(euk), nchar(una), nchar(no18S)))

pdf("output/Seq_length_kingdoms_cutoff.pdf", width= 10, height= 4, useDingbats = F)
ggplot(hist_all, aes(x=Length, color=Kingdom, fill=Kingdom)) +
  geom_histogram(binwidth = 2, position="dodge", alpha=0.2)  +
  geom_vline(xintercept=116, color="blue",linetype="dashed", size=0.2) +
  #geom_vline(xintercept=154, color="blue",linetype="dashed", size=0.2) +
  labs(title="ASVs length distribution - cutoff at 116 bp for the unassigned")
dev.off()

## gather the garbage and setting the length cutoff at 116bp for the unassigned
garbage <- unique(c(bac, arc, no18S, una[nchar(una) < 116]))
### how many garbage? 
message(length(garbage), " garbage ASVs / ", ncol(asv), " ASVs - ", round(length(garbage) / ncol(asv) * 100, 2), "%")
# 300513 garbage ASVs / 430212 ASVs - 69.85%

##########################################################################################
##########################################################################################
###### documenting the distribution bacteria / eukaryote across datasets before removing the garbage
## first aggreg by taxa --> see routines scripts (companion_functions.R)
agg_tax_all <- aggregByTaxa(asv, taxo = taxo_silva80[colnames(asv),], taxo_rank = 1)

## bind the dataset info before melting 
agg_tax_all_facet <- cbind(agg_tax_all$norm, dataset = comp$dataset, samples = comp$sampleID)
dat <- reshape2::melt(agg_tax_all_facet, measure.vars = colnames(agg_tax_all_facet)[1:4])
colnames(dat) <- c("Dataset", "Samples", "Taxa", "Proportion")
gg1 <- ggplot(data=dat, aes(x=Samples, y=Proportion, fill=Taxa)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
  facet_wrap(. ~ Dataset, scale="free", strip.position="bottom") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(title = "Abundance") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

agg_tax_all_bin_facet <- cbind(agg_tax_all$bin, dataset = comp$dataset, samples = comp$sampleID)
dat <- reshape2::melt(agg_tax_all_bin_facet, measure.vars = colnames(agg_tax_all_bin_facet)[1:4])
colnames(dat) <- c("Dataset", "Samples", "Taxa", "Proportion")
gg2 <- ggplot(data=dat, aes(x=Samples, y=Proportion, fill=Taxa)) +  
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
  facet_wrap(. ~ Dataset, scale="free", strip.position="bottom") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=1)) + labs(title = "Richness") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

#### plot the grid into pdf
pdf("output/SILVA_80_AbundRichness_rank1_allSamples.pdf", width = 12, height=7)
grid.arrange(gg1 + theme(legend.position="none"), gg2, widths = c(2.7), heights = c(1, 1))
dev.off()

##########################################################################################
##########################################################################################
###### filtering out prokaryotes and garbage from putative eukaryotes-only ASVs

asv_raw <- asv                ## keeping the raw in mem
# now remove the FIRST garbage
asv <- asv[,!colnames(asv) %in% garbage]  ## asv is the new table to analyse 

### anything loosely matching (best match on PR2) prokaryotes or organelle is also removed
table(sapply(strsplit(taxo02[colnames(asv),"taxon"], ";"), `[`, 1))
## Archaea  Bacteria Eukaryota Organelle 
## 398      9127    277072    375 

### second garbage
tmp <- taxo02[colnames(asv),]
bac2 <- rownames(tmp[grep("Bacteria", tmp$taxon),]) 
arc2 <- rownames(tmp[grep("Archaea", tmp$taxon),])  
org  <- rownames(tmp[grep("Organelle", tmp$taxon),])  

## and keep only ASVs starting with GTCG in the first 4 nt (highly conserved in euks)
WITHmotif <- colnames(asv)[grep("GTCG", substr(colnames(asv), 1, 4))]
NOmotif <-  colnames(asv)[!colnames(asv) %in% WITHmotif]

## collecting second garbage
garbage2 <- unique(c(bac2, arc2, org, NOmotif))

# now remove the garbage
asv <- asv[,!colnames(asv) %in% garbage2]  ## asv is the new table to analyse 

### all garbage represent 
message(round(100 - sum(asv) * 100 / sum(asv_raw), 2), "% of the dataset")   ## 13,12% 

## table features
dim(asv) ## 1944 242'485 
sum(asv) ## 1'958'402'792 

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

## subset per biome
euphotic_all <- subset(asv, comp$biome == "Euphotic")
aphotic_all  <- subset(asv, comp$biome == "Aphotic")
sediment_all <- subset(asv, comp$biome == "Sediment")

#### get ASVs sequences of each biomes 
asvs_euph <- colnames(euphotic_all[,colSums(euphotic_all) > 0])
asvs_aph <- colnames(aphotic_all[,colSums(aphotic_all) > 0])
asvs_sed <- colnames(sediment_all[,colSums(sediment_all) > 0])

###### and filtering plankton from sediment (benthic) / plankton ones in sediment (sinking plankton)
asvs_plankton <- unique(c(asvs_euph, asvs_aph))
asvs_sed_only <- asvs_sed[!asvs_sed %in% asvs_plankton]
asvs_plan_sed <- asvs_sed[asvs_sed %in% asvs_plankton]

##### ANDY / JAN / TRISTAN --> manual curation of metazoan in sinking plankton (benthic species with planktonic larval stages...)
# ## get planktonic metazoan and export file for manual curation 
# tt <- taxo85[asvs_plan_sed,]
# metaz <- grep("Metazoa", tt$taxon)
# metaz_plan_sed <- tt[metaz,]
# metaz_plan_sed$total_abundance_sediment <- colSums(sediment_all[,rownames(metaz_plan_sed)]) 
# metaz_plan_sed$total_abundance_euphotic <- colSums(euphotic_all[,rownames(metaz_plan_sed)]) 
# metaz_plan_sed$total_abundance_aphotic <- colSums(aphotic_all[,rownames(metaz_plan_sed)]) 
# metaz_plan_sed$sequences <- rownames(metaz_plan_sed)
# rownames(metaz_plan_sed) <- 1:nrow(metaz_plan_sed)
# write.table(metaz_plan_sed, "Planktonic_metazoa_in_sediment.tsv", quote = F, row.names = T, sep="\t")

## MANUAL CURATION OF METAZOAN to reclassify them as benthic / planktonic
# import manual curation results
curate_metazoan <- read.table("manual_curation_metazoa_plank_bent.txt", header = T, sep = "\t")
rownames(curate_metazoan) <- curate_metazoan$sequences
table(curate_metazoan$bent_plank)
# Benthic Contaminants   Planktonic  Terrestrial 
# 224     3              302         17 

add_asvs_sed_only <- rownames(subset(curate_metazoan, curate_metazoan$bent_plank == "Benthic"))
contaminants <- rownames(subset(curate_metazoan, curate_metazoan$bent_plank == "Contaminants" | curate_metazoan$bent_plank == "Terrestrial"))

## now remake the ASVs vectors 
# we remove contaminants from plankton and curated benthic
asvs_plankton <- asvs_plankton[!asvs_plankton %in% contaminants]
asvs_euph     <- asvs_euph[!asvs_euph %in% contaminants]
asvs_aph      <- asvs_aph[!asvs_aph %in% contaminants]

# we add the curated benthic animals to the benthic only
asvs_sed_only <- unique(c(add_asvs_sed_only, asvs_sed_only))
asvs_sed_only <- asvs_sed_only[!asvs_sed_only %in% contaminants]
# remove the curated benthic asvs from the plankton in sediment
asvs_plan_sed <- asvs_plan_sed[!asvs_plan_sed %in% add_asvs_sed_only]
asvs_plan_sed <- asvs_plan_sed[!asvs_plan_sed %in% contaminants]

## remove the contaminants from the datasets
asv <- asv[,!colnames(asv) %in% contaminants]
euphotic_all <- euphotic_all[,!colnames(euphotic_all) %in% contaminants]
aphotic_all  <- aphotic_all[,!colnames(aphotic_all) %in% contaminants]
sediment_all <- sediment_all[,!colnames(sediment_all) %in% contaminants]

# new dimensions
dim(asv) # 1944 242'465
sum(asv) # 1'958'390'988

### aggregate DNA/RNA per extraction for the deep sea dataset ### see below script for figure S12 documenting DNA vs RNA patterns 
asv_ds_agg <- aggregate(asv_ds, by=list(comp_ds$pseudoreplicate), FUN=sum) # <-- **very long** 
rownames(asv_ds_agg) <- asv_ds_agg[,1]
asv_ds_agg <- asv_ds_agg[,2:ncol(asv_ds_agg)]
saveRDS(asv_ds_agg, "output/asv_ds_agg_pseudoreplicates.rds")

## making the comp 
comp_ds_agg <- c()
for (i in rownames(asv_ds_agg)) {
  tmp <- subset(comp_ds, comp_ds$pseudoreplicate == i)
  if (nrow(tmp) > 1) {
    comp_ds_agg <- rbind(comp_ds_agg, tmp[1,])
  } else {
    comp_ds_agg <- rbind(comp_ds_agg, tmp)
  }
}
## check 
identical(as.character(comp_ds_agg[,"pseudoreplicate"]), rownames(asv_ds_agg)) # TRUE
rownames(comp_ds_agg) <- as.character(comp_ds_agg$pseudoreplicate)

## remerge the deep sea dataset to the rest 
asv_others <- subset(asv, comp$dataset != "deep_sea")
comp_other <- comp[rownames(asv_others),]

## and merge 
asv_agg <- rbind(asv_ds_agg, asv_others)
comp_agg <- rbind(comp_ds_agg, comp_other)
# check 
identical(rownames(comp_agg), rownames(asv_agg)) # TRUE

# new dimensions
dim(asv_agg) # 1685 242'465
sum(asv_agg) # 1'958'390'988

# export ASVs vectors 
saveRDS(asvs_plankton, "output/asvs_plankton.rds")
saveRDS(asvs_euph, "output/asvs_euph.rds")
saveRDS(asvs_aph, "output/asvs_aph.rds")
saveRDS(asvs_sed_only, "output/asvs_sed_only.rds")
saveRDS(add_asvs_sed_only, "output/add_asvs_sed_only.rds")
saveRDS(asvs_plan_sed, "output/asvs_plan_sed.rds")

# and tables
## subset asv_agg per biome (but keep all ASVs columns)
euphotic_all <- subset(asv_agg, comp_agg$biome == "Euphotic")
aphotic_all  <- subset(asv_agg, comp_agg$biome == "Aphotic")
sediment_all <- subset(asv_agg, comp_agg$biome == "Sediment")
# export
saveRDS(asv, "output/asv_table_euks_clean.rds")
saveRDS(asv_agg, "output/asv_table_euks_clean_agg.rds")
# mapped metadata
saveRDS(comp_agg, "output/comp_table_euks_clean.rds")
# per biome (with all ASVs of the full dataset)
saveRDS(euphotic_all, "output/asv_table_euphotic_all.rds")
saveRDS(aphotic_all, "output/asv_table_aphotic_all.rds")
saveRDS(sediment_all, "output/asv_table_sediment_all.rds")


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

#### Gathering ASVs meaningfull taxogroups, functional annotations and sizes

### get the curated meaningful taxo super-groups list 
taxo_groups <- read.table("Ref-V9-63400-taxogroups_noQuote.tsv", header = T, sep = "\t")
rownames(taxo_groups) <- taxo_groups$name

#### and redo the consensus between the 3 candidates 
taxo_group_asvs <- taxo85[colnames(asv),]
taxo_group_asvs$taxo_group <- "consensus"

# loop over them 
for (i in 1:nrow(taxo_group_asvs)) {
  # get the candidates
  tmp <- unlist(strsplit(as.character(taxo_group_asvs[i,"reference.ids"]), ";"))
  # if more than 1
  if (length(tmp) > 1) {
    # get the taxo groups for those
    tmp_groups <- data.frame(taxo_groups[tmp,c("supergroup", "taxogroup1", "taxogroup2")], stringsAsFactors = F)
    ## if same taxogroup, ok, otherwise, upper test one
    if (length(unique(tmp_groups$taxogroup2)) == 1) {
      taxo_group_asvs$taxo_group[i] <- paste(as.character(tmp_groups[1,1]), as.character(tmp_groups[1,2]), as.character(tmp_groups[1,3]), sep = ";")
    } else if (length(unique(tmp_groups$taxogroup1)) == 1) {
      taxo_group_asvs$taxo_group[i] <- paste(as.character(tmp_groups[1,1]), as.character(tmp_groups[1,2]), "Unassigned_3", sep = ";")
    } else if (length(unique(tmp_groups$supergroup)) == 1) {
      taxo_group_asvs$taxo_group[i] <- paste(as.character(tmp_groups[1,1]), "Unassigned_2", "Unassigned_3", sep = ";")
    } else {
      taxo_group_asvs$taxo_group[i] <- "Unassigned_1;Unassigned_2;Unassigned_3"
    }
  } else if (length(tmp) == 1) {
    tmp_groups <- data.frame(taxo_groups[tmp,c("supergroup", "taxogroup1", "taxogroup2")], stringsAsFactors = F)
    taxo_group_asvs$taxo_group[i] <- paste(as.character(tmp_groups[1,1]), as.character(tmp_groups[1,2]), as.character(tmp_groups[1,3]), sep = ";")
  } else {
    #message("Unassigned, line ", i)
    taxo_group_asvs$taxo_group[i] <- "Unassigned_1;Unassigned_2;Unassigned_3"
  }
  if (i %% 10000 == 0) message(i)
}

table(taxo_group_asvs$taxo_group)
## 4 NA;NA;NA  <- must be because we had non consensual candidates for the first rank (Euk + Bac ?)
taxo_group_asvs[rownames(subset(taxo_group_asvs, taxo_group_asvs$taxo_group == "NA;NA;NA")),"taxo_group"] <- "Unassigned_1;Unassigned_2;Unassigned_3"

### functional annotations for each planktonic ASV
funct_db <- read.table("refv9_corresp_v20171106_with_functions_R.tsv", header=TRUE, sep="\t")
taxo_annot <- taxo[asvs_plankton,] ## <--- using the 95 % similarity annotation to be more confident 
taxo_annot85 <- taxo85[asvs_plankton,] ## <--- but will add the annotation at 85% similarity (for taxo group that are based on that)
rownames(funct_db) <- funct_db$name

out <- c()
for (i in 1:nrow(taxo_annot)) {
  tmp <- taxo_annot[i,"reference.ids"]
  tmp <- strsplit(as.character(tmp), ";")[[1]]
  annot <- c()
  if (length(tmp) > 0) {
    if (!is.na(tmp)) {
      tmp_annot <- funct_db[tmp, c("chloroplast", "symb_small", "symb_host", "silicification", "calcification", "strontification")]
      for (j in 1:ncol(tmp_annot)) {
        if (length(unique(tmp_annot[,j])) == 1 & as.character(unique(tmp_annot[,j])) != "") {
          annot <- c(annot, as.character(tmp_annot[1,j]))
        } else {
          annot <- c(annot, NA)
        }
      } 
    } else  { annot <- rep(NA, 6) } 
  } else  {
    annot <- rep(NA, 6)
  } 
  out <- rbind(out, annot)
  if (i %% 5000 == 0) message(i)
}

colnames(out) <- c("chloroplast", "symb_small", "symb_host", "silicification", "calcification", "strontification")
### planktonic ASVs functional annotations
asvs_funct <- cbind(taxo_annot, out, taxo_groups = taxo_group_asvs[asvs_plankton,"taxo_group"])

### make a vector of funational groups
asvs_funct$functional_groups <- c()
asvs_funct[unique(c(grep("no", asvs_funct$symb_host), grep("no", asvs_funct$chloroplast))),"functional_groups"] <- "Other heterotrophs"
asvs_funct[grep("yes", asvs_funct$chloroplast),"functional_groups"] <- "Phototrophs"
asvs_funct[grep("photo", asvs_funct$symb_host),"functional_groups"] <- "Photosymbionts"
asvs_funct[grep("Metazoa", asvs_funct$taxo_groups),"functional_groups"] <- "Other metazoa"  ## <-- less stringent to ascribe copepods or other metazoa, taxo groups made at 85% similarity
asvs_funct[grep("Copepoda", asvs_funct$taxo_groups),"functional_groups"] <- "Copepoda"      ## <-- less stringent to ascribe copepods or other metazoa, taxo groups made at 85% similarity
asvs_funct[grep("parasite", asvs_funct$symb_small),"functional_groups"] <- "Parasite"      ## Overwrite if previously ascribed to another group
asvs_funct$functional_groups[is.na(asvs_funct$functional_groups)] <- "Unknown"
#### 
table(asvs_funct$functional_groups)
#Copepoda    Other heterotrophs  Other metazoa    Parasite     Photosymbionts     Phototrophs       Unknown 
#7346        21932               9104             4808         3099               8493              88190

#### and map all the infos to the taxo_group_asvs file that contains all ASVs
taxo_group_asvs$functional_groups <- c()
taxo_group_asvs[asvs_plankton,"functional_groups"] <- asvs_funct$functional_groups

## sinking plankton?
taxo_group_asvs$Sinking <- ifelse(rownames(taxo_group_asvs) %in% asvs_plan_sed, "yes", "no")

## inferring planktonic ASVs sizes
asv_size <- subset(asv_agg[,asvs_plankton], comp_agg[rownames(asv_agg),"lower_size_.µm"] != 0 & comp_agg[rownames(asv_agg),"biome"] == "Euphotic")
comp_size <- comp_agg[rownames(asv_size),]
## total sum normalization
asv_size <- decostand(asv_size, method="total")

asvs_sizes <- cbind(ASVs = asvs_plankton, Size = NA)
rownames(asvs_sizes) <- asvs_plankton

## weighted average
for (i in asvs_plankton) {
  tmp <- sum(asv_size[,i] * comp_size$lower_size_.µm) / sum(asv_size[,i])
  #if (tmp == "NaN") message(sum(otu_size[,i]), comp_size$lower_size_µm, sum(otu_size[,i])) ### the ones that are not computable are the aphotic ASVs for which no size fraction... 
  asvs_sizes[i,"Size"] <- tmp
  #if (i %% 1000 == 0) message(i)
}

#### and map all the sizes 
taxo_group_asvs$sizes <- c()
taxo_group_asvs[asvs_plankton,"sizes"] <- asvs_sizes[,"Size"]

### export the taxo_group_asvs file 
saveRDS(taxo_group_asvs, "output/asvs_taxo_functional_annotations.rds")

## check other scripts for data anlaysis and figures
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

### Figure S12
taxo_group_asvs <- readRDS("output/asvs_taxo_functional_annotations.rds")

### Focus on deep sea dataset DNA / RNA diversity and composition patterns
asv_ds <- subset(asv, comp$dataset == "deep_sea")
asv_ds <- subset(asv_ds, rowSums(asv_ds) >= 10000)
asv_ds <- asv_ds[,colSums(asv_ds)>0]
comp_ds<- comp[rownames(asv_ds),]

## Venn between DNA / RNA
df <- data.frame(DNA = colSums(asv_ds[comp_ds$molecule == "DNA",]),
                 RNA = colSums(asv_ds[comp_ds$molecule == "RNA",]),
                 Source = ifelse(colnames(asv_ds) %in% asvs_plan_sed, "Plankton", "Benthos"))

ven <- limma::vennCounts(df[,1:2])
pdf("output/Figure_S12A_venn_diagram.pdf")
limma::vennDiagram(ven, lwd=0.9)
dev.off()

## What proportion of reads represent shared ASVs between DNA and RNA
df_bin <- df[,1:2]
df_bin[df_bin > 0] <- 1
asvs_dnarna <- rownames(df_bin[rowSums(df_bin) == 2,])
message("Common ASVs in DNA and RNA represent ", round(sum(df[asvs_dnarna,1:2]) *100 / sum(df[,1:2]), 2), " % of the total reads") ## 89.82%

asvs_dnaonly <- rownames(df_bin[df_bin[,"DNA"] == 1,])
asvs_dnaonly <- asvs_dnaonly[!asvs_dnaonly %in% asvs_dnarna]
message("DNA-only ASVs represent ", round(sum(df[asvs_dnaonly,1:2]) *100 / sum(df[,1:2]), 2), " % of the total reads") ## 6.33%

asvs_rnaonly <- rownames(df_bin[df_bin[,"RNA"] == 1,])
asvs_rnaonly <- asvs_rnaonly[!asvs_rnaonly %in% asvs_dnarna]
message("DNA-only ASVs represent ", round(sum(df[asvs_rnaonly,1:2]) *100 / sum(df[,1:2]), 2), " % of the total reads") ## 3.85%

## Plankton vs Benthic in DNA/RNA
df_norm <- df
df_norm <- aggregate(df_norm[,c(1:2)], by = list(as.character(df_norm$Source)), FUN = sum)
colnames(df_norm)[1] <- "Source"
df_norm_melt <- reshape2::melt(df_norm, measure.vars = colnames(df_norm)[2:3])

# get the %
p_title <- paste0("Plankton in DNA:", round(df_norm_melt[2,3] *100 / sum(df_norm_melt[1:2,3]), 1), "%\nPlankton in RNA:", round(df_norm_melt[4,3] *100 / sum(df_norm_melt[3:4,3]), 1), "%") 
p_source <- ggplot(data=df_norm_melt, aes(x=variable, y=value, fill=Source)) +  
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
  labs(title = p_title) + coord_flip()
ggsave("output/Figure_S12B_DNA_RNA_source.pdf", p_source, width = 5, height = 2)

## are compositional patterns different? 
asv_ds_norm <- cssNorm(asv_ds, samples_as_rows = T)
asv_ds_bray <- vegan::vegdist(asv_ds_norm, distance = "bray")

## Permanova: is interaction region x molecule significant? yes 
permanova <- vegan::adonis(asv_ds_bray ~ comp_ds$watling * comp_ds$molecule, permutations = 999, parallel = 8)

## PCOA 
pcoa_df <- ape::pcoa(asv_ds_bray)
PCOA <- pcoa_df$vectors[,1:2]
var_axis_pcoa <- round(pcoa_df$values$Relative_eig * 100, 2)[1:2]

comp_ds <- cbind(comp_ds, PCOA)

pcoa_dnarna <- ggplot(comp_ds, aes(x = Axis.1, y = Axis.2, color = watling, shape = molecule)) +
  geom_point() + 
  labs(x= paste0(var_axis_pcoa[1],"%"),  y = paste0(var_axis_pcoa[2], "%"), 
       title = paste0("Provinces - R2: ", round(permanova$aov.tab[1,"R2"], 3), " p < 0.001*** \n", 
                      "Mol. - R2: ", round(permanova$aov.tab[2,"R2"], 3), " p < 0.001*** \n", 
                      "Provinces x Mol. - R2: ", round(permanova$aov.tab[3,"R2"], 3), " p < 0.001***"))


ggsave("output/Figure_S12C_PCOA_DNA_RNA_watling.pdf", pcoa_dnarna, width = 7, height = 6, device = cairo_pdf)


# and treemap 
## Aggreg by groups
tmp_taxo <- taxo_group_asvs[rownames(df),]

tree_df_benthic <- aggregate(df[df$Source == "Benthos",1:2], by = list(tmp_taxo[rownames(df[df$Source == "Benthos",1:2]),"taxo_group"]), FUN = sum)
tree_df_benthic$Group.1 <- gsub("NA;NA;NA", "Unassigned;Unassigned;Unassigned", tree_df_benthic$Group.1)

tree_df_plankton <- aggregate(df[df$Source == "Plankton",1:2], by = list(tmp_taxo[rownames(df[df$Source == "Plankton",1:2]),"taxo_group"]), FUN = sum)
tree_df_plankton$Group.1 <- gsub("NA;NA;NA", "Unassigned;Unassigned;Unassigned", tree_df_plankton$Group.1)

## getting the structure back
taxo_benthic <- data.frame(rank1 = sapply(strsplit(tree_df_benthic$Group.1, ";"), `[`, 1),
                           rank2 = sapply(strsplit(tree_df_benthic$Group.1, ";"), `[`, 2), 
                           rank3 = sapply(strsplit(tree_df_benthic$Group.1, ";"), `[`, 3), stringsAsFactors=FALSE)
taxo_benthic <- cbind(taxo_benthic, tree_df_benthic)

taxo_plankton <- data.frame(rank1 = sapply(strsplit(tree_df_plankton$Group.1, ";"), `[`, 1),
                            rank2 = sapply(strsplit(tree_df_plankton$Group.1, ";"), `[`, 2), 
                            rank3 = sapply(strsplit(tree_df_plankton$Group.1, ";"), `[`, 3), stringsAsFactors=FALSE)
taxo_plankton <- cbind(taxo_plankton, tree_df_plankton)

## focus on groups
focus <- c("Metazoa", "Rhizaria", "Alveolata", "Stramenopiles", "Euglenozoa", "Fungi", "Amoebozoa")
taxo_benthic$rank1 <- ifelse(taxo_benthic$rank1 %in% focus, taxo_benthic$rank1, "Others")
taxo_plankton$rank1 <- ifelse(taxo_plankton$rank1 %in% focus, taxo_plankton$rank1, "Others")

# melt
taxo_benthic <- reshape2::melt(taxo_benthic, idvar = c("DNA", "RNA"))
taxo_plankton <- reshape2::melt(taxo_plankton, idvar = c("DNA", "RNA"))


## treemap 
library(treemapify)
p_taxo_benthic <- ggplot(data=taxo_benthic, aes(area = value, fill = rank1, subgroup = rank1, subgroup2 = rank2, subgroup3 = rank3)) +
  geom_treemap(layout = "squarified", start = "topleft") +
  geom_treemap_subgroup_border(start = "topleft", colour = "black", size = 1.5) +
  geom_treemap_subgroup2_border(start = "topleft", colour = "black", size = 1.5) +
  geom_treemap_subgroup3_border(start = "topleft", colour = "white", size = 0.5) +
  #geom_treemap_subgroup_text(start = "topleft", place = "topleft", reflow = T, alpha = 1, colour =
  #                             "black", fontface = "bold", min.size = 4, size = 13) +
  geom_treemap_subgroup2_text(start = "topleft", place = "topleft", grow = F, reflow = T, alpha = 1, colour =
                                "black", fontface = "bold", min.size = 2, size = 8) + 
  geom_treemap_subgroup3_text(start = "topleft", place = "center", grow = F, reflow = T, alpha = 1, colour =
                                "black", fontface = "italic", min.size = 2, size = 7) +
  scale_fill_brewer(palette="Set3") + labs(title = "Benthic reads proportion") +  
  facet_wrap(~ variable, ncol = 1) 

p_taxo_plankton <- ggplot(data=taxo_plankton, aes(area = value, fill = rank1, subgroup = rank1, subgroup2 = rank2, subgroup3 = rank3)) +
  geom_treemap(layout = "squarified", start = "topleft") +
  geom_treemap_subgroup_border(start = "topleft", colour = "black", size = 1.5) +
  geom_treemap_subgroup2_border(start = "topleft", colour = "black", size = 1.5) +
  geom_treemap_subgroup3_border(start = "topleft", colour = "white", size = 0.5) +
  #geom_treemap_subgroup_text(start = "topleft", place = "topleft", reflow = T, alpha = 1, colour =
  #                             "black", fontface = "bold", min.size = 4, size = 13) +
  geom_treemap_subgroup2_text(start = "topleft", place = "topleft", grow = F, reflow = T, alpha = 1, colour =
                                "black", fontface = "bold", min.size = 2, size = 8) + 
  geom_treemap_subgroup3_text(start = "topleft", place = "center", grow = F, reflow = T, alpha = 1, colour =
                                "black", fontface = "italic", min.size = 2, size = 7) +
  scale_fill_brewer(palette="Set3") + labs(title = "Plankton reads proportion in sediment") +  
  facet_wrap(~ variable, ncol = 1) 

p_tree <- p_taxo_plankton / p_taxo_benthic + plot_layout(guides = "collect")
ggsave("output/Figure_S12B_DNA_RNA_treemap.pdf", p_tree, width = 10, height = 8)








