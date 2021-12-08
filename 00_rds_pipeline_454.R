### baobab pipeline
path <- '/path/to/paired_ends_fastqs/'
setwd(path)

## Cores to use
CORES <- 16
library('doMC')
registerDoMC(cores = CORES)
library(ShortRead)
library(dada2)
library(Biostrings)
library(tictoc)

### some functions
source("/home/cordiert/Rfunctions/fastqUtils.R")

## primers sequences in 5'-3' orientation 
FWD <- "GTACACACCGCCC"  ### <- for the 454 raw data from Pawlowski et al., 2011, the fwd primers is incomplete
REV <- "CCTTCYGCAGGTTCACCTAC"

## trimming primers from the reads with custom R function
fastqq <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))
for (i in 1:length(fastqq)) fastqTrimBothPrimer(FWD, REV, fastqR1In = fastqq[i])

## path for the fastq files
fnFs <- sort(list.files(path, pattern="trim.fastq.gz", full.names = TRUE))
## adjust to your need
sample.names <- paste0(sapply(strsplit(basename(fnFs), "_"), `[`,1))

### create directories
path.tmp1 <- file.path(path, "tmp1")
if(!dir.exists(path.tmp1)) dir.create(path.tmp1)
path.tmp2 <- file.path(path, "tmp2")
if(!dir.exists(path.tmp2)) dir.create(path.tmp2)
path.rds <- file.path(path, "output_rds")
if(!dir.exists(path.rds)) dir.create(path.rds)

####################################################################################
### starting the workflow
filtFs <- file.path(path.tmp1, basename(fnFs))

# basic filtering with default DADA2 settings
filtering <- filterAndTrim(fwd=fnFs, filt=filtFs, multithread=TRUE, verbose=TRUE)

# collecting summary stats
readsCounts <- array(NA, c(nrow(filtering), 6))
colnames(readsCounts) <- c("reads.in", "reads.out", "orient", "cutadapt", "primer.free", "length")
rownames(readsCounts) <- sample.names
readsCounts[,1] <- filtering[,1]
readsCounts[,2] <- filtering[,2]

# keep only samples with >100 reads
keep <- filtering[,"reads.out"] > 100

# paths list update
filtFs_kept <- filtFs[keep]

## the remaining reads that would contain primers are artefacts (allowing 1 mismatch)
fnFs.clean <- file.path(path.tmp2, basename(filtFs_kept))

## FWD and REV are passed to the function but the RC of both is performed within 
par_primers <- foreach(i = 1:length(fnFs.clean)) %dopar% fastqFilterPrimersMatchs(fwdPrimer = FWD, revPrimer = REV, 
                                                                                  fastqR1In = filtFs_kept[i],  
                                                                                  fastqR1Out = fnFs.clean[i], maxMismatch = 1, sampleName = sample.names[i])
## collect stats
for (i in 1:length(par_primers)) readsCounts[par_primers[[i]]$sampleID,"primer.free"] <- par_primers[[i]]$readsCount

## remove the fastq in the tmp1
system("rm tmp1/*.fastq.gz")

## some may be too short now..
fnFs.length <- file.path(path.tmp1, basename(fnFs.clean))
#fnRs.length <- file.path(path.tmp1, basename(fnRs.cut))

par_length <- foreach(i = 1:length(fnFs.length)) %dopar% fastqFilterWidth(fastqR1In = fnFs.clean[i],   
                                                                          fastqR1Out= fnFs.length[i], sampleName = sample.names[i])
## collect stats
for (i in 1:length(par_length)) readsCounts[par_length[[i]]$sampleID,"length"] <- par_length[[i]]$readsCount

## export the summary stats of filtering steps
write.table(readsCounts, file = "readsCounts_stats.tsv", quote = F, sep="\t", row.names = F, fileEncoding = "UTF-8")


###############
# now we have clean, oriented and primer trimmed fastq

# collect summary stats
# store last error estimates after the 10 rounds
lastE <- array(NA, c(length(fnFs.length), 4))
colnames(lastE) <- c("sample_id", "lastE_FWD", "lastE_REV", "duration")


for (i in seq_along(fnFs.length))
{
  tic()
  sampleName <- strsplit(basename(fnFs.length[[i]]), "_")[[1]][1]
  drpF <- derepFastq(fnFs.length[[i]])
  # error model
  errF <- learnErrors(drpF, multithread=CORES, randomize=TRUE)
  ## save the error model for both fwd and rev
  saveRDS(errF, file.path(path.rds, paste0(sampleName, "_errF.rds")))
  ## and the plot
  pdf(file.path(path.rds, paste0(sampleName, "_plot_errF.pdf")))
  print(plotErrors(errF, nominalQ = T))
  dev.off()
  ## denoising with advised parameters for 454 --> https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
  ddF <- dada(drpF, err=errF, selfConsist=F, multithread=CORES, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
  # export the merger file
  saveRDS(ddF, file.path(path.rds, paste0(sampleName, "_merger.rds")))
  ## get the last errors estimates if not converged
  lastE[i,"sample_id"] <- sampleName
  lastE[i,"lastE_FWD"] <- tail(dada2:::checkConvergence(errF),1)
  tt <- toc()
  lastE[i,"duration"] <- tt$toc - tt$tic
}

### export the stats
write.table(lastE, file = "dada2_stats.tsv", quote = F, sep="\t", row.names = F, fileEncoding = "UTF-8")



