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
FWD <- "TTGTACACACCGCCC"  
REV <- "CCTTCYGCAGGTTCACCTAC" 

### cutadapt path (or if it is in the PATH, like here) 
cutadapt <- "cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R

### first path for the fastq files
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2", full.names = TRUE))
# adjust to your need
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`,1)

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
filtRs <- file.path(path.tmp1, basename(fnRs))

# basic filtering with default DADA2 settings
filtering <- filterAndTrim(fwd=fnFs, rev=fnRs, filt=filtFs, filt.rev=filtRs, 
                           multithread=TRUE, verbose=TRUE, orient.fwd=FWD)

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
filtRs_kept <- filtRs[keep]

# cutadapt for trimming the fwd primer and rev primer (if already trimmed, comment these lines), removing reads for which FWD not in R1 or/and REV not in R2
fnFs.cut <- file.path(path.tmp2, basename(filtFs_kept))
fnRs.cut <- file.path(path.tmp2, basename(filtRs_kept))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
foreach(i = 1:length(fnFs.cut)) %dopar% system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "--discard-untrimmed", "-pair-filter=both", # -n required to remove FWD and REV from reads
                                                                   "-o", fnFs.cut[i], "-p", fnRs.cut[i],  # "-e", 0.4, # output files
                                                                   filtFs_kept[i], filtRs_kept[i])) # input files

# count reads after cutadapt
par_countCuta <- foreach(i = 1:length(fnFs.cut)) %dopar% fastqCountReads(fastqR1In = fnFs.cut[i], fastqR2In = fnRs.cut[i],
                                                                         sampleName = sample.names[i])
# collect stats
for (i in 1:length(par_countCuta)) readsCounts[par_countCuta[[i]]$sampleID,"cutadapt"] <- par_countCuta[[i]]$readsCount

## remove the fastq in the tmp2
system("rm tmp1/*.gz")

## the remaining reads that would contain primers are artefacts (allowing 1 mismatch)
fnFs.clean <- file.path(path.tmp1, basename(filtFs_kept))
fnRs.clean <- file.path(path.tmp1, basename(filtRs_kept))

## FWD and REV are passed to the function but the search of RC of both is performed within fastqFilterPrimersMatchs function
par_primers <- foreach(i = 1:length(fnFs.clean)) %dopar% fastqFilterPrimersMatchs(fwdPrimer = FWD, revPrimer = REV, 
                                                                                  fastqR1In = fnFs.cut[i], fastqR2In = fnRs.cut[i], 
                                                                                  fastqR1Out = fnFs.clean[i], fastqR2Out = fnRs.clean[i], 
                                                                                  maxMismatch = 1, sampleName = sample.names[i])

## collect stats
for (i in 1:length(par_primers)) readsCounts[par_primers[[i]]$sampleID,"primer.free"] <- par_primers[[i]]$readsCount

## remove the fastq in the tmp1
system("rm tmp2/*.gz")

## some reads are too short...
fnFs.length <- file.path(path.tmp2, basename(fnFs))
fnRs.length <- file.path(path.tmp2, basename(fnRs))

## removing reads <= 20 bp
par_length <- foreach(i = 1:length(fnFs.length)) %dopar% fastqFilterWidth(fastqR1In = fnFs.clean[i],  fastqR2In = fnRs.clean[i], 
                                                                          fastqR1Out= fnFs.length[i], fastqR2Out= fnRs.length[i], 20, 
                                                                          sampleName = sample.names[i])
## collect stats
for (i in 1:length(par_length)) readsCounts[par_length[[i]]$sampleID,"length"] <- par_length[[i]]$readsCount


## remove the fastq in the tmp2
system("rm tmp1/*.gz")

## export the summary stats of filtering steps
readsCounts <- cbind(sample_ID = rownames(readsCounts), readsCounts)
write.table(readsCounts, file = "readsCounts_stats.tsv", quote = F, sep="\t", row.names = F, fileEncoding = "UTF-8")

###############
# now we have clean, oriented and primer trimmed fastqs

# collect summary stats
# store last error estimates after the 10 rounds
lastE <- array(NA, c(length(fnFs.length), 4))
colnames(lastE) <- c("sample_id", "lastE_FWD", "lastE_REV", "duration")


for (i in seq_along(fnFs.length))
{
  tic()
  sampleName <- strsplit(basename(fnFs.length[[i]]), "_")[[1]][1]
  # error model
  errF <- learnErrors(fnFs.length[[i]], multithread=CORES, randomize=TRUE)
  errR <- learnErrors(fnRs.length[[i]], multithread=CORES, randomize=TRUE)
  ## save the error model for both fwd and rev
  saveRDS(errF, file.path(path.rds, paste0(sampleName, "_errF.rds")))
  saveRDS(errR, file.path(path.rds, paste0(sampleName, "_errR.rds")))
  # derep now
  drpF <- derepFastq(fnFs.length[[i]])
  drpR <- derepFastq(fnRs.length[[i]])

  ## export the plots
  pdf(file.path(path.rds, paste0(sampleName, "_plot_errF.pdf")))
  print(plotErrors(errF, nominalQ = T))
  dev.off()
  pdf(file.path(path.rds, paste0(sampleName, "_plot_errR.pdf")))
  print(plotErrors(errR, nominalQ = T))
  dev.off()
  ## denoising
  ddF <- dada(drpF, err=errF, selfConsist=F, multithread=CORES)
  ddR <- dada(drpR, err=errR, selfConsist=F, multithread=CORES)
  # trimOverhang is important here as we have fully overlapping reads
  merger <- mergePairs(ddF, drpF, ddR, drpR, trimOverhang=TRUE, verbose=TRUE)
  # export the merger file
  saveRDS(merger, file.path(path.rds, paste0(sampleName, "_merger.rds")))
  ## get the last errors estimates if not converged
  lastE[i,"sample_id"] <- sampleName
  lastE[i,"lastE_FWD"] <- tail(dada2:::checkConvergence(errF),1)
  lastE[i,"lastE_REV"] <- tail(dada2:::checkConvergence(errR),1)
  tt <- toc()
  lastE[i,"duration"] <- tt$toc - tt$tic
}

### export the summary stats
write.table(lastE, file = "dada2_stats.tsv", quote = F, sep="\t", row.names = F, fileEncoding = "UTF-8")



