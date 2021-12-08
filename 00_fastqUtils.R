
## Utilities for fastq files

################## removing reads that still match any of the primers
fastqFilterPrimersMatchs <- function(fwdPrimer, revPrimer, fastqR1In, fastqR2In = NULL, fastqR1Out, fastqR2Out = NULL, sampleName = NULL, maxMismatch = 0, splitHeader = " ") {
  require(ShortRead)
  require(Biostrings)
  require(dada2)
  # reading both fastq 
  tmpR1 <- readFastq(fastqR1In)
  if (!is.null(fastqR2In)) tmpR2 <- readFastq(fastqR2In)
  # sample name
  if (!is.null(sampleName)) sample.name <- sampleName else sample.name <- strsplit(basename(fastqR1In), ".fastq")[[1]][1]
  
  if (!is.null(fastqR2In)) 
  {
    # IDs are supposed to be matching all along ! 
    if (!identical(sapply(strsplit(as.vector(id(tmpR1)), splitHeader), `[`, 1), sapply(strsplit(as.vector(id(tmpR2)), splitHeader), `[`, 1))) 
    {
      message("not the same order in the fastq...")
      break
    }
  }
  # make reverseComplement of both primers
  fwd <- fwdPrimer
  fwdRC <- dada2:::rc(fwdPrimer)
  rev <- revPrimer
  revRC <- dada2:::rc(rev)
  
  for (p in c(fwd, fwdRC, rev, revRC))
  {
    # check fwd primer in all R1 reads
    idx <- vmatchPattern(p, sread(tmpR1), fixed=FALSE, max.mismatch=maxMismatch)
    # get the indexing
    xx <- elementNROWS(idx)
    tmpR1 <- tmpR1[which(xx==0)]
    if (!is.null(fastqR2In)) tmpR2 <- tmpR2[which(xx==0)]
  }
  if (!is.null(fastqR2In)) 
  {
    for (p in c(fwd, fwdRC, rev, revRC))
    {
      # check fwd primer in all R2 reads
      idy <- vmatchPattern(p, sread(tmpR2), fixed=FALSE, max.mismatch=maxMismatch)
      # get the indexing
      yy <- elementNROWS(idy)
      tmpR1 <- tmpR1[which(yy==0)]
      tmpR2 <- tmpR2[which(yy==0)]
    }
  }
  # export the reads that did not match the primers
  writeFastq(tmpR1, fastqR1Out, "a")
  if (!is.null(fastqR2In)) writeFastq(tmpR2, fastqR2Out, "a")
  return(list("sampleID" = sample.name, "readsCount" = length(tmpR1)))
}

############################################################################################################
############################################################################################################
################## filter by length fastq
fastqFilterWidth <- function(fastqR1In, fastqR2In = NULL, fastqR1Out, fastqR2Out = NULL, widthFilter = 20, sampleName = NULL, splitHeader = " ") {
  require(ShortRead)
  require(Biostrings)
  # reading both fastq 
  tmpR1 <- readFastq(fastqR1In)
  if (!is.null(fastqR2In)) tmpR2 <- readFastq(fastqR2In)
  # sample name
  if (!is.null(sampleName)) sample.name <- sampleName else sample.name <- strsplit(basename(fastqR1In), ".fastq")[[1]][1]
  #sample.name <- paste0(strsplit(basename(fastqR1In), "_")[[1]][5], "_", strsplit(basename(fastqR1In), "_")[[1]][6]) ## deep-sea samples names
  if (!is.null(fastqR2In)) 
  {
    # IDs are supposed to be matching all along ! 
    if (!identical(sapply(strsplit(as.vector(id(tmpR1)), splitHeader), `[`, 1), sapply(strsplit(as.vector(id(tmpR2)), splitHeader), `[`, 1))) 
    {
      message("not the same order in the fastq...check header of the fastq to split correctly?")
      break
    }
  }
  # Filter R1 & R2 with width of R1
  keepFWD <- width(tmpR1) >= widthFilter
  tmpR1 <- tmpR1[keepFWD]
  if (!is.null(fastqR2In)) tmpR2 <- tmpR2[keepFWD]
  # Filter R1 & R2 with width of R2
  if (!is.null(fastqR2In)) 
  {
    keepREV <- width(tmpR2) >= widthFilter
    tmpR1 <- tmpR1[keepREV]
    tmpR2 <- tmpR2[keepREV]
  }
  # export the reads that did not match the primers
  writeFastq(tmpR1, fastqR1Out, "a")
  if (!is.null(fastqR2In)) writeFastq(tmpR2, fastqR2Out, "a")
  return(list("sampleID" = sample.name, "readsCount" = length(tmpR1)))
}

############################################################################################################
############################################################################################################
################## Count reads
fastqCountReads <- function(fastqR1In, fastqR2In = NULL, sampleName = NULL, splitHeader = " ") {
  require(ShortRead)
  require(Biostrings)
  # reading both fastq 
  tmpR1 <- readFastq(fastqR1In)
  if (!is.null(fastqR2In)) tmpR2 <- readFastq(fastqR2In)
  # sample name
  if (!is.null(sampleName)) sample.name <- sampleName else sample.name <- strsplit(basename(fastqR1In), ".fastq")[[1]][1]
  # IDs are supposed to be matching all along ! 
  if (!is.null(fastqR2In))
  {
    if (!identical(sapply(strsplit(as.vector(id(tmpR1)), splitHeader), `[`, 1), sapply(strsplit(as.vector(id(tmpR2)), splitHeader), `[`, 1))) 
    {
      message("not the same order in the fastq...check header of the fastq to split correctly?")
      break
    }
  }
  # retunr reads count
  return(list("sampleID" = sample.name, "readsCount" = length(tmpR1)))
}

