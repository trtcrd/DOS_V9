## companion R functions

### data summary for ggplot
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

### encoding a vector into numeric (even if factor...)
varNumEncode <- function(variable) {
  uni <- unique(as.character(variable))
  uni_vec <- as.character(variable)
  num_vec <- c(1:length(variable))
  # now go along
  for (i in 1:length(variable))
  {
    if (!is.na(uni_vec[i])) {
      num_vec[i] <- which(uni_vec[i] == uni)
    } else {
      message(paste0("some empty values encoded into ", length(uni)+1))
      num_vec[i] <- length(uni)+1
    }
  }
  ### now for relevant pch
  for (i in 1:length(num_vec)) 
  {
    if (num_vec[i] > 25) num_vec[i] <- num_vec[i] + 7
  }
  return(num_vec)
}

### Cumulative Sum Scaling
cssNorm <- function(count_table, samples_as_rows) 
{
  require(metagenomeSeq)
  if (samples_as_rows == TRUE) count_table <- t(count_table)
  data.metagenomeSeq = newMRexperiment(count_table, # sample as rows
                                       featureData=NULL, libSize=NULL, normFactors=NULL) #using filtered data 
  p = cumNormStatFast(data.metagenomeSeq) # default is 0.5
  data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
  data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
  if (samples_as_rows) data.CSS <- t(data.CSS)
  return(data.CSS)
}



### aggreg taxa data from otu table
aggregByTaxa <- function(otu_table, taxo, taxo_rank, rare_threshold = NULL, remove_unass = F,
                         useDataTable = F, functional = F, sort_abund = T) {
  
  require(data.table)
  # if the column "taxon" is a factor, convert it into character
  if (is.factor(taxo[,"taxon"])) taxo[,"taxon"] <- as.character(taxo[,"taxon"])
  
  ## test congruence
  if (!identical(colnames(otu_table), rownames(taxo))) {
    message("not the same order between columns in table and rows in taxo... please fix it.")
    break
  }
  ### if remove unassiged
  if (remove_unass) 
  {
    unass <- grep("unassigned", taxo[,"taxon"], ignore.case = T)
    otu_table <- otu_table[,-unass]
    taxo <-taxo[-unass,]
  }
  
  ###
  # extract the "taxo_rank" assignment
  if (taxo_rank == "last")
  {
    taxa <- unlist(sapply(strsplit(as.character(taxo$taxon), ";"), function(x)x[length(x)]))
    taxa[is.na(taxa)] <- "zz_Unassigned"
    taxa[taxa == "unassigned"] <- "zz_Unassigned"
  } else {
    ## split at the tax rank
    taxa <- sapply(strsplit(as.character(taxo$taxon), ";"), `[`,taxo_rank)
    taxa[is.na(taxa)] <- "zz_Unassigned"
    taxa[taxa == "unassigned"] <- "zz_Unassigned"
  }
  # creating a table for sum by levels of tax rank
  dat <- as.data.frame(array(NA, c(nrow(otu_table), length(levels(as.factor(taxa))))))
  rownames(dat) <- rownames(otu_table)
  ll <- levels(as.factor(taxa))
  ## if functional we keep long complicated names of levels... 
  if (!functional) for (i in 1:length(ll)) ll[i] <- strsplit(ll[i], split="(", fixed=TRUE)[[1]][1]    # need to remove parenthesis if any
  colnames(dat) <- as.character(ll)
  ## binary
  dat_bin <- dat
  otu_table_bin <- otu_table
  otu_table_bin[otu_table_bin != 0] <- 1
  
  # summing the abundance by tax_rank
  for (i in ll)
  {
    ### !!! don't get it why sometimes, fixed =T, useBytes = T is required... sometimes not.
    tt     <- otu_table[,grep(i, taxa,  fixed =T, useBytes = T)]
    tt_bin <- otu_table_bin[,grep(i, taxa,  fixed =T, useBytes = T)]
    if (!is.null(ncol(tt))) { 
      dat[,i] <- rowSums(tt) 
      dat_bin[,i] <- rowSums(tt_bin) 
    } else { 
      dat[,i] <- tt
      dat_bin[,i] <- tt_bin
    }
  }
  
  # normalizing the otu table for % 
  dat_raw <- dat
  otu_norm <- colSums(dat)                  ### for the rare_threshold cutoff
  otu_norm <- otu_norm / sum(otu_norm) *100 ### for the rare_threshold cutoff
  dat <- dat / rowSums(dat) *100
  dat[dat == "NaN"] <- 0
  
  ## remove rare families groupings them in "Others"
  if (!is.null(rare_threshold) & ncol(dat) > 1) 
  {
    dat_others <-  dat[,otu_norm <= rare_threshold]
    dat_ <- dat[,otu_norm > rare_threshold]
    dat_ <- cbind(dat_, Others = rowSums(dat_others))
    # for the raw counts
    dat_others_raw <- dat_raw[,otu_norm <= rare_threshold]
    dat_raw_ <- dat_raw[,otu_norm > rare_threshold]
    dat_raw_ <- cbind(dat_raw_, Others = rowSums(dat_others_raw))
    # for the bin dat
    dat_bin_others <- dat_bin[,otu_norm <= rare_threshold]
    
    dat_bin_ <- dat_bin[,otu_norm > rare_threshold]
    dat_bin_ <- cbind(dat_bin_, Others = rowSums(dat_bin_others))
    
    # swapping unassigned and others
    if (remove_unass)  {
      if (sort_abund) {
        dat <- dat_[,colnames(dat_)[order(colSums(dat_[,1:ncol(dat_)]), decreasing = T)]]
        dat_raw <- dat_raw_[,colnames(dat_raw_)[order(colSums(dat_raw_[,1:ncol(dat_raw_)]), decreasing = T)]]
        # same sorting for binary
        dat_bin <- dat_bin_[,colnames(dat_)[order(colSums(dat_[,1:ncol(dat_)]), decreasing = T)]] 
      } else {
        dat <- dat_
        dat_raw <- dat_raw_
        dat_bin <- dat_bin_
      }
    }
    if (!remove_unass) {
      if (sort_abund) {
        dat <- dat_[,c(colnames(dat_)[order(colSums(dat_[,1:(ncol(dat_)-2)]), decreasing = T)], "Others", "zz_Unassigned")]
        dat_raw <- dat_raw_[,c(colnames(dat_)[order(colSums(dat_[,1:(ncol(dat_)-2)]), decreasing = T)], "Others", "zz_Unassigned")]
        # same sorting for binary
        dat_bin <- dat_bin_[,c(colnames(dat_)[order(colSums(dat_[,1:(ncol(dat_)-2)]), decreasing = T)], "Others", "zz_Unassigned")] 
      } else {
        dat <- dat_
        dat_raw <- dat_raw_
        dat_bin <- dat_bin_
      }
      colnames(dat) <- gsub("Others", paste0("Others < ", rare_threshold, "%"), colnames(dat))
      colnames(dat_raw) <- gsub("Others", paste0("Others < ", rare_threshold, "%"), colnames(dat_raw))
      colnames(dat_bin) <- gsub("Others", paste0("Others < ", rare_threshold, "%"), colnames(dat_bin))
    }
  } else {
    ## reorder by abundance (keeping the zz_unassigned as last)
    if (remove_unass)  {
      if (sort_abund) {
        dat <- dat[,colnames(dat)[order(colSums(dat[,1:ncol(dat)]), decreasing = T)]]
        dat_raw <- dat_raw[,colnames(dat)[order(colSums(dat[,1:ncol(dat)]), decreasing = T)]]
        dat_bin <- dat_bin[,colnames(dat)[order(colSums(dat[,1:ncol(dat)]), decreasing = T)]]
      } 
    }
    if (!remove_unass) {
      if (sort_abund) {
        #if (colnames(dat) %in% "zz_Unassigned") {
        dat <- dat[,c(colnames(dat)[order(colSums(dat[,1:(ncol(dat)-1)]), decreasing = T)], "zz_Unassigned")]
        dat_raw <- dat_raw[,c(colnames(dat)[order(colSums(dat[,1:(ncol(dat)-1)]), decreasing = T)], "zz_Unassigned")]
        dat_bin <- dat_bin[,c(colnames(dat)[order(colSums(dat[,1:(ncol(dat)-1)]), decreasing = T)], "zz_Unassigned")]
        #} else {message("no unassigned...")}
      }
    }
  }
  ## finally rename unassigned
  if (!remove_unass) {
    colnames(dat) <- gsub("zz_Unassigned", "Unassigned", colnames(dat))
    colnames(dat_raw) <- gsub("zz_Unassigned", "Unassigned", colnames(dat_raw))
    colnames(dat_bin) <- gsub("zz_Unassigned", "Unassigned", colnames(dat_bin))
  }
  
  # preparing the output
  output <- list("norm" = dat,
                 "raw" = dat_raw,
                 "bin" = dat_bin,
                 "tax_rank" = taxo_rank,
                 "rare_threshold" = rare_threshold)
  
  return(output)
}



ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}


### distance-decay parameters if log-linear regression form (see Soininen et al., 2007, Ecography)
decayparam <- function(sim, geo, init_dist = 1) ## default init_dist = 1 (in km if geo is in km)
{
  reg <- lm(sim~geo)
  # estimate
  init_sim <- (exp(init_dist*reg$coefficients[[2]]+reg$coefficients[[1]]))*100
  halv_dist <- init_dist-log(2)/reg$coefficients[[2]]
  slope <- reg$coefficients[[2]]
  # lower confidence interval
  init_sim_low <- (exp(init_dist*confint(reg)[2,1]+confint(reg)[1,1]))*100
  halv_dist_low <- init_dist-log(2)/confint(reg)[2,1]
  slope_low <- confint(reg)[2,1]
  # upper confidence interval
  init_sim_up <- (exp(init_dist*confint(reg)[2,2]+confint(reg)[1,2]))*100
  halv_dist_up <- ifelse (confint(reg)[2,2]<0, init_dist-log(2)/confint(reg)[2,2], "+inf")
  slope_up <- confint(reg)[2,2]
  list_return <- list (
    "InitSim(%)"=paste0(round(init_sim,2)," (", round(init_sim_low,2),",", round(init_sim_up,2),")"), 
    "Slope"=paste0(round(slope,8), " (", round(slope_low,8),",", round(slope_up,8),")"),
    "HalvDist(km)"=paste0(round(halv_dist,0), " (",round(halv_dist_low,0),",",ifelse(halv_dist_up=="+inf","+inf",round(halv_dist_up,0)),")"))
  return (list_return)
}
