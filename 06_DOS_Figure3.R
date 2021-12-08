##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Figure 3 S8 S9 S10 - Table S4 S5
##################################################################################


## setting path **to be adjusted*
setwd("~/path/to/be/adjusted/")

path.fig <- paste0(getwd(), "/Figure3")
if(!dir.exists(path.fig)) dir.create(path.fig)

## custom functions **to be adjusted to the companion_functions.R path**
source("~/path/to/be/adjusted/companionFunctions.R")

## libraries
library('doMC')
registerDoMC(cores = 16)
library(vegan)
library(ape)
library(ggplot2)
library(ggpubr)
library(cowplot)



###############################################################################################################
###############################################################################################################
###############################################################################################################
###### Import required data

asv <- readRDS("asv_table_euks_clean_agg.rds")
comp <- readRDS("comp_table_euks_clean.rds")
## typo on Arctic in the metadata :/
comp$basin <- gsub("Artic Ocean", "Arctic Ocean", comp$basin)
## Taxo file with functional and size annotations 
taxo <- readRDS("asvs_taxo_functional_annotations.rds")

# check mapping
identical(rownames(asv), rownames(comp)) # TRUE
identical(colnames(asv), rownames(taxo)) # TRUE

### ASVs vectors
asvs_euph <- readRDS("asvs_euph.rds")
asvs_aph  <- readRDS("asvs_aph.rds")
asvs_plankton  <- readRDS("asvs_plankton.rds")
asvs_sed_only <- readRDS("asvs_sed_only.rds")
asvs_plan_sed <- readRDS("asvs_plan_sed.rds")

###############################################################################################################
###############################################################################################################
###############################################################################################################
## Benthic alpha diversity

# Focus on abyssal plains and remove coastal stations: Lie (California gulf) and MDW_ST38 (Gibraltar) 
asv_benthic <- subset(asv[,asvs_sed_only], comp$biome == "Sediment" & comp$area != "Mediterranean Sea" & comp$dataset != "Lie" & comp$station != "MDW_ST38")
comp_benthic <- comp[rownames(asv_benthic),]
comp_benthic$abs_latitude <- abs(comp_benthic$latitude)


## GAMs for benthic alpha diversity
groups_check <- c("All", "Nematoda","Foraminifera", "Platyhelminthes","Ciliophora", "polychaetes",
                  "Mollusca") #,"Echinodermata","Malacostraca")

yvar <- c("normRich", "Shannon")
xvar <- c("latitude", "primprod", "POC_export", "POC_seafloor")
ynames <- c("Norm. richness", "Shannon div.")
xnames <- c("latitude", "primary prod.", "POC surface", "POC seafloor")

gp_all <- list()
cpt_all <- 1

for (g in groups_check) {
  if (g == "All") {
    asv_benthic_group <- asv_benthic
    comp_benthic_group <- comp_benthic
  } else {
    # Only taxo groups
    tmp_tax <- taxo[asvs_sed_only,]
    tmp_tax <- tmp_tax[grep(g, tmp_tax$taxo_group),]
    asv_benthic_group <- subset(asv_benthic[,rownames(tmp_tax)], rowSums(asv_benthic[,rownames(tmp_tax)]) >= 100) ## remove samples with too few reads
    comp_benthic_group <- comp_benthic[rownames(asv_benthic_group),]
  }
  # Diversity
  sha_ben <- vegan::diversity(asv_benthic_group, index = "shannon")
  comp_benthic_group$Shannon <- sha_ben
  comp_benthic_group$normRich <- specnumber(rrarefy(asv_benthic_group, min(rowSums(asv_benthic_group))))
  comp_benthic_tmp <- comp_benthic_group[,c(yvar, xvar)]
  
  gp_list <- list()
  cpt <- 1
  
  for (j in xvar) {
    for (i in yvar) {
      gam_mod <- mgcv::gam(comp_benthic_tmp[,i] ~ s(comp_benthic_tmp[,j], k=3))
      summ_gam <- summary(gam_mod)
      ex_dev <- round(summ_gam$dev.expl * 100, 1)
      p <- summ_gam$s.table[,'p-value']
      if(p < 0.05) {
        if(p >= 0.01) sig <- "*"
        if(p >= 0.001) sig <- "**"
        if(p < 0.001) sig <- "***"
      } else {sig <- "ns" }
      p1 <- ggplot(comp_benthic_tmp, aes_string(x=j, y = i)) + 
        geom_point(color = c("red", "green")[which(i == yvar)]) + 
        stat_smooth(method = "gam", formula = y ~ s(x, k=3), size = 1, alpha = 0.2) +
        theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) + 
        labs(x = xnames[which(j == xvar)], y = ynames[which(i == yvar)], title = paste0("Exp. dev. = ", ex_dev,"%", sig)) +
        theme(text = element_text(size=8))
      ## add lengend on top 
      if (g == "All") p1 <- annotate_figure(p1, top = text_grob(ynames[which(i == yvar)], face = "bold", size = 10))
      gp_list[[cpt]] <- p1
      cpt = cpt +1
    }
  }
  
  pLat <- ggarrange(gp_list[[1]], gp_list[[2]], ncol = 2, nrow = 1)
  pProd <- ggarrange(gp_list[[3]], gp_list[[4]], ncol = 2, nrow = 1)
  pPOCs <- ggarrange(gp_list[[5]], gp_list[[6]], ncol = 2, nrow = 1)
  pPOCb <- ggarrange(gp_list[[7]], gp_list[[8]], ncol = 2, nrow = 1)
  
  if (g == "All") {
    g <- "All benthic ASVs"
    pLat  <- annotate_figure(pLat, top = text_grob(xnames[1], face = "bold", size = 12))
    pProd <- annotate_figure(pProd, top = text_grob(xnames[2], face = "bold", size = 12))
    pPOCs <- annotate_figure(pPOCs, top = text_grob(xnames[3], face = "bold", size = 12))
    pPOCb <- annotate_figure(pPOCb, top = text_grob(xnames[4], face = "bold", size = 12))
  } 
  
  gp_all[[cpt_all]] <- ggarrange(pLat, pProd, pPOCs, pPOCb, ncol = 4, nrow = 1)
  gp_all[[cpt_all]] <- annotate_figure(gp_all[[cpt_all]], left = text_grob(g, face = "bold", size = 10, rot = 90))
  cpt_all <- cpt_all + 1
  
}

p_export <- ggarrange(plotlist =  gp_all, ncol = 1, nrow = 7, heights = c(1,rep(0.75,6)))
ggsave("Figure3/Figure_S8_GAMs_benthic_richness_diversity.pdf", p_export, width = 13, height = 9.5)


###############################################################################################################
###############################################################################################################
###############################################################################################################
## Benthic beta-diversity - Figure 3

### CSS normalization and Bray-Curtis (with all benthic ASVs)
asv_benthic_css  <- cssNorm(asv_benthic, samples_as_rows = TRUE) 
vegdist_benthic  <- vegdist(asv_benthic_css, method="bray")

## PCoA benthic communities 
pcoa_benthic <- pcoa(vegdist_benthic)
PCOA <- pcoa_benthic$vectors[,1:2]
var_axis_pcoa <- round(pcoa_benthic$values$Relative_eig * 100, 2)[1:2]

## envfit
envFit <- envfit(ord=PCOA, env = comp_benthic[,c('sb_nitrate', 'sb_silicate','sb_o2dissolve', 'sb_temp', 'sb_salinity', 'POC_seafloor', 'POC_export', 'primprod')], permutations = 999, na.rm = T)

pdf("Figure3/Figure_3A_PCOA_benthic.pdf", width = 10, height=6, useDingbats = F)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(PCOA, col=varNumEncode(comp_benthic$watling) +1, pch=varNumEncode(comp_benthic$watling) +1, cex = 1,
     main="PCoA bray CSS norm", xlab = paste0("PC1: ", var_axis_pcoa[1], "%"), ylab = paste0("PC2: ", var_axis_pcoa[2], "%"))
ordisurf(PCOA, comp_benthic$abs_latitude,  add = T, col="black", labcex = 0.8, lwd.cl = 0.2, family = "gaussian")
leg <- unique(as.character(comp_benthic$watling))
legend("topright", inset=c(-0.5,0), leg, col=unique(varNumEncode(comp_benthic$watling)+1), 
       pch=unique(varNumEncode(comp_benthic$watling)+1), box.lty=0)
plot(envFit, p.max = 0.05, col = "red")
dev.off()

## adonis - watling provinces
adonis(vegdist_benthic ~ comp_benthic$watling, permutations = 999, parallel = 8)
## adonis - absolute latitude 
adonis(vegdist_benthic ~ comp_benthic$abs_latitude, strata = comp_benthic$basin, permutations = 999, parallel = 8)

## contrained ordination with environmental variables selection
met_std <- subset(comp_benthic, comp_benthic$POC_seafloor != "")
veg_std <- as.dist(as.matrix(vegdist_benthic)[rownames(met_std),rownames(met_std)])
met_std <- decostand(met_std[,c('sb_nitrate', 'sb_silicate','sb_o2dissolve', 'sb_temp', 'sb_salinity', 'POC_seafloor', 'POC_export', 'primprod')], method = "standardize", na.rm = T)
mod <- dbrda(veg_std ~ ., data = met_std)
anova(mod) ## yes
rdaBase <- dbrda(veg_std ~ 1, data=met_std)
rdaTout <- dbrda(veg_std ~ ., data=met_std)
rdaSel <- ordiR2step(rdaBase,scope=formula(rdaTout))

## get the proportion of variance explained
rdsSel_tab <- as.data.frame(rdaSel$anova)
rdsSel_tab$Inertia    <- c(rdaSel$CCA$eig, sum(rdaSel$CCA$eig))
rdsSel_tab$Proportion <- rdsSel_tab$Inertia *100 / rdaSel$tot.chi
write.table(rdsSel_tab, "Figure3/Table_S4_dbRDA_benthic_explained_variance.tsv", quote = F, sep = "\t", dec = ",")


##### Distance decay 

## pairwise shared ASVs
dd_matrix <- vegdist(asv_benthic, method="jaccard", binary = T)
dd <- (1 - as.matrix(dd_matrix)) *100 ## % of shared ASVs
identical(rownames(as.matrix(dd)), rownames(comp_benthic)) # TRUE
## geographic distance matrix
geo_dist <- comp_benthic[rownames(as.matrix(dd)),c("sampleID", "latitude", "longitude")]
colnames(geo_dist) <- c("index", "lat", "lon")
geo_dist <- GeoDistanceInMetresMatrix(geo_dist) ## <-- see companion function
geo_dist <- as.matrix(geo_dist)/1000 ## in km
diag(geo_dist) <- NA
dimnames(geo_dist) <- list(rownames(comp_benthic), rownames(comp_benthic))
## Melt 
melt_geo <- reshape::melt(geo_dist)
colnames(melt_geo)[3] <- "geo_dist"
# Vegdist
melt_bray <- reshape::melt(dd)
colnames(melt_bray)[3] <- "bray"

melted <- cbind(melt_bray, distance = melt_geo$geo_dist)
# for each basin 
melted_basin <- c()
for (i in unique(comp_benthic$basin)) {
  tmp_cmp_eup <- subset(comp_benthic, comp_benthic$basin == i)
  tmp_geo_eup <- geo_dist[rownames(tmp_cmp_eup), rownames(tmp_cmp_eup)]
  tmp_veg_eup <- dd[rownames(tmp_cmp_eup), rownames(tmp_cmp_eup)]
  tmp_geo_melted <- reshape::melt(tmp_geo_eup)
  colnames(tmp_geo_melted)[3] <- "geo_dist"
  tmp_veg_melted <- reshape::melt(tmp_veg_eup)
  colnames(tmp_veg_melted)[3] <- "bray"
  melted_basin <- rbind(melted_basin, cbind(tmp_veg_melted, distance = tmp_geo_melted$geo_dist))
}

melted <- melted_basin
melted <- subset(melted, melted$distance <= 5000)

## bins 
dist_bins <- c()
for (i in 1:nrow(melted)) {
  if (!is.na(melted$distance[i])) {
    if (melted$distance[i] >= 0 & melted$distance[i] <= 1) dist_bins <- c(dist_bins, "1")
    if (melted$distance[i] > 1 & melted$distance[i] <= 5) dist_bins <- c(dist_bins, "5")
    if (melted$distance[i] > 5 & melted$distance[i] <= 10) dist_bins <- c(dist_bins, "10")
    if (melted$distance[i] >= 10 & melted$distance[i] <= 100) dist_bins <- c(dist_bins, "100")
    if (melted$distance[i] > 100 & melted$distance[i] <= 1000) dist_bins <- c(dist_bins, "1000")
    if (melted$distance[i] > 1000 & melted$distance[i] <= 2000) dist_bins <- c(dist_bins, "2000")
    if (melted$distance[i] > 2000 & melted$distance[i] <= 3000) dist_bins <- c(dist_bins, "3000")
    if (melted$distance[i] > 3000 & melted$distance[i] <= 4000) dist_bins <- c(dist_bins, "4000")
    if (melted$distance[i] > 4000 & melted$distance[i] <= 5000) dist_bins <- c(dist_bins, "5000")
  } else { dist_bins <- c(dist_bins, "NA") }
  if (i %% 5000 == 0) message(i)
}

melted_bins <- cbind(melted, dist_bins) 

p <- ggplot(melted_bins, aes(x=dist_bins, y=bray)) + 
  geom_violin() + xlim("1", "5","10", "100", "1000", "2000","3000","4000", "5000") +
  stat_summary(fun=mean, aes(group=1), geom="line", colour="darkgoldenrod1", size=1) +
  theme(legend.position="none") + theme_classic() + ylab("% of shared benthic ASVs") + xlab("Distance (km)") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="darkgoldenrod1")

pdf("Figure3/Figure_3B_dist_decay_shared_benthic_asvs.pdf", width = 3.2, height = 2, useDingbats = F)
print(p)
dev.off()

### Betadisper as function of sampling scales

#### beta disper by sampling spatial scale
bd <- vegdist_benthic
cp <- comp_benthic

Group_basin <- paste0(cp$basin)
Group_area <- paste0(cp$basin, "/",cp$watling)
Group_station <- paste0(cp$basin, "/",cp$watling, "/", cp$station)
Group_deploy <- paste0(cp$basin, "/",cp$watling, "/", cp$station, "/", cp$deployment)
Group_core <- paste0(cp$basin, "/",cp$watling, "/", cp$station, "/", cp$deployment, "/", cp$core)

out <- array(NA, c(nrow(cp), 5))
colnames(out) <- c("Basin", "Province", "Station", "Deployment", "Core")

bas <- betadisper(bd, group = Group_basin, type = "median")
out[,"Basin"] <- bas$distances
are <- betadisper(bd, group = Group_area, type = "median")
out[,"Province"] <- are$distances
sta <- betadisper(bd, group = Group_station, type = "median")
out[,"Station"] <- sta$distances
dep <- betadisper(bd, group = Group_deploy, type = "median")
out[,"Deployment"] <- dep$distances
cor <- betadisper(bd, group = Group_core, type = "median")
out[,"Core"] <- cor$distances

out[out == 0] <- NA


### ggplot
df <- reshape::melt(out)
p <- ggplot(df, aes(x=X2, y=value)) + 
  geom_violin() + xlim("Core", "Deployment", "Station","Province", "Basin") +
  stat_summary(fun=mean, aes(group=1), geom="line", colour="darkgoldenrod1", size=1) +
  theme(legend.position="none") + theme_classic() + 
  xlab("Sampling scale") + ylab("Variation in beta-diversity\n(distance to group centroid)") +
  stat_summary(fun.data=data_summary, geom="pointrange", color="darkgoldenrod1")

pdf("Figure3/Figure_3C_Betadisper_benthic_sampling_scales.pdf", width = 3.3, height = 2.5, useDingbats = F)
print(p)
dev.off()



###############################################################################################################
###############################################################################################################
###############################################################################################################


### Check distance decay parameters for selected benthic groups
groups_check <- c("All", "Nematoda","Foraminifera", "Platyhelminthes","Ciliophora", "polychaetes",
                  "Mollusca","Amoebozoa")

## prepare the output
dd_output <- array(NA, c(length(groups_check), 4))
rownames(dd_output) <- groups_check
colnames(dd_output) <- c("Correlation with distance", "Initial similarity (%)", "Slope", "Halving distance (km)")

for (g in groups_check) {
  
  if (g == "All") {
    tmp_dd <- asv_benthic
    tmp_comp_dd <- comp_benthic
  } else {
    tmp_g <- taxo[colnames(asv_benthic),]
    tmp_g <- tmp_g[grep(g, tmp_g$taxo_group),]
    tmp_dd <- asv_benthic[,rownames(tmp_g)]
    tmp_dd <- subset(tmp_dd, rowSums(tmp_dd) > 200)
    tmp_comp_dd <- comp_benthic[rownames(tmp_dd),]
  }
  ## average Sørenson distance matrix 
  #veg_dist <- avgdist(otu_benthic, sample = min(rowSums(otu_benthic)), iterations = 10, dmethod="bray", binary = T) <-- too slow
  vegList <- foreach(i = 1:10) %dopar% assign(paste0("vegList", i), vegdist(rrarefy(tmp_dd, min(rowSums(tmp_dd))), method= "bray", binary = T)) # <-- **ram-hungry** 
  gc()
  ## now average the matrices
  vegList_avg <- Reduce("+",vegList)/length(vegList)
  ## melt
  vegList_avg_vec <- as.vector(vegList_avg)
  ## geo dist
  geo_dist <- tmp_comp_dd[rownames(as.matrix(vegList_avg)),c("sampleID", "latitude", "longitude")]
  colnames(geo_dist) <- c("index", "lat", "lon")
  geo_dist <- GeoDistanceInMetresMatrix(geo_dist)
  ## melt
  geo_dist_vec <- as.vector(as.dist(as.matrix(geo_dist)/1000)) ## in km
  
  ## mantel 
  mant <- mantel(vegList_avg, geo_dist)
  if (mant$signif < 0.05) {
    if (mant$signif >= 0.01) sig <- "*"
    if (mant$signif >  0.001) sig <- "**"
    if (mant$signif <= 0.001) sig <- "***"
  } else {sig <- "ns" }
  ## paste the result
  dd_output[g,"Correlation with distance"] <- paste0(round(mant$statistic, 3),sig)
  
  ## dd 
  melted <- data.frame(distance = geo_dist_vec, bray = vegList_avg_vec)
  dd_params <- decayparam(log(1.01-melted$bray), melted$distance) ## <-- see companion function
  ## paste the result
  dd_output[g,"Initial similarity (%)"] <- dd_params$`InitSim(%)`
  dd_output[g,"Slope"] <- dd_params$Slope
  dd_output[g,"Halving distance (km)"] <- dd_params$`HalvDist(km)`
  message(g)
}

write.table(dd_output, "Figure3/Table_S5_DD_parameters_benthic.tsv", quote = F, row.names = T, sep = "\t")


#### compare with plankton
dd_output_plank <- array(NA, c(2, 4))
rownames(dd_output_plank) <- c("Euphotic", "Aphotic")
colnames(dd_output_plank) <- c("Correlation with distance", "Initial similarity (%)", "Slope", "Halving distance (km)")

## Euphotic pico nano size fractions
comp_euph_pico_nano <- subset(comp, comp$biome == "Euphotic" & comp$latitude != "" & (comp$size_class == "nano" | comp$size_class == "pico")) 
asv_euph_pico_nano <- asv[rownames(comp_euph_pico_nano),asvs_euph]
## Aphotic pico nano
comp_aph_pico_nano <- subset(comp, comp$biome == "Aphotic"& comp$latitude != ""  & (comp$size_class == "nano" | comp$size_class == "pico")) 
asv_aph_pico_nano <- asv[rownames(comp_aph_pico_nano),asvs_aph]


for (g in rownames(dd_output_plank)) {
  if (g == "Euphotic") {
    tmp_dd      <- asv_euph_pico_nano
    tmp_comp_dd <- comp_euph_pico_nano
  }
  if (g == "Aphotic") {
    tmp_dd      <- asv_aph_pico_nano
    tmp_comp_dd <- comp_aph_pico_nano
  }
  vegList <- foreach(i = 1:10) %dopar% assign(paste0("vegList", i), vegdist(rrarefy(tmp_dd, min(rowSums(tmp_dd))), method= "bray", binary = T))
  gc()
  ## now average the matrices
  vegList_avg <- Reduce("+",vegList)/length(vegList)
  ## melt
  vegList_avg_vec <- as.vector(vegList_avg)
  ## geo dist
  geo_dist <- tmp_comp_dd[rownames(as.matrix(vegList_avg)),c("sampleID", "latitude", "longitude")]
  colnames(geo_dist) <- c("index", "lat", "lon")
  geo_dist <- GeoDistanceInMetresMatrix(geo_dist)
  ## melt
  geo_dist_vec <- as.vector(as.dist(as.matrix(geo_dist)/1000)) ## in km
  
  ## mantel 
  mant <- vegan::mantel(vegList_avg, geo_dist)
  if (mant$signif < 0.05) {
    if (mant$signif >= 0.01) sig <- "*"
    if (mant$signif >  0.001) sig <- "**"
    if (mant$signif <= 0.001) sig <- "***"
  } else {sig <- "ns" }
  ## paste the result
  dd_output_plank[g,"Correlation with distance"] <- paste0(round(mant$statistic, 3),sig)
  
  ## dd 
  melted <- data.frame(distance = geo_dist_vec, bray = vegList_avg_vec)
  dd_params <- decayparam(log(1.01-melted$bray), melted$distance)
  ## paste the result
  dd_output_plank[g,"Initial similarity (%)"] <- dd_params$`InitSim(%)`
  dd_output_plank[g,"Slope"] <- dd_params$Slope
  dd_output_plank[g,"Halving distance (km)"] <- dd_params$`HalvDist(km)`
}

write.table(dd_output_plank, "Figure3/Table_S5_DD_parameters_pelagic.tsv", quote = F, row.names = T, sep = "\t")



###############################################################################################################
###############################################################################################################
###############################################################################################################


#### compare Sørensen distance decay between Euphotic Aphotic (pico nano fractions) and benthic

# benthic 
vegList <- foreach(i = 1:10) %dopar% assign(paste0("vegList", i), vegdist(rrarefy(asv_benthic, min(rowSums(asv_benthic))), method= "bray", binary = T))
gc()
## now average the matrices
vegList_avg <- Reduce("+",vegList)/length(vegList)
vegList_avg_vec_ben <- as.vector(vegList_avg)

## geo dist
geo_dist <- comp_benthic[rownames(as.matrix(vegList_avg)),c("sampleID", "latitude", "longitude")]
colnames(geo_dist) <- c("index", "lat", "lon")
geo_dist <- GeoDistanceInMetresMatrix(geo_dist)
## melt
geo_dist_vec_ben <- as.vector(as.dist(as.matrix(geo_dist)/1000)) ## in km

# euphotic
vegList_eu <- foreach(i = 1:10) %dopar% assign(paste0("vegList", i), vegdist(rrarefy(asv_euph_pico_nano, min(rowSums(asv_euph_pico_nano))), method= "bray", binary = T))
gc()
## now average the matrices
vegList_avg <- Reduce("+",vegList_eu)/length(vegList_eu)
vegList_avg_vec_euph <- as.vector(vegList_avg)

## geo dist
geo_dist <- comp_euph_pico_nano[rownames(as.matrix(vegList_avg)),c("sampleID", "latitude", "longitude")]
colnames(geo_dist) <- c("index", "lat", "lon")
geo_dist <- GeoDistanceInMetresMatrix(geo_dist)
## melt
geo_dist_vec_euph <- as.vector(as.dist(as.matrix(geo_dist)/1000)) ## in km

# aphotic
vegList_ap <- foreach(i = 1:10) %dopar% assign(paste0("vegList", i), vegdist(rrarefy(asv_aph_pico_nano, min(rowSums(asv_aph_pico_nano))), method= "bray", binary = T))
gc()
## now average the matrices
vegList_avg <- Reduce("+",vegList_ap)/length(vegList_ap)
vegList_avg_vec_aph <- as.vector(vegList_avg)

## geo dist
geo_dist <- comp_aph_pico_nano[rownames(as.matrix(vegList_avg)),c("sampleID", "latitude", "longitude")]
colnames(geo_dist) <- c("index", "lat", "lon")
geo_dist <- GeoDistanceInMetresMatrix(geo_dist)
## melt
geo_dist_vec_aph <- as.vector(as.dist(as.matrix(geo_dist)/1000)) ## in km

## concatenate
dd_all <- data.frame(rbind(cbind(Biome = "Benthic", Sor = vegList_avg_vec_ben, Dist = geo_dist_vec_ben),
                           cbind(Biome = "Euphotic", Sor = vegList_avg_vec_euph, Dist = geo_dist_vec_euph),
                           cbind(Biome = "Aphotic", Sor = vegList_avg_vec_aph, Dist = geo_dist_vec_aph)), stringsAsFactors = F)
dd_all$Sor <- as.numeric(dd_all$Sor)
dd_all$Dist <- as.numeric(dd_all$Dist)
dd_all$Biome <- factor(dd_all$Biome, levels=c("Euphotic", "Aphotic", "Benthic") )


p <- ggplot(dd_all, aes(x=Dist+1, y=1.01-Sor, color = Biome)) + 
  geom_point(size = c(0.001,0.001,0.05)[as.numeric(as.factor(dd_all$Biome))],
             color = c(alpha("deepskyblue3", 0.05), alpha("blue4", 0.6), alpha("orange3",0.1))[as.numeric(as.factor(dd_all$Biome))]) + 
  stat_smooth(method="lm", se = F) + 
  scale_color_manual(values=c(Benthic= "orange3", Euphotic= "deepskyblue3", Aphotic= "blue4")) + 
  scale_y_log10() + theme_classic() + 
  xlab("Distance separation (km)") + ylab("log community similarity (Sørensen)")

ggsave("Figure3/Figure_S9_DistanceDecay.png", p, width = 7, height= 6)



###############################################################################################################
###############################################################################################################
###############################################################################################################
## Neutral models (Sloan et al., 2006) 

## Neutral distribution model 
## See function here ## https://github.com/Russel88/MicEco/blob/master/R/neutral.fit.R
source("https://raw.githubusercontent.com/Russel88/MicEco/master/R/neutral.fit.R")
library(bbmle) # required package
library(Hmisc) # required package

### agg at the station level **long**
euph_sta <- aggregate(asv[comp$biome == "Euphotic",asvs_euph], by = list(comp$station[comp$biome == "Euphotic"]), FUN = sum) 
aph_sta  <- aggregate(asv[comp$biome == "Aphotic",asvs_aph], by = list(comp$station[comp$biome == "Aphotic"]), FUN = sum) 
ben_sta  <- aggregate(asv[comp$biome == "Sediment",asvs_sed_only], by = list(comp$station[comp$biome == "Sediment"]), FUN = sum) 

rownames(euph_sta) <- euph_sta$Group.1
euph_sta_ <- euph_sta[,2:ncol(euph_sta)]
rownames(aph_sta) <- aph_sta$Group.1
aph_sta_ <- aph_sta[,2:ncol(aph_sta)]
rownames(ben_sta) <- ben_sta$Group.1
ben_sta_ <- ben_sta[,2:ncol(ben_sta)]

saveRDS(euph_sta_, "asv_table_euph_station.rds")
saveRDS(aph_sta_, "asv_table_aph_station.rds")
saveRDS(ben_sta_, "asv_table_ben_station.rds")

## euphotic
euph_sloan <- neutral.fit(euph_sta_[,colSums(euph_sta_) > 0])
euph_sloan_df <- euph_sloan[[2]]

neut <- c()
for (i in 1:nrow(euph_sloan_df)) {
  if (euph_sloan_df[i,"freq"] < euph_sloan_df[i,"Lower"]) neut <- c(neut, "Below") 
  if (euph_sloan_df[i,"freq"] >= euph_sloan_df[i,"Lower"] & euph_sloan_df[i,"freq"] <= euph_sloan_df[i,"Upper"]) neut <- c(neut, "Neutral") 
  if (euph_sloan_df[i,"freq"] > euph_sloan_df[i,"Upper"]) neut <- c(neut, "Above") 
}

euph_sloan_df$Neutral <- neut
euph_sloan_df$taxo_groups <- taxo[rownames(euph_sloan_df),"taxo_group"]

distrib_euph <- ggplot(data=euph_sloan_df, aes(x=log(p), y=freq)) + 
  geom_point(size = 0.2, aes(colour = Neutral)) +
  geom_line(aes(y=freq.pred)) +
  geom_line(aes(y=Lower), linetype="dotted") +
  geom_line(aes(y=Upper), linetype="dotted") + 
  xlab("log (mean relative abundance)") + ylab("Occurrence frequency") + labs(title=paste0("Euphotic - R²=", round(euph_sloan[[1]][,"gRsqr"], 2)))

assign_tax_groups_check <- c("Copepoda", "Dinophyceae", "Collodaria", "Acantharea", "Diatomeae", "Eupelagonemidae", "Fungi", "Hydrozoa", "Haptophyta")
plot_list = list()

## ggplot default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggcol <- gg_color_hue(3)

for (i in assign_tax_groups_check) {
  tmp_taxo <- euph_sloan_df[grep(i, euph_sloan_df$taxo_groups),]
  ## cumulative relative abundances 
  pp <- ggplot(tmp_taxo, aes(x = Neutral, fill = Neutral)) +
    geom_bar() + xlim("Above", "Neutral", "Below") +
    labs(title=paste0("Euphotic ", i), x="", y = "#ASVs") +
    scale_fill_manual(values=c(Above= ggcol[1], Neutral= ggcol[3], Below= ggcol[2])) +
    theme(legend.position="none", axis.text = element_text(size = rel(0.7)), title = element_text(size = rel(0.7)))
  plot_list[[i]] <- pp
}


p_euph <- ggarrange(plotlist = plot_list, nrow = 3, ncol =3)



## Aphotic
aph_sloan <- neutral.fit(aph_sta_[,colSums(aph_sta_) > 0])
aph_sloan_df <- aph_sloan[[2]]

neut <- c()
for (i in 1:nrow(aph_sloan_df)) {
  if (aph_sloan_df[i,"freq"] < aph_sloan_df[i,"Lower"]) neut <- c(neut, "Below") 
  if (aph_sloan_df[i,"freq"] >= aph_sloan_df[i,"Lower"] & aph_sloan_df[i,"freq"] <= aph_sloan_df[i,"Upper"]) neut <- c(neut, "Neutral") 
  if (aph_sloan_df[i,"freq"] > aph_sloan_df[i,"Upper"]) neut <- c(neut, "Above") 
}

aph_sloan_df$Neutral <- neut
aph_sloan_df$taxo_groups <- taxo[rownames(aph_sloan_df),"taxo_group"]

distrib_aph <- ggplot(data=aph_sloan_df, aes(x=log(p), y=freq)) + 
  geom_point(size = 0.2, aes(colour = Neutral)) +
  #geom_smooth(method = "loess", se = T)
  geom_line(aes(y=freq.pred)) +
  geom_line(aes(y=Lower), linetype="dotted") +
  geom_line(aes(y=Upper), linetype="dotted") + 
  xlab("log (mean relative abundance)") + ylab("Occurrence frequency") + labs(title=paste0("Aphotic - R²=", round(aph_sloan[[1]][,"gRsqr"], 2)))

assign_tax_groups_check <- c("Copepoda", "Dinophyceae", "Collodaria", "Acantharea", "Diatomeae", "Eupelagonemidae", "Fungi", "Hydrozoa", "Haptophyta")
plot_list = list()

for (i in assign_tax_groups_check) {
  tmp_taxo <- aph_sloan_df[grep(i, aph_sloan_df$taxo_groups),]
  ## cumulative relative abundances 
  pp <- ggplot(tmp_taxo, aes(x = Neutral, fill = Neutral)) +
    geom_bar() + xlim("Above", "Neutral", "Below") +
    labs(title=paste0("Aphotic ", i), x="", y = "#ASVs") +
    scale_fill_manual(values=c(Above= ggcol[1], Neutral= ggcol[3], Below= ggcol[2])) +
    theme(legend.position="none", axis.text = element_text(size = rel(0.7)), title = element_text(size = rel(0.7)))
  plot_list[[i]] <- pp
}

p_aph <- ggarrange(plotlist = plot_list, nrow = 3, ncol =3)


## Benthic
ben_sloan <- neutral.fit(ben_sta_[,colSums(ben_sta_) > 0])
ben_sloan_df <- ben_sloan[[2]]

neut <- c()
for (i in 1:nrow(ben_sloan_df)) {
  if (ben_sloan_df[i,"freq"] < ben_sloan_df[i,"Lower"]) neut <- c(neut, "Below") 
  if (ben_sloan_df[i,"freq"] >= ben_sloan_df[i,"Lower"] & ben_sloan_df[i,"freq"] <= ben_sloan_df[i,"Upper"]) neut <- c(neut, "Neutral") 
  if (ben_sloan_df[i,"freq"] > ben_sloan_df[i,"Upper"]) neut <- c(neut, "Above") 
}

ben_sloan_df$Neutral <- neut
ben_sloan_df$taxo_groups <- taxo[rownames(ben_sloan_df),"taxo_group"]

distrib_ben <- ggplot(data=ben_sloan_df, aes(x=log(p), y=freq)) + 
  geom_point(size = 0.2, aes(colour = Neutral)) +
  #geom_smooth(method = "loess", se = T)
  geom_line(aes(y=freq.pred)) +
  geom_line(aes(y=Lower), linetype="dotted") +
  geom_line(aes(y=Upper), linetype="dotted") + 
  xlab("log (mean relative abundance)") + ylab("Occurrence frequency") + labs(title=paste0("Benthic - R²=", round(ben_sloan[[1]][,"gRsqr"], 2)))

assign_tax_groups_check <- c("Copepoda", "Dinophyceae", "polychaetes", "Nemertea", "Nematoda", "Platyhelminthes", "Bivalvia", "Foraminifera","Gromia", "Ascetosporea", "Ciliophora", "Porifera")
plot_list = list()

for (i in assign_tax_groups_check) {
  tmp_taxo <- ben_sloan_df[grep(i, ben_sloan_df$taxo_groups),]
  ## cumulative relative abundances 
  pp <- ggplot(tmp_taxo, aes(x = Neutral, fill = Neutral)) +
    geom_bar() + xlim("Above", "Neutral", "Below") +
    labs(title=paste0("Benthic ", i), x="", y = "#ASVs") +
    scale_fill_manual(values=c(Above= ggcol[1], Neutral= ggcol[3], Below= ggcol[2])) +
    theme(legend.position="none", axis.text = element_text(size = rel(0.7)), title = element_text(size = rel(0.7)))
  plot_list[[i]] <- pp
}

p_ben <- ggarrange(plotlist = plot_list, nrow = 4, ncol =3)

#### insets 
inset.ben <- ggplot(ben_sloan_df, aes(x = Neutral, fill = Neutral)) +
  geom_bar() + xlim("Above", "Neutral", "Below") + 
  labs(title="#ASVs") + theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size=8), axis.title.x=element_blank(), axis.title.y=element_blank())

inset.euph <- ggplot(euph_sloan_df, aes(x = Neutral, fill = Neutral)) +
  geom_bar() + xlim("Above", "Neutral", "Below") + 
  labs(title="#ASVs") + theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size=8), axis.title.x=element_blank(), axis.title.y=element_blank())

inset.aph <- ggplot(aph_sloan_df, aes(x = Neutral, fill = Neutral)) +
  geom_bar() + xlim("Above", "Neutral", "Below") + 
  labs(title="#ASVs") + theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size=8), axis.title.x=element_blank(), axis.title.y=element_blank())

distrib_ben.inset <- ggdraw() +
  draw_plot(distrib_ben + theme(legend.position="none")) +
  draw_plot(inset.ben + theme(legend.position="none"), x = 0.16, y = .6, width = .3, height = .3)

distrib_euph.inset <- ggdraw() +
  draw_plot(distrib_euph + theme(legend.position="none")) +
  draw_plot(inset.euph + theme(legend.position="none"), x = 0.16, y = .6, width = .3, height = .3)

distrib_aph.inset <- ggdraw() +
  draw_plot(distrib_aph + theme(legend.position="none")) +
  draw_plot(inset.aph + theme(legend.position="none"), x = 0.16, y = .6, width = .3, height = .3)


### ggarrange
p_sloan <- ggarrange(distrib_euph.inset, distrib_aph.inset, distrib_ben.inset, p_euph, p_aph, p_ben, nrow = 2, ncol = 3)
ggsave("Figure3/Figure_S10_Sloan_neutrals_models.pdf", p_sloan, width= 16, height = 10)
ggsave("Figure3/Figure_S10_Sloan_neutrals_models.png", p_sloan, width= 16, height = 10, dpi = 300)

###############################################################################################################
###############################################################################################################
###############################################################################################################


