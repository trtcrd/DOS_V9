##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Figure 4, 5 - Figure S11 - Table S6
##################################################################################


## setting path **to be adjusted*
setwd("~/path/to/be/adjusted/")

path.fig4 <- paste0(getwd(), "/Figure4")
path.fig5 <- paste0(getwd(), "/Figure5")
if(!dir.exists(path.fig4)) dir.create(path.fig4)
if(!dir.exists(path.fig5)) dir.create(path.fig5)

## custom functions **to be adjusted to the companion_functions.R path**
source("~/path/to/be/adjusted/companionFunctions.R")

## libraries
library('doMC')
registerDoMC(cores = 16)
library(ggplot2)
library(ggpubr)
library(vegan)
library(mgcv)
library(ranger)
library(mixOmics)

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

# concatenated into vectors
df_abund <- readRDS("DF_biomes_abundance.rds")
df_rich  <- readRDS("DF_biomes_richness.rds")

# sediment asv table
asv_sediment <- readRDS("asv_table_sediment_all.rds")
comp_sediment <- comp[rownames(asv_sediment),]

# asvs vectors 
asvs_euph <- readRDS("asvs_euph.rds")
asvs_aph <- readRDS("asvs_aph.rds")
asvs_sed_only <- readRDS("asvs_sed_only.rds")
asvs_plan_sed <- readRDS("asvs_plan_sed.rds")

# check mapping
identical(rownames(asv_sediment), rownames(comp_sediment)) # TRUE
identical(colnames(asv), rownames(taxo)) # TRUE
identical(colnames(df_abund), rownames(taxo)) # TRUE
identical(colnames(df_rich), rownames(taxo))  # TRUE

###############################################################################################################
###############################################################################################################
###############################################################################################################

## Figure 4

#### Attributes of ASVs between sinking / not sinking
# make df to plot 
sinking_dof <- data.frame(t(df_abund), 
                          Pelagic = colSums(df_abund[c("Euphotic", "Aphotic"),]),
                          Sinking = taxo$Sinking,
                          Taxo_group = taxo$taxo_group,
                          Size = taxo$sizes,
                          Funct_group = taxo$functional_groups, stringsAsFactors = F)
                          
### Abundance
p_abund_sink <- ggplot(sinking_dof, aes(x = Sinking, y = Pelagic, fill = Sinking)) +
  scale_x_discrete(limits = c("yes", "no")) + 
  geom_violin() + 
  scale_y_continuous(trans = 'log10') + 
  labs(y = "Pelagic ASVs abundance (#Reads)", x = "Sinking to the deep ocean floor") + 
  stat_summary(fun.data=data_summary, geom="pointrange", color="blue") + 
  stat_compare_means(comparisons = list( c("yes", "no") ), label = "p.signif") +
  theme_classic() + guides(fill = FALSE)

### Size distrib
p_size_sink <- ggplot(sinking_dof, aes(x = Sinking, y = as.numeric(Size), fill = Sinking)) +
  scale_x_discrete(limits = c("yes", "no")) + 
  geom_violin() + scale_y_continuous(trans = 'log10') + 
  labs(y = "Pelagic ASVs size distribution (μm)", x = "Sinking to the deep ocean floor") + 
  stat_summary(fun.data=data_summary, geom="pointrange", color="blue") + 
  stat_compare_means(comparisons = list( c("yes", "no") ), label = "p.signif") +
  theme_classic() + guides(fill = FALSE)

### Functional groups
## aggregete pelagic per functional groups sinking / not sinking
sinking_dof$Pelagic_binary <- ifelse(sinking_dof$Pelagic > 0, 1, 0)
agg_func <- aggregate(sinking_dof$Pelagic_binary, by = list(sinking_dof$Funct_group, sinking_dof$Sinking), FUN = sum)
colnames(agg_func) <- c("Funct", "Sinking", "Reads")
agg_func$proportion <- NA
agg_func[agg_func$Sinking == "yes", "proportion"] <- agg_func[agg_func$Sinking == "yes", "Reads"] *100 / sum(agg_func[agg_func$Sinking == "yes", "Reads"])
agg_func[agg_func$Sinking == "no", "proportion"] <- agg_func[agg_func$Sinking == "no", "Reads"] *100 / sum(agg_func[agg_func$Sinking == "no", "Reads"])
agg_func$Funct <- factor(agg_func$Funct, levels=rev(c("Parasite", "Phototrophs", "Photosymbionts", "Other heterotrophs", "Copepoda","Other metazoa","Unknown")))
agg_func$Sinking <- factor(agg_func$Sinking, levels=c("yes", "no"))

p_funct <- ggplot(agg_func, aes(x = Sinking, y = proportion, fill = Funct)) + 
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette="Set3") + 
  ylab("Proportion of ASVs") + xlab("Sinking to the deep ocean floor") +
  theme_classic()


#### Attributes of sinking plankton along latitudinal gradient in the sediment 
# data

asv_plank_sed <- asv_sediment[,asvs_plan_sed]

## aggreg
taxoLat_agg <- aggregate(t(asv_plank_sed), by = list(sapply(strsplit(taxo[asvs_plan_sed,"taxo_group"], ";"), `[`, 3)), FUN = sum)
rownames(taxoLat_agg) <- taxoLat_agg[,1]
taxoLat_agg <- taxoLat_agg[,2:ncol(taxoLat_agg)]
taxoLat_agg <- t(taxoLat_agg)
## rare ones
keep <- colSums(taxoLat_agg) / sum(taxoLat_agg) > 0.02
taxoLat_agg_ <- cbind(taxoLat_agg[,keep], Others = rowSums(taxoLat_agg[,keep == F])) 
## latitudes bins
taxoLat_agg_lat <- data.frame(lat60 = colMeans(subset(taxoLat_agg_, comp_sediment$latitude > 60) *100 / rowSums(subset(taxoLat_agg_, comp_sediment$latitude > 60))),
                              lat30 = colMeans(subset(taxoLat_agg_, comp_sediment$latitude <= 60 & comp_sediment$latitude > 30) *100 / rowSums(subset(taxoLat_agg_, comp_sediment$latitude <= 60 & comp_sediment$latitude > 30))),
                              lat0  = colMeans(subset(taxoLat_agg_, comp_sediment$latitude <= 30 & comp_sediment$latitude > 0) *100 / rowSums(subset(taxoLat_agg_, comp_sediment$latitude <= 30 & comp_sediment$latitude > 0))),
                              lat40  = colMeans(subset(taxoLat_agg_, comp_sediment$latitude <= 0 & comp_sediment$latitude > -40) *100 / rowSums(subset(taxoLat_agg_, comp_sediment$latitude <= 0 & comp_sediment$latitude > -40))),
                              latSouth  = colMeans(subset(taxoLat_agg_, comp_sediment$latitude < -40) *100 / rowSums(subset(taxoLat_agg_, comp_sediment$latitude < -40))))

taxoLat_agg_lat <- taxoLat_agg_lat[rownames(taxoLat_agg_lat)[order(rowSums(taxoLat_agg_lat), decreasing = T)],]

rownames(taxoLat_agg_lat) <- gsub("Unassigned_3", "Unassigned", rownames(taxoLat_agg_lat))
rownames(taxoLat_agg_lat) <- gsub("Others", "Others < 2%", rownames(taxoLat_agg_lat))

taxoLat_agg_lat_gr <- taxoLat_agg_lat[!rownames(taxoLat_agg_lat) %in% c("Unassigned", "Others < 2%"),] 
taxoLat_agg_lat <- rbind(taxoLat_agg_lat_gr, taxoLat_agg_lat[c("Others < 2%", "Unassigned"),])

dftax <- reshape::melt(t(taxoLat_agg_lat), id.vars = rownames(taxoLat_agg_lat))
colnames(dftax) <- c("latitude", "taxa", "prop") 

col <- colorRampPalette(c("dodgerblue", "gold2", "darkred", "darkolivegreen1", 
                          "firebrick1", "darkseagreen1", "bisque3", "maroon3", 
                          "blue", "green", "black", "grey"), 
                        bias=1, interpolate = "linear")(nrow(taxoLat_agg_lat))

dftax$taxa <- factor(dftax$taxa, levels=rev(rownames(taxoLat_agg_lat)))

p_tax <- ggplot(dftax, aes(x = latitude, y = prop, fill = taxa)) + 
  geom_bar(stat = "identity") + 
  xlim(rev(c("lat60","lat30","lat0", "lat40","latSouth"))) +
  scale_fill_manual(values = rev(col)) +
  ylab("Relative abundance") + xlab("") +
  theme_classic() + coord_flip()


### Size 
asvs_lat_size_lat_60 <- subset(asv_plank_sed, comp_sediment$latitude > 60)[,colSums(subset(asv_plank_sed, comp_sediment$latitude > 60)) > 0]
asvs_lat_size_lat_30 <- subset(asv_plank_sed, comp_sediment$latitude <= 60 & comp_sediment$latitude > 30)[,colSums(subset(asv_plank_sed, comp_sediment$latitude <= 60 & comp_sediment$latitude > 30)) > 0]
asvs_lat_size_lat_0 <- subset(asv_plank_sed, comp_sediment$latitude <= 30 & comp_sediment$latitude > 0)[,colSums(subset(asv_plank_sed, comp_sediment$latitude <= 30 & comp_sediment$latitude > 0)) > 0]
asvs_lat_size_lat_40 <- subset(asv_plank_sed, comp_sediment$latitude <= 0 & comp_sediment$latitude > -40)[,colSums(subset(asv_plank_sed, comp_sediment$latitude <= 0 & comp_sediment$latitude > -40)) > 0]
asvs_lat_size_lat_South <- subset(asv_plank_sed, comp_sediment$latitude < -40)[,colSums(subset(asv_plank_sed, comp_sediment$latitude < -40)) > 0]

lat_size_dist <- data.frame(rbind(cbind(ASV_seq = colnames(asvs_lat_size_lat_60), Size = taxo[colnames(asvs_lat_size_lat_60),"sizes"], Bins = "lat60", Abund = colSums(asvs_lat_size_lat_60)),
                                  cbind(ASV_seq = colnames(asvs_lat_size_lat_30), Size = taxo[colnames(asvs_lat_size_lat_30),"sizes"], Bins = "lat30", Abund = colSums(asvs_lat_size_lat_30)),
                                  cbind(ASV_seq = colnames(asvs_lat_size_lat_0), Size = taxo[colnames(asvs_lat_size_lat_0),"sizes"], Bins = "lat0", Abund = colSums(asvs_lat_size_lat_0)),
                                  cbind(ASV_seq = colnames(asvs_lat_size_lat_40), Size = taxo[colnames(asvs_lat_size_lat_40),"sizes"], Bins = "lat40", Abund = colSums(asvs_lat_size_lat_40)),
                                  cbind(ASV_seq = colnames(asvs_lat_size_lat_South), Size = taxo[colnames(asvs_lat_size_lat_South),"sizes"], Bins = "latSouth", Abund = colSums(asvs_lat_size_lat_South))), stringsAsFactors = F)
lat_size_dist$Size <- as.numeric(lat_size_dist$Size)
lat_size_dist$Abund <- as.numeric(lat_size_dist$Abund)

## make size bins 
size_bins <- c()
for (i in 1:nrow(lat_size_dist)) {
  if (!is.na(lat_size_dist$Size[i])) {
    if (lat_size_dist$Size[i] >= 100) size_bins <- c(size_bins, ">100μm")
    if (lat_size_dist$Size[i] < 100 & lat_size_dist$Size[i] >= 20) size_bins <- c(size_bins, ">20μm")
    if (lat_size_dist$Size[i] < 20 & lat_size_dist$Size[i] >= 5) size_bins <- c(size_bins, ">5μm")
    if (lat_size_dist$Size[i] < 5) size_bins <- c(size_bins, "<5μm")
  } else { size_bins <- c(size_bins, NA) }
}
lat_size_dist$Size_bins <- size_bins
lat_size_dist$Size_bins <- factor(lat_size_dist$Size_bins, levels = c(">100μm", ">20μm", ">5μm", "<5μm"))

## aggregate abundance per size_bins and lat_bins
size_bins_agg <- aggregate(lat_size_dist$Abund, by = list(lat_size_dist$Size_bins, lat_size_dist$Bins), FUN = sum)
colnames(size_bins_agg) <- c("Size", "Latitude", "x")
# relative abundance at each lab bin
size_bins_agg$prop <- NA
size_bins_agg[size_bins_agg$Latitude == "lat0", "prop"] <- size_bins_agg[size_bins_agg$Latitude == "lat0", "x"] / sum(size_bins_agg[size_bins_agg$Latitude == "lat0", "x"])
size_bins_agg[size_bins_agg$Latitude == "lat30", "prop"] <- size_bins_agg[size_bins_agg$Latitude == "lat30", "x"] / sum(size_bins_agg[size_bins_agg$Latitude == "lat30", "x"])
size_bins_agg[size_bins_agg$Latitude == "lat40", "prop"] <- size_bins_agg[size_bins_agg$Latitude == "lat40", "x"] / sum(size_bins_agg[size_bins_agg$Latitude == "lat40", "x"])
size_bins_agg[size_bins_agg$Latitude == "lat60", "prop"] <- size_bins_agg[size_bins_agg$Latitude == "lat60", "x"] / sum(size_bins_agg[size_bins_agg$Latitude == "lat60", "x"])
size_bins_agg[size_bins_agg$Latitude == "latSouth", "prop"] <- size_bins_agg[size_bins_agg$Latitude == "latSouth", "x"] / sum(size_bins_agg[size_bins_agg$Latitude == "latSouth", "x"])

p_size_bins <- ggplot(size_bins_agg, aes(x = Latitude, y = prop, fill = Size)) +
  geom_bar(stat = "identity") +
  xlim(rev(c("lat60","lat30","lat0", "lat40","latSouth"))) + 
  theme_classic() + coord_flip() + ylab("Relative abundance")


### Functional variation latitude
funct_agg <- aggregate(t(asv_plank_sed), by = list(taxo[asvs_plan_sed,"functional_groups"]), FUN = sum)
rownames(funct_agg) <- funct_agg[,1]
funct_agg <- funct_agg[,2:ncol(funct_agg)] 
funct_agg <- t(funct_agg)
funct_lat <- data.frame(lat60 = colMeans(subset(funct_agg, comp_sediment$latitude > 60) *100 / rowSums(subset(funct_agg, comp_sediment$latitude > 60))),
                        lat30 = colMeans(subset(funct_agg, comp_sediment$latitude <= 60 & comp_sediment$latitude > 30) *100 / rowSums(subset(funct_agg, comp_sediment$latitude <= 60 & comp_sediment$latitude > 30))),
                        lat0  = colMeans(subset(funct_agg, comp_sediment$latitude <= 30 & comp_sediment$latitude > 0) *100 / rowSums(subset(funct_agg, comp_sediment$latitude <= 30 & comp_sediment$latitude > 0))),
                        lat40  = colMeans(subset(funct_agg, comp_sediment$latitude <= 0 & comp_sediment$latitude > -40) *100 / rowSums(subset(funct_agg, comp_sediment$latitude <= 0 & comp_sediment$latitude > -40))),
                        latSouth  = colMeans(subset(funct_agg, comp_sediment$latitude < -40) *100 / rowSums(subset(funct_agg, comp_sediment$latitude < -40))))


df <- reshape::melt(t(funct_lat), id.vars = rownames(funct_lat))
colnames(df) <- c("latitude", "funct", "prop") 

df$funct <- factor(df$funct, levels=rev(c("Parasite", "Phototrophs", "Photosymbionts", "Other heterotrophs", "Copepoda","Other metazoa","Unknown")))

p_funct_lat <- ggplot(df, aes(x = latitude, y = prop, fill = funct)) + 
  geom_bar(stat = "identity") + 
  xlim(rev(c("lat60","lat30","lat0", "lat40","latSouth"))) + 
  scale_fill_brewer(palette="Set3") + 
  ylab("Relative abundance") + xlab("") +
  theme_classic() + coord_flip()


##### ggarrange
p <- ggarrange(p_abund_sink, p_size_sink, p_funct, p_tax, p_size_bins, p_funct_lat, ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))

pdf("Figure4/Figure_4_functional_attributes_sinking_plankton.pdf", p, width = 12, height = 6.5, useDingbats = F)
print(p)
dev.off()


###############################################################################################################
###############################################################################################################
###############################################################################################################

##### Figure 5A

prop_plan_sed <- data.frame(cbind(Prop_reads = rowSums(asv_sediment[,asvs_plan_sed]) *100 / rowSums(asv_sediment), 
                                  latitude = comp_sediment$latitude, 
                                  depth = as.character(comp_sediment$depth),
                                  POC_seafloor = comp_sediment$POC_seafloor,
                                  POC_export = comp_sediment$POC_export), stringsAsFactors = F)

prop_plan_sed[] <- lapply(prop_plan_sed,as.numeric)
prop_plan_sed$station <- comp_sediment$station

### Relative abundance of sinking plankton in the sediment
gam_mod <- mgcv::gam(prop_plan_sed$Prop_reads ~ s(prop_plan_sed$latitude, k=3))
summ_gam <- summary(gam_mod)
ex_dev <- round(summ_gam$dev.expl * 100, 1)
p <- summ_gam$s.table[,'p-value']
if(p < 0.05) {
  if(p >= 0.01) sig <- "*"
  if(p >= 0.001) sig <- "**"
  if(p < 0.001) sig <- "***"
} else {sig <- "ns" }
summ_gam <- paste0("Expl. deviance: ", ex_dev,"%",sig)

p_gam <- ggplot(prop_plan_sed, aes(x=latitude, y=Prop_reads)) +
  geom_point(color = alpha("orange3", 0.3)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k=3), size = 1, alpha = 0.2) + 
  scale_x_continuous(breaks = c(-60,-30,0,30,30,60,90)) + 
  geom_vline(xintercept=0, linetype="dashed") + 
  theme(legend.position="none") + theme_classic() + ylab("Proportion of planktonic reads in the sediment") + xlab("Latitude") + ylim(c(0,100)) +
  annotate("text", x = min(prop_plan_sed$latitude), y = min(prop_plan_sed$Prop_reads), label = summ_gam, hjust = 0) +
  coord_flip()

ggsave("Figure5/Figure_5A_prop_plankton_lat.pdf", p_gam, width = 4, height = 5)


## removing coastal outliers
prop_poc <- subset(prop_plan_sed, prop_plan_sed$POC_export < 80)
prop_sea <- subset(prop_plan_sed, prop_plan_sed$POC_seafloor < 30) 

summary(lm(prop_poc$POC_export ~ prop_poc$Prop_reads))
summary(lm(prop_sea$POC_seafloor ~ prop_sea$Prop_reads))

ggplot(prop_poc, aes(x=POC_export, y=Prop_reads)) + 
  geom_point() + stat_smooth(method='lm')
ggplot(prop_sea, aes(x=POC_seafloor, y=Prop_reads)) + 
  geom_point() + stat_smooth(method='lm')


###### Figure 5B - predicting the POC export and POC seafloor from sinking plankton at the station level

## POC export 
## Aggregated plankton in sediment at the station level 
asv_plank_agg_station <- aggregate(asv_plank_sed, by = list(comp_sediment$station), FUN = sum)
rownames(asv_plank_agg_station) <- asv_plank_agg_station$Group.1
asv_plank_agg_station <- asv_plank_agg_station[,2:ncol(asv_plank_agg_station)]

## making the comp 
comp_plank_agg_station <- c()
for (i in 1:length(rownames(asv_plank_agg_station))) {
  tmp <- rownames(asv_plank_agg_station)[i]
  tmp_ <- subset(comp_sediment, comp_sediment$station == tmp)
  if (nrow(tmp_) > 1) {
    comp_plank_agg_station <- rbind(comp_plank_agg_station, tmp_[1,])
  } else {
    comp_plank_agg_station <- rbind(comp_plank_agg_station, tmp_)
  }
}
identical(as.character(comp_plank_agg_station$station), rownames(asv_plank_agg_station)) # TRUE
rownames(comp_plank_agg_station) <- rownames(asv_plank_agg_station)

# remove coastal stations to focus on abyssal plains
comp_plank_agg_station <- subset(comp_plank_agg_station, !comp_plank_agg_station$station %in% c("Lie_ST7","MDW_ST38", "EssNaut_PL06-07", "EssNaut_PL11-12", "CANHROV_CT") & !is.na(comp_plank_agg_station$POC_seafloor))


asv_plank_agg_station <- asv_plank_agg_station[rownames(comp_plank_agg_station),]
asv_plank_agg_station <- asv_plank_agg_station[,colSums(asv_plank_agg_station)>0]

## cssNorm
POC_export_ASVs <- cssNorm(asv_plank_agg_station, samples_as_rows = T)
POC_export_COMP <- comp_plank_agg_station

## station LOOCV
combined_exp_rf_pred <- c()
combined_exp_rf_obs <- c()


# for each station
for (j in 1:nrow(POC_export_COMP)) {
  ## subset 
  ASVtr <- data.frame(POC_export_ASVs[-j,])
  ASVte <- t(data.frame(POC_export_ASVs[j,]))
  
  POC_tr <- data.frame(POC_export_COMP[-j,])
  POC_te <- data.frame(POC_export_COMP[j,])
  
  ## training RF model
  set.seed(1234)
  rmod <- ranger(POC_tr$POC_export ~ ., data = ASVtr, mtry=floor(ncol(ASVtr) / 3), num.trees = 300, importance= "permutation", write.forest = T)
  ## predict
  preds <- predict(rmod, ASVte)
  combined_exp_rf_pred <- c(combined_exp_rf_pred, preds$predictions)
  combined_exp_rf_obs  <- c(combined_exp_rf_obs, POC_te$POC_export)
  message(paste0(j, " / ", nrow(POC_export_COMP)))
}

out_lm <- data.frame(Predicted = combined_exp_rf_pred,  Observed = combined_exp_rf_obs)
coef_exp <- summary(lm(out_lm$Predicted ~ out_lm$Observed)) 
p <- coef_exp$coefficients[,"Pr(>|t|)"][2]
if(p < 0.05) {
  if(p >= 0.01) sig <- "*"
  if(p >= 0.001) sig <- "**"
  if(p < 0.001) sig <- "***"
} else {sig <- "ns" }

coef_exp <- paste0("LOOCV R²=", round(coef_exp$adj.r.squared, 2), sig)

p_exp_sta <- ggplot(out_lm, aes(x = Observed, y = Predicted)) + 
  geom_point() + stat_smooth(method = "lm") +
  theme_classic() + xlab("Observed POC export") + ylab("Predicted POC export") + 
  annotate("text", x = min(out_lm$Observed), y = max(out_lm$Predicted), label = coef_exp, hjust = 0, size=3)


# POC seafloor
POC_floor_ASVs <- POC_export_ASVs
POC_floor_COMP <- POC_export_COMP

## station LOOCV
combined_flo_rf_pred <- c()
combined_flo_rf_obs <- c()


# for each station
for (j in 1:nrow(POC_floor_ASVs)) {
  ## subset 
  ASVtr <- data.frame(POC_floor_ASVs[-j,])
  ASVte <- t(data.frame(POC_floor_ASVs[j,]))
  
  POC_tr <- data.frame(POC_floor_COMP[-j,])
  POC_te <- data.frame(POC_floor_COMP[j,])
  
  ## training RF model
  set.seed(1234)
  rmod <- ranger(POC_tr$POC_seafloor ~ ., data = ASVtr, mtry=floor(ncol(ASVtr) / 3), num.trees = 300, importance= "permutation", write.forest = T)
  ## predict
  preds <- predict(rmod, ASVte)
  combined_flo_rf_pred <- c(combined_flo_rf_pred, preds$predictions)
  combined_flo_rf_obs  <- c(combined_flo_rf_obs, POC_te$POC_seafloor)
  message(paste0(j, " / ", nrow(POC_floor_COMP)))
}

out_lm <- data.frame(Predicted = combined_flo_rf_pred,  Observed = combined_flo_rf_obs)
coef_exp <- summary(lm(out_lm$Predicted ~ out_lm$Observed)) 
p <- coef_exp$coefficients[,"Pr(>|t|)"][2]
if(p < 0.05) {
  if(p >= 0.01) sig <- "*"
  if(p >= 0.001) sig <- "**"
  if(p < 0.001) sig <- "***"
} else {sig <- "ns" }

coef_exp <- paste0("LOOCV R²=", round(coef_exp$adj.r.squared, 2), sig)

p_floor_sta <- ggplot(out_lm, aes(x = Observed, y = Predicted)) + 
  geom_point() + stat_smooth(method = "lm") +
  theme_classic() + xlab("Observed POC seafloor") + ylab("Predicted POC seafloor") + 
  annotate("text", x = min(out_lm$Observed), y = max(out_lm$Predicted), label = coef_exp, hjust = 0, size=3)

p_poc <- ggarrange(p_exp_sta, p_floor_sta, ncol = 2)


pdf("Figure5/Figure_5B_poc_predictions.pdf", width = 4, height = 2, useDingbats = F)
print(p_poc)
dev.off()




########################## 
### sPLS analysis

POC_export_ASVs <- asv_plank_agg_station
POC_export_meta <- comp_plank_agg_station[,c("latitude", "POC_export","POC_seafloor", "primprod")]  #
POC_export_meta$latitude <- abs(POC_export_meta$latitude)
POC_export_meta <- decostand(POC_export_meta, method = "standardize")

ncomp = 15
result.spls <- spls((POC_export_ASVs +1), POC_export_meta, ncomp = ncomp, keepX = c(rep(ncomp, ncomp)), mode = 'regression', logratio='CLR')
cimplot <- cim(result.spls, comp = 1:ncomp, margins = c(7, 16))

## export table of correlation 
out_spls <- cbind(ASV_ID = paste0("ASV", taxo[rownames(cimplot$mat),"ASV_ID"]),
                  cimplot$mat,
                  taxo[rownames(cimplot$mat),c("functional_groups", "taxo_group","taxon","mean.similarity", "reference.ids")])

write.table(out_spls, "Figure4/Table_S6_sPLS_results.tsv", sep = "\t", quote = F, dec = ",")

## rename the ASVs in table to export the heatmap with formated taxo names
POC_names <- colnames(POC_export_ASVs)
POC_names <- paste0(taxo[POC_names,"ASV_ID"], " ",
                    sapply(strsplit(taxo[POC_names,"taxo_group"], ";"), `[`, 3), " ", 
                    unlist(sapply(strsplit(taxo[POC_names,"taxon"], ";"), function(x)x[length(x)])))
                    #ifelse(taxo[POC_names,"mean.similarity"] > 0.95, unlist(sapply(strsplit(taxo[POC_names,"taxon"], ";"), function(x)x[length(x)])), ""))
POC_names <- gsub("_X", "",POC_names)
POC_names <- gsub("[+]", " ",POC_names)
POC_names <- gsub("[+]", " ",POC_names)
POC_names <- gsub("Unassigned_3 unassigned", "Unassigned",POC_names)

colnames(POC_export_ASVs) <- POC_names
result.spls <- spls((POC_export_ASVs +1), POC_export_meta, ncomp = ncomp, keepX = c(rep(ncomp, ncomp)), mode = 'regression', logratio='CLR')
cimplot <- cim(result.spls, comp = 1:ncomp, margins = c(7, 16))

## and again for spotting the ASVs to zoom on
new_names <- c()
for (i in colnames(POC_export_ASVs)) {
  if (i %in% rownames(cimplot$mat)) {
    tmp <- cimplot$mat[i,]
    if (length(tmp) > 0) {
      if (tmp["POC_export"] > 0.3 & tmp["POC_seafloor"] < 0.3) {
        new_names <- c(new_names, paste0("_exp_", i))
      } else if (tmp["POC_export"] > 0.3 & tmp["POC_seafloor"] > 0.3) {
        new_names <- c(new_names, paste0("_exp_flo_", i))
      } else if (tmp["POC_export"] < 0.3 & tmp["POC_seafloor"] > 0.3) {
        new_names <- c(new_names, paste0("_flo_", i))
      } else new_names <- c(new_names, i)
    } 
  } else new_names <- c(new_names, i)
}

colnames(POC_export_ASVs) <- new_names
result.spls <- spls((POC_export_ASVs +1), POC_export_meta, ncomp = ncomp, keepX = c(rep(ncomp, ncomp)), mode = 'regression', logratio='CLR')

### Heatmap to be adjusted and zoom on ASVs above 0.3
pdf("Figure5/Figure_5C_sPLS_heatmap.pdf", width = 8, height = 28)
cim(result.spls, comp = 1:ncomp, margins = c(7, 16), mapping = "XY", cluster = "row", threshold = 0.1)
dev.off()



###############################################################################################################
###############################################################################################################
###############################################################################################################
## Figure S11 ranked relative abundance profiles across biomes (coarse taphonomy)

## Euphotic pico nano size fractions
comp_euph_pico_nano <- subset(comp, comp$biome == "Euphotic" & (comp$size_class == "nano" | comp$size_class == "pico")) 
asv_euph_pico_nano <- asv[rownames(comp_euph_pico_nano),]
## Aphotic pico nano
comp_aph_pico_nano <- subset(comp, comp$biome == "Aphotic" & (comp$size_class == "nano" | comp$size_class == "pico")) 
asv_aph_pico_nano <- asv[rownames(comp_aph_pico_nano),]

## Merge vectors with sediment
df_biome_norm <- data.frame(Euphotic = colMeans(decostand(asv_euph_pico_nano, method = "total")),
                            Aphotic  = colMeans(decostand(asv_aph_pico_nano, method = "total")),
                            Sediment = colMeans(decostand(asv_sediment, method = "total")))

#### we fetch all planktonic ASVs found in the sediment
tmp_taxo <- taxo[asvs_plan_sed,]  
## but get the pelagic ASVs within the pico nano fraction
asvs_euph_pico_nano <- colnames(asv_euph_pico_nano[,colSums(asv_euph_pico_nano)>0])
asvs_aph_pico_nano <- colnames(asv_aph_pico_nano[,colSums(asv_aph_pico_nano)>0])
## taxo for pelagic 
tmp_taxo_pel <- taxo[asvs_aph_pico_nano[asvs_aph_pico_nano %in% asvs_euph_pico_nano],]  


## taxo groups to be checked 
tax_sink <- c("All sinking ASVs", "Copepoda", "Dinophyceae", "Diatomeae", "Acantharea", "Collodaria", "Eupelagonemidae", "MALV-I","MALV-II",  
              "Hydrozoa", "Spirotrichea", "Spumellaria", "Chrysophyceae", "Globigerinacea") # "Prymnesiophyceae"
out_sink <- array(NA, c(length(tax_sink), 9))
colnames(out_sink)<- c("#ASVs_euph_ratio", "adjR_euph", "p_euph", "#ASVs_aph_ratio", "adjR_aph", "p_aph", "#ASVs_euph_aph_ratio", "adjR_euph_aph", "p_euph_aph")
rownames(out_sink)<- tax_sink

#### we fetch all planktonic ASVs found in the sediment
tmp_taxo <- taxo[asvs_plan_sed,]  
## but get the pelagic ASVs within the pico nano fraction
asvs_euph_pico_nano <- colnames(asv_euph_pico_nano[,colSums(asv_euph_pico_nano)>0])
asvs_aph_pico_nano <- colnames(asv_aph_pico_nano[,colSums(asv_aph_pico_nano)>0])
## taxo for pelagic 
tmp_taxo_pel <- taxo[asvs_aph_pico_nano[asvs_aph_pico_nano %in% asvs_euph_pico_nano],]  


out_plots <- list()
j <- 1
for (i in tax_sink) {
  ### if all, manually get all the asvs_plan_sed
  if (i == "All sinking ASVs") {
    tmp_taxo_ <- tmp_taxo
    tmp_taxo_pel_ <- tmp_taxo_pel
  } else if (i != "Globigerinacea") {  ### Globigerinacea is not in the taxo group column, but in the taxon column
    tmp_taxo_ <- tmp_taxo[grep(i, tmp_taxo$taxo_group),]
    tmp_taxo_pel_ <- tmp_taxo_pel[grep(i, tmp_taxo_pel$taxo_group),]
  } else {
    tmp_taxo_ <- tmp_taxo[grep(i, tmp_taxo$taxon),]
    tmp_taxo_pel_ <- tmp_taxo_pel[grep(i, tmp_taxo_pel$taxon),]
  }
  # get the sequences 
  seqs_euph <- rownames(tmp_taxo_)[rownames(tmp_taxo_) %in% asvs_euph_pico_nano]
  seqs_aph <- rownames(tmp_taxo_)[rownames(tmp_taxo_) %in% asvs_aph_pico_nano]
  seqs_euph_aph <- rownames(tmp_taxo_pel_)[rownames(tmp_taxo_pel_) %in% asvs_euph_pico_nano]
  
  ## ratio (compared to all ASVs of the biome)
  ## Euph -> Sed
  tmp_taxo_euph <- taxo[asvs_euph_pico_nano,]  
  if (i == "All sinking ASVs") {
    tmp_taxo_euph <- tmp_taxo_euph
  } else if (i != "Globigerinacea") {  
    tmp_taxo_euph <- tmp_taxo_euph[grep(i, tmp_taxo_euph$taxo_group),]
  } else {
    tmp_taxo_euph <- tmp_taxo_euph[grep(i, tmp_taxo_euph$taxon),]
  }
  ratio_sum_euph <- round(sum(asv_euph_pico_nano[,seqs_euph]) *100 / sum(asv_euph_pico_nano[,rownames(tmp_taxo_euph)]), 1)
  ## Aph -> Sed
  tmp_taxo_aph <- taxo[asvs_aph_pico_nano,]  
  if (i == "All sinking ASVs") {
    tmp_taxo_aph <- tmp_taxo_aph
  } else if (i != "Globigerinacea") {
    tmp_taxo_aph <- tmp_taxo_aph[grep(i, tmp_taxo_aph$taxo_group),]
  } else {
    tmp_taxo_aph <- tmp_taxo_aph[grep(i, tmp_taxo_aph$taxon),]
  }
  ratio_sum_aph <- round(sum(asv_aph_pico_nano[,seqs_aph]) *100 / sum(asv_aph_pico_nano[,rownames(tmp_taxo_aph)]), 1)
  ## Euph -> Aph
  tmp_taxo_EA <- taxo[asvs_euph_pico_nano,]  
  if (i == "All sinking ASVs") {
    tmp_taxo_EA <- tmp_taxo_EA
  } else if (i != "Globigerinacea") {
    tmp_taxo_EA <- tmp_taxo_EA[grep(i, tmp_taxo_EA$taxo_group),]
  } else {
    tmp_taxo_EA <- tmp_taxo_EA[grep(i, tmp_taxo_EA$taxon),]
  }
  ratio_sum_EA <- round(sum(asv_euph_pico_nano[,seqs_euph_aph]) *100 / sum(asv_euph_pico_nano[,rownames(tmp_taxo_EA)]), 1)
  
  out_sink[i,"#ASVs_euph_ratio"] <- paste0(length(seqs_euph), " sinking ASVs / ",nrow(tmp_taxo_euph))
  out_sink[i,"#ASVs_aph_ratio"] <- paste0(length(seqs_aph), " sinking ASVs / ",nrow(tmp_taxo_aph))
  out_sink[i,"#ASVs_euph_aph_ratio"] <- paste0(length(seqs_euph_aph), " sinking ASVs / ",nrow(tmp_taxo_EA))
  
  ## size of geom_point
  if (length(seqs_euph) >= 1000) size_points <- 0.2
  if (length(seqs_euph) < 1000) size_points <- 0.4
  if (length(seqs_euph) < 200) size_points <- 1
  
  # bind them 
  tplot_euph <- data.frame(cbind(Euphotic = df_biome_norm[seqs_euph,"Euphotic"], Sediment = df_biome_norm[seqs_euph,"Sediment"]))
  tplot_aph <- data.frame(cbind(Aphotic = df_biome_norm[seqs_aph,"Aphotic"], Sediment = df_biome_norm[seqs_aph,"Sediment"]))
  tplot_EA <- data.frame(cbind(Euphotic = df_biome_norm[seqs_euph_aph,"Euphotic"], Aphotic = df_biome_norm[seqs_euph_aph,"Aphotic"]))
  
  ## lm
  lmE <- summary(lm(log(tplot_euph$Euphotic) ~ log(tplot_euph$Sediment)))
  lmA <- summary(lm(log(tplot_aph$Aphotic) ~ log(tplot_aph$Sediment)))
  lmEA <- summary(lm(log(tplot_EA$Euphotic) ~ log(tplot_EA$Aphotic)))
  
  ## fetch the stats
  out_sink[i,"adjR_euph"] <- round(lmE$adj.r.squared, 3)
  out_sink[i,"p_euph"] <- lmE$coefficients[,"Pr(>|t|)"][2]
  if (as.numeric(out_sink[i,"p_euph"]) >= 0.05) sigE <- "ns"
  if (as.numeric(out_sink[i,"p_euph"]) < 0.05 & out_sink[i,"p_euph"] < 0.01) sigE <- "*"
  if (as.numeric(out_sink[i,"p_euph"]) < 0.01 & out_sink[i,"p_euph"] >= 0.001) sigE <- "**"
  if (as.numeric(out_sink[i,"p_euph"]) < 0.001) sigE <- "***"
  
  out_sink[i,"adjR_aph"] <- round(lmA$adj.r.squared, 3)
  out_sink[i,"p_aph"] <- lmA$coefficients[,"Pr(>|t|)"][2]
  if (as.numeric(out_sink[i,"p_aph"]) >= 0.05) sigA <- "ns"
  if (as.numeric(out_sink[i,"p_aph"]) < 0.05 & out_sink[i,"p_aph"] < 0.01) sigA <- "*"
  if (as.numeric(out_sink[i,"p_aph"]) < 0.01 & out_sink[i,"p_aph"] >= 0.001) sigA <- "**"
  if (as.numeric(out_sink[i,"p_aph"]) < 0.001) sigA <- "***"
  
  out_sink[i,"adjR_euph_aph"] <- round(lmEA$adj.r.squared, 3)
  out_sink[i,"p_euph_aph"] <- lmEA$coefficients[,"Pr(>|t|)"][2]
  if (as.numeric(out_sink[i,"p_euph_aph"]) >= 0.05) sigEA <- "ns"
  if (as.numeric(out_sink[i,"p_euph_aph"]) < 0.05 & out_sink[i,"p_euph_aph"] < 0.01) sigEA <- "*"
  if (as.numeric(out_sink[i,"p_euph_aph"]) < 0.01 & out_sink[i,"p_euph_aph"] >= 0.001) sigEA <- "**"
  if (as.numeric(out_sink[i,"p_euph_aph"]) < 0.001) sigEA <- "***"
  
  p1 <- ggplot(tplot_euph, aes(y = log(Sediment), x = log(Euphotic))) +
    geom_point(size = size_points, color = "orange3") +
    stat_smooth(method = "lm") + 
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(out_sink[i,"#ASVs_euph_ratio"], "\n",ratio_sum_euph,"% of reads, R²= ", out_sink[i,"adjR_euph"], sigE ))
  
  p2 <- ggplot(tplot_aph, aes(y = log(Sediment), x = log(Aphotic))) +
    geom_point(size = size_points, color = "orange3") +
    stat_smooth(method = "lm") + 
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(out_sink[i,"#ASVs_aph_ratio"], "\n",ratio_sum_aph,"% of reads, R²= ", out_sink[i,"adjR_aph"], sigA ))
  
  p3 <- ggplot(tplot_EA, aes(y = log(Aphotic), x = log(Euphotic))) +
    geom_point(size = size_points, color = "blue4") +
    stat_smooth(method = "lm") + 
    theme(plot.title = element_text(size = 10)) +
    labs(title = paste0(out_sink[i,"#ASVs_euph_aph_ratio"], "\n",ratio_sum_EA,"% of reads, R²= ", out_sink[i,"adjR_euph_aph"], sigEA ))
  
  
  p4 <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1) # + theme(plot.background = element_rect(fill = "transparent", color = "black", size =2))
  p4 <- annotate_figure(p4, top = text_grob(i, face = "bold", size = 13)) 
  
  out_plots[[j]] <- p4 
  j <- j+1
  
}


pdf("Figure4/Figure_S11_ranked_relative_abundances_LM_.pdf", width = 15, height = 18)
ggarrange(plotlist = out_plots, nrow=7, ncol=2)
dev.off()



###############################################################################################################
###############################################################################################################
###############################################################################################################






