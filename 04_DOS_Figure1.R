##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Figure 1 + S1 S2 S3; Table S3  
##################################################################################


## setting path **to be adjusted*
setwd("~/path/to/be/adjusted/")

path.fig <- paste0(getwd(), "/Figure1")
if(!dir.exists(path.fig)) dir.create(path.fig)

## custom functions **to be adjusted to the companion_functions.R path**
source("~/path/to/be/adjusted/companionFunctions.R")

## libraries
library('doMC')
registerDoMC(cores = 16)
library(vegan)
library(ggplot2)
library(ggpubr)
library(maps)
library(UpSetR)

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

### ASVs tables
euphotic_all <- readRDS("asv_table_euphotic_all.rds")
aphotic_all <- readRDS("asv_table_aphotic_all.rds")
sediment_all <- readRDS("asv_table_sediment_all.rds")

### ASVs vectors
asvs_euph <- readRDS("asvs_euph.rds")
asvs_aph  <- readRDS("asvs_aph.rds")
asvs_plankton  <- readRDS("asvs_plankton.rds")
asvs_sed_only <- readRDS("asvs_sed_only.rds")
add_asvs_sed_only <- readRDS("add_asvs_sed_only.rds")
asvs_plan_sed <- readRDS("asvs_plan_sed.rds")

###############################################################################################################
###############################################################################################################
###############################################################################################################
### Figure 1A: map

#### global sampling map
pdf(file = paste0("Figure1/Figure_1A_map.pdf"), width = 10, height=8, useDingbats = F)
map("world", border=NA, fill=T, col="grey", bg="white", ylim=c(-90, 90), mar=c(0,0,0,0))
points(comp[comp$biome == "Aphotic","longitude"],comp[comp$biome == "Aphotic","latitude"], bg=alpha("blue4", alpha = 1), pch=24, cex= .8)
points(comp[comp$biome == "Euphotic","longitude"],comp[comp$biome == "Euphotic","latitude"], bg=alpha("cadetblue1", alpha = 1), pch=21, cex= 0.65)
points(comp[comp$biome == "Sediment","longitude"],comp[comp$biome == "Sediment","latitude"], bg=alpha("darkgoldenrod1", alpha = 1), pch=23, cex= 1.2)
dev.off()

## inset depth distribution
p_inset <- ggplot(comp, aes(y=depth, fill=biome)) + 
  geom_histogram(binwidth = 400) +
  scale_fill_manual(values = c(Euphotic="cadetblue1", Aphotic="blue4", Sediment = "darkgoldenrod1")) +
  #scale_y_discrete(limits = seq(0,6400, by=400), trans = "reverse") + 
  scale_y_continuous(trans = "reverse", breaks = seq(0,6400, by=400)) +
  scale_x_discrete(limits = seq(0,1400, by=50), position = "top") + 
  theme_classic() + ylab("depth (m)") + xlab("# of samples") + 
  theme(legend.position = "none")

### inset to be adjusted within figure 1A...
ggsave("Figure1/Figure_1A_inset.pdf", p_inset, width = 8, height = 3)


### Figure 1B: pie charts + venn 
df <- rbind(Euphotic = colSums(euphotic_all),
            Aphotic  = colSums(aphotic_all),
            Sediment = colSums(sediment_all))

# turn to 0 the plankton ASVs that are obvious benthic animals... 
df["Euphotic",add_asvs_sed_only] <- 0
df["Aphotic",add_asvs_sed_only] <- 0

pie_df <- data.frame(Reads = rowSums(df), ASVs = specnumber(df), Biome = rownames(df))
pie_df$Biome <- factor(pie_df$Biome, levels = c("Euphotic", "Aphotic", "Sediment"))

## reads
p_reads <- ggplot(pie_df, aes(x="", y=Reads, fill=Biome)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=3*pi/2, direction = -1) +
  scale_fill_manual(values = c(Euphotic="cadetblue1", Aphotic="blue4", Sediment = "darkgoldenrod1")) + 
  theme_void() + theme(legend.position = "none")

## ASVs
p_asvs <- ggplot(pie_df, aes(x="", y=ASVs, fill=Biome)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=pi/2, direction = -1) +
  scale_fill_manual(values = c(Euphotic="cadetblue1", Aphotic="blue4", Sediment = "darkgoldenrod1")) + 
  theme_void() + theme(legend.position = "none")

ggsave("Figure1/Figure_1B_pieReads.pdf", p_reads, width = 3, height = 3)
ggsave("Figure1/Figure_1B_pieASVs.pdf", p_asvs, width = 3, height = 3)


## Venn (to be adjusted...)
ven <- limma::vennCounts(t(df))
pdf("Figure1/Figure_1B_venn_diagram.pdf")
limma::vennDiagram(ven, lwd=0.9)
dev.off()

### Table S3: aggregated ASVs richness and abondance per taxogroups
# add the benthic and sinking plankton
df <- rbind(df, 
            Benthic = df["Sediment",],
            Sinking_Plankton = df["Sediment",])
## muting plankton in sediment for benthic and benthic for sinking_plankton 
df["Benthic",asvs_plankton[!asvs_plankton %in% asvs_sed_only]] <- 0
df["Sinking_Plankton",asvs_sed_only] <- 0

df_rich <- df
df_rich[df_rich > 0] <- 1

df_rich_agg <- aggregate(t(df_rich), by= list(taxo[colnames(df_rich),"taxo_group"]), FUN = "sum")
df_abun_agg <- aggregate(t(df),      by= list(taxo[colnames(df_rich),"taxo_group"]), FUN = "sum")

### export tables 
write.table(df_rich_agg, "Figure1/Table_S3_aggreg_rawRichness.tsv", quote = F, row.names = F, sep="\t", dec = ",")
write.table(df_abun_agg, "Figure1/Table_S3_aggreg_rawAbundance.tsv", quote = F, row.names = F, sep="\t", dec = ",")

## Export df for Figure 2 and Figure S4
saveRDS(df[c("Euphotic", "Aphotic", "Benthic", "Sinking_Plankton"),], "DF_biomes_abundance.rds")
saveRDS(df_rich[c("Euphotic", "Aphotic", "Benthic", "Sinking_Plankton"),], "DF_biomes_richness.rds")

###############################################################################################################
### Figure S1: normalized richness and abundance -- upset plot
### rarefy 1000 times the 'common_per_biome' vectors and get the amount of reads and ASVs shared... 
ups_asvs <- array(NA, c(1000, 7))
ups_reads <- array(NA, c(1000, 7))
colnames(ups_asvs) <- colnames(ups_reads) <- c("Euphotic", "Aphotic", "Sediment", "Euphotic_Aphotic",
                                               "Euphotic_Sediment", "Aphotic_Sediment", "Euphotic_Aphotic_Sediment")
## % of biome reads within intersection
ups_perc <- array(NA, c(1000, 9))
colnames(ups_perc) <- c("EupAph_EUPH", "EupAph_APH",
                        "EupSed_EUPH", "EupSed_SED",
                        "AphSed_APH", "AphSed_SED",
                        "EuphAphSed_EUPH", "EuphAphSed_APH", "EuphAphSed_SED")

for (i in 1:1000) {
  tmp_eup <- rrarefy(df["Euphotic",], 1000000)
  tmp_aph <- rrarefy(df["Aphotic",], 1000000)
  tmp_sed <- rrarefy(df["Sediment",], 1000000)
  tmp_common <- rbind(tmp_eup, tmp_aph, tmp_sed)
  rownames(tmp_common) <- c("Euphotic", "Aphotic", "Sediment")
  tmp_common <- as.data.frame(t(tmp_common))
  ### ASVs
  ups_Euphotic = length(subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic == 0 & tmp_common$Sediment == 0)[,"Euphotic"])
  ups_Aphotic = length(subset(tmp_common, tmp_common$Euphotic == 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment == 0)[,"Aphotic"])
  ups_Sediment = length(subset(tmp_common, tmp_common$Euphotic == 0 & tmp_common$Aphotic == 0 & tmp_common$Sediment != 0)[,"Sediment"])
  ups_Euphotic_Aphotic = nrow(subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment == 0)[,c("Euphotic", "Aphotic")])
  ups_Euphotic_Sediment = nrow(subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic == 0 & tmp_common$Sediment != 0)[,c("Euphotic", "Sediment")])
  ups_Aphotic_Sediment = nrow(subset(tmp_common, tmp_common$Euphotic == 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment != 0)[,c("Aphotic", "Sediment")])
  ups_Euphotic_Aphotic_Sediment = nrow(subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment != 0)[,c("Euphotic", "Aphotic", "Sediment")])
  ups_asvs[i,] <- c(ups_Euphotic, ups_Aphotic, ups_Sediment, ups_Euphotic_Aphotic, 
                    ups_Euphotic_Sediment, ups_Aphotic_Sediment, ups_Euphotic_Aphotic_Sediment)
  ### Reads
  ups_Euphotic = sum(subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic == 0 & tmp_common$Sediment == 0)[,"Euphotic"])
  ups_Aphotic = sum(subset(tmp_common, tmp_common$Euphotic == 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment == 0)[,"Aphotic"])
  ups_Sediment = sum(subset(tmp_common, tmp_common$Euphotic == 0 & tmp_common$Aphotic == 0 & tmp_common$Sediment != 0)[,"Sediment"])
  ## both sum and perc 
  ups_Euphotic_Aphotic = subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment == 0)[,c("Euphotic", "Aphotic")]
  ups_Euphotic_Sediment = subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic == 0 & tmp_common$Sediment != 0)[,c("Euphotic", "Sediment")]
  ups_Aphotic_Sediment = subset(tmp_common, tmp_common$Euphotic == 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment != 0)[,c("Aphotic", "Sediment")]
  ups_Euphotic_Aphotic_Sediment = subset(tmp_common, tmp_common$Euphotic != 0 & tmp_common$Aphotic != 0 & tmp_common$Sediment != 0)[,c("Euphotic", "Aphotic", "Sediment")]
  ups_reads[i,] <- c(ups_Euphotic, ups_Aphotic, ups_Sediment, sum(ups_Euphotic_Aphotic), 
                     sum(ups_Euphotic_Sediment), sum(ups_Aphotic_Sediment), sum(ups_Euphotic_Aphotic_Sediment))
  ## now perc 
  EupAph_EUPH <- sum(ups_Euphotic_Aphotic[,"Euphotic"])
  EupAph_APH <- sum(ups_Euphotic_Aphotic[,"Aphotic"])
  EupSed_EUPH <- sum(ups_Euphotic_Sediment[,"Euphotic"])
  EupSed_SED <- sum(ups_Euphotic_Sediment[,"Sediment"])
  AphSed_APH <- sum(ups_Aphotic_Sediment[,"Aphotic"])
  AphSed_SED <- sum(ups_Aphotic_Sediment[,"Sediment"])
  EuphAphSed_EUPH <- sum(ups_Euphotic_Aphotic_Sediment[,"Euphotic"])
  EuphAphSed_APH <- sum(ups_Euphotic_Aphotic_Sediment[,"Aphotic"])
  EuphAphSed_SED <- sum(ups_Euphotic_Aphotic_Sediment[,"Sediment"])
  
  ups_perc[i,] <- c(EupAph_EUPH, EupAph_APH, EupSed_EUPH, EupSed_SED, 
                    AphSed_APH, AphSed_SED, EuphAphSed_EUPH, EuphAphSed_APH, EuphAphSed_SED)  
  
  if (i %% 100 == 0) message(i)
}

saveRDS(ups_reads, "Figure1/NORMups_reads_1million.rds")
saveRDS(ups_asvs, "Figure1/NORMups_asvs_1million.rds")
saveRDS(ups_perc, "Figure1/NORMups_perc_1million.rds")

## normalize
ups_reads_norm <- ups_reads / rowSums(ups_reads) * 100
ups_asvs_norm <- ups_asvs / rowSums(ups_asvs) * 100

## colmeans and sd
ups_summary <- array(NA, c(4,7))
colnames(ups_summary) <- colnames(ups_reads_norm)
rownames(ups_summary) <- c("ASVs average", "ASVs sd", "Reads average", "Reads sd")
ups_summary["ASVs average",] <- colMeans(ups_asvs_norm)
ups_summary["ASVs sd",] <- colSds(ups_asvs_norm)
ups_summary["Reads average",] <- colMeans(ups_reads_norm)
ups_summary["Reads sd",] <- colSds(ups_reads_norm)
# reverse order to match upset plot
ups_summary <- ups_summary[,c(7:1)]

# split
ups_means <- ups_summary[grep(" average", rownames(ups_summary)),]
ups_sd <- ups_summary[grep(" sd", rownames(ups_summary)),]

## make the df for ggplot
df <- data.frame(cbind(Inter = colnames(ups_summary), t(ups_summary)))
colnames(df) <- c("Inter", "mean %ASVs", "ASVs_sd", "mean %Reads", "Reads_sd")
for (i in 2:5) df[,i] <- as.numeric(as.character(df[,i]))
df <- reshape2::melt(df, measure.vars = colnames(df)[c(2,4)])
df$sd <- c(df$ASVs_sd[1:7], c(df$Reads_sd[8:14]))

### barplot intersections
pdf("Figure1/Figure_S1_top_panel.pdf", width=6, height=3)
ggplot(df, aes(x = Inter, y = value, fill = variable)) +
  geom_bar(stat="identity", position = position_dodge()) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("grey45","grey90"))  + 
  scale_x_discrete(limits=unique(df$Inter)) + 
  scale_y_continuous(breaks = seq(10, 70, by = 10)) +
  #geom_text(aes(label = floor(value), y= floor(value), vjust = -0.5), size = 3) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size = 0.2, linetype='dashed'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
dev.off()

### barplot sets for the lower left part of the upset plot
sets_summary <- array(NA, c(2,3))
colnames(sets_summary) <- c("Euphotic", "Aphotic", "Sediment")
rownames(sets_summary) <- c("ASVs average", "ASVs sd")
sets_summary["ASVs average","Euphotic"] <- sum(colMeans(ups_asvs[,c("Euphotic","Euphotic_Aphotic", "Euphotic_Sediment", "Euphotic_Aphotic_Sediment")]))
sets_summary["ASVs sd","Euphotic"] <- mean(colSds(ups_asvs[,c("Euphotic","Euphotic_Aphotic", "Euphotic_Sediment", "Euphotic_Aphotic_Sediment")]))
sets_summary["ASVs average","Aphotic"] <- sum(colMeans(ups_asvs[,c("Aphotic","Euphotic_Aphotic", "Aphotic_Sediment", "Euphotic_Aphotic_Sediment")]))
sets_summary["ASVs sd","Aphotic"] <- mean(colSds(ups_asvs[,c("Aphotic","Euphotic_Aphotic", "Aphotic_Sediment", "Euphotic_Aphotic_Sediment")]))
sets_summary["ASVs average","Sediment"] <- sum(colMeans(ups_asvs[,c("Sediment","Aphotic_Sediment", "Euphotic_Sediment", "Euphotic_Aphotic_Sediment")]))
sets_summary["ASVs sd","Sediment"] <- mean(colSds(ups_asvs[,c("Sediment","Aphotic_Sediment", "Euphotic_Sediment", "Euphotic_Aphotic_Sediment")]))


pdf("Figure1/Figure_S1_left_panel.pdf", width=1.5, height=1.5)
df <- data.frame(cbind(Set = colnames(sets_summary), t(sets_summary)))
for (i in 2:3) df[,i] <- as.numeric(as.character(df[,i]))
colors_biomes <- c("blue4", "cadetblue1","darkgoldenrod1")
ggplot(df, aes(x = Set, y = ASVs.average, fill = Set)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values = colors_biomes)  + 
  scale_y_reverse() +
  geom_errorbar(aes(ymin=ASVs.average-ASVs.sd, ymax=ASVs.average+ASVs.sd), width=.2,
                position=position_dodge(.9)) +
  scale_x_discrete(limits=rev(unique(df$Set))) + 
  geom_text(aes(label = floor(ASVs.average), y= floor(ASVs.average), vjust = -0.5), size = 3) +
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(colour = "black", size = 0.2),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none") + coord_flip()
dev.off()

#### Finally the percentage of reads from each biome at the intersections
# Across three biomes Euph / Aph / Sed 
eup_aph_sed <- ups_perc[,c("EuphAphSed_SED", "EuphAphSed_APH", "EuphAphSed_EUPH")]
eup_aph_sed <- eup_aph_sed / rowSums(eup_aph_sed) * 100
eup_aph_sed_mean <- colMeans(eup_aph_sed)
eup_aph_sed_sd <- colSds(eup_aph_sed)
names(eup_aph_sed_mean) <- names(eup_aph_sed_sd) <- c("Sediment", "Aphotic", "Euphotic")
pdf("Figure1/Figure_S1_piechart_Euph_Aph_Sed.pdf", width=3, height=3)
pie(eup_aph_sed_mean, labels = "", lty = 0, col = c(alpha("darkgoldenrod1", 1), alpha("blue4", 1), alpha("cadetblue1", 1)), 
    init.angle = 180)
dev.off()

# Euph / Aph
eup_aph <- ups_perc[,c("EupAph_APH", "EupAph_EUPH")]
eup_aph <- eup_aph / rowSums(eup_aph) * 100
eup_aph_mean <- colMeans(eup_aph)
eup_aph_sd <- colSds(eup_aph)
names(eup_aph_mean) <- names(eup_aph_sd) <- c("Aphotic", "Euphotic")
pdf("Figure1/Figure_S1_piechart_Euph_Aph.pdf", width=3, height=3)
pie(eup_aph_mean, labels = "", lty = 0, col = c(alpha("blue4", 1), alpha("cadetblue1", 1)), 
    init.angle = 180)
dev.off()

# Aph / Sed
sed_aph <- ups_perc[,c("AphSed_SED", "AphSed_APH")]
sed_aph <- sed_aph / rowSums(sed_aph) * 100
sed_aph_mean <- colMeans(sed_aph)
sed_aph_sd <- colSds(sed_aph)
names(sed_aph_mean) <- names(sed_aph_sd) <- c("Sediment", "Aphotic")
pdf("Figure1/Figure_S1_piechart_Aph_Sed.pdf", width=3, height=3)
pie(sed_aph_mean, labels = "", lty = 0, col = c(alpha("darkgoldenrod1", 1), alpha("blue4", 1)), 
    init.angle = 180)
dev.off()

# Euph / Sed 
sed_euph <- ups_perc[,c("EupSed_SED", "EupSed_EUPH")]
sed_euph <- sed_euph / rowSums(sed_euph) * 100
sed_euph_mean <- colMeans(sed_euph)
sed_euph_sd <- colSds(sed_euph)
names(sed_euph_mean) <- names(sed_euph_sd) <- c("Sediment", "Euphotic")
pdf("Figure1/Figure_S1_piechart_Euph_Sed.pdf", width=3, height=3)
pie(sed_euph_mean, labels = "", lty = 0, col = c(alpha("darkgoldenrod1", 1), alpha("cadetblue1", 1)), 
    init.angle = 180)
dev.off()

### And the upset plot --> raw number of reads of the last rarefying draw but we actually just want the main plot and add the normalized numbers...
input <- c(
  Euphotic = ups_Euphotic,
  Aphotic = ups_Aphotic,
  Sediment = ups_Sediment,
  "Euphotic&Aphotic" = sum(ups_Euphotic_Aphotic),
  "Euphotic&Sediment" = sum(ups_Euphotic_Sediment),
  "Aphotic&Sediment" = sum(ups_Aphotic_Sediment),
  "Euphotic&Aphotic&Sediment" = sum(ups_Euphotic_Aphotic_Sediment)
)
ups_plot <- fromExpression(input)

pdf("Figure1/Figure_S1_upset_plot.pdf", width=6, height=4, useDingbats = F)
upset(ups_plot, 
      nintersects = 7,
      nsets = 3, 
      intersections = rev(list(list("Euphotic","Aphotic","Sediment"), list("Euphotic","Aphotic"), list("Euphotic","Sediment"),
                               list("Aphotic","Sediment"),list("Euphotic"), list("Aphotic"), list("Sediment"))),
      order.by = "degree", keep.order = T, empty.intersections = "on",
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1)
dev.off()


###############################################################################################################
###############################################################################################################
###############################################################################################################
### Figure 1C: Accumulation curves sampling effort (planktonic pico nano size fractions) 

# only small plankton fractions to compare with the sediment (micro and meso fractions are overwhelmingly dominated by few copepods)
tmp <- comp[rownames(euphotic_all),]
tmp_euph <- subset(tmp, tmp$size_class == "nano" | tmp$size_class == "pico")
#tmp_euph <- subset(tmp, tmp$size_class != "meso" & tmp$size_class != "micro")
euphotic_nano_pico <- euphotic_all[rownames(tmp_euph),]
# aphotic
tmp <- comp[rownames(aphotic_all),]
tmp_aph <- subset(tmp, tmp$size_class == "nano" | tmp$size_class == "pico")
#tmp_aph <- subset(tmp, tmp$size_class != "meso" & tmp$size_class != "micro")
aphotic_nano_pico <- aphotic_all[rownames(tmp_aph),]

## specaccum 
sp_accum_benth_raref    <- specaccum(sediment_all[,asvs_sed_only], method = "random")  ### benthic
sp_accum_plan_sed_raref <- specaccum(sediment_all[,asvs_plan_sed], method = "random")  ### plankton in sediment
sp_accum_euph_raref     <- specaccum(euphotic_nano_pico, method = "random")   ### nano pico only
sp_accum_aph_raref      <- specaccum(aphotic_nano_pico, method = "random")    ### nano pico only

## Plot accumulation curves
pdf("Figure1/Figure_1C_ASVs_accumulation_curves_sampling_effort.pdf", width=3.5, height=4)
plot(sp_accum_euph_raref, col = "cadetblue1", ylab = "#ASVs (10³)", xlab = "#Samples", ci.type="poly", lwd=2, ylim= c(1,length(asvs_euph)+1000), ci.lty=0, ci.col=alpha("cadetblue1", alpha = 0.2), cex.axis=0.6)
lines(sp_accum_benth_raref, col = "darkgoldenrod1", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("darkgoldenrod1", alpha = 0.2))
lines(sp_accum_aph_raref, col = "blue4", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("blue4", alpha = 0.2))
lines(sp_accum_plan_sed_raref, col = "cornflowerblue", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("cornflowerblue", alpha = 0.2))
dev.off()


### Figure 1C inset: Inset Shannon
colors_biomes <- c("cadetblue1", "blue4","darkgoldenrod1")
shanon_euphotic_nano_pico <- vegan::diversity(euphotic_nano_pico[,asvs_euph], index = "shannon")
shanon_aphotic_nano_pico  <- vegan::diversity(aphotic_nano_pico[,asvs_aph], index = "shannon")
shanon_benthic            <- vegan::diversity(sediment_all[,asvs_sed_only], index = "shannon")

df <- data.frame(Shannon = c(shanon_euphotic_nano_pico,shanon_aphotic_nano_pico,shanon_benthic),
                 Biome = c(rep("Euphotic", length(shanon_euphotic_nano_pico)), rep("Aphotic", length(shanon_aphotic_nano_pico)), rep("Benthic", length(shanon_benthic)))) 

my_comparisons <- list( c("Euphotic", "Aphotic"), c("Euphotic", "Benthic"), c("Aphotic", "Benthic") )

p <- ggplot(df, aes(x=Biome, y=Shannon, fill = Biome)) + 
  geom_violin(trim=T, show.legend = FALSE) + xlim(as.character(unique(df$Biome))) +
  scale_fill_manual(breaks = c("Euphotic", "Aphotic", "Benthic"), values = col2hex(unique(colors_biomes))) +
  stat_summary(fun.data=data_summary, geom="pointrange", color="red") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(6.2, 7, 6.60), label = "p.signif") +
  labs(x="", y = "Shannon") + theme(legend.position = "none")


pdf("Figure1/Figure_1C_inset_Shannon_biomes.pdf", width = 1.8, height = 1.2, useDingbats = F)
print(p)
dev.off()

## Figure S2 left panel: effect of plankton size fractions on alpha-diversity
euphotic_all_meso <- subset(euphotic_all, comp[rownames(euphotic_all),"size_class"] == "meso")
euphotic_all_micro <- subset(euphotic_all, comp[rownames(euphotic_all),"size_class"] == "micro")
euphotic_all_nano <- subset(euphotic_all, comp[rownames(euphotic_all),"size_class"] == "nano")
euphotic_all_pico <- subset(euphotic_all, comp[rownames(euphotic_all),"size_class"] == "pico")
## speaccum
sp_accum_euph_meso <- specaccum(euphotic_all_meso, method = "random")
sp_accum_euph_micro <- specaccum(euphotic_all_micro, method = "random")
sp_accum_euph_nano <- specaccum(euphotic_all_nano, method = "random")
sp_accum_euph_pico <- specaccum(euphotic_all_pico, method = "random")

## no meso fraction in the aphotic and not enough micro samples (3...)
aphotic_all_nano <- subset(aphotic_all, comp[rownames(aphotic_all),"size_class"] == "nano")
aphotic_all_pico <- subset(aphotic_all, comp[rownames(aphotic_all),"size_class"] == "pico")
## speaccum
sp_accum_aph_nano <- specaccum(aphotic_all_nano, method = "random")
sp_accum_aph_pico <- specaccum(aphotic_all_pico, method = "random")

## Plot 
pdf("Figure1/Figure_S2_left_panel_accumulation_curves.pdf", width=3.5, height=4)
plot(sp_accum_euph_meso, col = "cadetblue4", ylab = "#ASVs (10³)", xlab = "#Samples", ci.type="poly", lwd=2, xlim= c(1,400), ylim= c(1,60000), ci.lty=0, ci.col=alpha("cadetblue4", alpha = 0.2))
lines(sp_accum_euph_micro, col = "cadetblue3", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("cadetblue3", alpha = 0.2))
lines(sp_accum_euph_nano, col = "cadetblue2", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("cadetblue2", alpha = 0.2))
lines(sp_accum_euph_pico, col = "cadetblue1", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("cadetblue1", alpha = 0.2))
lines(sp_accum_aph_nano, col = "blue3", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("blue3", alpha = 0.2))
lines(sp_accum_aph_pico, col = "blue", ci.type="poly", lwd=2, ci.lty=0, ci.col=alpha("blue", alpha = 0.2))
dev.off()




###############################################################################################################
###############################################################################################################
###############################################################################################################
### NMDS on Euphotic, Aphotic, Sediment --> All samples (all plankton size fractions)
# removing samples with too few reads
seq_depth_cutoff <- 1000
# and rares ASVs
abund_cutoff <- 100

##### For all samples
# filter the otu table for rare ASVs
asv_beta <- subset(asv, rowSums(asv) >= seq_depth_cutoff)
comp_beta <- comp[rownames(asv_beta),]
## remove rare ASVs now
asv_beta  <- asv_beta[,colSums(asv_beta) > abund_cutoff]

# CSS normalization and Bray-Curtis
asv_beta_css_norm <- cssNorm(asv_beta, samples_as_rows = TRUE) ## <-- see companion function
# vegdist **long**  
bray_all <- vegdist(asv_beta_css_norm, method="bray")  
saveRDS(bray_all, "Figure1/bray_all.rds")

### NMDS
nmds_three_biomes <-metaMDS(bray_all, distance = "bray", k = 2, binary=F, trymax = 50, autotransform =F,    
                            wascores = T, expand = T, trace = 1, plot = F, old.wa = FALSE, display = "sites") 


pdf("Figure1/nmds_3biomes.pdf", width = 10, height=6, useDingbats = F)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(nmds_three_biomes$points, col=colors_biomes[varNumEncode(comp_beta[,"biome"])], pch=1, cex = 1, main="NMDS - bray curtis")
ordisurf(nmds_three_biomes, comp_beta[,"depth"], add = T, col="black", labcex = 0.6, lwd.cl = 0.2, family = "gaussian", levels = c(200,1000,2000,3000,4000,5000))
ordisurf(nmds_three_biomes, abs(comp_beta[,"latitude"]), add = T, col="red", labcex = 0.6, lwd.cl = 0.4, family = "gaussian", levels = c(0,30,60))
leg <- unique(as.character(comp_beta[,"biome"]))
legend("topright", inset=c(-0.5,0), leg, col=colors_biomes[unique(varNumEncode(comp_beta[,"biome"]))], pch=16, box.lty=0)
legend("topleft", legend = paste("stress =", round(nmds_three_biomes$stress, 3)), box.lty=0, cex=.8)
dev.off()

## D2217 is an aphotic outlier ##

##### Figure 1D: NMDS with only pico nano plankton fractions
## Making the pelagic nano pico sediment dataset by removing the micro / meso fractions from pelagic samples
remove_micro_meso <- rownames(subset(comp_beta, comp_beta$size_class == "micro" |  comp_beta$size_class == "meso"))
asv_3biomes_pico_nano  <- asv_beta[!rownames(asv_beta) %in% remove_micro_meso,]
comp_3biomes_pico_nano <- comp[rownames(asv_3biomes_pico_nano),]
## we remove the samples for which no GPS and no depth is available
keep <- rownames(subset(comp_3biomes_pico_nano, !is.na(comp_3biomes_pico_nano$latitude) & !is.na(comp_3biomes_pico_nano$depth) & comp_3biomes_pico_nano$sampleID != "D2217"))
asv_3biomes_pico_nano  <- asv_3biomes_pico_nano[keep,]
comp_3biomes_pico_nano <- comp_3biomes_pico_nano[keep,]

# CSS normalization and Bray-Curtis
asv_beta_css_norm <- cssNorm(asv_3biomes_pico_nano, samples_as_rows = TRUE) ## <-- see companion function
# vegdist **long** 
bray_pico_nano <- vegdist(asv_beta_css_norm, method="bray")  
saveRDS(bray_pico_nano, "Figure1/bray_pico_nano.rds")

### NMDS
set.seed(1)
nmds_three_biomes_pico_nano <- metaMDS(bray_pico_nano, distance = "bray", k = 2, binary=F, trymax = 50, autotransform =F,    
                                       wascores = T, expand = T, trace = 1, plot = F, old.wa = F, display = "sites") 

colors_biomes <- c("darkgoldenrod1","cadetblue1", "blue4")
pdf("Figure1/Figure_1D_nmds_3biomes_pico_nano.pdf", width = 10, height=6, useDingbats = F)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(nmds_three_biomes_pico_nano$points, col=colors_biomes[varNumEncode(comp_3biomes_pico_nano[,"biome"])], pch=1, cex = 1, main="NMDS - bray curtis")
ordisurf(nmds_three_biomes_pico_nano, comp_3biomes_pico_nano[,"depth"], add = T, col="black", labcex = 0.6, lwd.cl = 0.2, family = "gaussian", levels = c(200,1000,2000,3000,4000,5000))
ordisurf(nmds_three_biomes_pico_nano, abs(comp_3biomes_pico_nano[,"latitude"]), add = T, col="red", labcex = 0.6, lwd.cl = 0.4, family = "gaussian", levels = c(0,30,60))
leg <- unique(as.character(comp_3biomes_pico_nano[,"biome"]))
legend("topright", inset=c(-0.5,0), leg, col=colors_biomes[unique(varNumEncode(comp_3biomes_pico_nano[,"biome"]))], pch=16, box.lty=0)
legend("topleft", legend = paste("stress =", round(nmds_three_biomes_pico_nano$stress, 3)), box.lty=0, cex=.8)
dev.off()

### Figure 1D inset: Beta dipersion
bray_pico_nano <- readRDS("Figure1/bray_pico_nano.rds")
comp_3biomes_pico_nano <- comp[rownames(as.matrix(bray_pico_nano)),]

beta <- betadisper(bray_pico_nano, group = comp_3biomes_pico_nano$biome, type = "median")

df <- as.data.frame(cbind(Beta = beta$distances, Biome = comp_3biomes_pico_nano$biome))
df$Beta <- as.numeric(as.character(df$Beta))

my_comparisons <- list( c("Euphotic", "Aphotic"), c("Euphotic", "Sediment"), c("Aphotic", "Sediment") )
p <- ggplot(df, aes(x=Biome, y=Beta, fill = Biome)) + 
  geom_violin(trim=T) + 
  xlim("Euphotic","Aphotic", "Sediment") + 
  scale_fill_manual(values = c("blue4", "cadetblue1", "darkgoldenrod1")) + 
  stat_summary(fun.data=data_summary, geom="pointrange", color="red") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.78, 0.82, 0.80), label = "p.signif") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())


pdf(width = 2, height = 1.3, file = "Figure1/Figure_1D_insetBetaDispersion.pdf", useDingbats = F)
print(p)
dev.off()

## adonis - biome
adonis(bray_pico_nano ~ comp_3biomes_pico_nano$biome, permutations = 999, parallel = 8)
## adonis - absolute latitude 
adonis(bray_pico_nano ~ abs(comp_3biomes_pico_nano$latitude), strata = comp_3biomes_pico_nano$biome, permutations = 999, parallel = 8)


### Figure S3: NMDS basins
pdf("Figure1/Figure_S3_nmds_basins_pico_nano.pdf", width = 10, height=6, useDingbats = F)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(nmds_three_biomes_pico_nano$points, col=varNumEncode(comp_3biomes_pico_nano[,"basin"]), pch=varNumEncode(comp_3biomes_pico_nano[,"biome"]), cex = 1, main="NMDS - bray curtis")
ordisurf(nmds_three_biomes_pico_nano, comp_3biomes_pico_nano[,"depth"], add = T, col="black", labcex = 0.6, lwd.cl = 0.2, family = "gaussian", levels = c(200,1000,2000,3000,4000,5000))
ordisurf(nmds_three_biomes_pico_nano, abs(comp_3biomes_pico_nano[,"latitude"]), add = T, col="red", labcex = 0.6, lwd.cl = 0.4, family = "gaussian", levels = c(0,30,60))
leg <- unique(as.character(comp_3biomes_pico_nano[,"biome"]))
leg2 <- unique(as.character(comp_3biomes_pico_nano[,"basin"]))
legend("topright", inset=c(-0.5,0), leg, col="black", pch=unique(varNumEncode(comp_3biomes_pico_nano[,"biome"])), box.lty=0)
legend("bottomright", inset=c(-0.5,0), leg2, col=unique(varNumEncode(comp_3biomes_pico_nano[,"basin"])), pch=16, box.lty=0)
legend("topleft", legend = paste("stress =", round(nmds_three_biomes_pico_nano$stress, 3)), box.lty=0, cex=.8)
dev.off()


## Figure S2 right panel: effect of plankton size fractions on beta-diversity
asv_pelagic <- subset(asv_beta, comp_beta$biome != "Sediment") 
comp_pelagic <- comp[rownames(asv_pelagic),]
## remove empty ASVs now
asv_pelagic  <- asv_pelagic[,colSums(asv_pelagic) > 0]

# CSS normalization and Bray-Curtis
asv_pelagic_css_norm <- cssNorm(asv_pelagic, samples_as_rows = TRUE) ## <-- see companion function
# vegdist **long**  
bray_pelagic <- vegdist(asv_pelagic_css_norm, method="bray")  
saveRDS(bray_pelagic, "bray_pelagic.rds")

### NMDS
set.seed(1)
nmds_pelagic <- metaMDS(bray_pelagic, distance = "bray", k = 2, binary=F, trymax = 50, autotransform =F,    
                        wascores = T, expand = T, trace = 1, plot = F, old.wa = F, display = "sites") 
### Plot
pdf("Figure1/Figure_S2_nmds_pelagic_size_fractions.pdf", width = 10, height=6, useDingbats = F)
par(mar=c(4, 4, 4, 24), xpd=TRUE)
plot(nmds_pelagic$points, col=varNumEncode(comp_pelagic[,"size_class"]), pch=varNumEncode(comp_pelagic[,"biome"]), cex = 1, lwd = 1.5, main="NMDS - bray curtis")
leg <- unique(as.character(comp_pelagic[,"biome"]))
leg2 <- unique(as.character(comp_pelagic[,"size_class"]))
legend("topright", inset=c(-0.5,0), leg, col="black", pch=unique(varNumEncode(comp_pelagic[,"biome"])), box.lty=0)
legend("bottomright", inset=c(-0.5,0), leg2, col=unique(varNumEncode(comp_pelagic[,"size_class"])), pch=16, box.lty=0)
legend("topleft", legend = paste("stress =", round(nmds_pelagic$stress, 3)), box.lty=0, cex=.8)
dev.off()


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################






