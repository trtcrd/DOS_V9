##################################################################################
## Patterns of eukaryotic diversity from the surface to the deep-ocean sediment
## --- 
## Figure 2 S4 S5 S6 S7
##################################################################################


## setting path **to be adjusted*
setwd("~/path/to/be/adjusted/")

path.fig <- paste0(getwd(), "/Figure2")
if(!dir.exists(path.fig)) dir.create(path.fig)

## libraries
library('doMC')
registerDoMC(cores = 16)
library(ggplot2)
library(ggpubr)
library(treemapify)
library(seqinr)

###############################################################################################################
###############################################################################################################
###############################################################################################################
###### Import required data

## Taxo file with functional and size annotations 
taxo <- readRDS("asvs_taxo_functional_annotations.rds")
## Aggregated vectors 
df_abund <- readRDS("DF_biomes_abundance.rds")
df_rich  <- readRDS("DF_biomes_richness.rds")
# same order... 
identical(rownames(taxo), colnames(df_abund)) ## TRUE

## aggregate per taxo groups
df_abund_agg <- aggregate(t(df_abund), by = list(taxo$taxo_group), FUN = sum)
df_rich_agg  <- aggregate(t(df_rich), by = list(taxo$taxo_group), FUN = sum)

## get the taxo path 
taxo_path <- data.frame(rank3 = sapply(strsplit(df_abund_agg$Group.1, ";"), `[`, 1),
                        rank4 = sapply(strsplit(df_abund_agg$Group.1, ";"), `[`, 2),
                        rank5 = sapply(strsplit(df_abund_agg$Group.1, ";"), `[`, 3), stringsAsFactors=FALSE)


## abundance
df_abund_agg <- data.frame(cbind(taxo_path, df_abund_agg[,2:ncol(df_abund_agg)]), stringsAsFactors = F)
df_abund_agg <- reshape2::melt(df_abund_agg, id.vars = c("rank3", "rank4", "rank5"))                         
df_abund_agg$value <- as.numeric(df_abund_agg$value)

## richness
df_rich_agg <- data.frame(cbind(taxo_path, df_rich_agg[,2:ncol(df_rich_agg)]), stringsAsFactors = F)
df_rich_agg <- reshape2::melt(df_rich_agg, id.vars = c("rank3", "rank4", "rank5"))                         
df_rich_agg$value <- as.numeric(df_rich_agg$value)

## focus on groups
focus <- c("Metazoa", "Rhizaria", "Alveolata", "Stramenopiles", "Euglenozoa", "Fungi", "Amoebozoa")
df_abund_agg$rank_col <- ifelse(df_abund_agg$rank3 %in% focus, df_abund_agg$rank3, "Others")
df_rich_agg$rank_col <- ifelse(df_rich_agg$rank3 %in% focus, df_rich_agg$rank3, "Others")
# ordering factors
df_abund_agg$variable <- factor(df_abund_agg$variable, levels = c("Euphotic", "Aphotic", "Sinking_Plankton","Benthic"))
df_rich_agg$variable <- factor(df_rich_agg$variable, levels = c("Euphotic", "Aphotic", "Sinking_Plankton","Benthic"))

#### Figure 2: treemap richness

pdf("Figure2/Figure_2_Treemap_richness_eupthotic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_rich_agg, df_rich_agg$variable == "Euphotic"), 
                 title = "Richness - Euphotic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


pdf("Figure2/Figure_2_Treemap_richness_apthotic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_rich_agg, df_rich_agg$variable == "Aphotic"), 
                 title = "Richness - Aphotic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


pdf("Figure2/Figure_2_Treemap_richness_sinking_plankton.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_rich_agg, df_rich_agg$variable == "Sinking_Plankton"), 
                 title = "Richness - Sinking plankton", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


pdf("Figure2/Figure_2_Treemap_richness_benthic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_rich_agg, df_rich_agg$variable == "Benthic"), 
                 title = "Abundance - Benthic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


#### Figure 2 B and C: 
taxo02 <- readRDS("taxo02_annotations.rds")
asvs_euph <- readRDS("asvs_euph.rds")
asvs_aph <- readRDS("asvs_aph.rds")
asvs_sed_only <- readRDS("asvs_sed_only.rds")

df_euph <- cbind(taxo02[asvs_euph,c("taxon", "mean.similarity")], Biome = "Euphotic", Abundance = df_abund["Euphotic",asvs_euph] *100 / sum(df_abund["Euphotic",asvs_euph]))
df_aph <- cbind(taxo02[asvs_aph,c("taxon", "mean.similarity")], Biome = "Aphotic", Abundance = df_abund["Aphotic",asvs_aph] *100 / sum(df_abund["Aphotic",asvs_aph]))
df_sed <- cbind(taxo02[asvs_sed_only,c("taxon", "mean.similarity")], Biome = "Benthic", Abundance = df_abund["Benthic",asvs_sed_only]  *100 / sum(df_abund["Benthic",asvs_sed_only]))
df_all <- data.frame(rbind(df_euph, df_aph, df_sed))
colnames(df_all) <- c("Taxon", "Similarity", "Biome", "Abundance") ## be carefull of the ASVs overlap between euph and aph, and it renames as unique ID in rownames... 
df_all$Similarity <- df_all$Similarity *100

df_euph <- cbind(taxo02[asvs_euph,c("taxon", "mean.similarity")], Biome = "Euphotic", Abundance = df_abund["Euphotic",asvs_euph] *100 / sum(df_abund["Euphotic",asvs_euph]))
df_aph <- cbind(taxo02[asvs_aph,c("taxon", "mean.similarity")], Biome = "Aphotic", Abundance = df_abund["Aphotic",asvs_aph] *100 / sum(df_abund["Aphotic",asvs_aph]))
df_sed <- cbind(taxo02[asvs_sed_only,c("taxon", "mean.similarity")], Biome = "Benthic", Abundance = df_abund["Benthic",asvs_sed_only]  *100 / sum(df_abund["Benthic",asvs_sed_only]))
df_all <- data.frame(rbind(df_euph, df_aph, df_sed))
colnames(df_all) <- c("Taxon", "Similarity", "Biome", "Abundance") ## be carefull of the ASVs overlap between euph and aph, and it renames as unique ID in rownames... 
df_all$Similarity <- df_all$Similarity *100


# plot distrib
dens_e <- density(df_euph$mean.similarity[!is.na(df_euph$mean.similarity)]*100)
max_euph <- dens_e$x[which(dens_e$y == max(dens_e$y))]
dens_a <- density(df_aph$mean.similarity[!is.na(df_aph$mean.similarity)]*100)
max_aph <- dens_a$x[which(dens_a$y == max(dens_a$y))]
dens_s <- density(df_sed$mean.similarity[!is.na(df_sed$mean.similarity)]*100)
max_sed <- dens_s$x[which(dens_s$y == max(dens_s$y))]

p_density <- ggplot(df_all, aes(x=Similarity, fill=Biome, color = Biome)) +
  geom_density(alpha = 0.2) + 
  scale_x_reverse() + 
  scale_fill_manual(values = c(Benthic = "darkgoldenrod1", Aphotic = "blue4", Euphotic = "cadetblue1")) +
  scale_color_manual(values = c(Benthic = "darkgoldenrod1", Aphotic = "blue4", Euphotic = "cadetblue1")) +
  geom_vline(xintercept=85, color="red",linetype="solid", size=0.4) + 
  annotate(geom = "curve", x = max_euph -6, y = max(dens_e$y) + 0.001, xend = max_euph, yend = max(dens_e$y), curvature = .1, arrow = arrow(length = unit(2, "mm"))) +
  annotate(geom = "text",  x = max_euph -6.5, y = max(dens_e$y) + 0.001, label = paste0(round(max_euph, 1), "%"), hjust = "left") +
  annotate(geom = "curve", x = max_aph -4, y = max(dens_a$y) + 0.001, xend = max_aph, yend = max(dens_a$y), curvature = .1, arrow = arrow(length = unit(2, "mm"))) +
  annotate(geom = "text",  x = max_aph -4.5, y = max(dens_a$y) + 0.001, label = paste0(round(max_aph, 1), "%"), hjust = "left") +
  annotate(geom = "curve", x = max_sed -4, y = max(dens_s$y) + 0.003, xend = max_sed, yend = max(dens_s$y), curvature = .3, arrow = arrow(length = unit(2, "mm"))) +
  annotate(geom = "text",  x = max_sed -4.5, y = max(dens_s$y) + 0.003, label = paste0(round(max_sed, 1), "%"), hjust = "left") +
  geom_vline(xintercept=max_sed, color="darkgoldenrod1", linetype="dotted", size=0.4) +
  geom_vline(xintercept=max_aph, color="blue4", linetype="dotted", size=0.4) +
  geom_vline(xintercept=max_euph, color="cadetblue1", linetype="dotted", size=0.4) +
  ylab("Density") + xlab("% similarity with best hit in PR²") +
  theme_classic() + theme(legend.position = "bottom")


## Cumulative relative abundance as function of similarity
tt_sed <- df_sed[order(df_sed$mean.similarity, decreasing = T),]
tt_sed$Cumulative <- cumsum(tt_sed$Abundance)
tt_aph <- df_aph[order(df_aph$mean.similarity, decreasing = T),]
tt_aph$Cumulative <- cumsum(tt_aph$Abundance)
tt_euph <- df_euph[order(df_euph$mean.similarity, decreasing = T),]
tt_euph$Cumulative <- cumsum(tt_euph$Abundance)
df_cum <- data.frame(rbind(tt_sed, tt_aph, tt_euph), stringsAsFactors = F)
#df_cum$Biome <- paste0("% reads ", df_cum$Biome)
df_cum$Line <- "%reads"

inter_euph_reads <- sum(subset(tt_euph, tt_euph$mean.similarity < 0.85)[,"Abundance"])
inter_aph_reads <- sum(subset(tt_aph, tt_aph$mean.similarity < 0.85)[,"Abundance"])
inter_sed_reads <- sum(subset(tt_sed, tt_sed$mean.similarity < 0.85)[,"Abundance"])

inter_euph_asvs <- nrow(subset(tt_euph, tt_euph$mean.similarity < 0.85)) *100 / nrow(tt_euph) 
inter_aph_asvs <- nrow(subset(tt_aph, tt_aph$mean.similarity < 0.85)) *100 / nrow(tt_aph) 
inter_sed_asvs <- nrow(subset(tt_sed, tt_sed$mean.similarity < 0.85)) *100 / nrow(tt_sed) 

## cumulative #ASVs
tt_sed$Cumulative <- cumsum(rep(1,length(tt_sed$Abundance))*100/length(tt_sed$Abundance))
tt_aph$Cumulative <- cumsum(rep(1,length(tt_aph$Abundance))*100/length(tt_aph$Abundance))
tt_euph$Cumulative <- cumsum(rep(1,length(tt_euph$Abundance))*100/length(tt_euph$Abundance))
df_cum_asvs <- data.frame(rbind(tt_sed, tt_aph, tt_euph), stringsAsFactors = F)
#df_cum_asvs$Biome <- paste0("% ASVs ", df_cum_asvs$Biome)
df_cum_asvs$Line <- "%ASVs"
df_cum <- data.frame(rbind(df_cum, df_cum_asvs), stringsAsFactors = F)
colnames(df_cum) <- c("Taxon", "Similarity", "Biome", "Abundance", "Cumulative", "Line")

## cumulative relative abundances 
p_cum <-ggplot(df_cum) +
  geom_line(alpha = 1,size=0.7, aes(x=Similarity *100, y=Cumulative, color=Biome, linetype = Line)) + 
  scale_color_manual(values = c(Benthic = "darkgoldenrod1", Aphotic = "blue4", Euphotic = "cadetblue1")) +
  geom_vline(xintercept=85, color="red",linetype="solid", size=0.4) + 
  scale_x_reverse() + 
  ylab("Cumulative proportion") + xlab("% similarity with best hit in PR²") +
  theme_classic() +  theme(legend.title=element_blank(), legend.position="bottom")



p_dens_cum <- ggarrange(p_density, p_cum, ncol = 2, labels = c("B","C"))
ggsave("Figure2/Figure_2_density.pdf", width = 8, height = 3)


###############################################################################################################
###############################################################################################################
###############################################################################################################


#### Figure S4: treemap abundances

pdf("Figure2/Figure_S4_Treemap_abundance_eupthotic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_abund_agg, df_abund_agg$variable == "Euphotic"), 
                 title = "Abundance - Euphotic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


pdf("Figure2/Figure_S4_Treemap_abundance_apthotic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_abund_agg, df_abund_agg$variable == "Aphotic"), 
                 title = "Abundance - Aphotic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


pdf("Figure2/Figure_S4_Treemap_abundance_sinking_plankton.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_abund_agg, df_abund_agg$variable == "Sinking_Plankton"), 
                 title = "Abundance - Sinking plankton", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()


pdf("Figure2/Figure_S4_Treemap_abundance_benthic.pdf", width = 8, height = 2) 
treemap::treemap(subset(df_abund_agg, df_abund_agg$variable == "Benthic"), 
                 title = "Abundance - Benthic", 
                 algorithm = "pivotSize", border.lwds = c(1,0.5,0.1), 
                 border.col = c("black", "black", "black"),
                 mapping = c(0,0,0),
                 index=c("rank3","rank4", "rank5"),
                 vSize="value",
                 vColor="rank_col",
                 type="index",
                 fontsize.labels=c(18,14,10),             # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                 fontcolor.labels=c("Grey", "orange","white"),    # Color of labels
                 fontface.labels=c(3,2,1),                # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                 bg.labels=c("transparent"),              # Background color of labels
                 align.labels=list(
                   c("left", "top"), 
                   c("left", "bottom"),
                   c("center", "center")
                 ),                                   # Where to place labels in the rectangle?
                 overlap.labels=0.8,                  # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                 inflate.labels=F,
                 force.print.labels = F) 
dev.off()

###############################################################################################################
###############################################################################################################
###############################################################################################################


### Figure S5

## combine plankton to compare with the benthic taxogroups
df_pel_ben <- data.frame(Pelagic = colSums(df_abund[c("Euphotic", "Aphotic"),]),
                         Benthic = df_abund["Benthic",])
  
## binary 
df_pel_ben[df_pel_ben > 0] <- 1

df_pel_ben <- aggregate(df_pel_ben, by = list(taxo$taxo_group), FUN = sum)

## get the taxo path 
taxo_path <- data.frame(rank3 = sapply(strsplit(df_pel_ben$Group.1, ";"), `[`, 1),
                        rank4 = sapply(strsplit(df_pel_ben$Group.1, ";"), `[`, 2),
                        rank5 = sapply(strsplit(df_pel_ben$Group.1, ";"), `[`, 3), stringsAsFactors=FALSE)

df_pel_ben <- data.frame(cbind(taxo_path, df_pel_ben[,2:ncol(df_pel_ben)]), stringsAsFactors = F)
df_pel_ben <- reshape2::melt(df_pel_ben, id.vars = c("rank3", "rank4", "rank5"))                         
df_pel_ben$value <- as.numeric(df_pel_ben$value)

## remove unassigned
df_pel_ben <- subset(df_pel_ben, df_pel_ben$rank3 != "Unassigned_1") 
df_pel_ben <- subset(df_pel_ben, df_pel_ben$rank5 != "Unassigned_3") 
df_pel_ben$rank3 <- as.character(df_pel_ben$rank3)

# focus on groups
focus <- c("Metazoa", "Rhizaria", "Alveolata", "Stramenopiles", "Euglenozoa", "Fungi", "Amoebozoa")
df_pel_ben$rank_col <- ifelse(df_pel_ben$rank3 %in% focus, df_pel_ben$rank3, "Others")
df_pel_ben$rank_col <- factor(df_pel_ben$rank_col, levels = c("Metazoa", "Rhizaria", "Alveolata", "Stramenopiles", "Euglenozoa", "Fungi", "Amoebozoa", "Others"))

#colors
col <- colorRampPalette(c("dodgerblue", "gold2", "darkred", "darkolivegreen1", 
                          "firebrick1", "darkseagreen1", "bisque3", "maroon3", 
                          "blue", "green", "black", "grey"), bias=1, interpolate = "linear")(8)

p <- ggplot(df_pel_ben, aes(area = value, fill = rank_col, subgroup = rank5))  +
  geom_treemap(layout = "squarified", start = "topleft") +
  geom_treemap_subgroup_border(start = "topleft", colour = "black", size = 1) +
  geom_treemap_subgroup_text(start = "topleft", place = "topleft", reflow = T, alpha = 1, colour =
                             "black", fontface = "bold", min.size = 4, size = 13) +
  scale_fill_manual(values = col) +
  facet_wrap(~ variable, ncol = 1)

pdf("Figure2/Figure_S5_treemap_taxogroup_richness_pelagic_benthic.pdf", width = 8, height = 6)
print(p)
dev.off()



###############################################################################################################
###############################################################################################################
###############################################################################################################

## Figure S6

## #OTUs as function of clustering cutoff --> OTUs made with vsearch
path.clustering <- paste0(getwd(), "/otus_vsearch")
if(!dir.exists(path.clustering)) dir.create(path.clustering)
# eupthotic
unass_euphotic <- taxo[asvs_euph,]
unass_euphotic <- subset(unass_euphotic, unass_euphotic$taxon == "unassigned")
write.fasta(as.list(rownames(unass_euphotic)), unass_euphotic$ASV_ID, paste0(path.clustering, "/unassigned_euphotic.fasta"))
# aphotic
unass_aphotic <- taxo[asvs_aph,]
unass_aphotic <- subset(unass_aphotic, unass_aphotic$taxon == "unassigned")
write.fasta(as.list(rownames(unass_aphotic)), unass_aphotic$ASV_ID, paste0(path.clustering, "/unassigned_aphotic.fasta"))
# benthic
unass_benthic <- taxo[asvs_sed_only,]
unass_benthic <- subset(unass_benthic, unass_benthic$taxon == "unassigned")
write.fasta(as.list(rownames(unass_benthic)), unass_benthic$ASV_ID, paste0(path.clustering, "/unassigned_benthic.fasta"))

## vsearch clustering (vsearch has to be in the PATH, or adjust the path below to the binary)
vsearch <- "vsearch"
system2(vsearch, args = "-h")

cutoff <- rev(c(.85,.86,.87,.88,.89,.90,.91,.92,.93,.94,.95,.96,.97, 1))
otus_stats <- data.frame(array(NA, c(length(cutoff), 4)))
colnames(otus_stats) <- c("cutoff", "benthic", "aphotic", "euphotic")
rownames(otus_stats) <- otus_stats$cutoff <- cutoff

for (i in cutoff) {
  if (i != 1) {
    system2(vsearch, args = c("--id", i, "--cluster_size", paste0(path.clustering, "/unassigned_euphotic.fasta"), "--uc", paste0(path.clustering,"/euph_", i, ".uc"), "--centroids", paste0(path.clustering,"/euph_repset_", i, ".fasta")))
    system2(vsearch, args = c("--id", i, "--cluster_size", paste0(path.clustering, "/unassigned_aphotic.fasta"), "--uc", paste0(path.clustering,"/aph_", i, ".uc"), "--centroids", paste0(path.clustering,"/aph_repset_", i, ".fasta")))
    system2(vsearch, args = c("--id", i, "--cluster_size", paste0(path.clustering, "/unassigned_benthic.fasta"), "--uc", paste0(path.clustering,"/bent_", i, ".uc"), "--centroids", paste0(path.clustering,"/bent_repset_", i, ".fasta")))
    # collect 
    otus_stats[as.character(i),"euphotic"] <- as.numeric(system2("grep", args = c("'>'", paste0(path.clustering, "/euph_repset_",i,".fasta"), "| wc -l"), stdout = TRUE))
    otus_stats[as.character(i),"aphotic"] <- as.numeric(system2("grep", args = c("'>'", paste0(path.clustering, "/aph_repset_",i,".fasta"), "| wc -l"), stdout = TRUE))
    otus_stats[as.character(i),"benthic"] <- as.numeric(system2("grep", args = c("'>'", paste0(path.clustering, "/bent_repset_",i,".fasta"), "| wc -l"), stdout = TRUE))
  } else {
    otus_stats[as.character(i),"euphotic"] <- as.numeric(system2("grep", args = c("'>'", paste0(path.clustering, "/unassigned_euphotic.fasta"), "| wc -l"), stdout = TRUE))
    otus_stats[as.character(i),"aphotic"] <- as.numeric(system2("grep", args = c("'>'", paste0(path.clustering, "/unassigned_aphotic.fasta"), "| wc -l"), stdout = TRUE))
    otus_stats[as.character(i),"benthic"] <- as.numeric(system2("grep", args = c("'>'", paste0(path.clustering, "/unassigned_benthic.fasta"), "| wc -l"), stdout = TRUE))
  }
}

otus_stats <- reshape2::melt(otus_stats, id.vars = 1)
otus_stats$cutoff <- otus_stats$cutoff *100
otus_stats$cutoff <- factor(otus_stats$cutoff, levels = unique(otus_stats$cutoff))


pdf("Figure2/Figure_S6_unassigned_ASVs_clusteredOTUs.pdf", width = 4.5, height=3.5, useDingbats = F)
ggplot(otus_stats, aes(x=cutoff, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())  + 
  scale_y_continuous(name="#OTUs", breaks=c(0,5000, 10000, 20000, 30000, 40000, 50000,60000)) + 
  scale_fill_manual(values=c(benthic="darkgoldenrod1", aphotic="blue4", euphotic="cadetblue1")) + 
  ylab("#OTUs / ASVs") + xlab("OTU clustering threshold") + theme_classic() +
  theme(legend.title=element_blank(), legend.position="bottom") 
dev.off()


###############################################################################################################
###############################################################################################################
###############################################################################################################

### Figure S7: Focus on selected benthic groups
groups_check <- c("polychaetes", "Nemertea", "Nematoda", "Foraminifera", "Ciliophora", "Flabellinia")

plot_list = list()

for (i in groups_check) {
  tmp_taxo <- taxo[asvs_sed_only,]
  tmp_taxo <- tmp_taxo[grep(i, tmp_taxo$taxo_group),]
  tmp_asv  <- df_abund["Benthic", rownames(tmp_taxo)]
  tmp <- cbind(tmp_taxo, Abundance = tmp_asv *100 / sum(tmp_asv))
  tmp_sorted_reads <- tmp[order(tmp$mean.similarity, decreasing = T),]
  tmp_sorted_reads$Cumulative <- cumsum(tmp_sorted_reads$Abundance)
  tmp_sorted_asvs <- tmp[order(tmp$mean.similarity, decreasing = T),]
  tmp_sorted_asvs$Cumulative <- cumsum(rep(1,length(tmp_sorted_asvs$Abundance))*100/length(tmp_sorted_asvs$Abundance))
  tmp_sorted_asvs$Line <- "%ASVs"
  tmp_sorted_reads$Line <- "%reads"
  tmp <- data.frame(rbind(tmp_sorted_asvs, tmp_sorted_reads), stringsAsFactors = F)
  ## plot
  if (i == "polychaetes") {
    pp <- ggplot(tmp) +
      geom_line(alpha = 1,size=0.7, aes(x=mean.similarity, y=Cumulative, linetype = Line)) +
      xlim(0.85,1) + labs(title=paste0("Benthic ", i), x="% similarity with best hit in PR²", y = "Cumulative proportion") +
      scale_x_reverse() + 
      theme(legend.title=element_blank(), legend.position=c(0.15, 0.9), legend.background = element_rect(fill = "transparent")) 
  } else {
    pp <- ggplot(tmp) +
      geom_line(alpha = 1,size=0.7, aes(x=mean.similarity, y=Cumulative, linetype = Line)) + 
      xlim(0.85,1) + labs(title=paste0("Benthic ", i), x="% similarity with best hit in PR²", y = "Cumulative proportion") +
      scale_x_reverse() + 
      theme(legend.title=element_blank(), legend.position="none")
  }
  plot_list[[i]] <- pp
}

p <- ggarrange(plotlist = plot_list, nrow = 2, ncol = 3)
ggsave("Figure2/Figure_S7_benthic_groups_distrib_bestMatch.pdf", p, width = 12, height=9)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################



