library("dplyr")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ggh4x")
library("tidyverse")
library("RColorBrewer")
library("cowplot")

#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Species_analysis/")

#Read in the counts
counts <- read.csv("Species_level_Counts.csv", header = TRUE, sep = ",")

speciesCounts.long <- long %>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count")
write.csv(speciesCounts.long, file = "speciesCounts.long.csv")

#Read in the speciesCounts.long with metadata
spec <- read.csv("speciesCounts.long.csv", header = TRUE, sep = ",")

#Read in the subset
subset <- read.csv("SignificantSpecies_NOT_belonging_to_SignificantFamilies.csv", header = FALSE, sep = ",")
colnames(subset) <- c("taxon")
subset_taxa <- subset$taxon

#Filter the spec df for just the subset
nonSig <- spec[spec$taxon %in% subset_taxa,]

nonSigSpecs <- unique(nonSig$taxon)


#Set working directory for output files
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Species_analysis/SignificantSpecies_in_NonsignificantFamilies/")
#Write a for loop to iterate through the significant families and make a timepoint plot
for(taxa in nonSigSpecs) {
  counts <- nonSig[nonSig$taxon %in% taxa,]
  
  #Set facet factors
  Age_order <- c("8 Weeks", "20 Weeks")
  genotype_order <- c("WT", "KO")
  
  #Reorder facet factors
  counts <- counts %>%
    mutate(Age = factor(Age, levels = Age_order))
  counts <- counts %>%
    mutate(Genotype = factor(Genotype, levels = genotype_order))
  
  #Plot
  timepoint <- ggplot(counts, aes(x = Genotype, y = count, fill = taxa)) +
    geom_boxplot(width = 0.8, outlier.shape = NA, outlier.color = "black") +
    #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
    facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    stat_compare_means(method = "wilcox", label = "p.format") +
    labs(x = "Strain", y = "Normalized Read Counts/Sample") +
    labs(title = taxa) +
    theme_bw() +
    theme(legend.position = "right")
  timepoint
  ggsave(paste(taxa, "timepoint_boxplot.pdf", sep = "_"), timepoint, width = 12, height = 8)
  
}


meta <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Data/Metadata.csv", header = TRUE, sep = ",")

#Get the significant families from the raw familyCounts data
specSub <- counts[counts$taxon %in% subset_taxa,]
rownames(specSub) <- specSub$taxon
specSub$taxon <- NULL

#RAW COUNTS
library("pheatmap")
specSub.mat <- as.matrix(specSub)
specSub.mat.t <- as.data.frame(t(specSub.mat))
specSub.mat.t$Sample_ID <- rownames(specSub.mat.t)
specSub.mat.t.meta <- merge(specSub.mat.t, meta, by = "Sample_ID", all = TRUE)
rownames(specSub.mat.t.meta) <- specSub.mat.t.meta$Sample_ID
specSub.mat.t.meta$Sample_ID <- NULL
# rownames(speciesCounts.filt.mat.meta) <- speciesCounts.filt.mat.meta$Sample_ID

Strain <- specSub.mat.t.meta$Strain
names(Strain) <- rownames(specSub.mat.t.meta)
Strain <- as.data.frame(Strain)
Strain$Genotype <- specSub.mat.t.meta$Genotype
Strain$Age <- specSub.mat.t.meta$Age

setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Family_Analysis/")
heatmap <- pheatmap(specSub.mat.t.meta[,1:27],
                    gaps_row = c(23,47,69), 
                    cluster_rows = FALSE,
                    cluster_cols = TRUE,
                    scale = "column",
                    fontsize = 6,
                    annotation_row = Strain,
                    main = "Significant Species in Non-significant Families")
ggsave("SigSpecies_with_nonsignificantFamilies_Heatmap.pdf", heatmap, width = 12, height = 8)


