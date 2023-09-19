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
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Family_Analysis/")

#Read in the raw family counts
familyCounts <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Family_Analysis/AJB6_Microbiome_FamilyLevel.csv", header = TRUE, sep = ",")

#Remove any unknown families
familyCounts <- familyCounts[familyCounts$taxon != "Unknown_family",]

#Set taxa names
names <- familyCounts$taxon
sums <- rowSums(familyCounts[,2:ncol(familyCounts)])
relabund <- familyCounts[,2:ncol(familyCounts)]/sums
rownames(relabund) <- names
relabund <- na.omit(relabund)

familyCounts.long <- familyCounts %>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count")
write.csv(familyCounts.long, file = "familyCounts.long.csv")

#Re-read in the familyCounts.long.csv with the metadata now included
long <- read.csv("familyCounts.long.csv", header = TRUE, sep = ",")

#Read in the significant family data
sigFam <- read.csv("significant_families.csv", header = TRUE, sep = ",")

#Get the top 10 most significant families
sigFam <- sigFam[order(sigFam$p.adj, decreasing = FALSE),]
sigFamNames <- sigFam$taxa
top10sig <- sigFam[1:10,]
top10names <- top10sig$taxa

#Subset for the top10sig taxa in the familyCounts.long data
sigCounts <- long[long$taxon %in% top10names,]
significant <- long[long$taxon %in% sigFamNames,]

Age_order <- c("8 Weeks", "20 Weeks")
genotype_order <- c("WT", "KO")
taxa_order <- c("Akkermansiaceae",
                "Atopobiaceae",
                "Cellulomonadaceae",
                "Hungateiclostridiaceae",
                "Bacteroidaceae",
                "Hymenobacteraceae",
                "Clostridiaceae",
                "Odoribacteraceae",
                "Rikenellaceae",
                "Muribaculaceae")
sigCounts <- sigCounts %>%
  mutate(Family = factor(taxon, levels = taxa_order))
sigCounts <- sigCounts %>%
  mutate(Age = factor(Age, levels = Age_order))
sigCounts <- sigCounts %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

#Make a timepoint Plot for the top10 significant families
timepoint <- ggplot(sigCounts, aes(x = Genotype, y = count, fill = Family)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = "Top 10 significant families") +
  theme_bw() +
  theme(legend.position = "right")
timepoint
ggsave("Top10_significant_Families_timepoint.pdf", timepoint, width = 12, height = 8)

Significant_families <- unique(sigFam$taxa)


#Set working directory for output files
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Family_Analysis/Statistically_Significant_Family_BoxPlots/")
#Write a for loop to iterate through the significant families and make a timepoint plot
for(taxa in Significant_families) {
  counts <- significant[significant$taxon %in% taxa,]
  
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

#Read in the metadata
meta <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Data/Metadata.csv", header = TRUE, sep = ",")

#Get the significant families from the raw familyCounts data
famCounts.filt <- familyCounts[familyCounts$taxon %in% Significant_families,]
rownames(famCounts.filt) <- famCounts.filt$taxon
famCounts.filt$taxon <- NULL

#RAW COUNTS
library("pheatmap")
famCounts.filt.mat <- as.matrix(famCounts.filt)
famCounts.filt.mat.t <- as.data.frame(t(famCounts.filt.mat))
famCounts.filt.mat.t$Sample_ID <- rownames(famCounts.filt.mat.t)
famCounts.filt.mat.t.meta <- merge(famCounts.filt.mat.t, meta, by = "Sample_ID", all = TRUE)
rownames(famCounts.filt.mat.t.meta) <- famCounts.filt.mat.t.meta$Sample_ID
famCounts.filt.mat.t.meta$Sample_ID <- NULL
# rownames(speciesCounts.filt.mat.meta) <- speciesCounts.filt.mat.meta$Sample_ID

Strain <- famCounts.filt.mat.t.meta$Strain
names(Strain) <- rownames(famCounts.filt.mat.t.meta)
Strain <- as.data.frame(Strain)
Strain$Genotype <- famCounts.filt.mat.t.meta$Genotype
Strain$Age <- famCounts.filt.mat.t.meta$Age

setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Family_Analysis/")
heatmap <- pheatmap(famCounts.filt.mat.t.meta[,1:27],
                    gaps_row = c(23,47,69), 
                    cluster_rows = FALSE,
                    cluster_cols = TRUE,
                    scale = "column",
                    fontsize = 6,
                    annotation_row = Strain,
                    main = "Significant Families")
ggsave("SigFamilies_Heatmap.pdf", heatmap, width = 12, height = 8)










