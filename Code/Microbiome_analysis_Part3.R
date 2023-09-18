#Futher analysis of the top significant species found in Microbiome_Analysis_Part2.R


#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/")

#Load libraries
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ggh4x")
library("tidyverse")
library("RColorBrewer")
library("broom")
library("tidyverse")
library("dplyr")

#Read in the top 5 significant species counts data
top5 <- read.csv("Top5_sigSpecies.csv", header = TRUE, sep = ",")

#Pivot longer the csv file
top5.long <- top5 %>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 
top5.long$Level <- "Species"
write.csv(top5.long, file = "Top5_sigSpecies_long.csv")

#Manually add in the metadata and re-read in the long format csv file
long.top5 <- read.csv("Top5_sigSpecies_long.csv", header = TRUE, sep = ",")

#Reorder factors
species_order <- c("Prevotella_marshii",
                 "Akkermansia_muciniphila",
                 "Faecalibaculum_unclassified",
                 "Turicibacter_unclassified",
                 "Akkermansia_unclassified")
                 

#Re-order
long.top5 <- long.top5 %>%
  mutate(taxon = factor(taxon, levels = species_order))

#Convert Age to a factor so 8 weeks comes before 20
Age_order <- c("8 Weeks",
               "20 Weeks")

#Re-order
long.top5 <- long.top5 %>%
  mutate(Age = factor(Age, levels = Age_order))

#Convert Phenotype order to factor and set desired order
Genotype_order <- c("WT", "KO")

#Re-order phenotype order
long.top5 <- long.top5 %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

#Plotting all top 5 together
Box <- ggplot(long.top5, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = 5, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "skyblue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  ylim(0,400) +
  labs(title = expression(italic("Significant hits 1-5"))) +
  theme_bw()
Box
ggsave("Top5_SigSpecies_timepoint.pdf", Box, width = 12, height = 8)

#P.marshii
marshii <- read.csv("P.marshii.csv", header = TRUE, sep = ",")

#Re-order
marshii <- marshii %>%
  mutate(Age = factor(Age, levels = Age_order))

#Re-order phenotype order
marshii <- marshii %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

#Plot
marshii <- ggplot(marshii, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Prevotella marshii"))) +
  theme_bw() +
  theme(legend.position = "none")
marshii
ggsave("Prevotella_marshii_timepoint.pdf", marshii, width = 12, height = 8)


#Read in Akkermansia muciniphila data
akk <- read.csv("Akkermansia_muciniphila.csv", header = TRUE, sep = ",")
#Re-order
akk <- akk %>%
  mutate(Age = factor(Age, levels = Age_order))

#Re-order phenotype order
akk <- akk %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

#Plot
akkPlot <- ggplot(akk, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Akkermansia muciniphila"))) +
  theme_bw() +
  theme(legend.position = "none")
akkPlot
ggsave("Akkermansia_muciniphila_timepoint.pdf", akkPlot, width = 12, height = 8)

#Read in Faecalibaculum unclassified
fae <- read.csv("Faecalibaculum_unclassified.csv", header = TRUE, sep = ",")
#Re-order
fae <- fae %>%
  mutate(Age = factor(Age, levels = Age_order))

#Re-order phenotype order
fae <- fae %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

#Plot
faePlot <- ggplot(fae, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Faecalibaculum unclassified"))) +
  theme_bw() +
  theme(legend.position = "none")
faePlot
ggsave("Faecalibaculum_unclassified_timepoint.pdf", faePlot, width = 12, height = 8)


#Read in Turibacter unclassified
turi <- read.csv("Turibacter_unclassified.csv", header = TRUE, sep = ",")
#Re-order
turi <- turi %>%
  mutate(Age = factor(Age, levels = Age_order))

#Re-order phenotype order
turi <- turi %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

#Plot
turiPlot <- ggplot(turi, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Turibacter unclassified"))) +
  theme_bw() +
  theme(legend.position = "none")
turiPlot
ggsave("Turibacter_unclassified_timepoint.pdf", turiPlot, width = 12, height = 8)

#Read in Akkermansia unclassified
akkU <- read.csv("Akkermansia_unclassified.csv", header = TRUE, sep = ",")
#Re-order
akkU <- akkU %>%
  mutate(Age = factor(Age, levels = Age_order))

#Re-order phenotype order
akkU <- akkU %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

#Plot
akkUPlot <- ggplot(akkU, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Akkermansia unclassified"))) +
  theme_bw() +
  theme(legend.position = "none")
akkUPlot
ggsave("Akkermansia_unclassified_timepoint.pdf", akkUPlot, width = 12, height = 8)


###Significant hits 6-10
top6.10 <- read.csv("top6-10_significantSpecies.csv", header = TRUE, sep = ",")

#Pivot longer the csv file
top6.10.long <- top6.10 %>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 
top6.10.long$Level <- "Species"
write.csv(top6.10.long, file = "Top6.10_sigSpecies_long.csv")

#Read in the meta-data included long format file
long6.10 <- read.csv("Top6.10_sigSpecies_long.csv", header = TRUE, sep = ",")

#Reorder factors
species_order <- c("Christensenella_timonensis",
                   "Cellulomonas_unclassified",
                   "Muribaculum_intestinale",
                   "Bacteroides_paurosaccharolyticus",
                   "Faecalibacterium_unclassified")


#Re-order
long6.10 <- long6.10 %>%
  mutate(taxon = factor(taxon, levels = species_order))

#Convert Age to a factor so 8 weeks comes before 20
Age_order <- c("8 Weeks",
               "20 Weeks")

#Re-order
long6.10 <- long6.10 %>%
  mutate(Age = factor(Age, levels = Age_order))

#Convert Phenotype order to factor and set desired order
Genotype_order <- c("WT", "KO")

#Re-order phenotype order
long6.10 <- long6.10 %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))


#Plotting all top 5 together
Box6.10 <- ggplot(long6.10, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = 5, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "skyblue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Significant hits 6-10"))) +
  theme_bw()
Box6.10
ggsave("Top6-10_SigSpecies_timepoint.pdf", Box6.10, width = 12, height = 8)

#Subset by C.timonensis
chris <- long6.10[long6.10$taxon == "Christensenella_timonensis", ]

#Plot
chrisPlot <- ggplot(chris, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Christensenella timonensis"))) +
  theme_bw() +
  theme(legend.position = "none")
chrisPlot
ggsave("Christensenella_timonensis_timepoint.pdf", chrisPlot, width = 12, height = 8)

#Subset by C.unclassified
cell <- long6.10[long6.10$taxon == "Cellulomonas_unclassified", ]

#Plot
cellPlot <- ggplot(cell, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Cellulomonas unclassified"))) +
  theme_bw() +
  theme(legend.position = "none")
cellPlot
ggsave("Cellulomonas_unclassified_timepoint.pdf", cellPlot, width = 12, height = 8)

#Subset by M.intestinale
muri <- long6.10[long6.10$taxon == "Muribaculum_intestinale", ]

#Plot
muriPlot <- ggplot(muri, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Muribaculum intestinale"))) +
  theme_bw() +
  theme(legend.position = "none")
muriPlot
ggsave("Muribaculum_intestinale_timepoint.pdf", muriPlot, width = 12, height = 8)

#Subset by B.paurosaccharolyticus
bact <- long6.10[long6.10$taxon == "Bacteroides_paurosaccharolyticus", ]

#Plot
bactPlot <- ggplot(bact, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Bacteroides paurosaccharolyticus"))) +
  theme_bw() +
  theme(legend.position = "none")
bactPlot
ggsave("Bacteroides_paurosaccharolyticus_timepoint.pdf", bactPlot, width = 12, height = 8)

#Subset by F. unclassified
faec <- long6.10[long6.10$taxon == "Faecalibacterium_unclassified", ]

#Plot
faecPlot <- ggplot(faec, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = expression(italic("Faecalibacterium unclassified"))) +
  theme_bw() +
  theme(legend.position = "none")
faecPlot
ggsave("Faecalibacterium_unclassified_timepoint.pdf", faecPlot, width = 12, height = 8)


#Show the interesting ones side by side
one <- rbind(akk, fae)
two <- rbind(one, turi)
three <- rbind(two, cell)
four <- rbind(three, muri)
exciting <- rbind(four, bact)

#Reorder factors
species_order <- c("Akkermansia_muciniphila",
                   "Faecalibaculum_unclassified",
                   "Turicibacter_unclassified",
                   "Cellulomonas_unclassified",
                   "Muribaculum_intestinale",
                   "Bacteroides_paurosaccharolyticus")


#Re-order
exciting <- exciting %>%
  mutate(taxon = factor(taxon, levels = species_order))

#Convert Age to a factor so 8 weeks comes before 20
Age_order <- c("8 Weeks",
               "20 Weeks")

#Re-order
exciting <- exciting %>%
  mutate(Age = factor(Age, levels = Age_order))

#Convert Phenotype order to factor and set desired order
Genotype_order <- c("WT", "KO")

#Re-order phenotype order
exciting <- exciting %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))


#Plotting all top 5 together
excitingPlot <- ggplot(exciting, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "skyblue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  ylim(0,400) +
  labs(title = expression(italic("Significant Species"))) +
  theme_bw() +
  ylim(0,400)
excitingPlot
ggsave("Good_looking_hits_from_top10_mostSignificantSpecies_timepoint.pdf", excitingPlot, width = 12, height = 8)


noBact <- exciting[exciting$taxon != "Bacteroides_paurosaccharolyticus", ]
excitingPlot <- ggplot(noBact, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "skyblue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  ylim(0,400) +
  labs(title = expression(italic("Significant Species"))) +
  theme_bw() +
  ylim(0,400)
excitingPlot
ggsave("Good_looking_hits_from_top10_mostSignificantSpecies_timepoint.pdf", excitingPlot, width = 12, height = 8)

#Read species counts and significant species
speciesCounts <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Species_analysis/Species_level_Counts.csv")
sigSpecs <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Species_analysis/significant_species.csv")
sigSpecsNames <- sigSpecs$taxa
speciesCounts <- speciesCounts[speciesCounts$taxon != "Unknown_genus_unclassified",]
speciesCounts <- speciesCounts[speciesCounts$taxon != "unknown_unclassified",]

#Take only the counts from taxa that are significant
speciesCounts.filt <- speciesCounts[speciesCounts$taxon %in% sigSpecsNames, ]
rownames(speciesCounts.filt) <- speciesCounts.filt$taxon
speciesCounts.filt$taxon <- NULL

#Read in the metadata
meta <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Metadata.csv", header = TRUE, sep = ",")

#RAW COUNTS
library("pheatmap")
speciesCounts.filt.mat <- as.matrix(speciesCounts.filt)
speciesCounts.filt.mat <- as.data.frame(t(speciesCounts.filt.mat))
speciesCounts.filt.mat$Sample_ID <- rownames(speciesCounts.filt.mat)
speciesCounts.filt.mat.meta <- merge(speciesCounts.filt.mat, meta, by = "Sample_ID", all = TRUE)
rownames(speciesCounts.filt.mat.meta) <- speciesCounts.filt.mat.meta$Sample_ID
speciesCounts.filt.mat.meta$Sample_ID <- NULL
# rownames(speciesCounts.filt.mat.meta) <- speciesCounts.filt.mat.meta$Sample_ID

Strain <- speciesCounts.filt.mat.meta$Strain
names(Strain) <- rownames(speciesCounts.filt.mat.meta)
Strain <- as.data.frame(Strain)
Strain$Genotype <- speciesCounts.filt.mat.meta$Genotype


heatmap <- pheatmap(speciesCounts.filt.mat.meta[,1:109],
                    gaps_row = 23, 
                    cluster_rows = FALSE,
                    cluster_cols = TRUE,
                    scale = "column",
                    fontsize = 6,
                    annotation_row = Strain,
                    main = "Significant Species")
ggsave("testing.pdf", heatmap)
              
#RELATIVE ABUNDANCE
relabund$taxon <- rownames(relabund)
RB.filt <- relabund[relabund$taxon %in% sigSpecs$taxa, ]

rownames(RB.filt) <- RB.filt$taxon
RB.filt$taxon <- NULL
RB.filt.t <- t(RB.filt)
RB.filt.t <- as.data.frame(RB.filt.t)
RB.filt.t$Sample_ID <- rownames(RB.filt.t)
RB <- merge(RB.filt.t, meta, by = "Sample_ID", all = TRUE)
rownames(RB) <- RB$Sample_ID
RB$Sample_ID <- NULL


Strain <- RB$Genotype
names(Strain) <- rownames(RB)
Strain <- as.data.frame(Strain)
Strain$Genotype <- RB$Strain

RB <- as.data.frame(RB)
write.csv(RB, file = "RelAbundHeatmap_Species.csv")

library(pheatmap)
RB <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/RelAbundHeatmap_Species.csv", header = TRUE, sep = ",")
rownames(RB) <- RB$X
HM <- pheatmap(RB[,2:110],
               gaps_row = c(23,47,69), 
               cluster_rows = FALSE,
               cluster_cols = TRUE,
               scale = "column",
               fontsize = 4,
               annotation_row = Strain,
               legend_width = 1,
               legend_height = 2)
ggsave("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/RelAbund_sigSpeces_Heatmap.pdf", HM, width = 12, height = 8)



              



















