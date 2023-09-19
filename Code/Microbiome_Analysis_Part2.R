#The purpose of this script is to visualize the order level of the AJ/B6 Microbiome
#I want to find the order level relative abundances and plot the top 10

#Load libraries
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

#Calculate the rowSums of the family counts
sums <- rowSums(familyCounts[,2:ncol(familyCounts)])
relabund <- familyCounts[,2:ncol(familyCounts)]/sums
rownames(relabund) <- names
relabund <- na.omit(relabund)

#Pivot longer the data frame
relabund$taxa <- rownames(relabund)
fam.relabund.long <- relabund %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "count") 
fam.relabund.long$Level <- "Family"
write.csv(fam.relabund.long, file = "family.relative.abundances.long.csv")

#Get the mean frequency
mean_freqs <- fam.relabund.long %>%
  group_by(Sample_ID, taxa) %>%
  summarise(famCount = sum(count),
            .groups = 'drop') %>%
  group_by(Sample_ID) %>%
  summarise(famFreq = famCount / sum(famCount),
            Family = taxa,
            .groups = 'drop') %>%
  group_by(Family) %>%
  summarise(mean = mean(famFreq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs)
write.csv(mean_freqs, file = "Mean_Frequency_relAbundance_FamilyLevel.csv")


#Let's take the top 10 taxa from the mean frequency list
top20Fams <- mean_freqs$Family[1:20]
familyCounts.filt <- familyCounts[familyCounts$taxon %in% top20Fams,]

#Calculate the relative abundance from this filtered dataframe
names.filt <- familyCounts.filt$taxon
sums.filt <- rowSums(familyCounts.filt[,2:ncol(familyCounts.filt)])
relabund.filt <- familyCounts.filt[,2:ncol(familyCounts.filt)]/sums.filt
rownames(relabund.filt) <- names.filt
relabund.filt <- na.omit(relabund.filt)

#Pivot longer the filtered data frame
relabund.filt$taxa <- rownames(relabund.filt)
fam.relabund.long.filt <- relabund.filt %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "count") 
fam.relabund.long.filt$Level <- "Family"
write.csv(fam.relabund.long.filt, file = "top20_families_based_on_meanRelAbund_frequency.csv")


#Read in the top 20 families data
top20 <- read.csv("top20_families_based_on_meanRelAbund_frequency.csv", header = TRUE, sep = ",")

#Get the mean frequency
mean_freqs.filt <- fam.relabund.long.filt %>%
  group_by(Sample_ID, taxa) %>%
  summarise(famCount = sum(count),
            .groups = 'drop') %>%
  group_by(Sample_ID) %>%
  summarise(famFreq = famCount / sum(famCount),
            Family = taxa,
            .groups = 'drop') %>%
  group_by(Family) %>%
  summarise(mean = mean(famFreq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs.filt)

#Convert phyla to factors with a specific order based on the order obtained above (mean frequency)
Family_order <- c("Ruminococcaceae",
                 "Muribaculaceae",
                 "Erysipelotrichaceae",
                 "Bacteroidaceae",
                 "Rikenellaceae",
                 "Spiroplasmataceae",
                 "Thermoanaerobacteraceae",
                 "Chitinophagaceae",
                 "Deferribacteraceae",
                 "Aurantimonadaceae",
                 "Lactobacillales_unclassified",
                 "Acholeplasmataceae",
                 "Sterolibacteriaceae",
                 "Desulfomicrobiaceae",
                 "Nannocystaceae",
                 "Clostridiales_Family_XVII._Incertae_Sedis",
                 "Nocardioidaceae",
                 "Polyangiaceae",
                 "Ectothiorhodospiraceae",
                 "Desulfuromonadales_unclassified")
Age_order <- c("8 Weeks", "20 Weeks")
genotype_order <- c("WT", "KO")

#Added in phenotype and strain columns in excel
top10Families <- read.csv("top10_families_based_on_meanRelAbund_frequency.csv", header = TRUE, sep = ",")
write.csv(top10Families, file = "top10_families_based_on_meanRelAbund_frequency.csv")


Family_order <- c("Aurantimonadaceae",
                  "Deferribacteraceae",
                  "Acholeplasmataceae",
                  "Nocardioidaceae",
                  "Ectothiorhodospiraceae",
                  "Desulfomicrobiaceae",
                  "Polyangiaceae",
                  "Clostridiales_Family_XVII._Incertae_Sedis",
                  "Bacteroidaceae",
                  "Desulfuromonadales_unclassified")

#Re-order
top10Families <- top10Families %>%
  mutate(taxa = factor(taxa, levels = Family_order))
top10Families <- top10Families %>%
  mutate(AGE = factor(AGE, levels = Age_order))
top10Families <- top10Families %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

#Re-order
top20 <- top20 %>%
  mutate(taxa = factor(taxa, levels = Family_order))
top20 <- top20 %>%
  mutate(AGE = factor(AGE, levels = Age_order))
top20 <- top20 %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

#Plot
top10Family_barplot <- ggplot(top10Families, (aes(x = Sample_ID, y = count))) +
  geom_bar(aes(fill = taxa), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  facet_nested_wrap(~Strain + AGE + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
          scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Family Level Relative Abundances")
top10Family_barplot
ggsave("top10Family_relAbund_barplot.pdf", top10Family_barplot, width = 12, height = 8)

#Timepoint Plots
#Need to get the families of interest subsetted from the familyLong counts data
top10Fams <- mean_freqs$Family[1:10]
familyCounts.filt <- familyCounts[familyCounts$taxon %in% top10Fams,]


familyCounts.filt.long <- familyCounts.filt %>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 
familyCounts.filt.long$Level <- "Family"
write.csv(familyCounts.filt.long, file = "top10_families_long.csv")

#Read in the csv
top10TP <- read.csv("top10_families_long.csv", header = TRUE, sep = ",")
top10TP <- top10TP %>%
  mutate(taxa = factor(taxon, levels = Family_order))

#Subset for Acholesplasmatacea
def <- top10TP[top10TP$taxon == "Deferribacteraceae",]
defplot <- ggplot(Ach, aes(x = Genotype, y = count, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = "Deferribacteraceae") +
  theme_bw() +
  theme(legend.position = "none")
defplot


timepoint <- ggplot(top10TP, aes(x = Genotype, y = LogCounts, fill = taxon)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, outlier.color = "black") +
  #geom_point(position = position_jitter(width = 0.09), alpha = 0.5, color = "blue") +
  facet_nested_wrap(~ Strain + Age, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  labs(x = "Strain", y = "Normalized Read Counts/Sample") +
  labs(title = "Top 10 Families") +
  theme_bw() +
  theme(legend.position = "right")
timepoint

#Plot
top20Family_barplot <- ggplot(top20, (aes(x = Sample_ID, y = count))) +
  geom_bar(aes(fill = taxa), stat = "identity", position = "fill", width = 1) +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  facet_nested_wrap(~Strain + AGE + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Top 20 Family Level Relative Abundances")
top20Family_barplot
ggsave("top20Family_relAbund_barplot.pdf", top20Family_barplot, width = 12, height = 8)


famCounts.long <- familyCounts%>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 


#Calculate alpha diversity metrics
richness <- function(x){
  sum(x > 0)
}

shannon <- function(x){
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
}


shared <- famCounts.long %>%
  group_by(Sample_ID) %>%
  summarize(sobs = richness(count),
            Shannon = shannon(count),
            Simpson = simpson(count))
write.csv(shared, file = "alphaDiversity_Metrics_FamilyLevel.csv")


#Read in the edited alphaDiversity_metric_FamilyLevel.csv
alphaDiv <- read.csv("alphaDiversity_Metrics_FamilyLevel.csv", header = TRUE, sep = ",")

alphaDiv <- alphaDiv %>%
  group_by(Genotype)

alphaDiv <- alphaDiv %>%
  mutate(Age = factor(Age, levels = Age_order))
alphaDiv <- alphaDiv %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

#See if the data is parametric
ShanHist <- hist(alphaDiv$Shannon)
SimpHist <- hist(alphaDiv$Simpson)
shapiro.test(alphaDiv$Shannon) #If P-value is > 0.05, we assume normality
shapiro.test(alphaDiv$Simpson)

#Plot distribution of Shannon Diversity
Shannons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Shannon, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  ylim(2.0,3.1) +
  labs(title = "Family level Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "bottom")
Shannons_boxplot
ggsave("Shannon_FamilyLevel_diversity.pdf", Shannons_boxplot, width = 12, height = 8)

#Plot distribution of Simpson Diversity
Simpsons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Simpson, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  labs(title = "Family level Simpson Diversity") +
  ylim(0,0.31) +
  theme_bw() +
  theme(legend.position = "bottom")
Simpsons_boxplot
ggsave("Simpsons_FamilyLevel_diversity.pdf", Simpsons_boxplot, width = 12, height = 8)

#Testing for significance of individual taxa
familyCounts.comp <- fam.relabund.long %>%
  group_by(Sample_ID) %>%
  group_by(Sample_ID) %>%
  ungroup() %>%
  inner_join(., alphaDiv, by = "Sample_ID")

library(broom)
#Get the significant families
sig_families <- familyCounts.comp %>%
  nest(data = -taxa) %>%
  mutate(test = map(.x = data, ~aov(count~Strain + Genotype, data = .x) %>%
                      tidy)) %>%
  unnest(test) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  filter(p.adj < 0.05) %>%
  select(taxa, p.adj)
write.csv(sig_families, file = "significant_families.csv")

sigFamNames <- sig_families$taxa
sigFamilies.df <- familyCounts.comp[familyCounts.comp$taxa %in% sigFamNames,]

#Plot the relative abundances of the statistically significant families
Significant_FamBarplot <- ggplot(sigFamilies.df, aes(x=Sample_ID, y = count)) +
  geom_bar(aes(fill = taxa), stat = "identity", position = "fill", width = 1) +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  facet_nested_wrap(~Strain + Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  guides(fill = guide_legend(ncol = 2))
Significant_FamBarplot
ggsave("StatisticallySignificant_familiesRelAbund.pdf", Significant_FamBarplot, width = 12, height = 8)

  
#############################
setwd("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/")

#Read in the raw family counts
speciesCounts <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Species_analysis/Species_level_Counts.csv", header = TRUE, sep = ",")

#Remove any unknown families
speciesCounts <- speciesCounts[speciesCounts$taxon != "Unknown_genus_unclassified",]
speciesCounts <- speciesCounts[speciesCounts$taxon != "unknown_unclassified",]

#Set taxa names
names <- speciesCounts$taxon

#Calculate the rowSums of the family counts
sums <- rowSums(speciesCounts[,2:ncol(speciesCounts)])
relabund <- speciesCounts[,2:ncol(speciesCounts)]/sums
rownames(relabund) <- names
relabund <- na.omit(relabund)

#Pivot longer the data frame
relabund$taxa <- rownames(relabund)
spec.relabund.long <- relabund %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "count") 
spec.relabund.long$Level <- "Species"
write.csv(spec.relabund.long, file = "species.relative.abundances.long.csv")

#Get the mean frequency
mean_freqs <- spec.relabund.long %>%
  group_by(Sample_ID, taxa) %>%
  summarise(specCount = sum(count),
            .groups = 'drop') %>%
  group_by(Sample_ID) %>%
  summarise(specFreq = specCount / sum(specCount),
            Species = taxa,
            .groups = 'drop') %>%
  group_by(Species) %>%
  summarise(mean = mean(specFreq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs)
write.csv(mean_freqs, file = "mean_Species_frequency.csv")


#Let's take the top 10 taxa from the mean frequency list
top10spec <- mean_freqs$Species[1:10]
specCounts.filt <- speciesCounts[speciesCounts$taxon %in% top10spec,]

#Calculate the relative abundance from this filtered dataframe
spec.names.filt <- specCounts.filt$taxon
spec.sums.filt <- rowSums(specCounts.filt[,2:ncol(specCounts.filt)])
spec.relabund.filt <- specCounts.filt[,2:ncol(specCounts.filt)]/spec.sums.filt
rownames(spec.relabund.filt) <- spec.names.filt
spec.relabund.filt <- na.omit(spec.relabund.filt)

#Pivot longer the filtered data frame
spec.relabund.filt$taxa <- rownames(spec.relabund.filt)
spec.relabund.long.filt <- spec.relabund.filt %>%
  pivot_longer(-taxa, names_to = "Sample_ID", values_to = "count") 
spec.relabund.long.filt$Level <- "Species"
write.csv(spec.relabund.long.filt, file = "top10_species_based_on_meanRelAbund_frequency.csv")

#Added in phenotype and strain columns in excel
top10species <- read.csv("Species_analysis/top10_species_based_on_meanRelAbund_frequency.csv", header = TRUE, sep = ",")

#Convert phyla to factors with a specific order based on the order obtained above (mean frequency)
species_order <- c("Alistipes sp. Marseille-P5997",
                   "Thioalkalivibrio unclassified",
                   "Nocardioidaceae unclassified",
                   "Phocea unclassified",
                   "Martelella unclassified",
                   "Cytophaga unclassified",
                   "Desulforegula unclassified",
                   "Coprococcus comes",
                   "Thermaerobacter unclassified",
                   "Clostridium sp. MF28")
Age_order <- c("8 Weeks", "20 Weeks")
genotype_order <- c("WT", "KO")


#Re-order
top10species <- top10species %>%
  mutate(taxa = factor(taxa, levels = species_order))
top10species <- top10species %>%
  mutate(Age = factor(Age, levels = Age_order))
top10species <- top10species %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

#Plot
top10species_barplot <- ggplot(top10species, (aes(x = Sample_ID, y = count))) +
  geom_bar(aes(fill = taxa), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  facet_nested_wrap(~Strain + Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  labs(title = "Species Level Relative Abundances")
top10species_barplot
ggsave("top10Species_relAbund_barplot.pdf", top10species_barplot, width = 12, height = 8)

specCounts.long <- speciesCounts%>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 


#Calculate alpha diversity metrics
richness <- function(x){
  sum(x > 0)
}

shannon <- function(x){
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
}

simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
}


shared <- specCounts.long %>%
  group_by(Sample_ID) %>%
  summarize(sobs = richness(count),
            Shannon = shannon(count),
            Simpson = simpson(count))
write.csv(shared, file = "alphaDiversity_Metrics_speciesLevel.csv")

#Read in the edited alphaDiversity_metric_FamilyLevel.csv
spec.alphaDiv <- read.csv("Species_analysis/alphaDiversity_Metrics_speciesLevel.csv", header = TRUE, sep = ",")

spec.alphaDiv <- spec.alphaDiv %>%
  group_by(Genotype)

spec.alphaDiv <- spec.alphaDiv %>%
  mutate(Age = factor(Age, levels = Age_order))
spec.alphaDiv <- spec.alphaDiv %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

library("stats")
#See if the data is parametric
spec.ShanHist <- hist(spec.alphaDiv$Shannon)
spec.SimpHist <- hist(spec.alphaDiv$Simpson)
shapiro.test(alphaDiv$Shannon) #If P-value is > 0.05, we assume normality
shapiro.test(alphaDiv$Simpson) #NOT NORMAL

#Plot distribution of Shannon Diversity
spec.Shannons_boxplot <- ggplot(spec.alphaDiv, aes(x = Strain, y = Shannon, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  #ylim(2.0,3.1) +
  labs(title = "Species level Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "bottom")
spec.Shannons_boxplot
ggsave("Shannon_SpeciesLevel_diversity.pdf", spec.Shannons_boxplot, width = 12, height = 8)

#Plot distribution of Simpson Diversity
spec.Simpsons_boxplot <- ggplot(spec.alphaDiv, aes(x = Strain, y = Simpson, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  labs(title = "Species level Simpson Diversity") +
  #ylim(0,0.31) +
  theme_bw() +
  theme(legend.position = "bottom")
Simpsons_boxplot
ggsave("Simpsons_SpeciesLevel_diversity.pdf", spec.Simpsons_boxplot, width = 12, height = 8)

#Testing for significance of individual taxa
specCounts.comp <- spec.relabund.long %>%
  group_by(Sample_ID) %>%
  group_by(Sample_ID) %>%
  ungroup() %>%
  inner_join(., spec.alphaDiv, by = "Sample_ID")

library(broom)
#Get the significant families
sig_specs <- specCounts.comp %>%
  nest(data = -taxa) %>%
  mutate(test = map(.x = data, ~aov(count~Strain + Genotype, data = .x) %>%
                      tidy)) %>%
  unnest(test) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  filter(p.adj < 0.05) %>%
  select(taxa, p.adj)
sig_specs <- sig_specs[,c(1,9)]
write.csv(sig_specs, file = "significant_species.csv")

sig.specNames <- sig_specs$taxa
sig_specs.df <- specCounts.comp[specCounts.comp$taxa %in% sig.specNames,]

#Plot the relative abundances of the statistically significant species
Significant_specBarplot <- ggplot(sig_specs.df, aes(x=Sample_ID, y = count)) +
  geom_bar(aes(fill = taxa), stat = "identity", position = "fill", width = 1) +
  theme(axis.text.x = element_text(angle = 90, size = 4.5),
        strip.text = element_text(face = "bold")) +
  facet_nested_wrap(~Strain + Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  scale_y_continuous(name = "Relative Abundance",
                     labels = scales::percent) +
  theme(strip.background = element_rect(color = "black", fill = "lightgray"),
        panel.spacing = unit(0.2, "lines")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
  guides(fill = guide_legend(ncol = 4))
Significant_specBarplot
ggsave("StatisticallySignificant_speciesRelAbund.pdf", Significant_specBarplot, width = 32, height = 8)









  










