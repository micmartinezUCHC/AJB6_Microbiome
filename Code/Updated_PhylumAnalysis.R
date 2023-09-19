#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/Phylum_Level/")

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

#Read phylum csv
phy <- read.csv("/Users/mikemartinez/Desktop/AJB6_KnockOut/Results/Microbiome/Phylum_Level/PhylumData_Normalized2000.csv",
                header = TRUE, sep = ",")

#Pivot longer the csv file
phy.long <- phy %>%
  pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 
phy.long$Level <- "Phylum"
write.csv(phy.long, file = "newPhylum.long.csv")

#Read in metadated-edited phylum data 
phylum <- read.csv("newPhylum.long.csv", header = TRUE, sep = ",")

#Convert Age to a factor so 8 weeks comes before 20
Age_order <- c("8 Weeks",
               "20 Weeks")

#Re-order
phylum <- phylum %>%
  mutate(Age = factor(Age, levels = Age_order))

#Convert Phenotype order to factor and set desired order
Genotype_order <- c("WT", "KO")

#Re-order phenotype order
phylum <- phylum %>%
  mutate(Phenotype = factor(Genotype, levels = Genotype_order))

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


shared <- phylum %>%
  group_by(Sample_ID) %>%
  summarize(sobs = richness(count),
            Shannon = shannon(count),
            Simpson = simpson(count))
write.csv(shared, file = "alphaDiversity_Metrics_phylumLevel.csv")

#Read in the edited alphaDiversity_metric_FamilyLevel.csv
alphaDiv <- read.csv("alphaDiversity_Metrics_phylumLevel.csv", header = TRUE, sep = ",")

alphaDiv <- alphaDiv %>%
  group_by(Genotype)

alphaDiv <- alphaDiv %>%
  mutate(Age = factor(Age, levels = Age_order))
alphaDiv <- alphaDiv %>%
  mutate(Genotype = factor(Genotype, levels = Genotype_order))

library("stats")
#See if the data is parametric
spec.ShanHist <- hist(alphaDiv$Shannon)
spec.SimpHist <- hist(alphaDiv$Simpson)
shapiro.test(alphaDiv$Shannon) #If P-value is > 0.05, we assume normality
shapiro.test(alphaDiv$Simpson) #NOT NORMAL

#Plot distribution of Shannon Diversity
phy.Shannons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Shannon, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  #ylim(2.0,3.1) +
  labs(title = "Phylum level Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "bottom")
phy.Shannons_boxplot
ggsave("Shannon_PhylumLevel_diversity.pdf", phy.Shannons_boxplot, width = 12, height = 8)

#Plot distribution of Simpson Diversity
Simpsons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Simpson, fill = Strain)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 2) +
  geom_point(size = 0.3, position = "jitter") +
  facet_nested_wrap(~Age + Genotype, nrow = 1, scale = "free_x", 
                    strip.position = "top") +
  stat_compare_means(paired = FALSE, label = "p.format") +
  labs(title = "Phylum level Simpson Diversity") +
  #ylim(0,0.31) +
  theme_bw() +
  theme(legend.position = "bottom")
Simpsons_boxplot
ggsave("Simpsons_PhylumLevel_diversity.pdf", Simpsons_boxplot, width = 12, height = 8)

#Read in phylum long data and alpha div data
long <- read.csv("newPhylum.long.csv", header = TRUE, sep = ",")
alpha <- read.csv("alphaDiversity_Metrics_PhylumLevel.csv", header = TRUE, sep = ",")

#Testing for significance of individual taxa
comp <- long %>%
  group_by(Sample_ID) %>%
  group_by(Sample_ID) %>%
  ungroup() %>%
  inner_join(., alpha, by = "Sample_ID")

library(broom)
#Get the significant families
sig_phyla <- comp %>%
  nest(data = -taxon) %>%
  mutate(test = map(.x = data, ~aov(count~Strain.x + Genotype.x, data = .x) %>%
                      tidy)) %>%
  unnest(test) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  filter(p.adj < 0.05) %>%
  select(taxon, p.adj)
write.csv(sig_phyla, file = "significant_phyla.csv")

#Get the mean frequency
mean_freqs.filt <- long %>%
  group_by(Sample_ID, taxon) %>%
  summarise(phyCount = sum(count),
            .groups = 'drop') %>%
  group_by(Sample_ID) %>%
  summarise(phyFreq = phyCount / sum(phyCount),
            Phylum = taxon,
            .groups = 'drop') %>%
  group_by(Phylum) %>%
  summarise(mean = mean(phyFreq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs.filt)
















