library("ggplot2")
library("dplyr")
library("tidyverse")

#Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/")

#Read in the main data and significant species
shoreline <- read.csv("Data/all_shoreline_data_copy.csv", header = TRUE, sep = ",")
sigSpecs <- read.csv("Species_analysis/significant_species.csv", header = TRUE, sep = ",")

specNames <- sigSpecs$taxa
specLines <- shoreline[shoreline$taxon %in% specNames,]
species <- specLines$taxon
names(species) <- specLines$rankID
species <- as.data.frame(species)
species$RankID <- rownames(species)
species$rankID <- sub("^(\\d+\\.\\d+\\.\\d+\\.\\d+\\.\\d+\\.\\d+).*", "\\1", species$RankID)

familyRankIDs <- species$rankID
species$RankID <- NULL

families <- shoreline[shoreline$taxlevel == 5, ]
families$X <- NULL
family <- families[,c(2,5)]
familyNames <- species$rankID
familySubset <- family[family$rankID %in% familyNames,]

familyNames <- unique(familySubset$rankID)
names(familyNames) <- familySubset$taxon

merged <- merge(family, species, by = "rankID", all.y = FALSE)
uniquemerged <- !duplicated(merged)
merged <- subset(merged, uniquemerged)

#Prepare a file
sigSpecs$species <- sigSpecs$taxa
final <- merge(merged, sigSpecs, by = "species", all = TRUE)
final$X <- NULL
colnames(final) <- c("Species", "Family_rankID", "Family", "taxa", "p_adj")
final$taxa <- NULL
final <- final[,c(1,3,2,4)]
write.csv(final, file = "SignificantSpecies_and_their_Families.csv")





