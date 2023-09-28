#Modularizing the AJB6 Microbiome code
library(dplyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(roxygen2)

#Set working directory
setwd("/Users/mikemartinez/Desktop/uBiome/")


#Read in the microbiome data
raw <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Data/all_shoreline_data_copy.csv",
                header = TRUE, sep = ',')
raw$X <- NULL

#Get a list of all the unique tax levels (exclusing the first one because we don't want the root)
tax_levels <- unique(raw$taxlevel)[-1]
names(tax_levels) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
tax_levels <- as.data.frame(tax_levels)
tax_levels$taxonomy <- rownames(tax_levels)

#Iterate through the raw data for each tax level, collating a dataframe of JUST that tax level
#Creae an empty list to hold the 8 data frames
taxonomy_dfs <- list()

for (i in 1:nrow(tax_levels)) {
  print(tax_levels$taxonomy[i])
  print(i)
  temp <- raw[raw$taxlevel == i,]
  taxonomy_dfs[[i]] <- temp
}

#Read in the metadata
meta <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Data/Metadata.csv",
                 header = TRUE, sep = ",")


#####----------Declare functions----------#####
#'Function to calculate relative abundances
#'@counts a dataframe where column 1 is taxa names at a given taxonomic level followed by raw counts for all samples
#'@relabund return value is a data frame containing the relative abundances of each taxa provided in the input df
#Function to calculate relative abundance
relAbund <- function(counts) {
  sums <- rowSums(counts[,2:ncol(counts)]) + 0.01
  relabund <- counts[,2:ncol(counts)]/sums
  return(relabund)
}

#'Function to calculate top 12 most frequent taxa
#'@relabund a relative abundance data frame that has been pivoted to long format
#'@freq a data frame where column 1 is the taxa name and column 2 is the mean frequency of that taxa
meanFreq <- function(relabund) {
  freq <- relabund %>%
    group_by(Sample_ID, taxon) %>%
    summarise(Count = sum(RelAbund),
              .groups = 'drop') %>%
    group_by(Sample_ID) %>%
    summarise(Freq = Count / sum(Count),
              Taxa = taxon,
              .groups = 'drop') %>%
    group_by(Taxa) %>%
    summarise(mean = mean(Freq),
              .groups = 'drop') %>%
    arrange(desc(mean))
  
  return(as.data.frame(freq))
}

#'Function for relative abundance barplots
#'@x a data frame of relative abundances only for the most frequent taxa
#'@tax_level the taxonomic level you are plotting
#'@taxa_order the desired order of the taxa you are plotting (decreasing mean frequency, as a factor)
barplot <- function(x, tax_level, taxa_order) {
  Age_order <- c("8 Weeks", "20 Weeks")
  phenotype_order <- c("WT", "KO")
  
  x <- x %>%
    mutate(Age = factor(Age, levels = Age_order))
  x <- x %>%
    mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
  x <- x %>%
    mutate(taxon = factor(taxon, levels = taxa_order))
  
  barplot <- ggplot(x, (aes(x = Sample_ID, y = RelAbund))) +
    geom_bar(aes(fill = taxon), stat = "identity", position = "fill", width = 1) +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90, size = 4.5),
          strip.text = element_text(face = "bold")) +
    facet_nested_wrap(~Strain + Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    scale_y_continuous(name = "Relative Abundance",
                       labels = scales::percent) +
    theme(strip.background = element_rect(color = "black", fill = "lightgray"),
          panel.spacing = unit(0.2, "lines")) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
    labs(title = paste(tax_levels$taxonomy[long_df], "Level Relative Abundances", sep = " "))
  ggsave(paste(tax_level, "RelAbund_Barplot.pdf", sep = "_"), barplot, width = 12, height = 8)
}

#'Function to test for significance at the taxa level
#'@x a relative abundance data frame
#'@significance return value of significant taxa at a BH-corrected pvalue of 0.05
significance <- function(x){
  significant <- x %>%
    nest(data = -taxon) %>%
    mutate(test = map(.x = data, ~aov(RelAbund~Strain + Phenotype, data = .x) %>%
                        tidy)) %>%
    unnest(test) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
    filter(p.adj < 0.05) %>%
    select(taxon, p.adj)
  return(significant)
}


#Initialize an empty list to hold the metadata-including long data frames
long_dfs <- list()
relative_abundance <- list()

#----------Iterate over the taxonomy_dfs list and pivot longer, merge metadata, and store to list
for(df in 1:length(taxonomy_dfs)) {
  
  write.csv(taxonomy_dfs[[df]], file = paste(tax_levels$taxonomy[df], "Counts.csv", sep = "_"))
  counts <- taxonomy_dfs[[df]][,5:96] #Only the counts columns
  names <- counts$taxon #Taxon names to append as a new column in the abundances df
  
  #Run the relative abundance function and pivot longer
  abundances <- relAbund(counts) #Run function
  abundances$taxon <- names #Append names to new column
  
  #Pivot longer and merge metadata
  abund.long <- abundances %>%
    pivot_longer(-taxon, names_to = "Sample_ID", values_to = "RelAbund")
  abund.long.merged <- merge(abund.long, meta, by = "Sample_ID", all.y = TRUE)
  relative_abundance[[df]] <- abund.long.merged
  
  #Pivot longer the data frame, only keep taxonomy column and non-control sample columns
  long <- taxonomy_dfs[[df]][,5:96] %>%
    pivot_longer(-taxon, names_to = "Sample_ID", values_to = "count") 
  long.merged <- merge(long, meta, by = "Sample_ID", all.y = TRUE)
  long_dfs[[df]] <- long.merged
}

#----------Iterate over the list including metadata-including long dataframes and write as a csv
for (long_df in 1:length(long_dfs)) {
  write.csv(long_dfs[[long_df]], file = paste(tax_levels$taxonomy[long_df], "Long_meta.csv", sep = "_"))
}



#----------Iterate through the list and get the relative abundance dataframes to plot relative abundance barplots
#####-----Declare Functions-----#####
#'Function to order the significant taxa from smallest padj value to largest
#'@sigTaxa an input dataframe of the significant taxa
#'@top12 a return value of the ordered top 12 taxa
topSig <- function(sigTaxa) {
  #Order the significant taxa from lowest p.adj value to highest
  sigTaxaOrdered <- sigTaxa[order(sigTaxa$p.adj, decreasing = FALSE),]
  
  #Take the top 12
  sigTaxaOrderedTop <- sigTaxaOrdered[1:12,]
  top12 <- as.vector(sigTaxaOrderedTop$taxon)
  return(top12)
}

#'Function to subset the significant taxa from the relative abundance data frames
#'@sigTaxaOrderedTop an input vector of the top 12 significant taxa names
#'@relAbundDF a data frame of unfiltered relative abundances in long format
#'@significantSubset a return value of a long pivoted data frame only containing relative abundances for the top 12 significant taxa
subsetSig <- function(sigTaxaOrderedTop, relAbundDF) {
  significantSubset <- relAbundDF[relAbundDF$taxon %in% sigTaxaOrderedTop,]
  return(significantSubset)
}

#'Function to plot the top 12 most significant taxa and their relative abundances
#'@x a dataframe of the subsetted taxa and their relative abundances
#'@tax_level the level of taxonomy being plotted
plotSig <- function(x, tax_level) {
  Age_order <- c("8 Weeks", "20 Weeks")
  phenotype_order <- c("WT", "KO")
  
  x <- x %>%
    mutate(Age = factor(Age, levels = Age_order))
  x <- x %>%
    mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
  
  barplot <- ggplot(x, (aes(x = Sample_ID, y = RelAbund))) +
    geom_bar(aes(fill = taxon), stat = "identity", position = "fill", width = 1) +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90, size = 4.5),
          strip.text = element_text(face = "bold")) +
    facet_nested_wrap(~Strain + Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    scale_y_continuous(name = "Relative Abundance",
                       labels = scales::percent) +
    theme(strip.background = element_rect(color = "black", fill = "lightgray"),
          panel.spacing = unit(0.2, "lines")) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
    labs(title = paste(tax_levels$taxonomy[long_df], "Significant Taxa", sep = " "))
  ggsave(paste(tax_level, "Significant_Barplot.pdf", sep = "_"), barplot, width = 12, height = 8)
}

#####-----Run loop-----#####
for (long_df in 2:length(relative_abundance)) {
  write.csv(relative_abundance[[long_df]], file = paste(tax_levels$taxonomy[long_df], "RelAbund.csv", sep = "_"))
  
  #Get the relative abundance data frame for the given taxonomic level and remove instances of "Unknown"
  rawRelativeAbundance <- relative_abundance[[long_df]]
  rawRelativeAbundance <- rawRelativeAbundance[!grepl("unknown", rawRelativeAbundance$taxon, ignore.case = TRUE), ]
  abundances <- rawRelativeAbundance

  #Test for significance of every taxa
  significant_taxa <- significance(abundances)
  write.csv(significant_taxa, file = paste(tax_levels$taxonomy[long_df], "SignificantTaxa.csv", sep = "_"))
  
  #Run the frequencies function and take the top 12 frequencines (often times, 1 of them is "Unknown")
  frequencies <- as.data.frame(meanFreq(abundances))
  top_freq_taxa <- head(frequencies,12)
  taxa_keep <- top_freq_taxa$Taxa
  
  #Filter the abundance df only for the top taxa
  freq_abundances <- abundances[abundances$taxon %in% taxa_keep,]
  print(length(unique(freq_abundances$taxon)))
  
  #Plot: supply a df, tax level, and factor order for the legend
  barplot(freq_abundances, tax_levels$taxonomy[long_df], taxa_keep)
  
  #If there are more than 12 significant hits, run the functions
  if (nrow(significant_taxa) > 12) {

    #Order, subset, and plot the top 12 most significant taxa
    top12SignificantTaxa <- plotSig(subsetSig(topSig(significant_taxa), relative_abundance[[long_df]]),tax_levels$taxonomy[long_df])

  } else {
    next
  }
}








#Pie chart for F/B ratio between groups
#Read in the phylum long relative abundance data
phylum <- read.csv("Phylum_RelAbund.csv", header = TRUE, sep = ",")
taxa_subset <- c("Bacteroidetes", "Firmicutes")
BFRatio <- phylum[phylum$taxon %in% taxa_subset,]

BFRatio$group <- ifelse(BFRatio$Strain == "AJ" & BFRatio$Phenotype == "WT", "AJ_WT",
                        ifelse(BFRatio$Strain == "AJ" & BFRatio$Phenotype == "KO", "AJ_KO",
                               ifelse(BFRatio$Strain == "B6" & BFRatio$Phenotype == "WT", "B6_WT", "B6_KO")))

unique_groups <- unique(BFRatio$group)

pieChart_list <- list()

for (group in unique_groups) {
  group_data <- BFRatio[BFRatio$group == group, ]

    pie_chart <- ggplot(group_data, aes(x = "", y = RelAbund, fill = taxon)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      labs(title = group) +
      theme_void() +
      theme(legend.position = "none")
    ggsave(paste(group, "BF_PieChart.pdf", sep = "_"), pie_chart)
    pieChart_list[[group]] <- pie_chart
}  

FBRatio <- cowplot::plot_grid(pieChart_list[["AJ_WT"]], pieChart_list[["AJ_KO"]], pieChart_list[["B6_WT"]], pieChart_list[["B6_KO"]])
FBRatio
ggsave("FB_ratio.pdf", FBRatio, width = 12, height = 8)




