#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Microbiome analysis pipeline: A/J B6 mPGES-1 KO analysis
#  Mike Martinez
#  Rosenberg Lab, University of Connecticut Health Center
#  September 28th, 2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----------Load libraries
library(dplyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(roxygen2)

#----------Set working directory
setwd("/Users/mikemartinez/Desktop/AJB6_Microbiome/uBiome/")
parentDir <- file.path("/Users/mikemartinez/Desktop/AJB6_Microbiome/uBiome/")

#-----------Read in the microbiome data
raw <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Data/all_shoreline_data_copy.csv",
                header = TRUE, sep = ',')
raw$X <- NULL

#----------Read in the metadata
meta <- read.csv("/Users/mikemartinez/Desktop/AJB6_Microbiome/Data/Metadata.csv",
                 header = TRUE, sep = ",")

#-----------Get a list of all the unique tax levels (exclusing the first one because we don't want the root)
tax_levels <- unique(raw$taxlevel)[-1]
names(tax_levels) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
tax_levels <- as.data.frame(tax_levels)
tax_levels$taxonomy <- rownames(tax_levels)
taxonomic_levels <- tax_levels$taxonomy

#----------Initialize empty lists
taxonomy_dfs <- list() #Holds the raw counts for each taxonomic level
long_dfs <- list() #Holds the counts for each taxonomic level pivoted in long format
relative_abundance <- list() #Holds the relative abundances for each taxonomic level in long format

#---------------------------------------------#
###############################################
#####----------Declare functions----------#####
###############################################
#---------------------------------------------#
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
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 12)) +
    facet_nested_wrap(~Strain + Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    scale_y_continuous(name = "Relative Abundance",
                       labels = scales::percent) +
    theme(strip.background = element_rect(color = "black", fill = "lightgray"),
          panel.spacing = unit(0.2, "lines")) +
    theme(legend.text = element_text(size = 12)) +
    theme(plot.title = element_text(size = 16)) +
    theme(text = element_text(family = "Helvetica")) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
    labs(x = NULL,
         title = paste(tax_levels$taxonomy[long_df], "Level Relative Abundances", sep = " "))
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
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 12)) +
    facet_nested_wrap(~Strain + Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    scale_y_continuous(name = "Relative Abundance",
                       labels = scales::percent) +
    theme(strip.background = element_rect(color = "black", fill = "lightgray"),
          panel.spacing = unit(0.2, "lines")) +
    theme(legend.text = element_text(size = 12)) +
    theme(plot.title = element_text(size = 16)) +
    theme(text = element_text(family = "Helvetica")) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = NA, color ="black") +
    labs(x = NULL, 
         title = paste(tax_levels$taxonomy[long_df], "Significant Taxa", sep = " "))
  ggsave(paste(tax_level, "Significant_Barplot.pdf", sep = "_"), barplot, width = 12, height = 8)
}


#'Function to calculate alpha diversity metrics
#'@x is a dataframe of raw counts in the long format
richness <- function(x){
  sum(x > 0)
}

#'Function to calculate Shannon alpha diversity metric
#'@x is a dataframe of raw counts in the long format
shannon <- function(x){
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
}

#'Function to calculate Simpson alpha diversity metric
#'@x is a dataframe of raw counts in the long format
simpson <- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
}

plotAlpha <- function(x, tax_level) {
  
  alphaDiv <- x %>%
    group_by(Phenotype)
  
  Age_order <- c("8 Weeks", "20 Weeks")
  phenotype_order <- c("WT", "KO")
  
  x <- x %>%
    mutate(Age = factor(Age, levels = Age_order))
  x <- x %>%
    mutate(Phenotype = factor(Phenotype, levels = phenotype_order))
  
  Shannons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Shannon, fill = Strain)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 2) +
    geom_point(size = 0.3, position = "jitter") +
    facet_nested_wrap(~Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    stat_compare_means(paired = FALSE, label = "p.format") +
    labs(title = paste(tax_levels$taxonomy[long_df], "Shannon Diversity", sep = " ")) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste(tax_level, "Shannon.pdf", sep = "_"), Shannons_boxplot, width = 12, height = 8)
  
  
  Simpsons_boxplot <- ggplot(alphaDiv, aes(x = Strain, y = Simpson, fill = Strain)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 2) +
    geom_point(size = 0.3, position = "jitter") +
    facet_nested_wrap(~Age + Phenotype, nrow = 1, scale = "free_x", 
                      strip.position = "top") +
    stat_compare_means(paired = FALSE, label = "p.format") +
    labs(title = paste(tax_levels$taxonomy[long_df], "Simpson Diversity", sep = " ")) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste(tax_level, "Simpson.pdf", sep = "_"), Simpsons_boxplot, width = 12, height = 8)
  
}



#---------------------------------------------#
###############################################
#----------       Code Blocks       ----------#
###############################################
#---------------------------------------------#

#----------Create a vector of output directories for each taxonomic level
outDirs <- list()
for (tax_level in seq_along(taxonomic_levels)){
  output_dir <- file.path("/Users/mikemartinez/Desktop/AJB6_Microbiome/uBiome", taxonomic_levels[tax_level])
  print(output_dir)
  outDirs[[tax_level]] <- output_dir
}

#----------For each taxonomic level, subset the raw counts to only include the i-th taxonomic level
for (i in 1:nrow(tax_levels)) {
  #Create directory for each taxonomic level
  dir.create(outDirs[[i]])
  temp <- raw[raw$taxlevel == i,]
  taxonomy_dfs[[i]] <- temp
  write.csv(temp, file = paste(outDirs[[i]], paste(tax_levels$taxonomy[i], "Counts.csv", sep = "_"), sep = "/"))
  
  #Reset working directory to parent directory
  setwd(parentDir)
}

#----------Iterate over the taxonomy_dfs list and pivot longer, merge metadata, and store to list
for(df in 1:length(taxonomy_dfs)) {
  
  #Set working directory
  setwd(outDirs[[df]])
  
  #Get just the taxa names and counts columns
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
  write.csv(long.merged, file = paste(outDirs[[df]], paste(tax_levels$taxonomy[df], "long.meta.csv", sep = "_"), sep = "/"))
  
  
  #Calculate alpha diversity metrics using the alpha diversity functions
  alpha <- long.merged %>%
    group_by(Sample_ID) %>%
    summarize(sobs = richness(count),
              Shannon = shannon(count),
              Simpson = simpson(count))
  
  #Merge alpha diversity metrics to counts and metadata and write to a csv file
  merged.alpha <- merge(meta, alpha, by = "Sample_ID", all.y = TRUE)
  write.csv(merged.alpha, file = paste(outDirs[[df]], paste(tax_levels$taxonomy[df], "AlphaDiversity.csv", sep = "_"), sep = "/"))
  
  #Plot alpha diversity plots using function
  plotAlpha(merged.alpha, tax_levels$taxonomy[df])
  
  #Reset parent directory
  setwd(parentDir)
  
}

#----------Iterate over the list including metadata-including long dataframes and write as a csv
for (long_df in 1:length(long_dfs)) {
  
  #Set working directory for specified taxonomic level and write the Long.Meta data to a csv file
  setwd(outDirs[[long_df]])
  write.csv(long_dfs[[long_df]], file = paste(outDirs[[long_df]], paste(tax_levels$taxonomy[long_df], "Long.Meta.csv", sep = "_"), sep = "/"))
  
  #Reset parent directory
  setwd(parentDir)
}

#----------Iterate through the list and get the relative abundance dataframes to plot relative abundance barplots
for (long_df in 2:length(relative_abundance)) {
  
  #Set working directory for specified taxonomic level and write the RelAbund data to a csv file
  setwd(outDirs[[long_df]])
  write.csv(relative_abundance[[long_df]], file = paste(outDirs[[long_df]], paste(tax_levels$taxonomy[long_df], "RelAbund.csv", sep = "_"), sep = "/"))
  
  #Get the relative abundance data frame for the given taxonomic level and remove instances of "Unknown"
  rawRelativeAbundance <- relative_abundance[[long_df]]
  rawRelativeAbundance <- rawRelativeAbundance[!grepl("unknown", rawRelativeAbundance$taxon, ignore.case = TRUE), ]
  abundances <- rawRelativeAbundance
  
  #Test for significance of every taxa
  significant_taxa <- significance(abundances)
  write.csv(significant_taxa, file = paste(outDirs[[long_df]], paste(tax_levels$taxonomy[long_df], "SignificantTaxa.csv", sep = "_"), sep = "/"))
  
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
    plotSig(subsetSig(topSig(significant_taxa), relative_abundance[[long_df]]),tax_levels$taxonomy[long_df])
    
    #Reset parent directory
    setwd(parentDir)
    
  } else {
    setwd(parentDir)
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