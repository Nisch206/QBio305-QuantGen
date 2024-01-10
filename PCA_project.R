# Principle Component Analysis (PCA)


# load all required packages

library(PopGenome) 
library(dplyr)
library(ggplot2)
library (readr)
library(tibble)
library(vcfR)
library(adegenet)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library(ggrepel)
library(gplots)
library(StAMPP)
library(RColorBrewer)
library(plotly)

# Set the R working directory to the location where you have stored your indexed 'input.vcf' file 
# and the 'genomic_positions.bed' file." 

setwd("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/pca")

# If you want to check number columns and start and stop positions then you have to read 
# your "input.vcf" file using "read.vcfR" function from "vcfR" package

vcf_file <- read.vcfR("group_3_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")

# Convert VCF to genind object
genind_vcf <- vcfR2genind(vcf_file)


# Scale genind object for PCA
genind_vcf_scaled = scaleGen(genind_vcf, NA.method = "mean")

# Perform PCA
pca <- dudi.pca(genind_vcf_scaled, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 10)

# Check PCA dimensions
axis_all = pca$eig * 100 / sum(pca$eig)
barplot(axis_all[1:10], main = "PCA eigenvalues")

str(pca)
summary(pca)
#ind <- pca$rank
pca_axis <- pca$li
pca_eigenval <- pca$eig[1:10]
str(pca_axis)

# set names
ind_names<- read.table("group_3_accession_names.txt", header= FALSE)
str(ind_names)

# Add a new column "ind" to pca_axis using the values from ind
pca_axis$ind <- ind_names$V1

# You can also directly provide a list of names
population_labels <- c("SWE-S", "SWE-S", "SWE-N", "SWE-N", "SWE-N", "SWE-N", "SWE-S", "SWE-N", "SWE-N", "SWE-S",
                       "SWE-S", "SWE-S", "SWE-S", "SWE-S", "SWE-N", "SWE-N", "SWE-N", "SWE-N", "SWE-N", "SWE-S",
                       "SWE-N", "SWE-S", "SWE-N", "IT-N", "IT-N", "SWE-S", "SWE-S", "SWE-N", "SWE-S", "SWE-S",
                       "SWE-N", "SWE-N", "SWE-N", "SWE-N", "SWE-S", "SWE-S", "SWE-S", "SWE-S", "SWE-N", "SWE-N",
                       "SWE-S", "SWE-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S",
                       "IT-S", "IT-S", "IT-N", "IT-N", "IT-N", "IT-N", "IT-N", "IT-N", "IT-N", "IT-N", "IT-N", "IT-N",
                       "IT-N", "IT-N", "IT-N", "IT-N", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-S", "IT-N",
                       "IT-N", "IT-N", "IT-N", "IT-S", "IT-S", "IT-S") 

# Add a new column named "Population" to your data frame
pca_axis$Population <- population_labels

pop<-pca_axis$Population

# remake data.frame
pca_2 <- as.tibble(data.frame(pca_axis, population_labels))

n <- length(pca_eigenval)
n # use this number PC=1:n

# first convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = pca_eigenval/sum(pca_eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

str(pca_2)

# plot pca PC1 and PC2
b <- ggplot(pca_2, aes(Axis1, Axis2, col = pop)) + 
  geom_point(size = 4) +
  scale_colour_manual(values = c("red", "orange", "blue", "lightblue")) + #use same color and same sequence for admixture
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ggtitle("Arabidopsis thaliana Accessions from Italy (IT) and Sweden (SWE)") +
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        plot.margin = margin(20, 20, 20, 20))   # Adjust margin for title

# Print the plot
print(b)

# plot pca PC1 and PC3

b <- ggplot(pca_2, aes(Axis1, Axis3, col = pop)) + 
  geom_point(size = 4) +
  scale_colour_manual(values = c("red", "orange", "blue", "lightblue")) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  ggtitle("Arabidopsis thaliana Accessions from Italy (IT) and Sweden (SWE)") +
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        plot.margin = margin(20, 20, 20, 20))   # Adjust margin for title

# Print the plot
print(b)

#3D plot

fig <- plot_ly(pca_2, x = ~Axis1, y = ~Axis2, z = ~Axis3, color =pop, colors = c("red", "orange", "blue", "lightblue")) %>%
  add_markers(size = 12)

fig <- fig %>%
  layout(
    title = "Arabidopsis thaliana Accessions from Italy (ITA) and Sweden(SWE)",
    scene = list(bgcolor = "white")
  )

fig










