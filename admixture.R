# Admixture Analysis

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
library(stringr)
library(ggplot2)
library(dplyr)

#set working directory
setwd("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/admixture")

cv <- read.table("cross_validation.txt")
cv
# Analyze the cross-validation results Then, add a K-cluster column indicating the number of K you test and select only two columns of interest, CV and K.

cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- cv[, c("V4", "K")]
CV
CV$K <- as.numeric(sub("K=", "", CV$K))

# Rename your two columns CV and K-cluster
colnames(CV) <- c("CV","K")

# Do a graph showing the cross validation results. Then select the optimal number of clusters regarding :
# the lowest cross validation error
# when the cross-validation error decrease the most

graph_title="Cross-Validation plot"
x_title="K"
y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV))
graph_1+geom_line()+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))

#Save the graph
ggsave("Admixture_cross-validation.pdf",width=7,height=5,dpi=600)
dev.off()

# One can modify the color palette to include the colors yellow2, blue, orange, and red2 and its shades 
# using the palette function and the colors parameter. You can define a custom palette
# with these colors. Here's how you can do it:
# Define a custom color palette
my_colors <- c("red", "orange", "blue", "lightblue")

# Set the custom palette
palette(my_colors)

#load admxiture function
source("admixFun.R") 

#list Q files and sort for K
files <- list.files("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/admixture", full = TRUE, pattern = "Q")
files <- files[order(as.numeric(sub('.*a.thaliana(\\d+)\\.Q$', '\\1', files)))]

#population file
pop <- scan("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/admixture/pop.txt",what="df",na="")
table(pop) 

# possible K
Kall <- 1:10

## read Qs
allQ <- list()
for(K in Kall)
  allQ[[K]]<-t(read.table(files[K-min(Kall)+1]))

#exact ordering (slow for large K)
plotMulti(allQ,Kall=2:3,reorder=1,pop,fast=T,lwd=1,lty=1)

#make smaller line and change type to solid line
plotMulti(allQ,Kall=2:10,reorder=1,pop,fast=T,lwd=1,lty=1)

plotMulti(allQ,Kall=2:10,reorder=1,pop,fast=T,lwd=0) 
