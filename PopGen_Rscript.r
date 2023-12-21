
#########################################################
# Packages
########################################################
options(scipen=999) # remove scientific notation of very low integers

library(qqman)
library(ggplot2)
library(car)
library(qtl)    # load the qtl library for linkage analysis

#########################################################
# Data
#########################################################

dta.main <- read.csv("dta_ANOVA_QTLanalysis.csv")
map <- read.csv("input_map.csv")

#########################################################
# Anova
#########################################################

# Using indexing
dta_subset <- dta.main[, 1:4]
dta.anl <- dta_subset[, 1:3]
dta.anl$Environment <- "Düsseldorf"
head(dta.anl, n=14)

length(unique(dta.anl$Genotype))
length(unique(dta.anl$GroupReplicate))
length(unique(dta.anl$Environment))
nrow(dta.anl)

anova.model <- aov(Height ~ Genotype, data = dta.anl)
resid.model <- anova.model$residuals

length(resid.model)

ggplot(data = data.frame(1:length(resid.model), resid.model), 
       aes(x = resid.model)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.model), 
                            sd = sd(resid.model)),
                color="darkred")

# plot qq-plot
qqPlot(resid.model, envelope = F)

plot(resid.model)

# removing outliers
dta.anl.clean <- dta.anl[c(-40,-14,-182,-126,-11,-94,-160,-211,-169,-150),]

anova.model.clean <- aov(Height ~ Genotype, data = dta.anl.clean)
resid.model.clean <- anova.model.clean$residuals

length(resid.model)

ggplot(data = data.frame(1:length(resid.model.clean), resid.model.clean), 
       aes(x = resid.model.clean)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.model.clean), 
                            sd = sd(resid.model.clean)),
                color="darkred")

# plot qq-plot
qqPlot(resid.model.clean, envelope = F)

plot(resid.model.clean)

# final anova with removed outliers
anova.model.final <- aov(Height ~ Genotype, data = dta.anl.clean)
summary(anova.model.final)


anova.model <- aov(Height ~ Genotype, data = dta.anl)
summary(anova.model)

##############################################################
# QTL - Analysis
##############################################################

# loop through all markers and conduct ANOVAs for each one of them
lis.out <- list()
for(m in 5:ncol(dta.main)){
  # extract the first two columns (group + spike length) and the m-th column (=marker)
  loop.dta <- dta.main[,c(2,3,m)]
  # save the name of the current marker
  loop.marker.name <- colnames(loop.dta)[3]
  # in the model in the aov command, we will have to define the independent 
  # variable name. since the marker name changes in each loop, we need to 
  # to change the column names here to have the same marker name in each loop
  colnames(loop.dta) <- c("group","trait", "allele")
  # conduct one-way ANOVA
  loop.aov <- aov(trait ~ allele + group, loop.dta)
  # extract the allele's p-value
  loop.pval.allele <- summary(loop.aov)[[1]][1,5]
  # extract the marker's genetic position from the genetic map
  loop.map <- map[which(map$marker == loop.marker.name),]
  # create the output data frame
  loop.df.out <- data.frame(loop.map, pval.allele=loop.pval.allele)
  # save the output data frame in the list
  lis.out[[m]] <- loop.df.out
}
# combine the loop's output data frames into a single data frame
res.aov <- do.call("rbind", lis.out)

# have a look at the output data frame
head(res.aov)
dim(res.aov)

# create a data frame that follows the requirements of the 
# manhattan command of the qqman library (run as is)
dta.plot <- data.frame(CHR = as.numeric(gsub("H","",res.aov$chr)),
                       BP = res.aov$pos,
                       SNP = 1:nrow(res.aov),
                       pval.marker = res.aov$pval.allele)

# calculate the negative logarithm of the Bonferroni corrected significance threshold
# choosing the threshold 0.0001 to consider only the one significant peak as visual in the manhatten plot
sig.threshold.BonfCorrected <- -log10(0.0001/nrow(dta.plot))

# plot genome-wide p-values
# for markers
manhattan(dta.plot, genomewideline = sig.threshold.BonfCorrected,
          suggestiveline = F, logp=T, p="pval.marker", type="l", 
          lwd=3, ylab="-log10(p-values)", main="Marker effect")

res.aov

# extracting only significant markers due to chosen threshold 0.0001
sig.marker <- res.aov[res.aov$pval.allele < (0.0001/nrow(dta.plot)), ]
sig.marker
dim(sig.marker)



#####################################################
#Linkage mapping
#####################################################

# read in linkage map and qualitative phenotypes as denoted in marker annotation
owb <- read.cross("csv", ".", "owb_linkage_map_qualt_phenotypes_WS2324.csv", genotypes=c("a","b"),
                  alleles=c("a", "b"), crosstype="dh")
# complete linkage map
plotMap(owb, main="", show.marker.names=T)

# To map the locus, calculate all pairwise recombination frequencies and LOD scores
sixrowed <- tryallpositions(owb, "2vs6row", error.prob=0)

# showing best linkage to a marker from each chromosome
summary(sixrowed)
plot(sixrowed)

# moving locus to best marker position
owb <- movemarker(owb, "2vs6row", "2H", 94.8)
#update map
plotMap(owb, main="", show.marker.names=T)

#Is there any known locus/gene described in literature for your trait and chromosome location?


