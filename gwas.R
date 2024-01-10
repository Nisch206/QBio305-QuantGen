

setwd("C:/Users/aless/OneDrive/Desktop/Project_pop&qGen/gwas")

library(qqman)

# Read GWAS output file created using gemma
gwas <- read.table("gwas_project_result.assoc.txt", header = TRUE)
head(gwas)

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas, chr="chr", bp="ps", snp="rs", p="p_wald" )

# modify colour of points
manhattan(x = gwas, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction)
pval_bonf = 0.05/dim(gwas)[[1]]


#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


