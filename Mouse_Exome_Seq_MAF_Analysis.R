### Mouse Exome Seq MAF Analysis ###
rm(list=ls())
setwd("/Users/benstansfield/Documents/Research/PhD/Ooi/Projects/DNA_Damage/Mouse_Exome_Seq/Data/")

# load libraries
library(maftools)
library(BSgenome.Mmusculus.UCSC.mm39)
library(NMF)

MMCID_26B <- read.maf(maf = "MMCID-26B.maf")
MMCID_30B <- read.maf(maf = "MMCID-30B.maf")
TMCK1_14B <- read.maf(maf = "TCMK1-14B.maf")
TMCK1_23B <- read.maf(maf = "TMCK1-23B.maf")

MMCID_maf <- merge_mafs(maf = list(MMCID_26B, MMCID_30B))
TMCK1_maf <- merge_mafs(maf = list(TMCK1_14B, TMCK1_23B))

pdf(file = "MMCID_Summary.pdf")
plotmafSummary(maf = MMCID_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf(file = "TMCK1_Summary.pdf")
plotmafSummary(maf = TMCK1_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


## adjust chromosome names so they match with the reference genome
MMCID_maf@data$Chromosome <- paste0("chr", MMCID_maf@data$Chromosome)
MMCID_maf@maf.silent$Chromosome <- paste0("chr", MMCID_maf@maf.silent$Chromosome)

TMCK1_maf@data$Chromosome <- paste0("chr", TMCK1_maf@data$Chromosome)
TMCK1_maf@maf.silent$Chromosome <- paste0("chr", TMCK1_maf@maf.silent$Chromosome)

MMCID_TNM <- trinucleotideMatrix(MMCID_maf, ref_genome = "BSgenome.Mmusculus.UCSC.mm39")
TMCK1_TNM <- trinucleotideMatrix(TMCK1_maf, ref_genome = "BSgenome.Mmusculus.UCSC.mm39")

# combine data from each sample and combine into total number of mutations
MMCID_TNM_combined <- (MMCID_TNM$nmf_matrix[1,] + MMCID_TNM$nmf_matrix[2,])
TMCK1_TNM_combined <- (TMCK1_TNM$nmf_matrix[1,] + TMCK1_TNM$nmf_matrix[2,])

# convert values to a percentage of total mutations
mmcid_total_mutations <- sum(MMCID_TNM_combined)
MMCID_TNM_percentages <- (MMCID_TNM_combined/mmcid_total_mutations)*100

tmck1_total_mutations <- sum(TMCK1_TNM_combined)
TMCK1_TNM_percentages <- (TMCK1_TNM_combined/tmck1_total_mutations)*100

colors <- c(rep("skyblue", 16), rep("black", 16), rep("red", 16), rep("grey", 16), rep("green", 16), rep("pink", 16))

# Create the bar plot
pdf(file = "MMCID_Mutational_Signature.pdf", height = 5, width = 10)
par(mar = c(9, 5, 4, 1) + 0.1, cex.axis = 0.25, mgp = c(3, 0.1, 0))
barplot(height = MMCID_TNM_percentages, col = colors, las = 2, border = NA)
mtext("% of base substitutions", side = 2, line = 3)
dev.off()

pdf(file = "TMCK1_Mutational_Signature.pdf", height = 5, width = 10)
par(mar = c(9, 5, 4, 1) + 0.1, cex.axis = 0.25, mgp = c(3, 0.1, 0))
barplot(height = TMCK1_TNM_percentages, col = colors, las = 2, border = NA)
mtext("% of base substitutions", side = 2, line = 3)
dev.off()