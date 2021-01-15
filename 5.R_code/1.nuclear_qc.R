# --------------------------- #
#
# Maerl Genomics Project
#
# Description:
# P.calcareum nuclear SNP QC analysis
# Explore and filter SNPs from VCF file
#
# https://github.com/Tom-Jenkins/maerl_genomics
#
# --------------------------- #

# Load packages
library(vcfR)
library(adegenet)
library(poppr)
library(hierfstat)
library(reshape2)
library(ggplot)
library(ggplotify)
library(ggpubr)
library(tidyverse)


# ----------------- #
#
# Data import and preparation
#
# ----------------- #

# Import vcf file
vcf = read.vcfR("Input_data/filter3.recode.vcf")

# Convert VCF to genind
maerl = vcfR2genind(vcf)
maerl

# Add population labels
indNames(maerl)
maerl$pop = factor(str_extract(indNames(maerl), "[aA-zZ]+"))
maerl$pop
maerl
summary(maerl$pop)


# ----------------- #
#
# Missing data check
#
# ----------------- #

# Missing data for loci
propTyped(maerl, by = "loc") %>% summary()

# Remove loci with missing data > X %
# maerl = missingno(maerl, type = "loci", cutoff = 0.05)
# maerl

# Missing data for individuals
indmiss = propTyped(maerl, by = "ind")

# Barplot
barplot(indmiss, ylim = c(0,1), ylab = "Complete genotypes (proportion)", las = 2)
indmiss[ which(indmiss < 0.90) ] # less than X % complete genotypes

# Remove individuals with missing data > X %
maerl = missingno(maerl, type = "geno", cutoff = 0.20)
maerl

# Check loci still polymorphic after filtering
isPoly(maerl) %>% summary

# Only keep polymorphic loci
poly_loci = names(which(isPoly(maerl) == TRUE))
maerl = maerl[loc = poly_loci]
maerl


# ----------------- #
#
# Heterozygosity
#
# ----------------- #

# Plot observed heterozygosity per SNP locus
maerl_summary = summary(maerl)
plot(maerl_summary$Hobs)

# Plot Hobs versus Hexp
plot(maerl_summary$Hobs, maerl_summary$Hexp, xlab = "Hobs", ylab = "Hexp")

# Test where H0 is Hexp = Hobs
bartlett.test(list(maerl_summary$Hexp, maerl_summary$Hobs))

# Keep loci with heterozygosity <= 0.60
loci_to_keep = names(maerl_summary$Hobs[ which(maerl_summary$Hobs <= 0.60) ])
maerl = maerl[loc=loci_to_keep]
maerl


# ----------------- #
#
# Clone analysis
#
# ----------------- #

# Calculate Prevosti's genetic distance
ind.dist = prevosti.dist(maerl)

# Plot heatmap of results
library(lattice)
heatmap1 = levelplot(as.matrix(ind.dist), scales = list(x=list(rot=90)),
                     xlab = "Site 1", ylab = "Site 2",
                     main = "Heatmap of pairwise Prevosti genetic distances")
heatmap1
heatmap1 = as.grob(heatmap1)

# Convert dist to dataframe
library(spaa)
dist.df = dist2list(ind.dist)

# Remove same-site comparisons
dist.df = dist.df[ which(dist.df$col != dist.df$row), ]

# Remove duplicated comparisons
dist.df = distinct(dist.df, value, .keep_all = TRUE)

# Comparisons with lowest pairwise distances
dist.df[order(dist.df$value), ] %>% head(n = 10)

# Histogram ggtheme
hist_theme = theme(
  axis.text = element_text(size = 11),
  axis.title = element_text(size = 12),
  plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
)

# Plot histogram of distances
hist1 = ggplot(dist.df, aes(x = value))+
  geom_histogram(fill = "grey60", colour = "black", binwidth = 0.005)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  xlab("Prevosti genetic distances")+
  ylab("Frequency")+
  ggtitle("Histogram of pairwise Prevosti genetic distances")+
  hist_theme
hist1

# Figure S1
figS1 = ggarrange(heatmap1,
                  hist1,
                  ncol = 2)
figS1
ggsave(plot = figS1, filename = "Figures/FigureS1.pdf", width = 16, height = 7)
ggsave(plot = figS1, filename = "Figures/FigureS1.png", width = 16, height = 7, dpi = 600)


# ----------------- #
#
# OutFLANK outlier test
#
# ----------------- #

# Remove Arm, Nor, and Roc sites from data set (low sample sizes)
maerl.outflank = popsub(maerl, blacklist = c("Arm","Nor","Roc"))
summary(maerl.outflank$pop)

# Check loci are still polymorphic
isPoly(maerl.outflank) %>% summary()

# Only keep polymorphic loci
poly_loci = names(which(isPoly(maerl.outflank) == TRUE))
maerl.outflank = maerl.outflank[loc = poly_loci]
maerl.outflank

# Run OutFLANK using dartR wrapper script
library(dartR)
outflnk = gl.outflank(maerl.outflank, qthreshold = 0.05)

# OutFLANK results
outflnk.df = outflnk$outflank$results

# Print number of outliers (TRUE)
outflnk.df$OutlierFlag %>% summary

# Remove duplicated rows for each SNP locus
rowsToRemove = seq(1, nrow(outflnk.df), by = 2)
outflnk.df = outflnk.df[-rowsToRemove, ]

# Convert Fsts <0 to zero
outflnk.df$FST[outflnk.df$FST < 0] = 0 
min(outflnk.df$FST, na.rm = TRUE)

# Italic labels
fstlab = expression(italic("F")[ST])
hetlab = expression(italic("H")[e])

# Plot He versus Fst
ggplot(data = outflnk.df)+
  geom_point(aes(x = He, y = FST))+
  ggtitle("OutFLANK outlier test")+
  xlab(hetlab)+
  ylab(fstlab)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))


# ----------------- #
#
# Export filtered data set
#
# ----------------- #

# Export filtered data set in Genind format
maerl
save(maerl, file = "Output_rdata_files/maerl_SNPs.RData")

# Export file in geno format for SNMF analysis
load("Output_rdata_files/maerl_SNPs.RData")
maerl
summary(maerl$pop)
maerl.snmf = popsub(maerl, blacklist = c("Arm","Nor","Roc")) # remove - low sample size
summary(maerl.snmf$pop) 
isPoly(maerl.snmf) %>% summary # Only keep polymorphic loci
maerl.snmf = maerl.snmf[loc = names(which(isPoly(maerl.snmf) == TRUE))] 
maerl.snmf

# Convert to matrix and encode genotypes as characters
maerl.geno = genind2df(maerl.snmf, usepop = FALSE)
maerl.geno = apply(maerl.geno, MARGIN = 2, FUN = as.character)

# Convert 00 to 2 (00 means two copies of the reference allele)
maerl.geno = gsub("00", "2", maerl.geno)

# Convert 11 to 0 (11 means zero copies of the reference allele)
maerl.geno = gsub("11", "0", maerl.geno)

# Convert 01 to 1 (01 or 10 means one copy of the reference allele)
maerl.geno = gsub("01", "1", maerl.geno)
maerl.geno = gsub("10", "1", maerl.geno)

# Convert NA to 9
maerl.geno[is.na(maerl.geno)] = 9

# Convert back to numeric
maerl.geno = apply(maerl.geno, MARGIN = 2, FUN = as.numeric)

# Export geno file
library(LEA)
write.geno(maerl.geno, output.file = "SNMF/maerl_SNPs.geno")

