# --------------------------- #
#
# Maerl Genomics Project
#
# Description:
# Heterozygosity, F statistics and PCA using filtered SNP data set
#
# https://github.com/Tom-Jenkins/maerl_genomics
# 
# --------------------------- #

# Load packages
library(adegenet)
library(poppr)
library(hierfstat)
library(ggpubr)
library(tidyverse)

# Load genind file
load("Output_rdata_files/maerl_SNPs.RData")
maerl
summary(maerl$pop)


# ----------------- #
#
# Visualise heterozygosity and FIS
#
# ----------------- #

# Calculate the number of private alleles per site
private_alleles(maerl) %>% apply(MARGIN = 1, FUN = sum)

# Calculate mean allelic richness per site
allelic.richness(genind2hierfstat(maerl))$Ar %>%
  apply(MARGIN = 2, FUN = mean) %>% 
  round(digits = 3)

# Calculate basic stats using hierfstat
# basic_stats = basic.stats(maerl, diploid = TRUE, digits = 3)
# save(basic_stats, file = "Output_rdata_files/hierfstat_basicstats.RData")
load("Output_rdata_files/hierfstat_basicstats.RData")

# Create tibble containing mean heterozygosity per SNP per site
head(basic_stats$Ho)
het.site = tibble(SNP_ID = rownames(basic_stats$Ho),
                  # Arm = basic_stats$Ho[, "Arm"],
                  Bor = basic_stats$Ho[, "Bor"],
                  Fal = basic_stats$Ho[, "Fal"],
                  Man = basic_stats$Ho[, "Man"],
                  Mor = basic_stats$Ho[, "Mor"],
                  # Nor = basic_stats$Ho[, "Nor"],
                  Ons = basic_stats$Ho[, "Ons"],
                  # Roc = basic_stats$Ho[, "Roc"],
                  Tre = basic_stats$Ho[, "Tre"],
                  Zar = basic_stats$Ho[, "Zar"]
)
het.site

# Mean and median observed heterozygosity per site
dplyr::select(het.site, -SNP_ID) %>% apply(MARGIN = 2, FUN = mean, na.rm = TRUE) %>% round(digits = 2)
dplyr::select(het.site, -SNP_ID) %>% apply(MARGIN = 2, FUN = median, na.rm = TRUE) %>% round(digits = 2)

# Convert tibble to long format
het.site = het.site %>%
  pivot_longer(cols = 2:ncol(het.site), names_to = "Site", values_to = "Hobs")
het.site

# Plot histogram for each site
het.site.plt = ggplot(het.site, aes(x = Hobs))+
  geom_histogram(fill = "grey40", colour = "black", binwidth = 0.15)+
  facet_wrap(~Site)+
  scale_y_continuous(limits = c(0, 10000))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.25))+
  xlab("Heterozygosity")+
  ylab("Frequency")+
  ggtitle("Mean observed heterozygosity per SNP per site")+
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    panel.grid.minor = element_blank()
  )
het.site.plt

# Create tibble containing mean FIS per SNP per site
head(basic_stats$Fis)
fis.site = tibble(SNP_ID = rownames(basic_stats$Fis),
                  Bor = basic_stats$Fis[, "Bor"],
                  Fal = basic_stats$Fis[, "Fal"],
                  Man = basic_stats$Fis[, "Man"],
                  Mor = basic_stats$Fis[, "Mor"],
                  Ons = basic_stats$Fis[, "Ons"],
                  Tre = basic_stats$Fis[, "Tre"],
                  Zar = basic_stats$Fis[, "Zar"]
)
fis.site

# Mean and median FIS per site
dplyr::select(fis.site, -SNP_ID) %>% apply(MARGIN = 2, FUN = mean, na.rm = TRUE) %>% round(digits = 2)
dplyr::select(fis.site, -SNP_ID) %>% apply(MARGIN = 2, FUN = median, na.rm = TRUE) %>% round(digits = 2)

# Convert tibble to long format
fis.site = fis.site %>%
  pivot_longer(cols = 2:ncol(fis.site), names_to = "Site", values_to = "FIS") %>% 
  drop_na
fis.site

# Plot histogram for each site
fis.site.plt = ggplot(fis.site, aes(x = FIS))+
  geom_histogram(binwidth = 0.20, fill = "grey40", colour = "black")+
  facet_wrap(~Site)+
  scale_y_continuous(limits = c(0, 5000))+
  xlab(expression(italic("F")[IS]))+
  ylab("Frequency")+
  ggtitle(expression("Mean" ~ italic("F")[IS] ~ "per SNP per site"))+
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    panel.grid.minor = element_blank()
  )
fis.site.plt

# Figure S2
figS2 = ggarrange(het.site.plt + labs(tag = "(a)"),
                  fis.site.plt + labs(tag = "(b)"),
                  nrow = 2)
figS2
# ggsave(plot = figS2, filename = "Figures/FigureS2.pdf", width = 8, height = 10)
# ggsave(plot = figS2, filename = "Figures/FigureS2.png", width = 8, height = 10, dpi = 600)


# ----------------- #
#
# Nuclear SNP PCA 
#
# ----------------- #

# Replace missing data with the mean allele frequencies
x = tab(maerl, NA.method = "mean")

# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, center = TRUE, nf = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
percent
barplot(percent, ylab = "Percent of genetic variance explained by eigenvectors",
        names.arg = round(percent, 2))

# Create a dataframe containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(maerl)

# Add a column with the population IDs
ind_coords$Pop = maerl$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Pop,
                     data = ind_coords,
                     FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Pop", suffix = c("",".cen"))
head(ind_coords)

# Colour definitions
levels(ind_coords$Pop)
# blue "#377eb8"
# orange "#ff7f00"
# blue-green "#009E73"
# purple "#984ea3"
# red "#e41a1c"
# dark red "#a50f15"
cols = data.frame("Bor" = "#e41a1c",
                  "Fal" = "#ff7f00",
                  "Man" = "#377eb8",
                  "Mor" = "#984ea3",
                  "Nor" = "yellow3",
                  "Ons" = "#a50f15",
                  "Roc" = "white",
                  "Arm" = "deeppink",
                  "Tre" = "#009E73",
                  "Zar" = "#377eb8",
                  check.names = FALSE)

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom ggplot2 theme
ggtheme = theme(legend.title = element_blank(),
                axis.text.y = element_text(colour="black", size=14),
                axis.text.x = element_text(colour="black", size=14),
                axis.title = element_text(colour="black", size=14),
                legend.position = "right",
                legend.text = element_text(size=15),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(0.7, "cm"),
                legend.box.spacing = unit(0, "cm"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                # title centered
                plot.title = element_text(hjust=0.5, size=15) 
)


# Scatter plot axis 1 vs. 2 [site labels only]
pca.plt1 = ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # labels
  # geom_label(aes(label = Ind, fill = Pop), size = 4, show.legend = FALSE)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Pop), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Pop), shape = 21, size = 4,  show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Pop, fill = Pop), size = 5, alpha = 0.9, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # title
  ggtitle("Principal components analysis (1 vs. 2)")+
  # custom theme
  ggtheme
pca.plt1

# Scatter plot axis 1 vs. 3 [site labels only]
ylab.axis3 = paste("Axis 3 (", format(round(percent[3], 1), nsmall=1)," %)", sep="")
pca.plt2 = ggplot(data = ind_coords, aes(x = Axis1, y = Axis3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # labels
  # geom_label(aes(label = Ind, fill = Pop), size = 4, show.legend = FALSE)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis3.cen, colour = Pop), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Pop), shape = 21, size = 4,  show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Pop, fill = Pop), size = 5, alpha = 0.9, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab.axis3)+
  # title
  ggtitle("Principal components analysis (1 vs. 3)")+
  # custom theme
  ggtheme
pca.plt2

# Figure 1
fig1 = ggarrange(pca.plt1, pca.plt2, labels = c("(a)","(b)"), ncol = 2)
fig1
# ggsave(plot = fig1, "./Figures/Figure1.png", width = 15, height = 7, dpi = 600)
# ggsave(plot = fig1, "./Figures/Figure1.pdf", width = 15, height = 7)



# ----------------- #
#
# Nuclear SNP WC Fst 
#
# ----------------- #

# Remove Arm, Nor and Roc for computing Fst
maerlfst = popsub(maerl, blacklist = c("Arm","Nor","Roc"))

# Compute Weir and Cockerham pairwise Fst
# fst = genet.dist(maerlfst, method = "WC84")
# save(fst, file = "Output_rdata_files/pairwiseWCfst.RData")
load("Output_rdata_files/pairwiseWCfst.RData")

# Convert dist to dataframe
mat.pairwise = as.matrix(fst)
ind = which( upper.tri(mat.pairwise), arr.ind = TRUE)
df.pairwise = data.frame( Site1 = dimnames(mat.pairwise)[[2]][ind[,2]],
                          Site2 = dimnames(mat.pairwise)[[1]][ind[,1]],
                          Fst = mat.pairwise[ ind ])

# Fst italic label
fstlab = expression(italic("F")[ST])

# Create heatmap of Fsts
Fst = ggplot(data = df.pairwise, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  scale_fill_gradient(name = fstlab, low = "white", high = "red")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.8,
                               title.position = "top", title.hjust = 0.5))+
  ggtitle(expression("Pairwise Weir & Cockerham (1984)" ~ italic("F")[ST]))+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=10, face = "bold"),
        axis.text.y = element_text(colour="black", size=10, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust=0.5, size=18, face="bold"))
Fst
# ggsave(plot = Fst, "Figures/FigureS3.png", width=8, height=7, dpi=600)
# ggsave(plot = Fst, "Figures/FigureS3.pdf", width=8, height=7)


# ----------------- #
#
# Fal Estuary allelic differences
#
# ----------------- #

# Subset Fal Estuary
fal = popsub(maerl, sublist = "Fal")

# Remove all missing data
fal = missingno(fal, type = "loci", cutoff = 0)
fal

# Check loci still polymorphic after filtering
isPoly(fal) %>% summary

# Only keep polymorphic loci
poly_loci = names(which(isPoly(fal) == TRUE))
fal = fal[loc = poly_loci]
fal

# Calculate number of allelic differences
fal.mat = poppr::diss.dist(fal, percent = FALSE, mat = TRUE)

# Convert to data.frame
ind = which( upper.tri(fal.mat), arr.ind = TRUE)
fal.df = data.frame( Site1 = dimnames(fal.mat)[[2]][ind[,2]],
                     Site2 = dimnames(fal.mat)[[1]][ind[,1]],
                     Allelic_diff = fal.mat[ ind ])
mean(fal.df$Allelic_diff)
min(fal.df$Allelic_diff)
max(fal.df$Allelic_diff)
