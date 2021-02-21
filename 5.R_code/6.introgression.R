#----------------------------
#
# Maerl Genomics Project
#
# Description:
# Perform PCA and introgression analysis
# P. calcareum vs. P. purpureum vs. L. corallioides
#
# https://github.com/Tom-Jenkins/maerl_genomics
#
#----------------------------

# Load packages
library(vcfR)
library(adegenet)
library(poppr)
library(hierfstat)
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(tidyverse)


# ----------------- #
#
# Data import and preparation
#
# ----------------- #

# Convert vcf file to genind object
vcf = read.vcfR("Input_data/introgression.filt.recode.vcf")
intro.snps = vcfR2genind(vcf)
intro.snps
rm(vcf)

# Individual labels
indNames(intro.snps)

# Add population labels
intro.snps$pop = factor(str_extract(indNames(intro.snps), "[aA-zZ]+"))
summary(intro.snps$pop)
intro.snps


# ----------------- #
#
# Missing data check
#
# ----------------- #

# Missing data for loci
locmiss = propTyped(intro.snps, by = "loc")
summary(locmiss)

# Missing data for individuals
indmiss = propTyped(intro.snps, by = "ind")

# Plot using barplot
barplot(indmiss, ylim = c(0,1), ylab = "Complete genotypes (proportion)", las = 2)
indmiss[ which(indmiss < 0.90) ] # less than 90% complete genotypes

# Remove individuals with missing data > X %
intro.snps = missingno(intro.snps, type = "geno", cutoff = 0.20)
intro.snps


# ----------------- #
#
# Polymorphic check
#
# ----------------- #

# Check for monomorphic loci
isPoly(intro.snps) %>% summary()

# Only keep polymorphic loci
poly_loci = names(which(isPoly(intro.snps) == TRUE))
intro.snps = intro.snps[loc = poly_loci]
intro.snps

# How many SNPs had no missing data?
missingno(intro.snps, type = "loci", cutoff = 0)

# Only keep SNPs were no missing data
# intro.snps = missingno(intro.snps, type = "loci", cutoff = 0)
# intro.snps


# ----------------- #
#
# PCA
#
# ----------------- #

# Replace missing data with the mean allele frequencies
x = tab(intro.snps, NA.method = "mean")

# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

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
ind_coords$Ind = indNames(intro.snps)

# Add a column with the population IDs
ind_coords$Pop = intro.snps$pop
head(ind_coords)

# Convert Sco to Ppur
ind_coords$Ind = gsub("Sco", "Ppur", ind_coords$Ind)
ind_coords$Pop = gsub("Sco", "Ppur", ind_coords$Pop)

# Colour definitions
levels(ind_coords$Pop)
# blue "#377eb8"
# orange "#ff7f00"
# blue-green "#009E73"
# purple "#984ea3"
# red "#e41a1c"
# dark red "#a50f15"
cols = data.frame("Fal" = "#ff7f00",
                  "Lcor" = "#d9d9d9",
                  "Man" = "#377eb8",
                  "Mil" = "#377eb8",
                  "Mor" = "#984ea3",
                  "Nor" = "yellow3",
                  "Roc" = "white",
                  "Ppur" = "#737373",
                  "Tre" = "#009E73",
                  "Zar" = "#377eb8",
                  check.names = FALSE)

# Change factor order
pop_order = c("Nor","Zar","Mil","Man","Fal","Mor","Tre","Roc","Ppur","Lcor")
ind_coords$Pop = factor(ind_coords$Pop,
                            levels = pop_order) 

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom ggplot2 theme
ggtheme = theme(legend.title = element_text(size = 13, face = "bold"),
                axis.text.y = element_text(colour="black", size=14),
                axis.text.x = element_text(colour="black", size=14),
                axis.title = element_text(colour="black", size=14),
                legend.position = "top",
                legend.text = element_text(size=13),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(0.7, "cm"),
                legend.box.spacing = unit(0, "cm"),
                plot.tag = element_text(size=14, face = "bold"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                # title centered
                plot.title = element_text(hjust=0.5, size=15) 
)

# Add species labels
species_labels = data.frame(labs = c("P. calcareum", "P. purpureum", "L. corallioides"),
                            X = c(-9.5, 25, 12),
                            Y = c( -2, -2, 44)
)
          
# Scatter plot axis 1 vs. 2
pca1 = ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # points
  geom_point(aes(fill = Pop), shape = 21, size = 4, alpha = 0.9, show.legend = TRUE)+
  # colouring
  scale_fill_manual("Site", values = cols)+
  # species labels
  geom_text(data = species_labels, aes(x = X, y = Y, label = labs), fontface = "italic",
            size = 5)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # title
  # ggtitle("Principal components analysis (1 vs. 2)")+
  # custom labels
  labs(x = xlab, y = ylab, tag = "(a)")+
  # custom theme
  ggtheme+
  guides(fill = guide_legend(nrow = 2, title.position = "left", title.hjust = 0.5))
pca1


# ----------------- #
#
# Snapclust introgression test
#
# ----------------- #

# Test with L. corallioides
lcor.intro.snps = popsub(intro.snps, blacklist = "Sco")
lcor.results = snapclust(lcor.intro.snps, k = 2, hybrids = TRUE)

# Plot with compoplot
compoplot(lcor.results, col.pal = hybridpal(), n.col = 2,
          lab = indNames(lcor.intro.snps), show.lab = TRUE, legend = TRUE)

# Prepare data.frame
lcor.df = lcor.results$proba %>% melt()
colnames(lcor.df) = c("ID","Group","Value")
colnames(lcor.df)

# Define levels of Group variable
lcor.df$Group = factor(lcor.df$Group, levels = c("A","B","0.5_A-0.5_B"))

# Plot results
lcor.bar = ggplot(data=lcor.df, aes(x=ID, y=Value, fill=Group))+
  geom_bar(stat="identity", width = 1, size = 0.25, show.legend = TRUE, colour="black")+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c("pink","#d9d9d9","black"),
                    labels = c(expression(italic("P. calcareum")),
                               expression(italic("L. corallioides")),
                               "Hybrid"))+
  ylab("Membership probability")+
  labs(tag = "(c)")+
  # xlab("Individual")+
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 90, colour = "black", size = 12, vjust = 0.5),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.tag = element_text(size=14, face = "bold"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 14))
lcor.bar

# Test with Scotland samples
sco.intro.snps = popsub(intro.snps, blacklist = "Lcor")
sco.results = snapclust(sco.intro.snps, k = 2, hybrids = TRUE)

# Plot with compoplot
compoplot(sco.results, col.pal = hybridpal(), n.col = 2,
          lab = indNames(sco.intro.snps), show.lab = TRUE, legend = TRUE)

# Prepare data.frame
ppur.df = sco.results$proba %>% melt()
colnames(ppur.df) = c("ID","Group","Value")
colnames(ppur.df)

# Define levels of Group variable
ppur.df$Group = factor(ppur.df$Group, levels = c("A","B","0.5_A-0.5_B"))

# Convert Sco to Ppur
ppur.df$ID = gsub("Sco", "Ppur", ppur.df$ID)
ppur.df$ID

# Plot results
ppur.bar = ggplot(data=ppur.df, aes(x=ID, y=Value, fill=Group))+
  geom_bar(stat="identity", width=1, size=0.25, show.legend = TRUE, colour="black")+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c("pink","#737373","black"),
                    labels = c(expression(italic("P. calcareum")),
                               expression(italic("P. purpureum")),
                               "Hybrid"))+
  ylab("Membership probability")+
  labs(tag = "(b)")+
  # xlab("Individual")+
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 90, colour = "black", size = 12, vjust = 0.5),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.tag = element_text(size=14, face = "bold"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 14))
ppur.bar

# Figure 4
library(ggpubr)
barplots = ggarrange(ppur.bar, lcor.bar, nrow = 2)
barplots
fig4 = ggarrange(pca1, barplots, ncol = 2)
fig4
ggsave(plot = fig4, "./Figures/Figure4.png", width = 15, height = 7, dpi = 600)
ggsave(plot = fig4, "./Figures/Figure4.pdf", width = 15, height = 7)


# ----------------- #
#
# D statistic introgression test
#
# ----------------- #

# D3 forumla from Hahn & Mark 2019 S Hibbins
# https://doi.org/10.1093/molbev/msz178
D3 = function(distAC, distBC){
  (distBC - distAC) / (distBC + distAC)
}

# A = Lcor
# B = Ppur
# C = Pcal

# Convert genind to data.frame objects
AC = popsub(intro.snps, sublist = c("Lcor","Fal","Man","Mor","Tre","Zar")) %>% genind2hierfstat
BC = popsub(intro.snps, sublist = c("Sco","Fal","Man","Mor","Tre","Zar")) %>% genind2hierfstat

# Change site ID to species ID
library(mgsub)
AC$pop = mgsub(as.character(AC$pop), c("Fal","Man","Mor","Tre","Zar"), rep("Pcal", times = 5))
BC$pop = mgsub(as.character(BC$pop), c("Fal","Man","Mor","Tre","Zar"), rep("Pcal", times = 5))

# Calculate genetic distances
dist_AC = genet.dist(AC, method = "WC84")
dist_BC = genet.dist(BC, method = "WC84")

# Calculate D statistic
D3(dist_AC, dist_BC)
