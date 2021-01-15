# --------------------------- #
#
# Maerl Genomics Project
#
# Description:
# Tree building analysis
#
# Data:
# 1. Mitogenomes FASTA CDS alignment
# 2. Plastomes FASTA CDS alignment
#
# https://github.com/Tom-Jenkins/maerl_genomics
#
# --------------------------- #

# Load packages
library(ape)
library(adegenet)
library(poppr)
library(ggplot2)


# ----------------- #
#
# Mitogenomes tree
#
# ----------------- #

# Import mitogenomes alignment
mito_aln = ape::read.FASTA("Input_data/Pcalcareum_mito24CDS_mafft.fasta", type = "DNA")
names(mito_aln)
names(mito_aln) = substr(names(mito_aln), 1, 5)
names(mito_aln)
mito_aln

# Compute Tamura and Nei 1993 genetic distances
mito_dist = ape::dist.dna(mito_aln, model = "TN93")

# Hierarchical clustering ("average" = UPGMA)
mito_hclust = hclust(mito_dist, method = "average")

# Plot quick dendrogram
plot(mito_hclust)

# Plot quick phylogram using ape
plot.phylo(as.phylo(mito_hclust), type = "phylogram")

# Extract dendrogram data from hclust results
library(ggdendro)
mito_dend = dendro_data(mito_hclust)
names(mito_dend)
head(mito_dend$segments)
head(mito_dend$labels)

# Create dataframe individuals and site labels
mito_df = data.frame(Ind = mito_dend$labels$label,
                     Site = substr(mito_dend$labels$label, 1, 3))
head(mito_df)
unique(mito_df$Site) ; length(unique(mito_df$Site))

# Make a colour palette
library(scales)
show_col(brewer.pal(11, "RdYlBu"))
cols = brewer.pal(11, "RdYlBu") ; cols

# Function to add a specific colour to dataframe depending on site
# red "#e41a1c"
# blue "#377eb8"
# purple "#984ea3"
# orange "#ff7f00"
# green "#4daf4a"
addcolours = function(x){
  # If pop label is present function will output the colour
  if(x=="Nor") y = "yellow3"
  if(x=="Zar") y = "#377eb8"
  if(x=="Mil") y = "#377eb8"
  if(x=="Man") y = "#377eb8"
  if(x=="Fal") y = "#ff7f00"
  if(x=="Mor") y = "#984ea3"
  if(x=="Tre") y = "#4daf4a"
  if(x=="Roc") y = "black"
  if(x=="Bor") y = "#e41a1c"
  if(x=="Ons") y = "#e41a1c"
  if(x=="Arm") y = "deeppink"
  return(y)
}
mito_df$tree_cols = sapply(mito_df$Site, addcolours)
head(mito_df)

# Italic label
mitolab = expression(italic("Phymatolithon calcareum") ~ "mitogenome CDS alignment")

# Plot custom dendrogram using ggplot2
tree1 = ggplot(mito_dend$segments, horiz = TRUE)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.8)+
  # geom_text(data = mito_dend$labels, aes(x = x, y = y, label = label),
  #           size = 5, fontface = "bold", hjust = -0.25, col = mito_df$tree_cols)+
  # coord_flip()+
  # scale_y_reverse(limits = c(max(mito_dend$segments$y)+0.00005,
  #                            min(mito_dend$segments$y)-0.00005
  #                            )
  #                 )+
  coord_cartesian(clip = "off")+
  geom_text(data = mito_dend$labels, aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", hjust = 1.2, vjust = 0.5, col = mito_df$tree_cols,
            angle = 90)+
  labs(title = mitolab,
       subtitle = "(32 variant sites)",
       tag = "(a)")+
  theme_dendro()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5,),
        plot.tag = element_text(size=18, face = "bold"),
        plot.margin = margin(t = 0, r = 0, b = 1, l = 0, unit = "cm"))
tree1
# ggsave("dendrogram_mitogenomes.png", width = 15, height = 10, dpi = 600)


# ----------------- #
#
# Plastomes tree
#
# ----------------- #

# Import plastomes alignment
plastid_aln = ape::read.FASTA("Input_data/Pcalcareum_plastid212CDS_mafft.fasta", type = "DNA")
names(plastid_aln)
names(plastid_aln) = substr(names(plastid_aln), 1, 5)
names(plastid_aln)
plastid_aln

# Compute Tamura and Nei 1993 genetic distances
plastid_dist = ape::dist.dna(plastid_aln, model = "TN93")

# Hierarchical clustering ("average" = UPGMA)
plastid_hclust = hclust(plastid_dist, method = "average")

# Plot dendrogram
plot(plastid_hclust)

# Plot dendrogram using ape
plot(as.phylo(plastid_hclust), type = "phylogram")

# Extract dendrogram data from hclust results
library(ggdendro)
plastid_dend = dendro_data(plastid_hclust)
names(plastid_dend)
head(plastid_dend$segments)
head(plastid_dend$labels)

# Create dataframe individuals and site labels
plastid_df = data.frame(Ind = plastid_dend$labels$label,
                     Site = substr(plastid_dend$labels$label, 1, 3))
head(plastid_df)
unique(plastid_df$Site) ; length(unique(plastid_df$Site))

# Function to add a specific colour to dataframe depending on site
addcolours = function(x){
  # If pop label is present function will output the colour
  if(x=="Nor") y = "yellow3"
  if(x=="Zar") y = "#377eb8"
  if(x=="Mil") y = "#377eb8"
  if(x=="Man") y = "#377eb8"
  if(x=="Fal") y = "#ff7f00"
  if(x=="Mor") y = "#984ea3"
  if(x=="Tre") y = "#4daf4a"
  if(x=="Roc") y = "black"
  if(x=="Bor") y = "#e41a1c"
  if(x=="Ons") y = "#e41a1c"
  if(x=="Arm") y = "deeppink"
  return(y)
}
plastid_df$tree_cols = sapply(plastid_df$Site, addcolours)
head(plastid_df)

# Italic label
plastidlab = expression(italic("Phymatolithon calcareum") ~ "plastome CDS alignment")

# Plot custom dendrogram using ggplot2
tree2 = ggplot(plastid_dend$segments, horiz = TRUE)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.8)+
  # geom_text(data = plastid_dend$labels, aes(x = x, y = y, label = label),
  #           size = 5, fontface = "bold", hjust = -0.25, col = plastid_df$tree_cols)+
  # coord_flip()+
  # scale_y_reverse(limits = c(max(plastid_dend$segments$y)+0.00005,
  #                            min(plastid_dend$segments$y)-0.00005
  #                            )
  #                 )+
  # scale_y_continuous(limits =  c(min(plastid_dend$segments$y)-0.0001,
  #                                max(plastid_dend$segments$y)))+
  coord_cartesian(clip = "off")+
  geom_text(data = plastid_dend$labels, aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", hjust = 1.2, vjust = 0.5, col = plastid_df$tree_cols,
            angle = 90)+
  labs(title = plastidlab,
       subtitle = "(417 variant sites)",
       tag = "(b)")+
  theme_dendro()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        plot.tag = element_text(size=18, face = "bold"),
        plot.margin = margin(t = 0, r = 0, b = 1, l = 0, unit = "cm"))
tree2
# ggsave("dendrogram_plastomes.png", width = 15, height = 10, dpi = 600)
# ggsave("dendrogram_plastomes.png", width = 10, height = 15, dpi = 600)


# ----------------- #
#
# Create composite figure
#
# ----------------- #

# Combine ggplots
library(ggpubr)
fig3 = ggarrange(tree1, tree2, nrow = 2)
annotate_figure(fig3,
                top = text_grob("Hierarchical clustering using Tamura and Nei (1993) genetic distances\n",
                                face = "bold", size = 20))
ggsave(plot = fig3, "./Figures/Figure3.png", width = 15, height = 10, dpi = 600)
ggsave(plot = fig3, "./Figures/Figure3.pdf", width = 15, height = 10)



# ----------------- #
#
# Species mitogenomes tree
#
# ----------------- #

# Load packages
library(ape)
library(adegenet)
library(poppr)
library(ggplot2)

# Import mitogenomes alignment
mito_aln = ape::read.FASTA("Input_data/Maerl_mito23CDS_mafft.fasta")
names(mito_aln)
names(mito_aln) = substr(names(mito_aln), 1, 5)
names(mito_aln)
mito_aln

# Compute Tamura and Nei 1993 genetic distances
mito_dist = ape::dist.dna(mito_aln, model = "TN93")

# Hierarchical clustering ("average" = UPGMA)
mito_hclust = hclust(mito_dist, method = "average")

# Plot quick dendrogram
plot(mito_hclust)

# Plot quick phylogram using ape
plot.phylo(as.phylo(mito_hclust), type = "phylogram")

# Extract dendrogram data from hclust results
library(ggdendro)
mito_dend = dendro_data(mito_hclust)
names(mito_dend)
head(mito_dend$segments)
head(mito_dend$labels)

# Create dataframe individuals and site labels
mito_df = data.frame(Ind = mito_dend$labels$label,
                     Site = substr(mito_dend$labels$label, 1, 3))
head(mito_df)
unique(mito_df$Site) ; length(unique(mito_df$Site))

# Make a colour palette
library(scales)
show_col(brewer.pal(11, "RdYlBu"))
cols = brewer.pal(11, "RdYlBu") ; cols

# Function to add a specific colour to dataframe depending on site
# red "#e41a1c"
# blue "#377eb8"
# purple "#984ea3"
# orange "#ff7f00"
# green "#4daf4a"
addcolours = function(x){
  # If pop label is present function will output the colour
  if(x=="Nor") y = "#c51b8a"
  if(x=="Zar") y = "#c51b8a"
  if(x=="Mil") y = "#c51b8a"
  if(x=="Man") y = "#c51b8a"
  if(x=="Fal") y = "#c51b8a"
  if(x=="Mor") y = "#c51b8a"
  if(x=="Tre") y = "#c51b8a"
  if(x=="Roc") y = "#c51b8a"
  if(x=="Bor") y = "#c51b8a"
  if(x=="Ons") y = "#c51b8a"
  if(x=="Arm") y = "#c51b8a"
  if(x=="Lco") y = "#d9d9d9"
  if(x=="Sco") y = "#737373"
  return(y)
}
mito_df$tree_cols = sapply(mito_df$Site, addcolours)
head(mito_df)

# Convert Sco to Ppur
mito_dend$labels$label = gsub("Sco", "Ppur", mito_dend$labels$label)
mito_dend$labels$label

# Italic label
mitolab = expression("Maerl species mitogenome CDS alignment")

# Plot custom dendrogram using ggplot2
tree3 = ggplot(mito_dend$segments, horiz = TRUE)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.8)+
  # geom_text(data = mito_dend$labels, aes(x = x, y = y, label = label),
  #           size = 5, fontface = "bold", hjust = -0.25, col = mito_df$tree_cols)+
  # coord_flip()+
  # scale_y_reverse(limits = c(max(mito_dend$segments$y)+0.00005,
  #                            min(mito_dend$segments$y)-0.00005
  #                            )
  #                 )+
  coord_cartesian(clip = "off")+
  geom_text(data = mito_dend$labels, aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", hjust = 1.2, vjust = 0.5, col = mito_df$tree_cols,
            angle = 90)+
  labs(title = mitolab,
       subtitle = "(60.9 % of sites were invariant)",
       tag = "(a)")+
  theme_dendro()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5,),
        plot.tag = element_text(size=18, face = "bold"),
        plot.margin = margin(t = 0, r = 0, b = 1.5, l = 0, unit = "cm"))
tree3
# ggsave("dendrogram_mitogenomes.png", width = 15, height = 10, dpi = 600)


# ----------------- #
#
# Plastomes tree
#
# ----------------- #

# Import plastomes alignment
plastid_aln = ape::read.FASTA("Input_data/Maerl_plastid212CDS_mafft.fasta", type = "DNA")
names(plastid_aln)
names(plastid_aln) = substr(names(plastid_aln), 1, 5)
names(plastid_aln)
plastid_aln

# Compute Tamura and Nei 1993 genetic distances
plastid_dist = ape::dist.dna(plastid_aln, model = "TN93")

# Hierarchical clustering ("average" = UPGMA)
plastid_hclust = hclust(plastid_dist, method = "average")

# Plot dendrogram
plot(plastid_hclust)

# Plot dendrogram using ape
plot(as.phylo(plastid_hclust), type = "phylogram")

# Extract dendrogram data from hclust results
library(ggdendro)
plastid_dend = dendro_data(plastid_hclust)
names(plastid_dend)
head(plastid_dend$segments)
head(plastid_dend$labels)

# Create dataframe individuals and site labels
plastid_df = data.frame(Ind = plastid_dend$labels$label,
                        Site = substr(plastid_dend$labels$label, 1, 3))
head(plastid_df)
unique(plastid_df$Site) ; length(unique(plastid_df$Site))

# Function to add a specific colour to dataframe depending on site
addcolours = function(x){
  # If pop label is present function will output the colour
  if(x=="Nor") y = "#c51b8a"
  if(x=="Zar") y = "#c51b8a"
  if(x=="Mil") y = "#c51b8a"
  if(x=="Man") y = "#c51b8a"
  if(x=="Fal") y = "#c51b8a"
  if(x=="Mor") y = "#c51b8a"
  if(x=="Tre") y = "#c51b8a"
  if(x=="Roc") y = "#c51b8a"
  if(x=="Bor") y = "#c51b8a"
  if(x=="Ons") y = "#c51b8a"
  if(x=="Arm") y = "#c51b8a"
  if(x=="Lco") y = "#d9d9d9"
  if(x=="Sco") y = "#737373"
  return(y)
}
plastid_df$tree_cols = sapply(plastid_df$Site, addcolours)
head(plastid_df)

# Convert Sco to Ppur
plastid_dend$labels$label = gsub("Sco", "Ppur", plastid_dend$labels$label)
plastid_dend$labels$label

# Italic label
plastidlab = expression("Maerl species plastome CDS alignment")

# Plot custom dendrogram using ggplot2
tree4 = ggplot(plastid_dend$segments, horiz = TRUE)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.8)+
  # geom_text(data = plastid_dend$labels, aes(x = x, y = y, label = label),
  #           size = 5, fontface = "bold", hjust = -0.25, col = plastid_df$tree_cols)+
  # coord_flip()+
  # scale_y_reverse(limits = c(max(plastid_dend$segments$y)+0.00005,
  #                            min(plastid_dend$segments$y)-0.00005
  #                            )
  #                 )+
  # scale_y_continuous(limits =  c(min(plastid_dend$segments$y)-0.0001,
  #                                max(plastid_dend$segments$y)))+
  coord_cartesian(clip = "off")+
  geom_text(data = plastid_dend$labels, aes(x = x, y = y, label = label),
            size = 5, fontface = "bold", hjust = 1.2, vjust = 0.5, col = plastid_df$tree_cols,
            angle = 90)+
  labs(title = plastidlab,
       subtitle = "(73.1 % of sites were invariant)",
       tag = "(b)")+
  theme_dendro()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        plot.tag = element_text(size=18, face = "bold"),
        plot.margin = margin(t = 0, r = 0, b = 1.5, l = 0, unit = "cm"))
tree4
# ggsave("dendrogram_plastomes.png", width = 15, height = 10, dpi = 600)
# ggsave("dendrogram_plastomes.png", width = 10, height = 15, dpi = 600)


# ----------------- #
#
# Create composite figure
#
# ----------------- #

# Combine ggplots
library(ggpubr)
figS6 = ggarrange(tree3, tree4, nrow = 2)
annotate_figure(figS6,
                top = text_grob("Hierarchical clustering using Tamura and Nei (1993) genetic distances\n",
                                face = "bold", size = 20))
figS6
ggsave(plot = figS6, "./Figures/FigureS6.png", width = 15, height = 10, dpi = 600)
ggsave(plot = figS6, "./Figures/FigureS6.pdf", width = 15, height = 10)
