# --------------------------- #
#
# Maerl Genomics Project
#
# Description:
# SNMF analysis
#
# https://github.com/Tom-Jenkins/maerl_genomics
#
# --------------------------- #

# Load packages
library(adegenet)
library(poppr)
library(LEA)
library(ggplot2)
library(RColorBrewer)
library(ggplotify)
library(reshape2)
library(scales)
library(ggpubr)

# Import Genind data
load("Output_rdata_files/maerl_SNPs_snmf.RData")
maerl = maerl.snmf
summary(maerl$pop)
nPop(maerl)


# ----------------- #
#
# Admixture analysis using SNMF
#
# ----------------- #

# Run snmf algorithm
# snmf1 = snmf("SNMF/maerl_SNPs.geno",
#              K = 1:nPop(maerl), # number of K ancestral populations to run
#              repetitions = 10, # ten repetitions for each K
#              entropy = TRUE, # calculate cross-entropy
#              project = "new",
#              seed = 1234)

# Load snmf project
snmf1 = load.snmfProject("SNMF/maerl_SNPs.snmfProject")

# Plot cross-entropy results to assess optimal number of K
# Smaller values of cross-entropy usually mean better runs
# A plateau usually represents the K that best fits the data
# png("./Figures/FigureS4a.png", width = 10, height = 7, unit = "in", res = 600)
plot(snmf1, col = "blue", cex = 1.5, pch = 19, cex.axis = 1, cex.lab = 1.5,
     main = "SNMF: Cross-entropy", cex.main = 1.5)
# dev.off()
# pdf("./Figures/FigureS4a.pdf", width = 10, height = 7)
# plot(snmf1, col = "blue", cex = 1.5, pch = 19, cex.axis = 1, cex.lab = 1.5,
#      main = "SNMF: Cross-entropy", cex.main = 1.5)
# dev.off()


# ----------------- #
#
# R function to plot each K
#
# ----------------- #

# Function to plot each K from snmf results
# Argument k = number of ancestral populations
# Argument colours = colour definitions of each K (default is ggplot2 palette)
plotk_func = function(k, colours = hue_pal()(k)){
  
  # Extract the cross-entropy of all runs where K = k
  ce = cross.entropy(snmf1, K = k)
  
  # Find the run with the lowest cross-entropy
  lowest.ce = which.min(ce)
  
  # Extract Q-matrix for the best run
  qmatrix = as.data.frame(Q(snmf1, K = k, run = lowest.ce))
  
  # Label column names of qmatrix
  cluster_names = c()
  for (i in 1:k){
    cluster_names[i] = paste("Cluster", i)
  }
  colnames(qmatrix) = cluster_names
  
  # Add individual IDs
  qmatrix$Ind = indNames(maerl)
  
  # Add site IDs
  qmatrix$Site = as.character(maerl$pop)
  
  # Convert dataframe to long format
  qlong = melt(qmatrix, id.vars=c("Ind","Site"), variable.name="Ancestry", value.name="Proportion")
  
  # Sort data.frame with ancestry order by overall proportion of each ancestry
  # qlong = qlong %>% arrange(Ind, desc(Proportion))
  
  # Change order of individuals by using the factor function
  ind.order = c(grep("Zar", indNames(maerl), value = T),
                grep("Man", indNames(maerl), value = T),
                grep("Fal", indNames(maerl), value = T),
                grep("Mor", indNames(maerl), value = T),
                grep("Tre", indNames(maerl), value = T),
                # grep("Roc", indNames(maerl), value = T),
                grep("Bor", indNames(maerl), value = T),
                grep("Ons", indNames(maerl), value = T))
  qlong$Ind = factor(qlong$Ind, levels = rev(ind.order))

  # Plot admixture barplot
  snmfbar = ggplot(data=qlong, aes(x=Ind, y=Proportion, fill=Ancestry))+
    geom_bar(stat = "identity", width=1, colour="black", size=0.25, show.legend=FALSE)+
    coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(values = as.character(colours))+
    ylab("Admixture proportion")+
    # xlab("Individual")+
    theme(axis.text.y = element_text(face = "bold", colour = "black", size = 10),
          axis.text.x = element_text(angle = 360, colour = "black", size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.margin = margin(t = 0, r = 20, b = 0, l = 0, unit = "pt")
          )
  
  # Return barplot
  return(snmfbar)
}

# Extract barplots for K2-7
snmfbar_K2 = plotk_func(k = 2)
snmfbar_K3 = plotk_func(k = 3)
snmfbar_K4 = plotk_func(k = 4)
snmfbar_K5 = plotk_func(k = 5)
snmfbar_K6 = plotk_func(k = 6)
snmfbar_K7 = plotk_func(k = 7)

# Combine barplots
snmfbars = ggarrange(
  snmfbar_K2 + labs(title = "K2") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  snmfbar_K3 + labs(title = "K3") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  snmfbar_K4 + labs(title = "K4") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  snmfbar_K5 + labs(title = "K5") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  snmfbar_K6 + labs(title = "K6"),
  snmfbar_K7 + labs(title = "K7"),
  nrow = 3, ncol = 2
  )

# Figure S4b
snmfbars = annotate_figure(snmfbars,
                           text_grob("SNMF admixture K2-7", face = "bold", size = 25))
snmfbars
# ggsave(plot = snmfbars, "Figures/FigureS4b.png", width = 15, height = 20, dpi = 600)
# ggsave(plot = snmfbars, "Figures/FigureS4b.pdf", width = 15, height = 20)


# ----------------- #
#
# Figure 2 barplot (K = 6)
#
# ----------------- #

# Colour definitions for colour argument
# blue "#377eb8"
# orange "#ff7f00"
# blue-green "#009E73"
# purple "#984ea3"
# red "#e41a1c"
# dark red "#a50f15"
cols = data.frame(
  "bor_cluster" = "#e41a1c",
  "tre_cluster" = "#009E73",
  "zar_cluster" = "#377eb8",
  "fal_cluster" = "#ff7f00",
  "Ons_cluster" = "#a50f15",
  "mor_cluster" = "#984ea3",
  check.names = FALSE
)

# Plot barplot
fig2_bar = plotk_func(k = 6, colours = cols)
fig2_bar


# ----------------- #
#
# Figure 2 pie charts (K = 6)
#
# ----------------- #

# Select K
k = 6

# Extract the cross-entropy of all runs where K = k
ce = cross.entropy(snmf1, K = k)

# Find the run with the lowest cross-entropy
lowest.ce = which.min(ce)

# Extract Q-matrix for the best run
qmatrix = as.data.frame(Q(snmf1, K = k, run = lowest.ce))

# Label column names of qmatrix
cluster_names = c()
for (i in 1:k){
  cluster_names[i] = paste("Cluster", i)
}
colnames(qmatrix) = cluster_names

# Add individual IDs
qmatrix$Ind = indNames(maerl)

# Add site IDs
qmatrix$Site = maerl$pop

# Calculate mean admixture proportions for each site
head(qmatrix)
clusters = grep("Cluster", names(qmatrix)) # indexes of cluster columns
avg_admix = aggregate(qmatrix[, clusters], list(qmatrix$Site), mean)

# Order alphabetically by site
avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
avg_admix

# Convert dataframe from wide to long format
avg_admix = melt(avg_admix, id.vars = "Group.1")
head(avg_admix)

# Define a function to plot pie charts using ggplot for each site
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, Group.1 == site),
         aes(x = "", y = value, fill = variable))+
    geom_bar(width = 1, alpha = 0.95, stat = "identity", colour = "black", size = 0.5, show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = as.character(cols))+
    theme_void()
}

# Test function on one site
pie_charts(avg_admix, site = "Man", cols = as.character(cols))

# Apply function to all sites using for loop
pies = list()
for (i in levels(avg_admix$Group.1)){
  pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = as.character(cols)) 
}


# ----------------- #
#
# Prepare basemap
#
# ----------------- #

# Load R packages
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(raster)

# Import csv file containing coordinates
coords = read.csv("Input_data/coordinates.csv")

# Subset coordinates
coords = coords[coords$Code %in% levels(avg_admix$Group.1), ]

# Order alphabetically by site
coords = coords[order(coords$Code), ] 
coords

# Check order matches coords order
as.character(avg_admix$Group.1) == as.character(coords$Code)

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-15, 5, 38, 56)
boundary

# Download a basemap
map.outlines = ne_countries(scale = "large")

# Crop to boundary and convert to dataframe
map.outlines = crop(map.outlines, y = boundary) %>% fortify()

# Plot basemap
basemap = ggplot()+
  geom_polygon(data=map.outlines, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", size=0.5)+
  coord_quickmap(xlim = c(-12,5), expand = FALSE)+
  ggsn::north(map.outlines, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outlines, dist = 100, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", border.size = 0.5,
                 location = "bottomleft", anchor = c(x = -10, y = 45),
                 st.bottom = FALSE, st.size = 4, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  geom_rect(aes(xmin=-5.8, xmax=-4, ymin=49.9, ymax=50.6),
            fill=NA, colour="#252525", linetype=1, size=0.6)+
  geom_segment(aes(x=-5.8, y=50.6, xend=-2, yend=54.4), colour="#252525", linetype=2)+
  geom_segment(aes(x=-4, y=49.9, xend=5, yend=51.6), colour="#252525", linetype=2)+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top"
  )
basemap


# ----------------- #
#
# Add pie charts to basemap
#
# ----------------- #

# Extract coordinates for each site
coord.list = list()
for (i in levels(avg_admix$Group.1)){
  coord.list[[i]] = c(subset(coords, Code == i)$Lon, subset(coords, Code == i)$Lat)
}
coord.list

# Define pie chart sizes
radius = c(1, 0.1, 0.1, 1, 1, 1, 1)

# Convert ggplot pie charts to annotation_custom layers
pies.ac = list()
for (i in 1:length(coord.list)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius[i],
                                   xmax = coord.list[[i]][[1]] + radius[i],
                                   ymin = coord.list[[i]][[2]] - radius[i],
                                   ymax = coord.list[[i]][[2]] + radius[i])
}

# Add layers to basemap
pie.map = basemap + pies.ac
pie.map
# ggsave("1.pie_charts_map.png", width = 8, height = 10, dpi = 300)


# ----------------- #
#
# Zoomed-in map of SW England
#
# ----------------- #

# Create a zoomed-in map of SW England
snmf.eng = ggplot() +
  geom_polygon(data=map.outlines, aes(x=long, y=lat, group=group), fill="grey",
               colour="black", size=0.5)+
  coord_quickmap(xlim=c(-5.8,-4.5), ylim=c(49.9,50.4), expand=F)+
  # coord_map(projection = "mercator", xlim=c(-5.8,-4.5), ylim=c(49.9,50.4))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text= element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin=grid::unit(c(0,0,-1,-1), "mm")
  )
snmf.eng

# Define pie chart sizes
radius = 0.15

# Extract only pies and coords for Fal and Man
coord.fal = coord.list[2:3]
pies.fal = pies[2:3]

# Convert ggplot pie charts to annotation_custom layers
fal.pies = list()
for (i in 1:length(coord.fal)){
  fal.pies[[i]] = annotation_custom(grob = ggplotGrob(pies.fal[[i]]),
                                   xmin = coord.fal[[i]][[1]] - radius,
                                   xmax = coord.fal[[i]][[1]] + radius,
                                   ymin = coord.fal[[i]][[2]] - radius,
                                   ymax = coord.fal[[i]][[2]] + radius)
}

# Add Falmouth and Manacles pies to SW England map
pie.eng = snmf.eng + fal.pies
pie.eng
# ggsave("2.pie_charts_Eng.png", width = 10, height = 8, dpi = 300)

# Add annotation labels
fal.pies[[1]]$aes_params
site_labs2 = data.frame(labels = c("Fal","Man"),
                       x = c(-4.877973, -5.177973),
                       y = c(50.31351, 50.01351))
site_labs2 = data.frame(labels = names(coord.fal),
                        x = c(coord.fal$Fal[1], coord.fal$Man[1]) + 0.22,
                        y = c(coord.fal$Fal[2], coord.fal$Man[2]))
pie.eng.lab = pie.eng +
  geom_label(data = site_labs2, aes(x = x, y = y, label = labels), 
             fontface = "bold", size = 5)
pie.eng.lab


# ----------------- #
#
# Create inset map
#
# ----------------- #

# Inset Falmouth zoomed-in map over main map
inset = pie.map +
  annotation_custom(grob = ggplotGrob(pie.eng.lab),
                    xmin=-2, xmax=5, ymin=51, ymax=55)
# inset
# ggsave("snmf_map", width = 10, height = 8, dpi = 300)

# Add labels to inset
site_labs1 = subset(coords, Code != "Fal" & Code != "Man")
site_labs1
fig2_map = inset +
  geom_label(data = site_labs1, aes(x = Lon + 1.4, y = Lat, label = Code),
             fontface = "bold", size = 5)
fig2_map

# Combine ggplots
fig2 = ggarrange(ncol = 2, widths = c(1,3),
                 fig2_bar + labs(tag = "(a)") + theme(plot.tag = element_text(size=15, face = "bold")),
                 fig2_map + labs(tag = "(b)") +
                    theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = -180, unit = "pt"),
                      plot.tag = element_text(size=15, face = "bold"))
                 )
fig2
# ggsave(plot = fig2, "./Figures/Figure2.png", width = 12.8, height = 9.97, dpi = 600)
# ggsave(plot = fig2, "./Figures/Figure2.pdf", width = 12.8, height = 9.97)
