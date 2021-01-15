#----------------------------
#
# Maerl Genomics Project
#
# Description:
# Plot maerl bed distribution map
# Plot Fal Estuary map
#
# https://github.com/Tom-Jenkins/maerl_genomics
#
#----------------------------

# Load packages
library(sf)
library(tidyverse)
library(ggspatial)


# ----------------- #
#
# Prepare OSPAR maerl bed records
#
# ----------------- #

# Import OSPAR shapefile
# Shapefile data downloaded from here:
# https://www.emodnet-seabedhabitats.eu/access-data/download-data/?linkid=ospar2018_poly
ospar = file.path("Input_data/OSPAR_declining_habitats_2018/Shapefiles/OSPARHabitats2018_points_PUBLIC.shp") %>%
  st_read()
st_crs(ospar)

# Transform coordinate system
ospar = st_transform(ospar, 4326)

# Extract maerl bed records
table(ospar$HabType)
maerlbeds = dplyr::filter(ospar, HabType == "Maerl beds")

# Filter out absences and uncertain records
# Only keep record ID, HabStatus, Certainty and geometry
maerlbeds = maerlbeds %>% 
  dplyr::filter(HabStatus == "Present" & Certainty == "Certain") %>%
  dplyr::select(RecordKey, HabStatus, Certainty, Longitude, Latitude)
maerlbeds


# ----------------- #
#
# Plot map of maerl bed records
#
# ----------------- #

# Import Europe basemap from rnaturalworld package
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
europe = ne_countries(scale = "large", returnclass = "sf") %>%
  dplyr::filter(continent == "Europe") %>% 
  select(name_en)

# Plot basemap with maerl beds
map1 = ggplot()+
  geom_sf(data = europe, fill = "grey")+
  geom_sf(data = maerlbeds, col = "pink")+
  coord_sf(xlim = c(-13, 10), ylim = c(36, 63), expand = F)
map1


# ----------------- #
#
# Prepare MPA data
#
# ----------------- #

# Import England boundary data
england = file.path("Input_data/MPA_boundary_data/england_ol_2001.shp") %>% st_read()
st_crs(england)

# Transform coordinate system to the EPSG (4326)
england = st_transform(england, 4326)
st_crs(england)
england

# Import Special Area of Conservation (SAC) boundary data
mcz = file.path("Input_data/MPA_boundary_data/Marine_Conservation_Zones_England.shp") %>% st_read()
st_crs(mcz)

# Subset The Manacles MCZ
mcz = filter(mcz, mcz_name == "The Manacles")
mcz

# Transform coordinate system to the EPSG (4326)
mcz = st_transform(mcz, 4326)
st_crs(mcz)
mcz

# Import Special Area of Conservation (SAC) boundary data
sac = file.path("Input_data/MPA_boundary_data/SACs_with_MC_Desig_Prop_20200602.shp") %>% st_read()
st_crs(sac)

# Subset Fal and Helford SAC boundary data
sac = filter(sac, SiteName == "Fal and Helford")
sac

# Check all CRS match
st_crs(england) == st_crs(mcz)
st_crs(england) == st_crs(sac)


# ----------------- #
#
# Plot map of Falmouth
#
# ----------------- #

# Point data as sfc object
fal_gps = st_point(c(-5.028, 50.164))
man_gps =  st_point(c(-5.05, 50.05))
points = st_sfc(fal_gps, man_gps, crs = 4326)

# Prepare attribute data and convert to sf object
attrib = data.frame(Site = c("Fal","Man"))
points = st_sf(attrib, geometry = points)

# Transform coordinate system to the EPSG (4326)
points = st_transform(points, 4326)
points

# Plot map of Falmouth with MPAs
falmouth = ggplot()+
  geom_sf(data = england, fill = "grey")+
  geom_sf(data = sac, aes(colour = SiteName), fill = "white")+
  geom_sf(data = mcz, aes(colour = mcz_name), fill = "white")+
  geom_sf(data = points, aes(fill = Site), size = 3, shape = 21, colour = "black")+
  coord_sf(xlim = c(-5.25,-4.9), ylim = c(50.02,50.28), expand = F)+
  annotation_scale(data = england, location = "br", width_hint = 0.2,
                   text_cex = 0.8)+
  annotation_north_arrow(data = england, location = "tr",
                         width = unit(0.6, "cm"), height = unit(0.6, "cm"))+
  scale_colour_manual("Marine reserves",
                      values = c("red","blue"),
                      labels = c("Fal and Helford SAC","The Manacles MCZ"))+
  scale_fill_manual("Sampling point",
                    values = c("#ff7f00","#377eb8"),
                    labels = c("Fal Estuary","The Manacles"))+
  theme(
    axis.text = element_text(colour = "black", size = 10),
    panel.grid = element_line(colour = "grey90", size = 0.5),
    panel.background = element_rect(fill = "lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.7, "cm"),
    legend.key = element_rect(fill = NA),
    legend.position = "right"
    )
falmouth
  

# ----------------- #
#
# Plot map of Europe
#
# ----------------- #

# Import basemap from rnaturalworld package
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
world = ne_countries(scale = "large", returnclass = "sf")

# Transform coordinate system to the EPSG (4326)
world = st_transform(world, 4326)
st_crs(world) == st_crs(england)

# Create a data.frame for country labels
countries = data.frame(Country = c("UK","Ireland","France"),
                       X = c(-1, -7.5, 2),
                       Y = c(52, 53.2, 49))

# Plot map of the British Isles and France
europe = ggplot()+
  geom_sf(data = world, fill = "grey")+
  coord_sf(xlim = c(-12,4.5), ylim = c(48,55), expand = F)+
  geom_rect(aes(xmin=-5.25, xmax=-4.9, ymin=50.02, ymax=50.28),
            fill = NA, colour = "black", linetype = 1, size = 1)+
  geom_text(data = countries, aes(X, Y, label = Country), size = 3)+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    plot.margin = grid::unit(c(0,0,-1,-1), "mm")
  )
# europe      


# ----------------- #
#
# Create inset map
#
# ----------------- #

# Inset Falmouth zoomed-in map over main map
inset = falmouth +
  annotation_custom(grob = ggplotGrob(europe),
                    xmin = -5.24, xmax = -5.08, ymin = 50.2, ymax = Inf)
ggsave(plot = inset, "./Figures/Figure5.png", width = 10, height = 7, dpi = 600)
ggsave(plot = inset, "./Figures/Figure5.pdf", width = 10, height = 7)

