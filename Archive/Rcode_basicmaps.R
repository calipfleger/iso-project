#### Map proxy locations ####
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtext)

# Define country borders and point types
data("wrld_simpl", package = "maptools")
proxy = rep(15, 6)

# Import lat/lon data and site names
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data")
coral_loc = read.csv("coral_proxy_locations.csv", head = T)

ggplot() +
  # Plot country borders
  geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
               fill = "grey", colour = "black", alpha = 0.2) +
  # Plot proxy locations
  geom_point(data = coral_loc, mapping = aes(x = lon, y = lat, 
                                             shape = archive), color = "blue", size = 7) +
  scale_shape_manual(values = proxy) + 
  # Label proxies
  geom_text(data = coral_loc, aes(x = lon, y = lat, label = archive),
            hjust = -0.3, cex = 4) +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Change map borders
  coord_equal(xlim = c(-180, 0), ylim = c(-90, 90)) +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())