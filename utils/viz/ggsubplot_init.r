library(ggplot2)
library(maps)
library(plyr)

# getbox by Heike Hoffman, 
# https://github.com/ggobi/paper-climate/blob/master/code/maps.r
getbox <- function (map, xlim, ylim) {
  # identify all regions involved
  small <- subset(map, (long > xlim[1]) & (long < xlim[2]) & (lat > ylim[1]) & (lat < ylim[2]))
  regions <- unique(small$region)
  small <- subset(map, region %in% regions)  
  
  # now shrink all nodes back to the bounding box
  small$long <- pmax(small$long, xlim[1])
  small$long <- pmin(small$long, xlim[2])
  small$lat <- pmax(small$lat, ylim[1])
  small$lat <- pmin(small$lat, ylim[2])
  
  # Remove slivvers
  small <- ddply(small, "group", function(df) {
    if (diff(range(df$long)) < 1e-6) return(NULL)
    if (diff(range(df$lat)) < 1e-6) return(NULL)
    df
  })
  
  small
}


## map layer
## adapted from map_nasa:
# https://github.com/ggobi/paper-climate/blob/master/code/maps.r

# assembling data
world <- map_data("world")
states <- map_data("state")
states$group <- max(world$group) + states$group
both <- rbind(world, states)
americas <- getbox(both, xlim = c(-115, -55), ylim = c(-21.1, 36.6))

# building americas layer
map_americas <- list(
  geom_polygon(aes(long, lat, group = group), data = americas, fill = "grey70", 
               colour = "grey60", inherit.aes = FALSE, show_guide = FALSE),
  scale_x_continuous("", breaks = NULL, expand = c(0.02, 0)),
  scale_y_continuous("", breaks = NULL, expand = c(0.02, 0)))


# building afghanistan layer
afghanistan <- getbox(world, c(60,75), c(28, 39))
map_afghanistan <- list(
  geom_polygon(aes(long, lat, group = group), data = afghanistan, 
               fill = "grey80", colour = "white", inherit.aes = FALSE, 
               show_guide = FALSE),
  scale_x_continuous("", breaks = NULL, expand = c(0.02, 0)),
  scale_y_continuous("", breaks = NULL, expand = c(0.02, 0)))