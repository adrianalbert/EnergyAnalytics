plotmatrix2 <- function (data, mapping = aes())
{
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
               x = data[, xcol], y = data[, ycol], data)
  }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))
  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
    data.frame(xvar = names(data)[i], yvar = names(data)[i], 
               x = data[, i])
  }))
  densities$xvar <- factor(densities$xvar, levels = names(data))
  densities$yvar <- factor(densities$yvar, levels = names(data))
  mapping <- defaults(mapping, aes_string(x = "x", y = "y"))
  class(mapping) <- "uneval"
  ggplot(all) + facet_grid(xvar ~ yvar, scales = "free") + 
    geom_point(mapping, na.rm = TRUE) + stat_density(aes(x = x, 
                                                         y = ..scaled.. * diff(range(x)) + min(x)), data = densities, 
                                                     position = "identity", colour = "grey20", geom = "line")
}


plotmatrix2(mtcars[,1:3],aes(colour = factor(cyl)))