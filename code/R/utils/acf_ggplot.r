# acf_ggplot.r
# 
# Produces an ACF plot using ggplot.
# Can compare theoretical vs empirical correlations.

acf_ggplot <- function(x, x.ref = NULL, conf.level = 0.95, 
                 max.lag = NULL, min.lag = 0, title = "") {
   ciline <- qnorm((1 - conf.level)/2)/sqrt(length(x))
   bacf   <- acf(x, plot = FALSE, lag.max = max.lag)
   if (!is.null(x.ref)) bacf.ref <- acf(x.ref, plot = FALSE, lag.max = max.lag)
   
   bacfdf    <- with(bacf, data.frame(lag, acf))
   bacfdf.ref<- with(bacf.ref, data.frame(lag, acf))
   bacfdf    <- merge(bacfdf, bacfdf.ref, 
                      by = 'lag', 
                      suffixes = c('', '.ref'))
   if (min.lag > 0) {
     bacfdf  <- bacfdf[-seq(1, min.lag), ]
   }
   
   significant <- (abs(bacfdf[, 2]) > abs(ciline))^2
   bacfdf <- cbind(bacfdf, significant)
   bacfdf.mlt = melt(bacfdf, measure.vars=c('acf', 'acf.ref'))
   q <- ggplot(bacfdf.mlt, 
               aes(x=lag, y=value, 
                   fill=factor(variable))) +                    
                 geom_bar(position="dodge", stat="identity")
   q <- q + theme_bw() +
    opts(panel.grid.major = theme_blank(),
         panel.grid.minor = theme_blank(),
         panel.background = theme_blank(),
         axis.title.y = theme_blank(),
         axis.title.x = theme_blank(),
         legend.position = c(0.9,0.6),
         panel.background = theme_rect(fill = "transparent",colour = NA),
         axis.ticks = theme_blank() ) + 
    opts( title=title, size=1)

   q <- q + geom_hline(yintercept = -ciline, color = "blue", size = 0.8, alpha=0.5)
   q <- q + geom_hline(yintercept = ciline, color = "blue", size = 0.8, alpha = 0.5)
   q <- q + geom_hline(yintercept = 0, color = "red", size = 0.9, alpha = 0.5)
   q <- q + scale_fill_hue(name = "Legend", 
                           breaks = c('acf', 'acf.ref'), labels = c("Empirical ACF", "Model ACF"))
   return(q)
}
