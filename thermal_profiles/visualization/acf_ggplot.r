# acf_ggplot.r
# 
# Produces an ACF plot using ggplot.
# Can compare theoretical vs empirical correlations.

acf_ggplot <- function(x, x.ref = NULL, conf.level = 0.95, 
                 max.lag = NULL, min.lag = 0, title = "", PACF = F) {
   ciline <- qnorm((1 - conf.level)/2)/sqrt(length(x))
   
   if (!PACF) {
     bacf      <- acf(x, plot = FALSE, lag.max = max.lag)
     labels    <- c("Empirical ACF", "Model ACF")
     bacf.ref  <- acf(x.ref, plot = FALSE, lag.max = max.lag)
     bacfdf    <- with(bacf, data.frame(lag, acf))
     bacfdf.ref<- with(bacf.ref, data.frame(lag, acf))
   } else {
     bacf      <- pacf(x, plot = FALSE, lag.max = max.lag)$acf[,1,1]
     labels    <- c("Empirical PACF", "Fit PACF")
     bacf.ref  <- pacf(x.ref, plot = FALSE, lag.max = max.lag)$acf[,1,1]
     bacfdf    <- data.frame(acf = bacf, lag = 1:length(bacf))       
     bacfdf.ref<- data.frame(acf = bacf.ref, lag = 1:length(bacf))       
   }   

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
    theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.text.y      = element_text(size=15),
         axis.text.x      = element_text(size=15),
         axis.title.y     = element_text(size=18),
         axis.title.x     = element_text(size=18),
         plot.title       = element_text(size=20),
         axis.title.y     = element_blank(),
         legend.text      = element_text(size=15),
         legend.position  = c(0.8,0.6),
         panel.background = element_rect(fill = "transparent",colour = NA),
         axis.ticks = element_blank() ) + 
    ggtitle( title ) + xlab('Lag Order')

   q <- q + geom_hline(yintercept = -ciline, color = "blue", size = 0.8, alpha=0.5)
   q <- q + geom_hline(yintercept = ciline, color = "blue", size = 0.8, alpha = 0.5)
   q <- q + geom_hline(yintercept = 0, color = "red", size = 0.9, alpha = 0.5)
   q <- q + scale_fill_hue(name = "Legend", 
                           breaks = c('acf', 'acf.ref'), labels = labels)
   return(q)
}

pacf_ggplot_single <- function(x, conf.level = 0.95, 
                              max.lag = NULL, min.lag = 0, title = "") {
  ciline <- qnorm((1 - conf.level)/2)/sqrt(length(x))
  bacf   <- pacf(x, plot = FALSE, lag.max = max.lag)$acf[,1,1]
  bacfdf <- data.frame(pacf = bacf, lag = 1:length(bacf))  
  
  if (min.lag > 0) {
    bacfdf  <- bacfdf[-seq(1, min.lag), ]
  }
  
  significant <- (abs(bacfdf[, 2]) > abs(ciline))^2
  bacfdf <- cbind(bacfdf, significant)
  bacfdf.mlt = melt(bacfdf, measure.vars=c('pacf'))
  q <- ggplot(bacfdf.mlt, 
              aes(x=lag, y=value, fill=factor(variable))) +                    
    geom_bar(position="dodge", stat="identity", width = 0.2)
  q <- q + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = c(0.9,0.6),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.ticks = element_blank() ) + 
    ggtitle( title )
  
  q <- q + geom_hline(yintercept = -ciline, color = "blue", size = 0.8, alpha=0.5)
  q <- q + geom_hline(yintercept = ciline, color = "blue", size = 0.8, alpha = 0.5)
  q <- q + geom_hline(yintercept = 0, color = "red", size = 0.9, alpha = 0.5)
  return(q)
}

