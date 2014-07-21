# plot_metrics.r
#
# Plot metrics to validate individual model performance. 
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

library(ggplot2)
library(reshape)
library(gplots)

# creates plot of multiple Gaussian densities on the same axis
plot_gaussian_distr = function(distr, type = NULL, state = NULL,
                               title = 'Distribution Comparison',
                               xlab = 'Response [kWh/deg F]', ylab = 'pdf') {
  if (is.null(type)) type = rep('Estimate', length(distr))
  if (is.null(state)) state = 1:length(distr)

  dat = lapply(1:length(distr), function(i) {
    d = rnorm(1000, mean = distr[[i]]['mu'], sd = distr[[i]]['sd'])
    data.frame(data = d, type = type[i], state = state[i])
  })
  dat = do.call('rbind', dat);
  
  p = ggplot(dat, aes(data, color = state)) + geom_density(aes(linetype = type), size=2)
  p = p + facet_wrap(~state, scales = 'free')
  p = p + scale_linetype_manual(values=setdiff(1:(length(unique(dat$type))+1), 2)) 
  p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),        
          legend.title     = element_text(size=18),    
          axis.ticks       = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20))
  p = p + ggtitle(title)  + xlab(xlab) + ylab(ylab)
    
  return(p)
}

# plots CDF/Risk functions per state
plot_error_state = function(P, support,
                            title = 'Error Probability',
                            xlab = 'Response [kWh/deg F]', ylab = 'P(error)') {
  
  dat = as.data.frame(t(P)); names(dat) = 1:ncol(dat); dat$support = support
  dat = melt(dat, id.vars = 'support'); dat$variable = as.factor(dat$variable)
  
  p = ggplot(dat, aes(support, value, color = variable)) 
  p = p + geom_line(size=2) + geom_point(aes(shape=variable), size=4)
  p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),        
          legend.title     = element_text(size=18),    
          axis.ticks       = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20))
  p = p + ggtitle(title)  + xlab(xlab) + ylab(ylab)
  
  return(p)
}

# plots noisy response curves with respect to a given variable
plot_noisy_curve = function(curves, 
                            title = 'Error Probability',
                            xlab = 'Temperature [deg F]', ylab = 'Response [kwh/deg F]') {
  
  labels = names(curves)
  if (is.null(labels)) labels = as.factor(1:length(curves))
  df = lapply(1:length(curves), function(t) {
    cur = curves[[t]]; 
    data.frame(mu = cur$mu, sd = cur$sd, var = cur$var, curve = labels[t])
  })
  df = do.call('rbind', df)
  
  p = ggplot(df, aes(var, mu, color = curve))
  p = p + geom_line(size = 1.5) + geom_point(aes(shape = curve), size = 3.5)
  p = p + geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.1)  
  p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),        
          legend.title     = element_text(size=18),    
          axis.ticks       = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20))
  p = p + ggtitle(title)  + xlab(xlab) + ylab(ylab)
  
  return(p)
}

# plot overlapping densities
plot_overlapping_densities = function(df, var = 'TemperatureF', by = 'Season', 
                                      title = NULL, xlab = var, ylab = 'pdf') {
  
  p  = ggplot(df, aes_string(x = var, color = by, fill = by)) + geom_density(size = 2, alpha = 0.5)
  p = p + theme_bw() + 
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),        
          legend.title     = element_text(size=18),    
          axis.ticks       = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20))
  p = p + ggtitle(title) + xlab(xlab) + ylab(ylab)
  

}

# plots noisy response curves with respect to a given variable
plot_prob_profile = function(P, var = NULL,
                            title = 'Error Probability',
                            xlab = 'Temperature [deg F]', ylab = 'P(state)') {
  
  if (is.null(var)) var = 1:nrow(P)
  colnames(P) = 1:ncol(P)
  P = as.data.frame(cbind(var, P)) 
  var_name = names(P)[1]
  df = melt(P, id.vars = var_name)
  
  p = ggplot(df, aes_string(x = var_name, y = 'value', color = 'variable', group = 'variable'))
  p = p + geom_line(size = 1.5) + geom_point(size = 3.5)
  p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),        
          legend.title     = element_text(size=18),    
          axis.ticks       = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20))
  p = p + ggtitle(title)  + xlab(xlab) + ylab(ylab)
  
  return(p)
}

# plots trellis plot (dotplot)
# TODO: still needs some work...
ggplot_dotplot = function(df, facet = NULL, by = NULL){
      
  p = ggplot(df, aes_string(x='value', y=by)) + geom_point(size=4)
  p = p + geom_line(aes(xmin = 0, xmax = value), size=1.5)
  p = p + facet_wrap(as.formula(paste("~", facet)), scales = 'free')
  p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_text(size = 18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=14),        
          legend.title     = element_text(size=14),  
          legend.position  = c(0.84, 0.8),
          axis.ticks       = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20))
  p = p + ggtitle('Avg. thermal response distribution: models and true')
  p = p + xlab('Response [kWh/deg F]') + ylab('pdf')

  return(p)
}
  

