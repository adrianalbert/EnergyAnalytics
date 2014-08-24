
library('ggplot2')

# create "square wave" schedules
define_schedule = function(eta, beta, tau = 24) {  
  u = rep(0, tau);
  t1 = max(eta + beta - tau, 0)
  u[max(eta,1):min(eta+beta, tau)] = 1
  if (t1>0) u[1:t1] = 1
  return(u)
} 

u = define_schedule

test_plot = 1
if (test_plot == 1) {

  df = rbind(data.frame(hour = 1:24, eta = 15, beta = 3, effort = u(15, 3)),
             data.frame(hour = 1:24, eta = 19, beta = 3, effort = u(19, 3)),
             data.frame(hour = 1:24, eta = 15, beta = 7, effort = u(15, 7)),
             data.frame(hour = 1:24, eta = 19, beta = 7, effort = u(19, 7)))             
  df$beta = as.factor(df$beta)
  df$eta  = as.factor(df$eta)
  
  PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'
  pdf(file = paste(PLOTS_PATH,'example-fixed-schedules.pdf', sep=''), width = 5, height=3)
  plt = ggplot(df, aes(y = effort, x = hour, color = eta, linetype = beta, shape = eta)) + 
    geom_line(size = 1.3) + geom_point(size = 3)
  plt = plt + theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18, angle = -20),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=14),
          legend.title      = element_text(size=14),
          legend.position  = c(0.25, 0.5),
          axis.ticks = element_blank()) + 
    ggtitle(paste("Effort schedules: examples")) + ylab('Effort') + xlab('Hour of day')
  print(plt)
  dev.off()
}