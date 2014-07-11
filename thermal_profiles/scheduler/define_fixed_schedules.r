
library('ggplot2')

define_schedule = function(gamma, eta, beta, tau = 24) {
  
  u = rep(0, tau);
  u[eta:min(eta+gamma, 24)] = beta / gamma
  return(u)
} 

u = define_schedule

test_plot = 1
if (test_plot == 1) {

  beta = 3
  df = rbind(data.frame(hour = 1:24, eta = 15, gamma = 3, effort = u(3,15, beta)),
             data.frame(hour = 1:24, eta = 18, gamma = 3, effort = u(3,18, beta)),
             data.frame(hour = 1:24, eta = 15, gamma = 5, effort = u(5,15, beta)),
             data.frame(hour = 1:24, eta = 18, gamma = 5, effort = u(5,18, beta)))             
  df$gamma = as.factor(df$gamma)
  df$eta = as.factor(df$eta)
  
  PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'
  pdf(file = paste(PLOTS_PATH,'example-fixed-schedules.pdf', sep=''), width = 5, height=3)
  plt = ggplot(df, aes(y = effort, x = hour, color = eta, linetype = gamma, shape = gamma)) + 
    geom_line(size = 1.5) + geom_point(size = 4)
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