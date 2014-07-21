# #########################################################################
# cluster_profiles.r
# -------------------------
#
# Clustering analysis for Bakersfield profiles data.
# 
# Adrian Albert
# Last modified: July 2014.
# #########################################################################

rm(list=ls())
options(error = recover)

#______________________________________________________________
# Prepare data

# load data
load('~/energy-data/bakersfield_profiles.RData')
library('pracma')
library('ggplot2')
library('reshape')

PLOTS_PATH = '~/Dropbox/OccupancyStates/plots/scheduling/'

source('../../clustering/kError.r', chdir = T)

cur_inputs = inputs#[1:400]

# format inputs
Abar = do.call('rbind', lapply(cur_inputs, function(l) as.numeric(l$a$mu[,1])))
Abar[which(Abar < 0)] = 0
rownames(Abar)= names(cur_inputs)
W = lapply(cur_inputs, function(l) diag(diag(l$a$covmat)))
# W = lapply(cur_inputs, function(l) l$a$covmat)

#______________________________________________________________
# Cluster profiles

res = mclapply(c(1:20)*2, 
               mc.cores = 5, 
               function(k) {
                 cat(paste('k=', k, '\n'))
                 done = F; att = 0
                 while (!done & att < 5) {
                   fit = try(kError(Abar, W, k, iter = 100))
                   done= class(fit) != 'try-error'
                   if (!done) {
                     print(paste('Error at K=', k, '... re-trying', att, '/5'))
                     att = att + 1
                   }
                 }
                 return(fit)
               })

obj = sapply(1:length(res), function(k) rev(res[[k]]$objective)[1])
cl = res[[10]]; K = 20

#______________________________________________________________
# Plots for clustering analysis: objective trace vs K

pdf(file = paste(PLOTS_PATH,'number-of-clusters.pdf', sep=''), width = 6, height=4)
plotyy(c(1:20)*2, obj/1000, col.y1 = 'black', col.y2 = 'darkred',
       c(2:20)*2, abs(diff(obj)) / obj[-1],
       lwd = 4, cex.main = 2, cex.lab = 1.5, type = 'b', cex.axis = 1.5,
       main = 'Choosing the number of segments', 
       ylab = 'Objective (x1000)', xlab = '# clusters K')
abline(v = 20, lwd = 4, lty = 3, col = 'darkgrey')
legend(17, 100, c('Obj(K)', '[Diff. Obj(K)]/Obj(K)'), fill = c('black', 'darkred'))
dev.off()

# based on the above analysis, a good number of clusters is K = 12

#______________________________________________________________
# Plots for clustering analysis: centers

# extract error diagonals
sd = as.data.frame(t(sapply(cl$errors, function(x) sqrt(abs(diag(x))))))
names(sd) = 1:ncol(sd)
sd$Segment = sprintf("%02d", 1:nrow(sd))
sd = melt(sd, id.vars = 'Segment')
names(sd)[3] = 'sd'

# extract centers
df = as.data.frame(cl$centers)                            
names(df) = 1:ncol(df)
df$Segment = sprintf("%02d", 1:nrow(df))
df = melt(df, id.vars = 'Segment')
names(df)[3] = 'mu'

# form plotting data
dfp = merge(df, sd, by = c('Segment', 'variable'))
dfp$variable = as.numeric(dfp$variable)

# percentages of membership              
tab = table(cl$assignment)
tab = round(tab / sum(tab), digits = 4)
tab = tab * 100
dfp$Segment = paste(dfp$Segment, ': ', tab[as.numeric(dfp$Segment)], '%', sep = '')

plt = ggplot(dfp, aes(y = mu, x = variable)) + 
  geom_point(size = 2.5) + geom_line(size=1.5)
plt = plt + geom_errorbar(aes(ymax = mu + sd, ymin = mu - sd))        
plt = plt + facet_wrap(~Segment, nrow = 4, scales = 'free')
plt = plt + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),
        legend.position  = "none",
        axis.ticks = element_blank()) + 
  ggtitle(paste("Segments")) + ylab('Response [kWh/F]') + xlab('Time of day')

pdf(file = paste(PLOTS_PATH, 'cluster_centers.pdf', sep=''), width = 14, height=6)
print(plt)
dev.off()

#______________________________________________________________
# Plots for clustering analysis: cluster membership

# extract cluster members
df = lapply(1:nrow(cl$centers), function(k) {
  x_k = Abar[which(cl$assignment == k), ]
  x_k = as.data.frame(x_k); 
  names(x_k) = 1:ncol(x_k)
  x_k$Segment = sprintf("%02d", k)
  x_k = melt(x_k, id.vars = 'Segment')
  names(x_k)[3] = 'mu'  
  return(x_k)
})
df = do.call('rbind', df)
df$variable = as.numeric(df$variable)
df$Segment = paste(df$Segment, ': ', tab[as.numeric(df$Segment)], '%', sep = '')

plt = ggplot(df, aes(y = mu, x = variable, group = variable)) + geom_boxplot(size=0.7)
plt = plt + facet_wrap(~Segment, nrow = 4, scales = 'free')
plt = plt + scale_x_continuous(breaks=c(1, 6, 12, 18, 24))
plt = plt + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),
        legend.position  = "none",
        axis.ticks = element_blank()) + 
  ggtitle(paste("Segments")) + ylab('Response [kWh/F]') + xlab('Time of day')

pdf(file = paste(PLOTS_PATH, 'cluster_members.pdf', sep=''), width = 16, height=6.5)
print(plt)
dev.off()

save.image('~/energy-data/bakersfield_clustering.RData')
