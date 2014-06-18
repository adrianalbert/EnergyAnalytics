# pecan_ranking_AC.r
#
# Correlate ground truth end-uses with model results. 
#
# Adrian Albert
# Last modified: June 2014.
# ---------------------------------------------------------

rm(list = ls())
options(error = recover)

# ------------------------------------------
# Initializations...
# ------------------------------------------

library('ggplot2')
library('fields')
library('xtable')
source('../../thermal_profiles/validator/metrics.r')
source('../../thermal_profiles/validator/plot_metrics.r')
source('../../utils/select_data.r')

DUMP_PATH = '~/energy-data/pecan_street/models/'
PLOT_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street/ranking/'
dir.create(PLOT_PATH)

# load user names
unames = read.csv('~/energy-data/pecan_street/metadata/user_names_ids.csv')
user_names = unames$name; names(user_names) = unames$ID

# load metrics
load(file = file.path(DUMP_PATH, 'metrics.RData'))
metrics_state = lapply(metrics, function(l) l$response_state)
metrics_state = do.call('rbind', metrics_state)
metrics_prob  = lapply(metrics, function(l) l$response_prob)
metrics_prob  = do.call('rbind', metrics_prob)
metrics_temp  = lapply(metrics, function(l) l$response_temp)
metrics_temp  = do.call('rbind', metrics_temp)

# ___________________________________________________________
# Plots distribution of response at different temperatures

df = subset(metrics_temp, resolution == '60min' & quantity == 'mu' & metric %in% c('a_hat_S', 'a_S', 's_SW'))
df = df[,c('metric', '60', '70', '80', '90')]
df$metric = droplevels(df$metric)
levels(df$metric) = c('Ground Truth', 'Model (summer)', 'Debiased (summer)')
df = melt(df, id.vars = c('metric'))
p = ggplot(df, aes(value, color = metric)) + geom_density(size=2)
p = p + facet_wrap(~variable, scales = 'free')
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

pdf(file=paste(PLOT_PATH, 'density_response_temp.pdf',sep='/'),width=8,height=4.5)
print(p)
dev.off()

# ----------------------------------------------
# Compute rankings of users by ground truth
# ----------------------------------------------

# ____________________________________
# Compute rankings @ 70F

# select first 10 users by a_hat_S
df   = subset(metrics_temp, metric == 'a_hat_S' & resolution == '60min' & quantity == 'mu', select = c('UID', '70'))
df = cbind(UserName = user_names[df$UID], df)
df$UID = NULL
df$UserName = droplevels(df$UserName)
names(df)[2] = 'value'
df = df[with(df, order(value, decreasing=T)), ]
rownames(df) = NULL
names(df)[2] = '$\\hat{a}$@70F'

# dump to latex table
tab <- xtable(df[1:10,])
digits(tab)[3] <- 4
print(tab,floating=FALSE, sanitize.colnames.function = function(x) x,
      file = paste(PLOT_PATH, 'example_ranking_temp_a_hat_S.tex',sep='/'))

# select first 10 users by a_SW
df   = subset(metrics_temp, metric == 's_SW' & resolution == '60min' & quantity == 'mu', select = c('UID', '70'))
df = cbind(UserName = user_names[df$UID], df)
df$UID = NULL
df$UserName = droplevels(df$UserName)
names(df)[2] = 'value'
df = df[with(df, order(value, decreasing=T)), ]
rownames(df) = NULL
names(df)[2] = '$a_{SW}$@70F'

# dump to latex table
tab <- xtable(df[1:10,])
digits(tab)[3] <- 4
print(tab,floating=FALSE, sanitize.colnames.function = function(x) x,
      file = paste(PLOT_PATH, 'example_ranking_temp_a_S.tex',sep='/'))

# ____________________________________
# Compute rankings @ 90F

# select first 10 users by a_hat_S
df   = subset(metrics_temp, metric == 'a_hat_S' & resolution == '60min' & quantity == 'mu', select = c('UID', '90'))
df = cbind(UserName = user_names[df$UID], df)
df$UID = NULL
df$UserName = droplevels(df$UserName)
names(df)[2] = 'value'
df = df[with(df, order(value, decreasing=T)), ]
rownames(df) = NULL
names(df)[2] = '$\\hat{a}$@90F'

# dump to latex table
tab <- xtable(df[1:10,])
digits(tab)[3] <- 4
print(tab,floating=FALSE, sanitize.colnames.function = function(x) x,
      file = paste(PLOT_PATH, 'example_ranking_temp_a_hat_S_90.tex',sep='/'))

# select first 10 users by a_SW
df   = subset(metrics_temp, metric == 's_SW' & resolution == '60min' & quantity == 'mu', select = c('UID', '90'))
df = cbind(UserName = user_names[df$UID], df)
df$UID = NULL
df$UserName = droplevels(df$UserName)
names(df)[2] = 'value'
df = df[with(df, order(value, decreasing=T)), ]
rownames(df) = NULL
names(df)[2] = '$a_{SW}$@90F'

# dump to latex table
tab <- xtable(df[1:10,])
digits(tab)[3] <- 4
print(tab,floating=FALSE, sanitize.colnames.function = function(x) x,
      file = paste(PLOT_PATH, 'example_ranking_temp_a_S_90.tex',sep='/'))

# -------------------------------------------------------
# Compare rankings by different measures
# -------------------------------------------------------
# _________________________________________________________________
# Compute standard correlations between rankings: avg. response
# compute Kendall tau, Pearson rho, along with standard errors of these rankings

rank_corr = lapply(50:91, function(t){
  result = lapply(c('15min', '60min'), function(r) {
    df = subset(metrics_temp, resolution == r & quantity == 'mu' & metric %in% c('a_hat_S', 'a_S', 's_SW'))
    df = df[,c('UID', 'metric', as.character(t))]
    a_hat_S = subset(df, metric == 'a_hat_S')[,as.character(t)]
    a_SW    = subset(df, metric == 's_SW')[,as.character(t)]
    a_S     = subset(df, metric == 'a_S')[,as.character(t)]
    
    rho.SW = cor.test(a_hat_S, a_SW, method = 'spearman', exact = T)
    tau.SW = cor.test(a_hat_S, a_SW, method = 'kendall', exact = T)
    rho.S  = cor.test(a_hat_S, a_S, method = 'spearman', exact = T)
    tau.S  = cor.test(a_hat_S, a_S, method = 'kendall', exact = T)
    m.SW = c(m = rho.SW$estimate, 
             p.rho = rho.SW$p.value,
             m = tau.SW$estimate,
             p.tau = tau.SW$p.value)
    m.S  = c(m = rho.S$estimate, 
             p.rho = rho.S$p.value,
             m = tau.S$estimate,
             p.tau = tau.S$p.value)
    df = cbind(Estimate = c('a_SW', 'a_S'), as.data.frame(rbind(m.SW, m.S)))
    return(df)
  })  
  result = as.data.frame(do.call('rbind', result))
  result = as.data.frame(cbind(Temperature = t, res = rep(c('15min', '60min'), each=2), result))
  return(result)
})
rank_corr = do.call('rbind', rank_corr)

# ______________________________________________________________________
# Compute Kendall tau between rankings: avg. response top q%

rank_corr_top = lapply(50:91, function(t){
  
  cat(paste('t=', t, '\n'))
  
  df = subset(metrics_temp, 
              quantity == 'mu' & resolution == '60min' & metric %in% c('a_hat_S', 'a_S', 's_SW'))
  df = df[,c('UID', 'metric', as.character(t))]
  a_hat_S = subset(df, metric == 'a_hat_S')[,as.character(t)]
  a_SW    = subset(df, metric == 's_SW')[,as.character(t)]
  a_S     = subset(df, metric == 'a_S')[,as.character(t)]
  result = lapply(seq(10, 90, length.out = 9), function(q) {
    n_q = trunc(length(a_hat_S) * q / 100)
    idx_SW  = order(a_SW, decreasing = T)[1:n_q]
    idx_S   = order(a_S, decreasing = T)[1:n_q]
    
    tau.SW = cor.test(a_hat_S[idx_SW], a_SW[idx_SW], method = 'kendall', exact = T)
    tau.S  = cor.test(a_hat_S[idx_S], a_S[idx_S], method = 'kendall', exact = T)
    m.SW = c(m = tau.SW$estimate,
             p.tau = tau.SW$p.value)
    m.S  = c(m = tau.S$estimate,
             p.tau = tau.S$p.value)
    df = cbind(Estimate = c('a_SW', 'a_S'), Top = q, as.data.frame(rbind(m.SW, m.S)))
  })

  result = as.data.frame(do.call('rbind', result))
  result = as.data.frame(cbind(Temperature = t, result))
  return(result)
})
rank_corr_top = do.call('rbind', rank_corr_top)

# ----------------------------------------------
# Plots
# ----------------------------------------------

# ________________________________________________________
# Plot ranking match: average response

# overall
df = melt(rank_corr[,-c(4,5,7)], id.vars = c('Temperature', 'res', 'Estimate'))
df$variable = 'Kendall tau'
p = ggplot(df, aes(Temperature, value, color = res))
p = p + geom_line(aes(linetype = Estimate), size=2) + geom_point(aes(shape = Estimate), size=4)
p = p + facet_wrap(~variable, scales = 'free')
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
p = p + ggtitle('Kendall tau rank correlation statistic') + 
  xlab('Temperature [deg F]') + ylab('tau')
pdf(file=paste(PLOT_PATH, 'ranking_compare_mu_temp.pdf',sep='/'),width=8,height=4.5)
print(p)
dev.off()

# top k%
df = rank_corr_top[,-5]
df$variable = 'Kendall tau (top k%)'
p = ggplot(df, aes(Temperature, m.tau, color = Estimate))
p = p + geom_line(aes(linetype = Estimate), size=2) + geom_point(aes(shape = Estimate), size=4)
p = p + facet_wrap(~Top, scales = 'free')
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
p = p + ggtitle('Kendall tau rank correlation statistic (top k%)') + 
  xlab('Temperature [deg F]') + ylab('tau')
pdf(file=paste(PLOT_PATH, 'ranking_compare_mu_temp_top_k.pdf',sep='/'),width=8,height=4.5)
print(p)
dev.off()


