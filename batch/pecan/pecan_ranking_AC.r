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
source('../../thermal_profiles/validator/drank.r')
source('../../utils/select_data.r')

DUMP_PATH = '~/energy-data/pecan_street/models/'
PLOT_PATH = '~/Dropbox/OccupancyStates/plots/pecan-street/ranking/'
dir.create(PLOT_PATH)

# load user names
unames = read.csv('~/energy-data/pecan_street/metadata/user_names_ids.csv')
user_names = unames$name; names(user_names) = unames$ID

# load weather data 
weather = read.csv('~/energy-data/pecan_street/weather//weather_hourly.csv')
tempF.grid = quantile(weather$TemperatureF, probs = seq(0, 1, length.out = 5))

# load metrics
load(file = file.path(DUMP_PATH, 'metrics.RData'))
metrics_state = lapply(metrics, function(l) l$response_state)
metrics_state = do.call('rbind', metrics_state)
metrics_prob  = lapply(metrics, function(l) l$response_prob)
metrics_prob  = do.call('rbind', metrics_prob)
metrics_temp  = lapply(metrics, function(l) l$response_temp)
metrics_temp  = do.call('rbind', metrics_temp)

# ----------------------------------------------
# Plots of basic metrics
# ----------------------------------------------

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

# __________________________________________________________________________
# Plots distribution of absolute relative error>10%, 30%, 100%, 50% at T=90

df = subset(metrics_prob, 
            resolution == '60min' & metric == 'P_rel.abs_SW')
df = df[,c('eta', '90')]
df$eta = round(df$eta * 100, digits = 0)
names(df) = c('e', 'T90F')
p = ggplot(df, aes(as.factor(e), T90F)) + geom_boxplot()
p = p + scale_x_discrete(labels = seq(0, 100, length.out = 20), breaks = seq(0, 100, length.out = 20))
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
p = p + ggtitle(bquote(P(abs(e)<eta) ~ ' distribution @ T=90F'))
p = p + xlab(bquote(eta ~ '[%]')) + ylab(bquote(P(abs(e)<eta)))

pdf(file=paste(PLOT_PATH, 'density_error_temp_90.pdf',sep='/'),width=8,height=4.5)
print(p)
dev.off()

# __________________________________________________________________________
# Plots distribution of absolute relative error>10%, 30%, 100%, 50% at T=70

df = subset(metrics_prob, 
            resolution == '60min' & metric == 'P_rel.abs_SW')
df = df[,c('eta', '70')]
df$eta = round(df$eta * 100, digits = 0)
names(df) = c('e', 'T70F')
p = ggplot(df, aes(as.factor(e), T70F)) + geom_boxplot()
p = p + scale_x_discrete(labels = seq(0, 100, length.out = 20), breaks = seq(0, 100, length.out = 20))
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
p = p + ggtitle(bquote(P(abs(e)<eta) ~ ' distribution @ T=70F'))
p = p + xlab(bquote(eta ~ '[%]')) + ylab(bquote(P(abs(e)<eta)))

pdf(file=paste(PLOT_PATH, 'density_error_temp_70.pdf',sep='/'),width=8,height=4.5)
print(p)
dev.off()

# __________________________________________________________________________
# Plots distribution of magnitude error at T = 90

df = subset(metrics_prob, 
            resolution == '60min' & metric == 'P_magnitude_SW')
df = df[,c('eta', '90')]
df$eta = round(df$eta, digits = 2)
names(df) = c('delta', 'T90F')
p = ggplot(df, aes(as.factor(delta), T90F)) + geom_boxplot()
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
p = p + ggtitle(bquote(P(delta<eta) ~ ' distribution @ T=90F'))
p = p + xlab(bquote(eta ~ '[%]')) + ylab(bquote(P(abs(e)<eta)))

pdf(file=paste(PLOT_PATH, 'density_error_temp_mag_90.pdf',sep='/'),width=8,height=4.5)
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
    a_hat_S = subset(df, metric == 'a_hat_S')[, as.character(t)]
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

# if pval > 0.05 reject H0 
idx.p = which(rank_corr_top$p.tau > 0.05)
if (length(idx.p)>0) rank_corr_top$m.tau[idx.p] = 0

# __________________________________________________________________________
# Compute rank distance d_rank for for temperature ranges for avg. response

names_ranges = c('Low', 'Medium-Low', 'Medium-High', 'High')
ranges = round(tempF.grid)
ranges[which(ranges>91)] = 91
d_rank_corr_top = lapply(1:(length(ranges)-1), function(r){
  
  cur_range = seq(ranges[r], ranges[r+1])
  r_name    = names_ranges[r]
  cat(paste('T=', paste(ranges[r:(r+1)], collapse = '...'), '\n'))
  
  df = subset(metrics_temp, 
              quantity == 'mu' & resolution == '60min' & metric %in% c('a_hat_S', 'a_S', 's_SW'))
  df = df[,c('UID', 'metric', as.character(cur_range))]  
  a_hat_S = t(as.matrix(subset(df, metric == 'a_hat_S')[,as.character(cur_range)]))
  a_SW    = t(as.matrix(subset(df, metric == 's_SW')[,as.character(cur_range)]))
  a_S     = t(as.matrix(subset(df, metric == 'a_S')[,as.character(cur_range)]))
  y       = colMeans(a_hat_S)
  idx.y   = order(y, decreasing = T)

  result = lapply(seq(10, 100, length.out = 10), function(q) {
    n_q = trunc(length(y) * q / 100)        
    X_SW= a_SW[,idx.y[1:n_q]]
    X_S = a_S[,idx.y[1:n_q]]
  
    res.SW = compute_drank(y[1:n_q], X_SW, B = 500); 
    res.S  = compute_drank(y[1:n_q], X_S, B = 500); 
    
    df = cbind(Estimate = c('a_SW', 'a_S'), Top = q, as.data.frame(rbind(res.SW, res.S)))
  })
  
  result = as.data.frame(do.call('rbind', result))
  result = as.data.frame(cbind(Range = r_name, result))
  return(result)
})
d_rank_corr_top = do.call('rbind', d_rank_corr_top)


# ________________________________________________________________________________
# Compute rank distance d_rank for for temperature ranges for P(error) (magnitude)

ranges = round(tempF.grid)
ranges[which(ranges>91)] = 91
df_prob = subset(metrics_prob, resolution == '60min' & 
              metric %in% c('P_magnitude_SW_own', 'P_magnitude_Shat_own', 'P_magnitude_S_own'))
eta_grid = c(-Inf, seq(-0.1, 0.1, length.out = 30))
uids = as.character(unique(metrics_prob$UID))
d_rank_corr_prob_mag = lapply(1:(length(ranges)-1), function(r){
   
  cur_range = seq(ranges[r], ranges[r+1])
  r_name    = names_ranges[r]
  cat(paste('T=', paste(ranges[r:(r+1)], collapse = '...'), '\n'))
  df = df_prob[,c('UID', 'eta', 'metric', as.character(cur_range))]  
  mat = array(0, c(3,length(uids), length(eta_grid), length(cur_range)))
  
  result = lapply(1:(length(eta_grid)-1), function(q) {
    a_hat_S = subset(df, eta >= eta_grid[q] & eta < eta_grid[q+1] & metric == 'P_magnitude_Shat_own')[,c('UID', 'eta', as.character(cur_range))]    
    if (nrow(a_hat_S) > 0) {
      names(a_hat_S)[-1] = paste('T', names(a_hat_S)[-1], sep = '_')
      a_hat_S$eta = eta_grid[q+1]
      a_hat_S = aggregate(as.formula(paste('cbind(', paste('T_', cur_range, collapse = ',', sep=''), ')~UID+eta')), data = a_hat_S, FUN = mean)
      ids      = which(uids %in% as.character(a_hat_S$UID))
      mat[1,ids,q,] = as.matrix(a_hat_S[,-c(1,2)])
    }
    a_S = subset(df, eta >= eta_grid[q] & eta < eta_grid[q+1] & metric == 'P_magnitude_S_own')[,c('UID', 'eta', as.character(cur_range))]    
    if (nrow(a_S) > 0) {
      names(a_S)[-1] = paste('T', names(a_S)[-1], sep = '_')
      a_S$eta = eta_grid[q+1]
      a_S = aggregate(as.formula(paste('cbind(', paste('T_', cur_range, collapse = ',', sep=''), ')~UID+eta')), data = a_S, FUN = mean)
      ids      = which(uids %in% as.character(a_S$UID))
      mat[2,ids,q,] = as.matrix(a_S[,-c(1,2)])
    }  
    a_SW = subset(df, eta >= eta_grid[q] & eta < eta_grid[q+1] & metric == 'P_magnitude_SW_own')[,c('UID', 'eta', as.character(cur_range))]    
    if (nrow(a_SW) > 0) {
      names(a_SW)[-1] = paste('T', names(a_SW)[-1], sep = '_')
      a_SW$eta = eta_grid[q+1]
      a_SW     = aggregate(as.formula(paste('cbind(', paste('T_', cur_range, collapse = ',', sep=''), ')~UID+eta')), data = a_SW, FUN = mean)
      ids      = which(uids %in% as.character(a_SW$UID))
      mat[3,ids,q,] = as.matrix(a_SW[,-c(1,2)])
    }
    
    y   = colMeans(mat[1,,q,])
    X_SW= mat[3,,q,]
    X_S = mat[2,,q,]
    
    res.SW = compute_drank(1-y, 1-X_SW, B = 500); 
    res.S  = compute_drank(1-y, 1-X_S, B = 500); 
    
    df = cbind(Estimate = c('a_SW', 'a_S'), eta = eta_grid[q+1], as.data.frame(rbind(res.SW, res.S)))
  })
  
  result = as.data.frame(do.call('rbind', result))
  result = as.data.frame(cbind(Range = r_name, result))
  return(result)
})
d_rank_corr_prob_mag = do.call('rbind', d_rank_corr_prob_mag)

# -------------------------------------------------------
# Compare aggregate statistics of rankings: model vs GT
# -------------------------------------------------------

# compute aggregate of top q% users by different metrics
compute_agg_stats = function(df, q = 100) {
  df.mu = cast(data = df[,-4], UID ~ metric)
  df.sd = cast(data = df[,-3], UID ~ metric)
  n_q   = trunc(nrow(df.mu) * q / 100)      
  l = list()
  for (i in setdiff(names(df.mu), 'UID')) {
    idx   = order(df.mu[,i], decreasing = T)
    mu.ag = sum(df.mu[idx[1:n_q], i])
    sd.ag = sqrt(sum(df.sd[idx[1:n_q], i]^2))
    l[[i]] = c(mu = mu.ag, sd = sd.ag)
  }
  l.tot = cbind(estimate = setdiff(names(df.mu), 'UID'), as.data.frame(do.call('rbind', l)))
  l.avg = l.tot; l.avg[,2] = l.avg[,2] / n_q; l.avg[,3] = l.avg[,3] / sqrt(n_q); 
  return(list(tot = l.tot, avg = l.avg))
}

# ________________________________________________________________________________
# Compare group aggregate flexibiliy by different estimates/temperature ranges

flex_agg = lapply(51:91, function(t){
  
  cat(paste('T=', t, '\n'))
  
  df = subset(metrics_temp, resolution == '60min' & metric %in% c('a_hat_S', 'a_S', 's_SW'))
  df = df[,c('UID', 'metric', 'quantity', as.character(t))]  
  df = cast(data = df, UID + metric ~ quantity)
  
  
  result = lapply(seq(10, 100, length.out = 10), function(q) {    
    res = compute_agg_stats(df, q = q)          
    ret = rbind(cbind(type = 'total', res$tot), cbind(type = 'average', res$avg))
    ret$Temperature = t
    ret$TopPercent = q
    return(ret)
  })
  
  result = as.data.frame(do.call('rbind', result))
  return(result)
})
flex_agg = do.call('rbind', flex_agg)

# ________________________________________________________________________________
# Compute metrics of validation for the top k% profiles

f_vec_list = function(a) {
  lapply(1:nrow(a), function(i) { z = c(a[['mu']][i], a[['sd']][i]); names(z) = c("mu", "sd"); z})
}
  
l = list()
for (q in 10*(1:10)) {
  a_hat_S = subset(flex_agg, type == 'average' & TopPercent == q & estimate == 'a_hat_S' & Temperature < 91,
                   select = c('mu', 'sd', 'Temperature'))
  a_S = subset(flex_agg, type == 'average' & TopPercent == q & estimate == 'a_S' & Temperature < 91,
                   select = c('mu', 'sd', 'Temperature'))
  a_SW = subset(flex_agg, type == 'average' & TopPercent == q & estimate == 's_SW' & Temperature < 91,
                   select = c('mu', 'sd', 'Temperature'))
  metr.temp.S  = compute_comparisons_states(f_vec_list(a_hat_S), f_vec_list(a_S))
  metr.temp.SW = compute_comparisons_states(f_vec_list(a_hat_S), f_vec_list(a_SW))
  
  l[[as.character(q)]] = list(S = metr.temp.S, SW = metr.temp.SW)
}

temp.grid = rev(rev(unique(flex_agg$Temperature))[-1])

# ----------------------------------------------
# Plots: ranking
# ----------------------------------------------

# ________________________________________________________
# Plot ranking match: average response

# overall
df = rank_corr[,-c(4,5,7)]
df$Range = names_ranges[sapply(df$Temperature, function(t) findInterval(t, tempF.grid, all.inside = T))]
df$Range = factor(df$Range, levels = names_ranges)
df$Temperature = NULL
p = ggplot(df, aes(x = Range, y = m.tau, color = Estimate))
# p = p + geom_line(aes(linetype = Estimate), size=2) + geom_point(aes(shape = Estimate), size=4)
p = p + geom_boxplot(size=2)
p = p + facet_wrap(~res, scales = 'free')
p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.title.x     = element_text(size = 18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18, angle = 30),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),        
        legend.title     = element_text(size=18),    
        axis.ticks       = element_blank() ) + 
  theme(plot.title=element_text(family="Times", face="bold", size=20))
p = p + ggtitle('Kendall tau rank correlation statistic') + 
  xlab('Temperature Ranges') + ylab('tau')
pdf(file=paste(PLOT_PATH, 'ranking_compare_mu_temp.pdf',sep='/'),width=10,height=4.5)
print(p)
dev.off()

# top k%: kendall tau
df = rank_corr_top[,-5] #subset(rank_corr_top[,-5], Estimate == 'a_SW'); df$Estimate = NULL
df$Range = names_ranges[sapply(df$Temperature, function(t) findInterval(t, tempF.grid, all.inside = T))]
df$Range = factor(df$Range, levels = names_ranges)
df$Temperature = NULL
p = ggplot(df, aes(Range, m.tau, color = Estimate))
p = p + facet_wrap(~Top, scales = 'fixed')
p = p + geom_boxplot(size=2)
p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.title.x     = element_text(size = 18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18, angle = 30),
        axis.title.y     = element_text(size=18),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),        
        legend.title     = element_text(size=18),    
        axis.ticks       = element_blank() ) + 
  theme(plot.title=element_text(family="Times", face="bold", size=20))
p = p + ggtitle('Kendall tau rank correlation statistic (top k%)') + 
  xlab('Temperature [deg F]') + ylab('tau')
pdf(file=paste(PLOT_PATH, 'ranking_compare_mu_temp_top_k.pdf',sep='/'),width=10,height=6)
print(p)
dev.off()

# which rankings may be similar?
d_zero = subset(d_rank_corr_top, p > 0.1)
rownames(d_zero) = NULL; names(d_zero)[1] = 'Temperature'
tab <- xtable(d_zero[,-4])
print(tab,floating=FALSE,
      file = paste(PLOT_PATH, 'performance_drank.tex',sep='/'))

# ________________________________________________________
# Plot ranking match: probability of response

df = d_rank_corr_prob_mag[,-4] #subset(rank_corr_top[,-5], Estimate == 'a_SW'); df$Estimate = NULL
df$eta = round(df$eta, digits = 3)
df = subset(df, p > 0.1 & Estimate == 'a_SW' & eta %in% c(0.003, 0.024, 0.038, 0.052, 0.072, 0.1))
p = ggplot(df, aes(Range, p, color = as.factor(eta), shape = as.factor(eta), group = as.factor(eta)))
p = p + geom_point(size=4) + geom_line(size=2)
# p = p + facet_wrap(~Ranges, scales = 'fixed')
p = p + theme_bw() + #scale_y_continuous(limits = c(0, 100))
  theme(panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.x     = element_text(size=18),
        axis.title.x     = element_text(size = 18),
        axis.text.y      = element_text(size=18), 
        axis.text.x      = element_text(size=18, angle = 30),
        axis.title.y     = element_text(size=20),
        axis.title.x     = element_text(size=18),
        plot.title       = element_text(size=20),            
        legend.text      = element_text(size=18),        
        legend.title     = element_text(size=18),    
        axis.ticks       = element_blank() ) + 
  theme(plot.title=element_text(family="Times", face="bold", size=20))
p = p + scale_color_discrete(name = 'Threshold')
p = p + scale_shape_discrete(name = 'Threshold')
p = p + ggtitle('p-value for rank distance hypothesis test (magnitude)') + 
  xlab('Temperature [deg F]') + ylab('p')
pdf(file=paste(PLOT_PATH, 'ranking_compare_prob_mag_temp.pdf',sep='/'),width=9.3,height=4)
print(p)
dev.off()

# ________________________________________________________
# Plot aggregate response profiles

df = subset(flex_agg, type == 'average' & TopPercent %in% c(10, 20, 30) & Temperature < 91)
p = ggplot(df, aes(Temperature, mu, color = estimate))
p = p + geom_line(size = 1.5) + geom_point(aes(shape = estimate), size = 3.5)
p = p + geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.1)  
p = p + facet_wrap(~TopPercent, scales = 'fixed')
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
        legend.position  = c(0.8,0.7),
        axis.ticks       = element_blank() ) + 
  theme(plot.title=element_text(family="Times", face="bold", size=20))
p = p + ggtitle('Avg. aggregate flexibility of top k%')  + xlab('Temperature [deg F]') + ylab('Avg. Response [kWh/deg F]')

pdf(file=paste(PLOT_PATH, 'avg_flex_group_top20.pdf',sep='/'),width=13,height=4)
print(p)
dev.off()

# ________________________________________________________
# Plot aggregate response profiles

#top 10%
pdf(file=paste(PLOT_PATH, 'avg_flex_group_top10_prob_profile.pdf',sep='/'),width=9.3,height=4)
print(image.plot(temp.grid[-(1:10)], l[['10']]$SW$support.rab, l[['10']]$SW$P.rel.abs[-(1:10),],
                 xlab = 'Temperature [deg F]', 
                 ylab = expression(eta), 
                 main = bquote(P(abs(e)<eta) ~ ': Aggregate error probability profile (top 10% users)')))
dev.off()

#top 30%
pdf(file=paste(PLOT_PATH, 'avg_flex_group_top30_prob_profile.pdf',sep='/'),width=9.3,height=4)
print(image.plot(temp.grid[-(1:10)], l[['30']]$SW$support.rab, l[['30']]$SW$P.rel.abs[-(1:10),],
                 xlab = 'Temperature [deg F]', 
                 ylab = expression(eta), 
                 main = bquote(P(abs(e)<eta) ~ ': Aggregate error probability profile (top 30% users)')))
dev.off()

