# plot_utils.r
#
# Plots for HMM/OLS analysis.
# 
# Adrian Albert
# Last modified: December 2012.
# -----------------------------------------------------------------------

library(igraph)
library('useful')
library('grid')
library('reshape2')
library('RColorBrewer')

# __________________________________________
# Plots time series HMM color-coded by state

plot_hmm_ts = function(hmm.means, hmm.sigma, states, timestamps, observed, y.lab = 'observed', title = 'HMM-ts') {
	# construct plotting data frames
	df.hmm = data.frame(Mean   = hmm.means, 
		          Sigma  = hmm.sigma, 
		          State  = states,      
			  Index  = 1:length(states),				                              
		          Time   = format(timestamps, "%a,%m/%d %H:00"))              
	df.obs = data.frame(Observed = observed, 
		          Time   = format(timestamps, "%a,%m/%d %H:00"),
			  Index  = 1:length(states),				                              
		          Sigma    = rep(0, length(observed)))

	df.mlt.obs = melt(df.obs, id.vars = c('Time', 'Sigma', 'Index'))
	df.mlt.hmm = melt(df.hmm, id.vars = c('Time', 'Sigma', 'State', 'Index'), measure.vars = c('Mean'))
	df.mlt.hmm$variable = paste(df.mlt.hmm$variable, df.mlt.hmm$State)
	df.mlt.hmm$State = NULL

	# plot current zoom-in
	pd <- position_dodge(.1)
	plt = ggplot(df.mlt.hmm, aes(x = Index, y = value)) + 
	geom_errorbar(aes(ymin = value-Sigma, ymax=value+Sigma), 
		      colour="black", width=.1, position=pd, alpha=0.2) +
	geom_line(position=pd) + geom_point(position=pd, aes(colour=variable)) + 
	geom_line(col='black', data = df.mlt.obs, alpha = 0.4) + 
	geom_point(col='black', data = df.mlt.obs, alpha = 0.4) + 
	scale_x_continuous(breaks = seq(1,length(df.mlt.hmm$Time), length.out=10), 
		           labels = format(df.mlt.hmm$Time[seq(1,length(df.mlt.hmm$Time), length.out=10)], 
		                           format = "%a,%m/%d %H:00")) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	       panel.grid.minor = element_blank(),
	       panel.background = element_blank(),
	       axis.title.x = element_blank(),
	       panel.background = element_rect(fill = "transparent",colour = NA),
	       axis.ticks = element_blank() ) + 
	ylab(y.lab) + xlab('Time') + 
	theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
	ggtitle( title )

	return(plt)
}

# ________________________________________
# Plot underlying MC of a HMM using ggplot

plot_HMM_MC = function(hmm.means, hmm.sigma, transition, title = 'HMM MC Structure') {
	nStates = length(hmm.means)
  	df = data.frame(Mean   = hmm.means,
  			Sigma  = hmm.sigma,
		        State  = 1:nStates)
	state.grid = expand.grid(State.i = as.factor(1:1:nStates), 
		               State.j = as.factor(1:1:nStates))
	trans = sapply(1:nrow(state.grid), function(s) {
	i = state.grid[s,1]
	j = state.grid[s,2]
	tr= paste(i,'->', j, sep='')
	ret = transition[1,tr]
	})
		    
	df.P = state.grid
	df.P$P = trans
	df.P$Mean.i  = as.numeric(as.character(df.P$State.i))
	df.P$Mean.j  = as.numeric(as.character(df.P$State.j))
	df.P$Sigma.i = as.numeric(hmm.sigma[df.P$State.i])
	df.P$Sigma.j = as.numeric(hmm.sigma[df.P$State.j])
	
	plt = ggplot(df) + 
	geom_point(aes(x = State, y = Sigma, color = State, size = Mean)) + 
	scale_size(range = c(2, 12)) + 
	geom_segment(aes(x = Mean.i, xend = Mean.j, 
		         y = Sigma.i, yend = Sigma.j, 
		         size=P, alpha=P), data=df.P) +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	     panel.grid.minor = element_blank(),
	     panel.background = element_blank(),
	#        axis.title.y = element_blank(),
	#        axis.title.x = element_blank(),
	     legend.position = "none",
	     panel.background = element_blank(),
	     axis.ticks = element_blank()) + 
	ylab('Sigma') + xlab('Mean') + 
	ggtitle(title)
	return(plt)
}

# _________________________________
# Confidence Interval plot for HMM

plot_HMM_CI = function(p.vals, states, timestamps, title = 'Confidence Levels') {

	df = data.frame( Prob = p.vals, State = as.factor(states), 
		       Weekday = weekdays(timestamps), Hour = hour(timestamps) )
	df.st = aggregate(Prob ~ State, FUN = mean, data = df)
	names(df.st) = c('State')
	df.mu = aggregate(Prob ~ Weekday + Hour, FUN = mean, data = df)
	names(df.mu) = c('Weekday', 'Hour', 'Mean')
	df.sd = aggregate(Prob ~ Weekday + Hour, FUN = sd, data = df)
	names(df.sd) = c('Weekday', 'Hour', 'Sigma')
	df.sd.mu = merge(df.sd, df.mu)

	pd <- position_dodge(.1)
	plt = ggplot(df.sd.mu, aes(x = Hour, y = Mean, color = Weekday)) + 
	geom_errorbar(aes(ymin = Mean-Sigma, ymax=Mean+Sigma), 
		      colour="black", width=.1, position=pd, alpha=0.2) +
	geom_line(position=pd) + geom_point(position=pd) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	     panel.grid.minor = element_blank(),
	     panel.background = element_blank(),
	     panel.background = element_rect(fill = "transparent",colour = NA),
	     axis.ticks = element_blank()) + 
	     ylab('Confidence') + xlab('Hour') + 
	ggtitle( title )

	return(plt)              
}

# ______________________________________
# Plot MC structure of HMM using igraph

plot_HMM_MC_NET = function(hmm.means, hmm.sigma, transition) {

      nStates = length(hmm.means)
      # organize MC information
      df = data.frame(Mean   = hmm.means,
                      Sigma  = hmm.sigma,
                      State  = nStates)              
      state.grid = expand.grid(State.i = 1:nStates, State.j = 1:nStates)
      trans = sapply(1:nrow(state.grid), function(s) {
        i = state.grid[s,1]
        j = state.grid[s,2]
        tr= paste(i,'->', j, sep='')
        ret = transition[1,tr]
      })
                
      # form igraph object
      g = graph.empty()              
      g = add.vertices(g, nrow(df), State=as.character(df$State), Mean =df$Mean, Sigma = df$Sigma)              
      g = add.edges(g, t(as.matrix(state.grid)), P = trans)                            
      
      # add graph plotting attributes
      scale <- function(v, a, b) {
        v <- v-min(v) ; v <- v/max(v) ; v <- v * (b-a) ; v+a
      }              
      pallete    = colorRampPalette(brewer.pal(9,"Blues"))(100)
      V(g)$color = 'grey'              
      V(g)$size  <- scale(V(g)$Mean, 10, 20)
      E(g)$color <- pallete[scale(E(g)$P, 1,100)]
      
      # plot
      return(g)             
}

# ________________________________
# Residual analysis for HMM

plot_HMM_res = function(hmm.residual, states, hmm.sigma, sw_test = NULL) {

      # densities of residuals
      df = data.frame(State = factor(states), Residual = hmm.residual)
      plt1 = ggplot(df) + 
        geom_density(aes(x = Residual, colour = State), size=2) + 
        theme(legend.position = c(0.9,0.7)) + 
        theme_bw() + facet_wrap( ~ State, nrow=1, scales = 'free') + 
        theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             panel.background = element_rect(fill = "transparent",colour = NA),
             axis.ticks = element_blank()) + 
             ylab('pdf') + xlab('Residual') + 
        ggtitle("HMM Residuals")

      # SW test on normality of residuals
      df      = data.frame(State = factor(states), Residual = hmm.residual)
      if (!is.null(sw_test)) levels(df$State) = paste(levels(df$State), round(sw_test, digits=5), sep=':')
      for (j in 1:length(unique(states))) {
        idx = which(df$Residual == j)
        df$Residual[idx] = df$Residual[idx] / hmm.sigma[j]^2
      }
      plt2 = ggplot(df, aes(sample = Residual)) + facet_wrap(~State) +  
        stat_qq(geom = "point", size = 2, position = "identity") + 
                #dparams = x@HMM$response$stdev) + 
        theme(legend.position = c(0.9,0.7)) + 
        theme_bw() + facet_wrap( ~ State, nrow=1) + 
        theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             panel.background = element_rect(fill = "transparent",colour = NA),
             axis.ticks = element_blank()) + 
             ylab('Sample Quantiles') + xlab('Theoretical Quantiles')
           grid.newpage() 

      # plot 
      grid.newpage() 
      pushViewport(viewport(layout = grid.layout(2, 1))) 
      print(plt1, vp = vplayout(1,1))
      print(plt2, vp = vplayout(2,1))
}

# _______________________________________
# Contributions plot: time series

plot_components_ts = function(df, timestamps, title = 'Covariate Contributions') {

	myTime   = format(timestamps, "%a,%m/%d %H:00")
	df$Index = 1:nrow(df)
	df.obs   = subset(df, select = c('fit', 'consumption', 'Index'))
	df$fit   = NULL
	df$consumption = NULL
	df.mlt   = melt(df, id.vars = c('Index'))
	lev.var  = levels(df.mlt$variable)
	lev.var.n= c(setdiff(lev.var, c(wthr_vars_all, 'mean')), intersect(lev.var, wthr_vars_all), 'mean')
	df.mlt$variable <- ordered( df.mlt$variable, levels = lev.var.n)

	df.mlt.p = df.mlt
	df.mlt.p$value[df.mlt.p$value < 0] = 0
	df.mlt.n = df.mlt
	df.mlt.n$value[df.mlt.n$value > 0] = 0
	
	p <- ggplot()   
	p <- p + geom_area(data = df.mlt.p, aes(x = Index, y = value, fill = variable)) 
	p <- p + geom_area(data = df.mlt.n, aes(x = Index, y = value, fill = variable)) 

	p <- p + geom_point(data = df.obs, aes(x = Index, y = fit), alpha = I(0.3), color='red') + 
		 geom_line(data = df.obs, aes(x = Index, y = fit), alpha = I(0.6), color='red') 
	p <- p + geom_point(data = df.obs, aes(x = Index, y = consumption), alpha = 0.3) + 
		 geom_line(data = df.obs, aes(x = Index, y = consumption), alpha = 0.6) 
	p <- p + scale_x_continuous(breaks = seq(1,length(myTime), length.out=10), 
                           	    labels = format(myTime[seq(1,length(myTime), length.out=10)], 
                                    format = "%a,%m/%d %H:00")) + 
        theme_bw() +
        theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.title.x = element_blank(),
               panel.background = element_rect(fill = "transparent",colour = NA),
               axis.ticks = element_blank() ) + 
        ylab("kWh [%]") + xlab('Time') + 
        theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
        ggtitle( title )

	return(p)
}

# ____________________________________________
# Dotplot tornado for covariate contributions 

plot_tornado = function(df, title = 'Components', nrow = 1) {

  df.mlt = melt(df, id = 'state')
  p <- ggplot(df.mlt, aes(x = 0, xend = value, y = variable, yend = variable, color = variable))
  p <- p + geom_point(aes(value, variable), size = 5) + geom_segment(size = 3) + ggtitle(title)
  p = p + facet_wrap(~state, nrow = nrow) + xlab('Component [%]')
  p = p + geom_vline(xintercept = 0)
  p = p + theme_bw() + 
        theme(panel.grid.major = element_blank(),
	       axis.text.y=element_text(size=15),
	       axis.text.x=element_text(size=15),
	       plot.title =element_text(size=20),
	       legend.position = 'none',
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               panel.background = element_rect(fill = "transparent",colour = NA),
               axis.ticks = element_blank() )

  return(p)
}

