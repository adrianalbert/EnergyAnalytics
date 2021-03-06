# #########################################################
#' plot_hmm_ts
#' Plots time series HMM color-coded by state
#'
#' @param hmm.means Fitted values (means) of HMM
#' @param hmm.sigma Standard deviation time series for HMM
#' @param states Time series of states for HMM
#' @param timestamps Time series of time stamps 
#' @param observed Time series of original data
#' @param y.lab Label for y-axis
#' @param title Plot title
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

library('reshape2')
library('ggplot2')
library('RColorBrewer')
library('grid')

plot_hmm_ts = function(hmm.means, hmm.sigma, states, timestamps, observed, y.lab = 'observed', title = 'HMM-ts') {
	# construct plotting data frames
	df.hmm = data.frame(Regime = hmm.means, 
		                  Sigma  = hmm.sigma, 
		                  State  = states,      
			                Index  = 1:length(states),				                              
		                  Time   = format(timestamps, "%a,%m/%d %H"))              
	df.obs = data.frame(Observed = observed, 
		                  Time     = format(timestamps, "%a,%m/%d %H"),
			                Index    = 1:length(states),				                              
		                  Sigma    = rep(0, length(observed)))

	df.mlt.obs = melt(df.obs, id.vars = c('Time', 'Sigma', 'Index'))
	df.mlt.hmm = melt(df.hmm, id.vars = c('Time', 'Sigma', 'State', 'Index'), measure.vars = c('Regime'))
	df.mlt.hmm$variable = paste(df.mlt.hmm$variable, df.mlt.hmm$State, sep=' ')
	df.mlt.hmm$State = NULL

  # add in error bars
	pd <- position_dodge(.1)
	plt = ggplot(df.mlt.hmm, aes(x = Index, y = value)) + 
	geom_errorbar(aes(ymin = value-Sigma, ymax=value+Sigma), 
		      colour="black", width = 1.4, alpha=0.4)
  
	# plot current zoom-in
	plt = plt + 
    geom_line(size = 1.2) + geom_point(aes(colour=variable), size = 3) + 
	  geom_line(col='black', data = df.mlt.obs, alpha = 0.4, size = 1.2) + 
	  geom_point(col='black', data = df.mlt.obs, alpha = 0.4, size = 3) + 
	  scale_x_continuous(breaks = seq(1,length(df.mlt.hmm$Time), length.out=5), 
		           labels = format(df.mlt.hmm$Time[seq(1,length(df.mlt.hmm$Time), length.out=5)], 
		                           format = "%a,%m/%d %H:00"))
	plt = plt + theme_bw() +
	theme(panel.grid.major  = element_blank(),
	       panel.grid.minor = element_blank(),
	       panel.background = element_blank(),
	       axis.title.x     = element_blank(),
	       axis.text.y      = element_text(size=20), 
	       axis.text.x      = element_text(size=20),
	       axis.title.y     = element_text(size=20),
	       axis.title.x     = element_text(size=20),
	       plot.title       = element_text(size=22),            
	      legend.text      = element_text(size=22),	      
	      legend.title     = element_blank(),	      
	      legend.position  = c(0.2,0.75),	      
	      axis.ticks = element_blank() ) + 
	ylab(y.lab) + xlab('Time') + ggtitle( title )

	return(plt)
}

# #########################################################
#' plot_hmm_coefs_ts
#' Plots time series HMM color-coded by state w/ coefficients for additional covars.
#'
#' @param covar_state Data frame containing covariates that influence the HMM.
#' @param states Time series of decoded states in the HMM.
#' @param timestamps Time series of time stamps 
#' @param observed Time series of original data
#' @param title Plot title
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################
plot_hmm_coefs_ts = function(covar_state, states, observed, timestamps, 
                             title = 'HMM-coefs-ts') {
  # construct plotting data frames
  df = data.frame(obs    = observed, 
                  State  = as.factor(states),      
                  Index  = 1:length(states),				                              
                  Time   = format(timestamps, "%a,%m/%d %H:00"))   

#   if (ncol(covar_state)>1) tmp = covar_state[states,] else {
#     tmp = data.frame(covar_state[states,])
#     names(tmp) = names(covar_state)
#   }
  tmp = covar_state
  df = cbind(df,tmp)
  
  print(summary(df))
  
  df.mlt = melt(df, id.vars = c('Time', 'State', 'Index'))
  
  pd <- position_dodge(.1)
  plt = ggplot(df.mlt, aes(x = Index, y = value)) + 
    geom_line(position=pd, size = 1.2, alpha = 0.5) + 
    geom_point(position=pd, aes(colour=State), size = 3) + 
    scale_x_continuous(breaks = seq(1,length(df.mlt$Time), length.out=5), 
                       labels = format(df.mlt$Time[seq(1,length(df.mlt$Time), length.out=5)], 
                                       format = "%a,%m/%d %H:00"))
  plt = plt + facet_wrap(~variable, ncol = 1, scale = 'free')
  plt = plt + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=18),
          axis.title.x     = element_blank(),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),	      
          axis.ticks = element_blank() ) + 
    xlab('Time') + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( title )
  
  return(plt)
}


# #########################################################
#' plot_HMM_MC
#' Visualize the underlying structure of Markov Chain in HMM.
#'
#' @param hmm.means Fitted values (means) of HMM
#' @param hmm.sigma Standard deviation time series for HMM
#' @param transition Transition matrix parameters for HMM
#' @param title Plot title
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

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
	scale_size(range = c(4, 16)) + 
	geom_segment(aes(x = Mean.i, xend = Mean.j, 
		         y = Sigma.i, yend = Sigma.j, 
		         size=P, alpha=P), data=df.P) +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	     panel.grid.minor = element_blank(),
	     panel.background = element_blank(),
       axis.text.y      = element_text(size=18), 
       axis.text.x      = element_text(size=18),
       axis.title.y     = element_text(size=18),
       axis.title.x     = element_text(size=18),
       plot.title       = element_text(size=20),            
 	     legend.position = "none",
	     panel.background = element_blank(),
	     axis.ticks = element_blank()) + 
	ylab('Sigma') + xlab('State') + 
	ggtitle(title)
	return(plt)
}

# #########################################################
#' plot_state_breakdown
#' Plot breakdown of tod/dow occupancy state profiles.
#'
#' @param states Time series of decoded states in the HMM.
#' @param timestamps Time series of time stamps 
#' @param state.type Labels for states (strings)
#' @param title Plot title
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_state_breakdown = function(states, timestamps, state.type = NULL, 
                                title = 'States Breakdown') {
  
  # create categories for time of day and day of week
  df = data.frame(State   = as.factor(states), 
                  Hour    = hour(timestamps), 
                  Day     = weekdays(timestamps), 
                  Weekend = isHoliday(timeDate(timestamps)))
  if (!is.null(state.type)) 
    df$State = droplevels(factor(df$State, labels = state.type))
  #  df$State = factor(df$State, labels = paste(levels(df$State), state.type, sep=':'))
  df$Slot = rep(NA, nrow(df))
  df$Slot[df$Hour <= 5 | df$Hour >= 23] = 'Night'
  df$Slot[df$Hour > 5 & df$Hour < 11] = 'Morning'
  df$Slot[df$Hour >= 11 & df$Hour < 19] = 'Afternoon'
  df$Slot[df$Hour >= 19 & df$Hour < 23] = 'Evening'
  df$Slot = as.factor(df$Slot)
  df$Weekend = factor(df$Weekend, labels = c('Weekday', 'Weekends'))
  df$Season= 'Summer'            
  df$Season[which(month(timestamps) %in% c(11,12,1,2,3,4))] = 'Winter'
  
  # construct plot
  p = ggplot(df, aes(x = factor(Hour), fill = State))
  p = p + geom_bar(width = 0.9) 
  p = p + facet_wrap( ~Season, ncol=1, scales = 'free')
  p = p + theme_bw() +
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=22),
          axis.title.x     = element_text(size = 22),
          axis.text.y      = element_text(size=22), 
          axis.text.x      = element_text(size=22, angle=45),
          axis.title.y     = element_text(size=22),
          axis.title.x     = element_text(size=22),
          plot.title       = element_text(size=24),            
          legend.text      = element_text(size=24),        
          legend.title      = element_text(size=22),        
          axis.ticks = element_blank() ) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle(title )  + xlab('Hour of Day')
  
  return(p)
}


# #########################################################
#' plot_state_heatmap
#' Plot heatmap of yearly occupancy state profiles.
#'
#' @param myMat Time series of (hourly) observations.
#' @param timestamps Time series of time stamps 
#' @param title Plot title
#' @return heatmap plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_state_heatmap = function(myMat, timestamps, title = 'State Heatmap') {
  
  # select appropriate color palette
  if (is.factor(myMat)) hmcol<-brewer.pal(nlevels(myMat),'Paired') else hmcol<-brewer.pal(10,'RdBu')
  
  # convert long profile to wide profile
  profile.mat = matrix(as.numeric(myMat), ncol = 24, byrow = T)
  
  # choose appropriate heatmap labels
  days        = unique(as.Date(timestamps))
  y.labels    = rep(NA, nrow(profile.mat))
  idx         = seq(1, nrow(profile.mat), by=30)
  y.labels[idx] = days[idx]
  y.labels    = as.Date(y.labels)
  x.labels    = c(paste(0:12, ':00 AM', sep=''), paste(1:11, ':00 PM', sep=''))
  
  
  heatmap(profile.mat, Rowv = NA, Colv=NA, 
          labRow = y.labels, 
          labCol = x.labels, 
          cexRow = 1.5, cexCol = 1.3,
          col = hmcol, main = title)
}

# #########################################################
#' plot_state_heatmap2
#' Plot heatmap of yearly state profiles (ggplot version).
#'
#' @param myMat Time series of (hourly) observations.
#' @param timestamps Time series of time stamps 
#' @param title Plot title
#' @return heatmap plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_state_heatmap2 = function(myMat, timestamps, title = 'State Heatmap') {
  
  # convert long profile to wide profile
  profile.mat = data.frame(State = as.factor(myMat), Hour = hour(timestamps), 
                           Date = as.Date(timestamps))
  
  plt = ggplot(profile.mat, aes(x = Hour, y = Date)) + 
    geom_tile(aes(fill = State), size = 2) + 
    facet_wrap(~State, ncol = 2)
  
  plt = plt + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y      = element_text(size=18), 
          axis.text.x      = element_text(size=18),
          axis.title.y     = element_text(size=18),
          axis.title.x     = element_text(size=18),
          plot.title       = element_text(size=20),            
          legend.text      = element_text(size=18),
          axis.ticks = element_blank()) + 
    ggtitle(title)
  return(plt)
}

# # _________________________________________________________
# # Plot underlying MC of a HMM using ggplot: covariates case
# 
# plot_HMM_MC_cov = function(P, pi, contrib, title = 'HMM Structure', P.thresh = 0.05) {
# 
#   nStates = nrow(P)
#   state.grid = expand.grid(State.i = as.factor(1:1:nStates), 
#                            State.j = as.factor(1:1:nStates))
#   hmm.means = contrib$mu
#   hmm.sigma = contrib$sigma
#   df.P      = state.grid
#   df.P$P    = as.numeric(t(P))
#   df.P$Mean.i = hmm.means[df.P$State.i]
#   df.P$Mean.j = hmm.means[df.P$State.j]
#   
#   df.S      = data.frame(State = 1:nStates, mu = hmm.means, sigma = hmm.sigma, Characteristic.Time = pi)
# 
#   # states characteristics
#   contrib <- contrib[,colSums(is.na(contrib))<nrow(contrib)]
#   df = melt(contrib, id.vars = c('mu', 'sigma2', 'state'))
#     
#   # add (mu, sigma) canvass
#   plt = ggplot(df) + 
#     geom_subplot(aes(x = state, y = mu, group = state,
#                      subplot = geom_bar(aes(x = variable, y = value, fill = variable), width = 0.95,
#                                         stat = 'identity', position = 'dodge')))
# 
#   # add in links
#   plt = plt +     
#     geom_segment(aes(y = Mean.i, yend = Mean.j, 
#                      x = State.i, xend = State.j, 
#                      alpha=P), data=df.P, size = 3) +
# #     scale_size(range = c(0.5, 5)) +     
#     theme_bw() +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank(),
#           axis.text.y      = element_text(size=18), 
#           axis.text.x      = element_text(size=18),
#           axis.title.y     = element_text(size=18),
#           axis.title.x     = element_text(size=18),
#           plot.title       = element_text(size=20),            
#           legend.text      = element_text(size=18),
#           axis.ticks = element_blank()) + 
#     ylab('Mean') + xlab('State') + 
#     ggtitle(title)
# 
#   # add state label
#   plt = plt + geom_point(data = df.S, aes(x = State, y = mu, size = Characteristic.Time, color = Characteristic.Time), show_guide = FALSE) + 
#     scale_size(range = c(4, 16)) + 
#     geom_text(data = df.S, aes(x = State, y = mu, label = State, size = Characteristic.Time), color = 'white', show_guide = FALSE) +
#     geom_text(data = df.S, aes(x = State, y = mu, label = paste('s=',round(sigma, digits=2))), hjust = 1, vjust = 4, show_guide = FALSE)
#     
#   return(plt)
# }

# #########################################################
#' plot_HMM_MC_NET
#' Plot MC structure of HMM using igraph.
#'
#' @param hmm.means Fitted values (means) of HMM
#' @param hmm.sigma Standard deviation time series for HMM
#' @param transition Transition matrix parameters for HMM
#' @return igraph plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

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

# #########################################################
#' plot_HMM_res
#' Residual analysis for HMM.
#'
#' @param hmm.residual Residuals time series for HMM.
#' @param states Time series of states for HMM
#' @param hmm.sigma Standard deviation time series for HMM
#' @param title Plot title
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_HMM_res = function(hmm.residual, states, hmm.sigma) {

      # densities of residuals
      df = data.frame(State = factor(states), Residual = hmm.residual)
      plt = ggplot(df) + 
        geom_density(aes(x = Residual, colour = State), size=2) + 
        theme(legend.position = c(0.9,0.7)) + 
        theme_bw() + facet_wrap( ~ State, nrow=1, scales = 'free') + 
        theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             axis.text.y      = element_text(size=18), 
             axis.text.x      = element_text(size=18),
             axis.title.y     = element_text(size=18),
             axis.title.x     = element_text(size=18),
             plot.title       = element_text(size=20),            
             legend.text      = element_text(size=18),
             axis.ticks = element_blank()) + 
             ylab('pdf') + xlab('Residual') + 
        ggtitle("HMM Residuals")

      print(plt)
}


# #########################################################
#' plot_components_ts
#' Plot emissions components (by covariates) of HMM fit. 
#'
#' @param df Data frame containing components calculations. 
#' @param timestamps Time stamps of time series. 
#' @param states (optional) Add scale to indicate states
#' @param covars (optional) Add panel with time series of responses to covariates in emissions. 
#' @param title Plot title. 
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_components_ts = function(df, timestamps, 
                              states = NULL, covars = NULL,
                              title = 'Covariate Contributions') {

  # ______________________________
  # Main components plot
  
	myTime   = format(timestamps, "%m/%d %H:00")
	df$Index = 1:nrow(df)
	df.obs   = subset(df, select = c('fit', 'obs', 'Index'))
  df$state = NULL
  df$fit   = NULL
  df$obs   = NULL
	df.mlt   = melt(df, id.vars = c('Index'))
  df.mlt$value = abs(df.mlt$value)
	lev.var  = levels(df.mlt$variable)
	lev.var.n= c('(Intercept)', setdiff(lev.var, '(Intercept)'))
	df.mlt$variable <- ordered( df.mlt$variable, levels = lev.var.n)

	df.mlt.p = df.mlt
	df.mlt.p$value[df.mlt.p$value < 0] = 0
	df.mlt.n = df.mlt
	df.mlt.n$value[df.mlt.n$value > 0] = 0
	
	p <- ggplot()   
# 	p <- p + geom_area(data = df.mlt.p, aes(x = Index, y = value, fill = variable)) # positive contributions
# 	p <- p + geom_area(data = df.mlt.n, aes(x = Index, y = value, fill = variable)) # negative contributions
  
  # add in states color bar
  if (!is.null(states)) {
    y_pos =  max(df.obs$fit)*1.5
    df_state = data.frame(Index = 1:length(states), State = as.factor(states), obs = y_pos)
  	p <- p + geom_point(data = df_state, aes(x = Index, y = obs, colour = State), size = 6, shape = 15)
    # p <- p + scale_color_grey()
  }	
  
  # add in fitted time series
	p <- p + geom_point(data = df.obs, aes(x = Index, y = fit), alpha = I(0.6), color='red', size = 2) + 
		 geom_line(data = df.obs, aes(x = Index, y = fit), alpha = I(0.7), color='red') 

  # add in reference time series
	p <- p + geom_point(data = df.obs, aes(x = Index, y = obs), alpha = 0.6, size = 2) + 
		 geom_line(data = df.obs, aes(x = Index, y = obs), alpha = 0.8) 
  
  # x axis time labels
	p <- p + scale_x_continuous(breaks = seq(1,length(myTime), length.out=5), 
                           	    labels = format(myTime[seq(1,length(myTime), length.out=5)], 
                                    format = "%a,%m/%d %H:00")) + 
        theme_bw() +
        theme(panel.grid.major  = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y      = element_text(size=18),
               axis.text.x      = element_text(size=18, angle = 30),
               axis.title.y     = element_text(size=18),
               axis.title.x     = element_blank(),
               plot.title       = element_text(size=20),
               legend.text      = element_text(size=18),              
               plot.margin      = unit(c(0,0,0,0), "cm"),
               axis.ticks       = element_blank() ) + 
        ylab("observation") + 
        theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
        ggtitle( title )

  # _____________________
  # Auxiliary covariates
  
    if (!is.null(covars)) {
    df2    = data.frame(Index = 1:length(myTime), state = as.factor(states))
    df2  = cbind(df2, covars)
    df2    = melt(df2, id.vars = c('Index', 'state'))    
#     theta  = atan(df2$value)
#     df2$dx = 1
#     df2$dy = df2$value    
#     
#     p2  = ggplot(df2) + 
#       geom_point(aes(x = Index, y = 0, color = state), size=2) 
#     p2 = p2 + geom_segment(aes(x = Index-dx/2, xend = Index + dx/2, 
#                                y = -dy/2, yend = dy/2, color = variable), 
#                            arrow = arrow(length = unit(0.1,"cm")), alpha = 0.7)
    p2  = ggplot(df2) + 
      geom_point(aes(x = Index, y = value, color = state), size=2) 
    p2 = p2 + geom_line(aes(x = Index, y = value))
    p2 = p2 + facet_wrap(~variable, ncol=1, scale = 'free') + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.text.x     = element_text(size=18),
            axis.text.y      = element_text(size=18),
            axis.text.x      = element_blank(),
            axis.title.y     = element_blank(),
            axis.title.x     = element_blank(),
            plot.title       = element_text(size=20),
            legend.text      = element_text(size=18),
            plot.margin      = unit(c(0,0,0,0), "cm"),
            axis.ticks       = element_blank() ) + 
      theme(plot.title=element_text(family="Times", face="bold", size=20))
    grid.newpage() 
    pushViewport(viewport(layout = grid.layout(2, 1, height = c(0.7, 0.3)))) 
    print(p,  vp = vplayout(1,1))
    print(p2, vp = vplayout(2,1))  
  } else print(p)
  
  return(NULL)
}

# #########################################################
#' plot_tornado
#' Dotplot tornado for covariate contributions.
#'
#' @param df Data frame containing components calculations. 
#' @param title Plot title
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_tornado = function(df, title = 'Components', nrow = 1) {

  df <- df[,colSums(is.na(df))<nrow(df)]  
  df.mlt = melt(df, id = 'state')
  p <- ggplot(df.mlt, aes(x = 0, xend = value, y = variable, yend = variable, color = variable))
  p <- p + geom_point(aes(value, variable), size = 5) + geom_segment(size = 3) + ggtitle(title)
  p = p + facet_wrap(~state, nrow = nrow) + xlab('Component [%]')
  p = p + geom_vline(xintercept = 0)
  p = p + theme_bw() + 
        theme(panel.grid.major = element_blank(),
	       axis.text.y=element_text(size=18),
	       axis.text.x=element_text(size=18),
         axis.title.y     = element_text(size=18),
         axis.title.x     = element_text(size=18),
         plot.title =element_text(size=20),
	       legend.position = 'none',
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.ticks = element_blank() )

  return(p)
}

# #########################################################
#' plot_dep_covar
#' Plots dependence of dependent variable on given covariate by state.
#'
#' @param x Covariate
#' @param y Response
#' @param state State of each response
#' @param title Plot title
#' @param x.lab Title of x-axis
#' @param y.lab Title of y-axis
#' @param separate Separate plot into panels for each state
#' @param highlight Assign each state a different symbol
#' @param markers Add in vertical markers at each of the elements in the vector. 
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_dep_covar = function(x, y, state, title = 'HMM-dep-covar', 
                          y.lab = '', x.lab = '', highlight = T, markers = c(), separate = F, ...) {
  
  df = data.frame(x = x, y = y, state = as.factor(state))
  
  # plot current zoom-in  
  plt = ggplot(df, aes(x = x, y = y))
  if (highlight) plt = plt + geom_point(aes(color = state, shape = state), size = 2) else plt = plt + geom_point()
  
  # add in vertical markers
  if (length(markers)>0 & highlight)
    plt = plt + geom_vline(data = data.frame(z = markers, labels = names(markers)), 
                           aes(xintercept = z), linetype = "longdash", size = 2)  
  #   plt = plt + stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) + scale_fill_gradient(limits=c(1e-3,1e-1))
  plt = plt + 
        scale_colour_hue(l=50) + # Use a slightly darker palette than normal
        geom_smooth(method=lm,   # Add linear regression lines
                    se=T,     # Don't add shaded confidence region
                    fullrange=F, # Extend regression lines
                    aes(group=state, color = state),alpha = 0.1, size = 0.5) 
  plt = plt + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y      = element_text(size=22),
          axis.text.x      = element_text(size=22),
          axis.title.y     = element_text(size=22),
          axis.title.x     = element_text(size=22),
          legend.title     = element_text(size=22),
          plot.title       = element_text(size=24),
          legend.text      = element_text(size=24),
          legend.position  = "none", #c(0.2,0.7),
          axis.ticks       = element_blank() ) + 
    ylab(y.lab) + xlab(x.lab) + ggtitle( title )
  nx = sqrt(length(unique(df$state)))
  if (separate) plt = plt + facet_wrap(~state, ncol = nx)
  return(plt)
}

# #########################################################
#' plot_tran_covar
#' Plots dependence of transition probability on given covariate.
#'
#' @param x_var Covariate grid.
#' @param dep List of distributions dependence on covariate for each state.
#' @param title Plot title
#' @param x.lab Title of x-axis
#' @param y.lab Title of y-axis
#' @param markers Add in vertical markers at each of the elements in the vector. 
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

plot_tran_covar = function(x_var, dep, 
                           type = paste('State', 1:length(dep), sep='.'),
                           title = 'HMM-dep-covar', 
                           y.lab = '', x.lab = '', markers = c()) {
  
  df       = as.data.frame(do.call('rbind', dep))
  names(df)= type
  df$From  = rep(1:length(dep), each = nrow(dep[[1]]))  
  covar    = names(x_var)
  df       = cbind(rep(x_var[,covar], length(dep)), df)
  names(df)[1] = covar
  df.mlt   = melt(df, id.vars = c(covar, 'From'))
  df.mlt$From = as.factor(df.mlt$From)
  
  # plot 
  plt = ggplot(df.mlt, aes_string(x = covar, y = 'value'))
  plt = plt + geom_point(aes(color = variable, shape = variable))
  nx = sqrt(length(unique(df$From)))  
  if (length(dep)>1) plt = plt + facet_wrap(~From, ncol = nx)
  
  # add in vertical markers
  if (length(markers)>0)
    plt = plt + geom_vline(data = data.frame(z = markers, labels = names(markers)), 
                           aes(xintercept = z), linetype = "longdash", size = 2)  
  # formatting
  plt = plt + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text.x     = element_text(size=22),
          axis.text.y      = element_text(size=22),
          axis.text.x      = element_text(size=22),
          axis.title.y     = element_text(size=22),
          axis.title.x     = element_text(size=22),
          plot.title       = element_text(size=24),
          legend.text      = element_text(size=24),
          legend.title     = element_text(size=22),
          legend.position  = "none", #c(0.15, 0.6),
          axis.ticks       = element_blank() ) + 
    ylab(y.lab) + xlab(x.lab) + 
    theme(plot.title=element_text(family="Times", face="bold", size=20)) + 
    ggtitle( title )
  return(plt)
}

# #########################################################
#' multiplot
#' Multiple plots on a page for ggplot objects. 
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param plotlist (optional) a list of ggplot objects
#' @param cols number of columns for page
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#' @return ggplot plot object
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
# #########################################################

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

hmap = function(data,colorMap=NULL,yvals=NA,xvals=NA,xlabs=NULL,ylabs=NULL,cex.axis=1,...) {
  n = dim(data)[1]
  m = dim(data)[2]
  # defailt values
  if(is.null(colorMap)) { colorMap = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)) }
  image(t(data),col=colorMap,axes=F,...)
  if(!is.null(xlabs)){ axis(1, at=(1:length(xlabs) - 1)/(length(xlabs)-1), labels=xlabs, cex.axis=cex.axis) }
  if(!is.null(ylabs)){ axis(2, at=(1:length(ylabs) - 1)/(length(ylabs)-1), labels=ylabs, cex.axis=cex.axis) }    
}
