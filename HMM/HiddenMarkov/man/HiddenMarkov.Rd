\name{HiddenMarkov}
\alias{HiddenMarkov}
\docType{package}

\title{Overview of Package HiddenMarkov}

\description{
In this topic we give an overview of the package.
}

\section{Classes of Hidden Markov Models Analysed}{
The classes of models currently fitted by the package are listed below. Each are defined within an object that contains the data, current parameter values, and other model characteristics.

\describe{
\item{Discrete Time Hidden Markov Model:}{
is described under the topic \code{\link{dthmm}}. This model can be simulated or fitted to data by defining the required model structure within an object of class \code{"\link{dthmm}"}.}
\item{Markov Modulated Generalised Linear Model:}{
is described under the topic \code{\link{mmglm1}}.}
\item{Markov Modulated Poisson Process:}{
is described under the topic \code{\link{mmpp}}. This model can be simulated or fitted to data by defining the required model structure within an object of class \code{"\link{mmpp}"}.}
}
}

\section{Main Tasks Performed by the Package}{
The main tasks performed by the package are listed below. These can be achieved by calling the appropriate generic function.
\describe{
\item{Simulation of HMMs:}{can be performed by the function \code{\link{simulate}}.}
\item{Parameter Estimation:}{can be performed by the functions \code{\link{BaumWelch}} (EM algorithm), or \code{\link{neglogLik}} together with \code{\link[stats]{nlm}} or \code{\link[stats]{optim}} (Newton type methods or grid searches).}
\item{Model Residuals:}{can be extracted with the function \code{\link{residuals}}.}
\item{Model Summary:}{can be extracted with the function \code{\link{summary}}.}
\item{Log-Likelihood:}{can be calculated with the function \code{\link{logLik}}.}
\item{Prediction of the Markov States:}{can be performed by the function \code{\link{Viterbi}}.}
}
All other functions in the package are called from within the above generic functions, and only need to be used if their output is specifically required. I have referred to some of these other functions as \dQuote{2nd level} functions, for example see the topic \code{\link{mmpp-2nd-level-functions}}.
}


\section{Organisation of Topics in the Package}{
\describe{
\item{Cited References:}{anywhere in the manual are only listed within this topic.}
\item{General Documentation:}{topics summarising general structure are indexed under the keyword \dQuote{documentation} in the Index.}
}
}


\section{Acknowledgement}{
Many of the functions contained in the package are based on those of Walter Zucchini (2005).
}

\references{
Baum, L.E.; Petrie, T.; Soules, G. & Weiss, N. (1970). A maximization technique occurring in the statistical analysis of probabilistic functions of Markov chains. \emph{Annals of Mathematical Statistics} \bold{41(1)}, 164--171. DOI: \url{http://dx.doi.org/10.1214/aoms/1177697196}

Charnes, A.; Frome, E.L. & Yu, P.L. (1976). The equivalence of generalized least squares and maximum likelihood estimates in the exponential family. \emph{JASA} \bold{71(353)}, 169--171. DOI: \url{http://dx.doi.org/10.2307/2285762}

Dempster, A.P.; Laird, N.M. & Rubin, D.B. (1977). Maximum likelihood from incomplete data via the EM algorithm (with discussion). \emph{J. Royal Statist. Society B} \bold{39(1)}, 1--38.

Elliott, R.J.; Aggoun, L. & Moore, J.B. (1994). \emph{Hidden Markov Models: Estimation and Control.} Springer-Verlag, New York.

Harte, D. (2006). \emph{Mathematical Background Notes for Package \dQuote{HiddenMarkov}.} Statistics Research Associates, Wellington. URL: \url{http://goo.gl/J6s7M}, file \dQuote{notes.pdf}.

Hartley, H.O. (1958). Maximum likelihood estimation from incomplete data. \emph{Biometrics} \bold{14(2)}, 174--194. DOI: \url{http://dx.doi.org/10.2307/2527783}

Klemm, A.; Lindemann, C. & Lohmann, M. (2003). Modeling IP traffic using the batch Markovian arrival process. \emph{Performance Evaluation} \bold{54(2)}, 149--173. DOI: \url{http://dx.doi.org/10.1016/S0166-5316(03)00067-1}

MacDonald, I.L. & Zucchini, W. (1997). \emph{Hidden Markov and Other Models for Discrete-valued Time Series.} Chapman and Hall/CRC, Boca Raton.

McCullagh, P. & Nelder, J.A. (1989). \emph{Generalized Linear Models (2nd Edition).} Chapman and Hall, London.

Rabiner, L.R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. \emph{Proceedings of the IEEE} \bold{77(2)}, 257--286. DOI: \url{http://dx.doi.org/10.1109/5.18626}.

Roberts, W.J.J.; Ephraim, Y. & Dieguez, E. (2006). On Ryden's EM algorithm for estimating MMPPs. \emph{IEEE Signal Processing Letters} \bold{13(6)}, 373--376. DOI: \url{http://dx.doi.org/10.1109/LSP.2006.871709}

Ryden, T. (1994). Parameter estimation for Markov modulated Poisson processes. \emph{Stochastic Models} \bold{10(4)}, 795--829. DOI: \url{http://dx.doi.org/10.1080/15326349408807323}

Ryden, T. (1996). An EM algorithm for estimation in Markov-modulated Poisson processes. \emph{Computational Statistics & Data Analysis} \bold{21(4)}, 431--447. DOI: \url{http://dx.doi.org/10.1016/0167-9473(95)00025-9}

Zucchini, W. (2005). \emph{Hidden Markov Models Short Course, 3--4 April 2005.} Macquarie University, Sydney.
}

\keyword{documentation}
