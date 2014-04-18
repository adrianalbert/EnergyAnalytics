\name{Change Log}
\alias{Changes}
\title{Changes Made to Package HiddenMarkov}
\description{
This page contains a listing of recent changes made to the package.
}

\details{
\enumerate{
\item Since we have included different classes of HMMs (see \code{\link{dthmm}}, \code{\link{mmglm0}} and \code{\link{mmpp}}), it is much tidier to use an object orientated approach. This ensures that the functions across all models follow a more consistent naming convention, and also the argument list for the different model functions are more simple and consistent (see \code{\link{HiddenMarkov}}). (14 Sep 2007)

\item The main tasks (model fitting, residuals, simulation, Viterbi, etc) can now be called by generic functions (see topic \code{\link{HiddenMarkov}}). The package documentation has been rearranged so that these generic functions contain the documentation for all model types (e.g. see \code{\link{BaumWelch}}). (14 Sep 2007)

\item There are a number of functions, still contained in the package, that are obsolete. This is either because they do not easily fit into the current naming convention used to implement the more object orientated approach, or their argument list is slightly complicated. These functions have been grouped in the topics \code{dthmm.obsolete} and \code{mmpp.obsolete}. (14 Sep 2007)

\item There are various \emph{second level} functions. For example, the model fitting is achieved by the generic \code{\link{BaumWelch}} function. However, this will call functions to do the E-step, M-step, forward and backward probabilities, and so on. At the moment, these \emph{second level} functions have not been modified into an object orientated approach. It is not clear at this point whether this would be advantageous. If one went down this route, then one would probably group all of the E-step functions (for all models) under the same topic. If not, then it may be best to group all second level functions for each model under the same topic (e.g. \code{\link{forwardback}}, \code{\link{probhmm}} and \code{\link{Estep}} would be grouped together, being the second level functions for the \code{\link{dthmm}} model). (14 Sep 2007)

\item The original function called \code{Viterbi} has been renamed to \code{\link{Viterbihmm}}, and \code{\link{Viterbi}} is now a generic function. (14 Sep 2007)

\item Programming code that uses old versions of the functions should still work with this revised version of the package. However, you will get warning messages stating that certain functions are deprecated, and suggesting a possible alternative. To get a quick overview of the programming style, have a look at the examples in topic \code{\link{dthmm}}. (09 Nov 2007)

\item \code{\link{forwardback}}: \code{for} loops replaced by Fortran code; much faster. The corresponding R code is still contained within the function in case the Fortran has incompatibility issues. (23 Nov 2007)

\item \code{\link{forwardback.mmpp}}: \code{for} loops replaced by Fortran code. The corresponding R code is still contained within the function in case the Fortran has incompatibility issues. (24 Nov 2007)

\item \code{\link{Estep.mmpp}}: \code{for} loops replaced by Fortran code. Cuts time considerably. These loops in R used far more time than the forward and backward equations. The corresponding R code is still contained within the function in case the Fortran has incompatibility issues. (27 Nov 2007)

\item \code{\link{forwardback.mmpp}}, \code{\link{forwardback}} and \code{\link{Estep.mmpp}}: argument \code{fortran} added. (3 Dec 2007)

\item \code{\link{forwardback}}, \code{\link{forwardback.mmpp}} and \code{\link{Estep.mmpp}}: inclusion of all variable sized arrays in the Fortran subroutine call to be compatible with non gfortran compilers (3 Dec 2007); more added for calls to Fortran subroutines \code{multi1} and \code{multi2}. (6 Dec 2007)

\item \code{\link{Estep.mmpp}}: error in Fortran code of loop 6; \code{j1=0} to \code{j1=1}. (5 Dec 2007)

\item \code{\link{BaumWelch.mmpp}}: \code{if (diff < 0) stop ...} to \code{if (diff < 0 & control$posdiff) stop ...}, consistent with \code{\link{BaumWelch.dthmm}}. (11 Dec 2007)

\item \code{\link{logLik.dthmm}}, \code{\link{logLik.mmglm0}}, \code{\link{logLik.mmpp}}: \code{for} loop replaced by Fortran code. (15 Feb 2008)

\item \code{\link{dthmm}}: argument \code{discrete} set automatically for known distributions, stops if not set for unknown distributions. (15 Feb 2008)

\item \code{\link{neglogLik}}, \code{\link{Pi2vector}}, \code{\link{vector2Pi}}, \code{\link{Q2vector}}, \code{\link{vector2Q}}: new functions providing an alternative means of calculating maximum likelihood parameter estimates. (18 Feb 2008)

\item \code{\link{dthmm}}: argument \code{nonstat} was not set correctly. (21 Jun 2008)

\item Hyperlinks on package vignettes page. (22 Jun 2008)

\item \code{\link{mmpp}}: argument \code{nonstat} was not set correctly. (23 Jun 2008)

\item The manual pages \code{\link{HiddenMarkov-dthmm-deprecated}} and \code{\link{HiddenMarkov-mmpp-deprecated}} have been given a keyword of \dQuote{internal}. This hides them from the listing of package functions. (3 Jul 2008)

\item All cited references are now only listed in the topic \code{\link{HiddenMarkov}}. (3 Jul 2008)

\item \code{\link{neglogLik}}: argument \code{updatep} has been renamed to \code{pmap}. (9 Jul 2008)

\item \code{\link{neglogLik}}: format of this function changed to be consistent with that in package \pkg{PtProcess}. Argument \code{p} renamed as \code{params}. (07 Aug 2008)

\item \code{\link{mmglm0}}: remove some LaTeX specific formatting to be compatible with \R 2.9.0. (26 Jan 2009)

\item \code{\link{Viterbi}}: Correct hyperlink to base function \code{which.max}. (10 Oct 2009)

\item Tidied HTML representation of equations in manual pages. (15 Dec 2009)

\item \code{mmglm}: Renamed to \code{\link{mmglm0}}, new version \code{\link{mmglm1}}. See manual page for more details. (5 Jan 2010)

\item \code{\link{mmglmlong1}}: new function for longitudinal data. (18 Jan 2010)

\item \code{\link{dthmm}}: clarify argument \code{distn} on manual page, and nature of parameter estimates when the Markov chain is stationary. (04 Feb 2010)

\item \code{\link{BaumWelch.mmglmlong1}}: new argument \code{tmpfile} added. (13 Feb 2010)

\item \code{\link{Viterbi}}: Methods added for objects of class \code{"\link{mmglm1}"} and \code{"\link{mmglmlong1}"}. (29 Jul 2010)

\item \code{\link{logLik}}: Method added for object of class \code{"\link{mmglmlong1}"}. (30 Jul 2010)

\item \code{\link{forwardback.dthmm}}, \code{\link{forwardback.mmpp}}: New argument \code{"fwd.only"}. (30 Jul 2010)

\item \code{\link{logLik.dthmm}}: Calls \code{\link{forwardback.dthmm}} to perform calculations. (30 Jul 2010)

\item \code{\link{logLik.mmpp}}: Calls \code{\link{forwardback.mmpp}} to perform calculations. (30 Jul 2010)

\item \code{\link{Viterbi}}: Now generates an error message when applied to objects of class \code{"\link{mmpp}"}. Method not currently available. (03 Aug 2010)

\item \code{"\link{Viterbi.mmglm1}"}: Fixed bug with number of Bernoulli trials specification when using a binomial family. (05 Aug 2010)

\item \code{\link{residuals.dthmm}}, \code{\link{probhmm}}: Modify for greater efficiency and generality to accommodate more general models. Arguments of \code{\link{probhmm}} have also been changed. (05 Aug 2010)

\item \code{residualshmm}: Made defunct, incompatible with revised \code{\link{probhmm}}. (05 Aug 2010)

\item \code{\link{residuals}}: Methods added for objects of class \code{"\link{mmglm1}"} and \code{"\link{mmglmlong1}"}. Generates an error message when applied to objects of class \code{"\link{mmpp}"}, currently no method available. (05 Aug 2010)

\item Add CITATION file. (24 Sep 2010)

\item \code{\link{makedistn}}: Change \code{eval(parse(text=paste(x, " list(log=log)))", sep="")))} to \code{eval(parse(text=paste(x, " list(log.p=log)))", sep="")))}. (19 Dec 2010)

\item \code{\link{pglm}}, \code{\link{pmmglm}}: Change all \code{log} arguments to \code{log.p}. (19 Dec 2010)

\item Revise examples in \code{/tests} directory. (02 May 2011)

\item Implement very basic NAMESPACE and remove file \code{/R/zzz.R}. (5 Nov 2011)

\item List functions explicitly in NAMESPACE. (19 Dec 2011)

\item \code{\link{mmglm}} and \code{\link{neglogLik}}: Restrict the number of iterations in examples on manual pages to minimise time during package checks. (19 Dec 2011)

\item \code{modify.func}: New function to allow the user to modify package functions in the NAMESPACE set-up. This function violates the CRAN policy as users are not supposed to change the NAMESPACE on such packages. Some examples where it is required to modify package functions will not work, for example, the second example in \code{\link{Mstep}}. (7 Mar 2012)

\item \code{modify.func}: Function removed. See the second example in \code{\link{Mstep}} for a work-around when package functions need to be modified. (14 Apr 2012)

\item \code{\link{Mstep}}: Revised documentation about distribution requirements and ways to include other distributions into the software framework. (14 Apr 2012)
}}


\section{Future Development}{
\enumerate{
\item The functions \code{\link{Viterbi}} and \code{\link{residuals}} need methods for objects of class \code{\link{mmpp}}.

\item A number of the original functions have names that are too general. For example \code{\link{forwardback}} calculates the forward-backward probabilities, but only for the model \code{\link{dthmm}}. The corresponding function for the \code{\link{mmpp}} model is \code{\link{forwardback.mmpp}}. It would be more consistent to attach to these original functions a \code{dthmm} suffix.

\item The demonstration examples are all for \code{\link{dthmm}}. Also need some for \code{\link{mmglm1}}, \code{\link{mmglmlong1}} and \code{\link{mmpp}}.
}}

\keyword{documentation}


