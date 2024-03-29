% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FTT_algorithmRed.R
\name{FTT_algorithmRed}
\alias{FTT_algorithmRed}
\title{Alternative algorithm to test for fixing taxon traceability}
\usage{
FTT_algorithmRed(data, verbose = F, c = 4, n)
}
\arguments{
\item{data}{Data.table as constructed by create input. All possible c-tuples given the taxon set with status information (c-tuple as input available, c-tuple not in input = unresolved). In addition, all c (c-1)-tuples possible by each c-tuple are listed.}

\item{verbose}{Logical parameter if message should be printed, Default: F}

\item{c}{Parameter indicating the size of the c-tuples. E.g. if c=4, one uploads quadruples to be tested for fixing taxon traceability. Default: 4}

\item{n}{Number of all taxa in the analysis}
}
\value{
The same data.table is returned, with updated status, fixing taxon & round of resolvement.
}
\description{
Testing the input data stepwise for fixing taxa and check if the whole set can be resolved or not. The algorithm is an adaption to the \emph{red} one, and it checks if any red c-tuple if can be solved by \emph{c} green ones.
}
\details{
See publication and master thesis
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[data.table]{copy}}
}
