\name{plot.ktsp}
\alias{plot.ktsp}
\title{
Graphical representation of ktsp objects
}
\description{
Show a graphical representation of ktsp objects based on the expressions of the genes present in the pairs constituting the k-TSP.
}
\usage{
  \method{plot}{ktsp}(x,...)
}

\arguments{
  \item{x}{A ktsp object.}
  \item{...}{Plotting arguments, not used}
}
\details{
This is the generic function for ktspplot.
}
\value{
No value is returned. Only a graph is ploted.
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene expression profiles from pairwise mRNA comparisons," 
  Statist. Appl. in Genetics and Molecular Biology, 3, 2004.

  A.C. Tan, D.Q. Naiman, L. Xu, R.L. Winslow, D. Geman, "Simple decision rules for classifying human cancers from gene expression
  profiles," Bioinformatics, 21: 3896-3904, 2005.

  J. Damond, supervised by S. Morgenthaler and S. Hosseinian, "Presentation and study of robustness for several methods to classify individuals based on their gene  
  expressions", Master thesis, Swiss Federal Institute of Technology Lausanne (Switzerland), 2011.
  
  J. Damond, S. Morgenthaler, S. Hosseinian, "The robustness of the TSP and the k-TSP and the computation of ROC curves", paper is submitted in
  Bioinformatics, December 2011.
 
  Jeffrey T. Leek <jtleek@jhu.edu> (). tspair: Top Scoring Pairs for
  Microarray Classification. R package version 1.10.0.
}
\author{
Julien Damond \email{julien.damond@gmail.com}
}

\seealso{
\code{\link{ktspcalc}}, \code{\link{kts.pair}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  ktsp <- kts.pair(dat,grp,3)
  plot(ktsp,select=1)
  plot(ktsp)
 }
}
\keyword{misc}
