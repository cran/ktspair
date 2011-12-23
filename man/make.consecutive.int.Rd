\name{make.consecutive.int}
\alias{make.consecutive.int}
\title{
Transform the group vector into a binary vector
}
\description{
This function transforms the vector corresponding to the groups of the observations into a binary vector.
}
\usage{
make.consecutive.int(y)
}
\arguments{
  \item{y}{A numeric or character vector.}
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
\value{
Return a vector containing only 0 and 1.
}

\author{
Julien Damond \email{julien.damond@gmail.com}
}

\seealso{
\code{\link{ktspcalc}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  make.consecutive.int(grp)
 }
}
\note{
This function is useful to link the functions kts.pair and ktspcalc (kts.pair accepts only a vector grp containing 0 and 1).
}
\keyword{misc}
