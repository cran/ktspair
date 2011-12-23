\name{ordertsp}
\alias{ordertsp}
\title{
Ordering of the pairs of genes
}
\description{
Determine the order (decreasing) of the pairs of genes represented by their scores.
}
\usage{
ordertsp(delta, gamma)
}
\arguments{
  \item{delta}{
  The score Delta of the pairs of genes
}
  \item{gamma}{
  The score Gamma of the pairs of genes
}
}
\value{
  The position in which each pair should appear to be ordered (based on delta and gamma).
}
\references{
  D. Geman, C. d'Avignon, D. Naiman and R. Winslow, "Classifying gene expression profiles from pairwise mRNA comparisons," 
  Statist. Appl. in Genetics and Molecular Biology, 3, 2004.

  A.C. Tan, D.Q. Naiman, L. Xu, R.L. Winslow, D. Geman, "Simple decision rules for classifying human cancers from gene expression
  profiles," Bioinformatics, 21: 3896-3904, 2005.

  Aik Choon Tan, Daniel Q. Naiman, Lei Xu, Raimond L. Winslow, and Donald Geman, "Simple decision rules for classifying human cancers 
  from gene expression profiles," Bioinformatics, 21:3896NAK3904, October 2005.

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
\code{\link{kts.pair}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}, \code{\link{ktspcalc2}}
}
\note{
This function is used in the function ktspcalc2() to order the selected pairs of genes.
}
\keyword{misc}
