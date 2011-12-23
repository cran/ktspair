\name{AUC.calc}
\alias{AUC.calc}
\title{
Compute the AUC
}
\description{
Compute the AUC for a given curve (the area under the curve).
}
\usage{
AUC.calc(at, values)
}
\arguments{
  \item{at}{
  The values of the curve on the x-axis
}
  \item{values}{
  The values of the curve on the y-axis
}
}
\value{
  No value is returned
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
\code{\link{kts.pair}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}, \code{\link{ROC.voting}}, \code{\link{ROC.offset}}, \code{\link{ROC.graphic.ktsp}}
}
\note{
This function is used in the function ROC.graphic() to compute the AUC for the ROC curves.
}
\keyword{misc}
