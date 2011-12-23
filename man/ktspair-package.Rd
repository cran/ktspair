\name{ktspair-package}
\alias{ktspair-package}
\alias{ktspair}
\docType{package}
\title{
Computation of the k-TSP
}
\description{
The package computes the k-TSP. This method uses pairs of genes to perform a classification which compares the relative ordering of the gene expressions within each profile. 
It ranks pairs of genes with respect to a score based on the sensibility and the specificity achieved by each pair.
It selects the k pairs that achieved the maximum score with the restriction that a gene can appear in at most one pair. The number of pairs of genes is computed through crossvalidation or can be chosen by the user.
This package also contains function to display graphical properties of the k-TSP, to summarize the performance (accuracy, sensitivity, specificity, ROC curve) of the k-TSP, to predict new individual, etc.
The functions deals also with "Expression Set" classes. The functions contained in this package are able to deal correctly with the presence of NAs in the dataset.
}
\details{
\tabular{ll}{
Package: \tab ktspair\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-04\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Julien Damond <julien.damond@gmail.com> 

Maintainer: Julien Damond <julien.damond@gmail.com> 
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
\keyword{misc}
\seealso{
\code{\link{kts.pair}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  ktsp <- ktspcalc(dat,grp,3)
  ktsp <- ktspcalc(eSet,1,3)
  ktsp <- ktspcalc(eSet,grp,3)
  ktsp
  plot(ktsp)
  summary(ktsp)
  predict(ktsp)
 }
}
