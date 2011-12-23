\name{predict.ktsp}
\alias{predict.ktsp}
\title{
Prediction using a ktsp object
}
\description{
This function is used to predict the group of an individual based on its gene expressions via a ktsp object.
}
\usage{
predict.ktsp(object, dat = NULL, select = NULL, display = TRUE,...)
}
\arguments{
  \item{object}{
A ktsp object
}
  \item{dat}{
A dataset from which predictions should be computed. It can be either under the form of a matrix or under the form of an eSet.
}
  \item{select}{
An integer to perform the prediction based on a single TSP present in the k-TSP.
}
  \item{display}{
  Allows the user to avoid the function ktspcalc() to print waring message over the loop.
}
  \item{...}{Plotting arguments, not used}
}
\details{
By default, this function computes predictions of the orginal dataset (the one used to compute the k-TSP) based on the k-TSP. It is possible to predict
new obsevations via the variable "dat". It can either be under the form of a matrix or of an eSet. The function searchs for the gene names of the ktsp and 
try to identify them in the rownames or featuresNames of the matrix of gene expressions "dat". If rownames or featureNames are not available, it uses the number of the 
row as name to perform the predictions. By default, the prediction is based on all the pairs present in the k-TSP and by using a voting system.
It is also possible to have predictions based on a single pair by using the variable "select". In order to obtain a well defined predicton, the number of pairs of genes
used has to be an odd number. If the number of pairs is even, the function will not consider the pair with the lowest score Delta in the calculation of the prediction.

}
\value{
A vector containing the class predictions of each individual based on the k-TSP is returned.
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
\code{\link{ktspcalc}}, \code{\link{ktspplot}},\code{\link{kts.pair}}, \code{\link{summary.ktsp}}
}
\section{warning}{
If rownames or featureNames of the gene expression matrix "dat" are not available, the number of the row will be used as name. To be correct the order
of the genes in the matrix "dat" must be the same as in the ktsp object (or equivalently the same as in the original dataset).
}

\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  ktsp <- kts.pair(dat,grp,3)
  predict(ktsp)
  predict(ktsp, select=1)
 }
}
\keyword{misc}
