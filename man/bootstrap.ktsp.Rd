\name{bootstrap.ktsp}
\alias{bootstrap.ktsp}
\title{
Bootstrap procedure for the k-TSP
}
\description{
  This function computes n bootstrap resamples of the original dataset, at each step it applies the k-TSP on the simulated dataset and records
  which pairs of genes were selected as well as the number k.
}
\usage{
bootstrap.ktsp(dat, grp, k = NULL, seed = NULL, n = 500, length = 40, display = FALSE, med = FALSE)
}
\arguments{
  \item{dat}{
  Can either be (a) a matrix of m lines (the gene expressions) and n columns (the observations) or (b) an eSet object.
}
  \item{grp}{
  Can either be (a) a character (or numeric) vector indicating the group of each observations or (b) an integer indicating the column of pData(dat) that represents the group of the observations.
}
  \item{k}{
  If the number of pair of genes should be the same at each step or if it should be computed every time (by crossvalidation).
}
  \item{seed}{
  If a seed should be set or not.
}
  \item{n}{
  The number of bootstrap resamples.
}
  \item{length}{
  To control the length of the list in the C code, see kts.pair() or ktspcalc() for more details.
}
  \item{display}{
  If warnings should be displayed or not.
}
  \item{med}{
  If the mean of the median between the two groups for each gene should be substracted to the dataset or not.
}
}
\value{
A bootstrap object with the following elements
  \item{score}{A n by 9 matrix containing the scores of the k-TSP computed on the bootstrap. May contain NA (if k<9).}
  \item{index}{A n by 18 matrix containing the index of the k-TSP computed on the bootstrap. May contain NA (if k<9).}
  \item{k_value}{A vector of length n containing the number of pairs for each k-TSP (k).}
  \item{k}{A booleen indicating if the value of k was computed at every step or chosen as a parameter by the user.}
  \item{n}{The number of boostrap resamples.}
  \item{genenames}{The names of the genes contained in the dataset. If none, a numerotation is used.}
  \item{ktsp}{A ktsp object computed on the original dataset.}
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
\code{\link{ktspcalc}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
   data(ktspdata)
   bootstrap1 <- bootstrap.ktsp(dat, grp, k=3, n=20)
   bootstrap2 <- bootstrap.ktsp(dat, grp, n=20, seed=1)
  }
}
\keyword{misc}
