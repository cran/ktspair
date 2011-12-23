\name{kts.pair}
\alias{kts.pair}
\title{
Calculation of the k top scoring pairs.
}
\description{
This function computes the k pairs of genes that achieved the maximum difference between the sensitivity and the specificity (in absolute value) between
two specific groups based on the comparison of the expressions of the two genes present in the pairs. The function ktspcalc() is the general function and uses
the function kts.pair once the dataset has been prepared for this function. The function ktspcalc() is also able to deal with eSets.
}
\usage{
kts.pair(dat, grp, k, display = TRUE, length = 40, med = FALSE)
}
\arguments{
  \item{dat}{
  A matrix of m lines (the gene expressions) and n columns (the observations).
}
  \item{grp}{
  A vector of 0 and 1 for the groups of the observations.
}
  \item{k}{
  The number of pairs of genes that the function will select.
}
  \item{display}{
  Allows the user to avoid the function kts.pair() to print warning message (mainly used in the function crossvalidation).
}
  \item{length}{
  This paramters allows the used to control the length of the list used in the C code.
}
  \item{med}{
  If the mean of the median between the two groups for each gene should be substracted to the dataset or not.
}
}
\details{
  This function only works with matrices and vector of group containing only 0 and 1. For a more general use (eSets and labels
  for the groups) see the function ktspcalc(). This classifier can only be used for classifications with two groups. 
  The k-TSP was introduced in Tan et al. (2005) and is an extension of the TSP, which was presented in Geman et al. (2004).
}
\value{
  A ktsp object with the following elements:
  \item{index}{A k by 2 matrix composed of genes where the ith row stands for the ith best pair of genes with the restriction that a gene can appear
  in only one pair. The pairs are seleceted with respect to the score Delta and Gamma (in case of ties), see Tan et al. (2005) for more details about the k-TSP.}
  \item{ktspcore}{A vector of size k containing the scores Delta achieved by each selected pair of genes. The score Delta is based on the
  sensitivity and the specificity of a pair, see Geman et al. (2004) for more details.}
  \item{grp}{The group for each observation in a binary form}
  \item{ktspdat}{The row i and the row i+k represents the expressions of the genes present in the ith pair.}
  \item{k}{The number of pairs of genes.}
  \item{labels}{The name of the two groups that were present in the original variable grp.}
  \item{rankscore}{The score Gamma achieved by each pair of genes, for more details on this score see Geman et al. (2004).}
  \item{accuracy}{A vector of the estimated percentage of correct prediction for the k-TSP with k=1,3,5,7,9.}
  \item{accuracy_k}{The estimated percentage of correct prediction of the k-TSP with the selected k.}
  \item{sensitivity}{A vector of the estimated sensitivity for the k-TSP with k=1,3,5,7,9.}
  \item{specificity}{A vector of the estimated specificity for the k-TSP with k=1,3,5,7,9.}
  \item{med}{If the mean of the medians within each group has been substracted to the dataset return the values of the mean of the median, return FALSE otherwise}
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

\section{note}{
The length sets to the list used in the C code (defined by the paramter length) has to be at least as big as k.
}

\seealso{
\code{\link{ktspcalc}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  ktsp <- kts.pair(dat,grp,3)
  ktsp$index
  ktsp$ktspscore

  ktsp2 <- kts.pair(dat, grp, 9)
  ktsp2
  ktsp2 <- kts.pair(dat, grp, 9, length=40)
  ktsp2
 }
}
\keyword{misc}
