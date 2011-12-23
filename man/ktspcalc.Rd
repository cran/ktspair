\name{ktspcalc}
\alias{ktspcalc}
\title{
Compute the k top scoring pairs based on a gene expression matrix and a group indicator vector.
}
\description{
This function computes the k pairs of genes that achieved the maximum difference between the sensitivity and the specificity (in absolute value) between two specific groups based on the comparison of the expressions of the two genes present in the pairs. 
}
\usage{
ktspcalc(dat, grp, k = NULL, cross = 5, performance = FALSE, healthy = NULL, seed = NULL, display = TRUE, length = 40, med = FALSE)
}
\arguments{
  \item{dat}{
  Can either be (a) a matrix of m lines (the gene expressions) and n columns (the observations) or (b) an eSet object.
}
  \item{grp}{
  Can either be (a) a character (or numeric) vector indicating the group of each observations or (b) an integer indicating the column of pData(dat) that represents the group of the observations.
}
  \item{k}{
  The number of pairs of genes used in the k-TSP (k). By default this parameter is computed through crossvalidation.
}
  \item{cross}{
  The number of fold that should be used for the crossvalidation estimation of the parameter k. By default the number of fold is 5.
}
  \item{performance}{
  An indicator if the performance of the model should be computed or not.
}
  \item{healthy}{
  This variable is used to determine which group will be considerer as the healthy group (specificity). Need to give the label of the group.
}
  \item{seed}{
  The crossvalidation and the computation of the performance of the model are based on a random partition of the dataset. The variable seed allows the user to fix a seed.
}
  \item{display}{
  Allows the user to avoid the function ktspcalc() to print warning message (mainly used in the function crossvalidation).
}
  \item{length}{
  This paramters allows the used to control the length of the list used in the C code.
}
  \item{med}{
  If the mean of the median between the two groups for each gene should be substracted to the dataset or not.
}
}
\details{
The original version of the k-TSP only works for two groups classification. It is possible to deal with multiclass classification by using trees or multiple steps methods
that reduce the problem to a combination of several two classes classifications, see Aik et al. (2005) for more details. The function computes 
the score Delta (see Geman et al. (2004) for more details about this score) for every possible pair of genes. This makes the required computational time to grow rapidly
in the number of genes. A pre-filtering step can be useful in some cases. This function is able to deal with NA present in the dataset. It considers only the patients
for which the gene expressions were measured and adapts the computation of the score Delta to the number of measures without NA.
This function also deals with "Expression Set" dataset. The group indicator can be replaced by the number of the column of pData(eSet) that contains the indicatior vectors of the group.
The user has the possibility to let the function compute a value for the parameter k. This value is computed through crossvalidation. The user can choose the number of
fold. The special case of the Leave-One-Out crossvalidation is when the number of fold is equal to the number of observations.
The user has the possibility to obtain the accuracy reached by the method via the variable performance. It will compute the accuracy using a partitioning of the dataset as
in the crossvalidation. If the crossvalidation had to be computed, the performance will be computed at the same time and will be based on the same partition. If a reference for a healthy group is given (via the variable healthy), the sensitivity and the specificity will be computed (estimated by crossvalidation).
The number of partition for the performance is also given by the variable cross. The user can also fix a seed via the variable seed.
The parameter length has to be used carefuly. Indeed, choosing a low value will speed up the computation but may produce results with
less pairs of genes than expected (k). In order to obtain k pairs of genes, one has to set a bigger value for this parameter. We note that, on some datasets, it is not possible to obtain k pairs of genes, this may be dued to the dataset itself and not the algorithm.
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

\section{note}{
The length sets to the list used in the C code (defined by the paramter length) has to be at least as big as k.
}

\section{Warning }{
If NAs are present in the dataset, the computation of the score Delta will be based only on observations for which the
measures of the two genes of the current comparison are not NA. This will reduce the number of observations used to compute the score Delta and
can produce lower accuracy of the estimation compared to the scores for others pairs.
} 

\seealso{
\code{\link{kts.pair}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  ktsp <- ktscalc(dat,grp,3)
  ktsp <- ktspcalc(eSet,grp,3)
  ktsp <- ktspcalc(eSet,1,3)
  ktsp
  ktsp$rankscore
  ktsp$accuracy_k
  ktsp$accuracy

  ktsp2 <- ktspcalc(dat, grp, 9)
  ktsp2
  ktsp2 <- ktspcalc(dat, grp, 9, length=40)
  ktsp2
 }
}
\keyword{misc}
