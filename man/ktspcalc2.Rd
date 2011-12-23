\name{ktspcalc2}
\alias{ktspcalc2}
\title{
The k-TSP with chosen pairs of genes
}
\description{
This function creates a ktsp object where the pairs of genes used have been selected by the user. 
}
\usage{
ktspcalc2(indice, dat, grp, length = 40, cross =5, display = display, med = FALSE, performance = FALSE, seed = NULL, healthy = NULL)
}
\arguments{
  \item{indice}{
  The pairs that should be used to construct the ktsp object
}
  \item{dat}{
  Can either be (a) a matrix of m lines (the gene expressions) and n columns (the observations) or (b) an eSet object.
}
  \item{grp}{
  Can either be (a) a character (or numeric) vector indicating the group of each observations or (b) an integer indicating the column of pData(dat) that represents the group of the observations.
}
  \item{length}{
  This paramters allows the used to control the length of the list used in the C code.
}
  \item{cross}{
  The number of fold that should be used in the crossvalidation.
}
  \item{display}{
  Allows the user to avoid the function ktspcalc() to print waring message over the loop.
}
  \item{med}{
  If the mean of the median between the two groups for each gene should be substracted to the dataset or not.
}
  \item{performance}{
  An indicator if the performance of the model should be computed or not.
}
  \item{seed}{
  The crossvalidation and the computation of the performance of the model are based on a random partition of the dataset. The variable seed allows the user to fix a seed.
}
  \item{healthy}{
  This variable is used to determine which group will be considerer as the healthy group (specificity). Need to give the label of the group in the same form as given in the variable grp.
}
}
\details{
This function would be useful when the user wants to select pairs of genes based on the results of the bootstrap. The pairs of genes will be 
then reordered based on the score Delta and Gamma from the best one to the worst one. When the genes are inserted (through the variable indice), the two first genes will represent the first pair, the two next genes the second pair, etc.
}
\value{
A ktsp object, see the function ktspcalc() for more details.
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
  ktsp2 <- ktspcalc2(c(100,200,300,400),dat, grp)
  ktsp2
 }
}
\keyword{misc}
