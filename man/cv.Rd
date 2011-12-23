\name{cv}
\alias{cv}
\title{
Crossvalidation for the parameter k
}
\description{
This function computes the value of k through crossvalidation. 
}
\usage{
cv(dat, grp, cross = 5, display = FALSE, length = 40, seed = NULL, med = FALSE, healthy = NULL)
}
\arguments{
  \item{dat}{
  Can either be (a) a matrix of m lines (the gene expressions) and n columns (the observations) or (b) an eSet object.
}
  \item{grp}{
  Can either be (a) a character (or numeric) vector indicating the group of each observations or (b) an integer indicating the column of pData(dat) that represents the group of the observations.
}
  \item{cross}{
  The number of fold that should be used in the crossvalidation.
}
  \item{display}{
  Allows the user to avoid the function ktspcalc() to print waring message over the loop.
}
  \item{length}{
  This paramters allows the used to control the length of the list used in the C code.
}
  \item{seed}{
  Allow the user to set a seed.
}
  \item{med}{
  If the mean of the median between the two groups for each gene should be substracted to the dataset or not.
}
  \item{healthy}{
  This variable is used to determine which group will be considerer as the healthy group (specificity). Need to give the label of the group.
}
}

\value{
  \item{k}{The selected value for k}
  \item{accuracy_k}{The estimated percentage of correct prediction achieved by the k-TSP with the selected k.}
  \item{accuracy}{A vector of the estimated percentage of correct prediction reached by the k-TSP with k = 1,3,5,7,9.}
  \item{sensitivity}{A vector of the estimated sensitivity reached by the k-TSP with k = 1,3,5,7,9.}
  \item{specificity}{A vector of the estimated specificity reached by the k-TSP with k = 1,3,5,7,9.}
}

\details{
This function computes the value of k through crossvalidation. The number of fold is given by the variable cross and by default is 5.
It also computes the percentage of correct prediction based on the same partition as for the crossvalidation for the possible values of k.
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
  cv <- cv(dat, grp,cross =10)
  ktsp <- ktspcalc(dat, grp, cv$k)
  ktsp
  cv
 }
}
\keyword{misc}
