\name{ROC.voting}
\alias{ROC.voting}
\title{
ROC curve based on the voting system for k-TSP method
}
\description{
This function creates a ROC object used by the function ROC.graphic.ktsp() to plot a ROC curve based on a cutoff on the majority of voting system.
}
\usage{
ROC.voting(dat, grp, n = 200, healthy = NULL, seed = NULL, para1 = 200, para2 = 1, length = 40, display = FALSE, med = FALSE)
}
\arguments{
  \item{dat}{
  Can either be (a) a matrix of m lines (the gene expressions) and n columns (the observations) or (b) an eSet object.
}
  \item{grp}{
  Can either be (a) a character (or numeric) vector indicating the group of each observations or (b) an integer indicating the column of pData(dat) that represents the group of the observations.
}
  \item{n}{
  The number of bootstrap resample that should be used to compute the ROC curve.
}
  \item{healthy}{
  This variable is used to determine which group will be considerer as the healthy group (specificity). Need to give the label of the group.
}
  \item{seed}{
  If a seed should be set or not.
}
  \item{para1}{
  This paramater is used to control the smoothness of the ROC curves computed on the bootstrap.
}
  \item{para2}{
  This paramater is also used to control the smoothness of the ROC curves computed on the bootstrap.
}
  \item{length}{
  This paramters allows the used to control the length of the list used in the C code.
}
  \item{display}{
  If warnings should be displayed or not.
}
  \item{med}{
  If the mean of the median between the two groups for each gene should be substracted to the dataset or not.
}
}
\details{
The parameters para1 and para2 are used to make the ROC curve smoother. The way of computing the ROC curve is based on a cutoff on the voting system. This procedure is described in Damond et al. (paper is submitted in Bioinformatic, December 2011).
}
\value{
A ROC object with the following elements
  \item{spec}{A n by 2(para1+para2)+1 matrix representing the sensitivities achieved among the boostraps.}
  \item{sens}{A n by 2(para1+para2)+1 matrix representing the specificities achieved among the boostraps.}
  \item{n}{The number of boostrap resample used.}
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
\code{\link{ktspcalc}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}, \code{\link{ROC.offset}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  roc <- ROC.voting(dat, grp, n=20, healthy="healthy")
}
}
\keyword{misc}
