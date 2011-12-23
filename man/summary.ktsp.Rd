\name{summary.ktsp}
\alias{summary.ktsp}
\title{
Summary of ktsp object
}
\description{
This function summarizes a ktsp object.
}
\usage{
summary.ktsp(object, select = NULL, printall = FALSE,...)
}
\arguments{
  \item{object}{
A ktsp object.
}
  \item{select}{
To obtain the summary for only a single TSP in the k-TSP.
}
  \item{printall}{
If the summary for the whole method k-TSP or if all the summary of the TSPs should be ploted in succession.
}
  \item{...}{Plotting arguments, not used}
}
\value{
No value is returned. Only a table is ploted.
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

\details{
It is possible to obtain the summary of one TSP present in the k-TSP via the variable select. To obtain the summary of all the TSPs in succession
one should use the variable printall with the value TRUE. By default, a summary of the whole method k-TSP will be ploted.
}

\seealso{
\code{\link{ktspcalc}}, \code{\link{ktspplot}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata) 
  ktsp <- kts.pair(dat,grp,3)
  summary.ktsp(ktsp, select=1)
  summary(ktsp, select=1)
  summary(ktsp, printall=FALSE)
  summary(ktsp, printall=TRUE)
 }
}
\note{
It has the same effect to write "summary.ktsp(x)" or "summary(x)" for a ktsp object x.
}
\keyword{misc}
