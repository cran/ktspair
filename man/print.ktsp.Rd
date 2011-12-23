\name{print.ktsp}
\alias{print.ktsp}

\title{
Print the results of the k-TSP
}
\description{
This function is used to print the results of a ktsp object.
}
\usage{
print.ktsp(x,...)
}

\arguments{
  \item{x}{
A ktsp object.
}
  \item{...}{Plotting arguments, not used}
}
\details{
It prints the results of the k-TSP. More precisely, the number of pairs used in the k-TSP, the gene present in the pairs and the score of each pair.
}
\value{
No value is returned. A table is ploted.
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
\code{\link{ktspcalc}}, \code{\link{kts.pair}},\code{\link{predict.ktsp}}, \code{\link{summary.ktsp}}
}
\examples{
  \dontrun{
  ## Load data
  data(ktspdata)
  ktsp <- ktspcalc(dat, grp, 3)
  ktsp.print(ktsp)
  ktsp
  }
}
\note{
Once the function ktsp.print is defined, writing "x" has the same effect as writing "print.ktsp(x)", where x is a ktsp object.
}
\keyword{misc}
