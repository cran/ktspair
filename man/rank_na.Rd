\name{rank_na}
\alias{rank_na}
\title{Rank the gene expression and Replace NA}
\description{
This function computes the ranks of the gene expression within each row (within one profile) and replaces the NAs present in the dataset by a given value.
}
\usage{
rank_na(dat, na = -1e+05)
}
\arguments{
  \item{dat}{
A gene expression dataset under the form of a matrix.
}
  \item{na}{
The value by which the NAs should be replaced, by default -1e+05.
}
}

\value{
A matrix of the same size as the matrix dat. The columns of this matrix contain the ranks of the gene expressions within each profile. In other words the ith row of the new
matrix contains the ranks of the gene expression of the ith column of the matrix dat. If NA are present in the matrix dat, they won't be take into account as computing the
ranks and will be replaced by the value na, which is, by default, -1e+05.
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
  rank_na(dat, 2000)
 }
}
\keyword{misc}
