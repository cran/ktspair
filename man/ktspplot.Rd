\name{ktspplot}
\alias{ktspplot}
\title{
Graphical representation of ktsp objects
}
\description{
Show a graphical representation of ktsp objects based on the expressions of the genes present in the pairs constituting the k-TSP.
}
\usage{
ktspplot(ktspobj, select = NULL)
}
\arguments{
  \item{ktspobj}{
  A ktsp object (computed through the functions ktspcalc or kts.pair).
  }
  \item{select}{
  If the number of TSPs in the k-TSP is greater than 1, one can choose the representation of a single TSP (select = 1,2,...,k).
  By default the graphical representation of every TSPs in succession is ploted.
  }
}
\details{
ktspplot() creates a two dimensional graph for a single TSP present in the k-TSP. The graph plots the expression of the first gene in the pair
versus the expression of the second gene in the pair of interest. The group of the observations are ploted in red and blue, the score of the 
current pair is shown on the top of the graph. The black line of slope 1 and intercept 0 represents the decision rule of the current TSP. 
If not only one TSP is selected, hitting return will plot the next TSP in the k-TSP.
}
\value{
No value is returned. Only a graph is ploted.
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
  ktsp <- ktspcalc(dat,grp,3)
  ktspplot(ktsp,select=1)
  ktspplot(ktsp)
 }
}
\keyword{misc}
