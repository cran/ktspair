\name{bootstrap.graphic.ktsp}
\alias{bootstrap.graphic.ktsp}
\title{
Graphical display of the bootstrap procedure
}
\description{
This function displays the results of the bootstrap obtained with the function bootstrap(). It plots the frequency of appearance of the single
genes as well as the pairs. If k was not set by the used an histogram of the values of k is also ploted.
}
\usage{
bootstrap.graphic.ktsp(bootstrap, para1 = 0, para2 = 0, title = NULL, mtext = NULL)
}
\arguments{
  \item{bootstrap}{
  A boostrap object obtained with the function bootstrap().
}
  \item{para1}{
  This parameter allows the user to control the lower bound of the frequency of appearance of the genes (the first plot).
}
  \item{para2}{
  This parameter allows the user to control the lower bound of the frequency of appearance of pairs of the genes (the second plot).  
}
  \item{title}{
  This is used to determine the title of the graph.
}
  \item{mtext}{
  Allow the user to specify a subtitle.
}
}
\details{
The pair(s) of genes with the red dot represent(s) the pair(s) obtained on the original dataset.
}
\value{
None, a graphic is ploted.
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
   bootstrap1 <- bootstrap(dat, grp, k=3, n=20)
   bootstrap2 <- bootstrap(dat, grp, n=30, seed=1)
   bootstrap.graphic.ktsp(bootstrap1)
   bootstrap.graphic.ktsp(bootstrap2)   
  }
}
\keyword{misc}
