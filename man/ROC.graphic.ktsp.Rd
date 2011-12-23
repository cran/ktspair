\name{ROC.graphic.ktsp}
\alias{ROC.graphic.ktsp}
\title{
Graphical display of the ROC curve
}
\description{
This functions displays the results of the ROC curve obtained through the function ROC().
}
\usage{
ROC.graphic.ktsp(roc, m = 1, boxplot = TRUE, AUC = TRUE, auc.x = 0.55, auc.y = 0.25, box.col = "yellow", line.col = "red", multiple.col = c("red", "orange", "blue"), maintitle=NULL, mtext = NULL, undertitle = NULL, graphic = TRUE)
}
\arguments{
  \item{roc}{
  A roc object obtained with the function ROC().
}
  \item{m}{
  Allows the user to determine the proportion (1/m) of boxplot that should be ploted.
}
  \item{boxplot}{
  A logical integer that determines if the boxplot should be ploted or not (if FALSE, the values are ploted as points).
}
  \item{AUC}{
  If the AUC should be computed or not.
}
  \item{auc.x}{
  If AUC=TRUE, the x coordinate where the AUC value will be ploted.
}
  \item{auc.y}{
  If AUC=TRUE, the y coordinate where the AUC value will be ploted.
}
  \item{box.col}{
  The color of the boxplots.
}
  \item{line.col}{
  The color of the line that represent the median of the sensitivity for several values of the specificity.
}
  \item{multiple.col}{
  If multiple cutoff were used to construct the roc object (in the function ROC()), the colors that should be used on the summary graph.
}
  \item{maintitle}{
  Allow the user to specify a title.
}
  \item{mtext}{
  Allow the user to specify a subtitle.
}
  \item{undertitle}{
  Allow the user to control the title display in the single graph display more accurately.
}
  \item{graphic}{
  If a grahic should be displayed or not.
}
}
\details{
In the function ROC() a bootstrap procedure was used to construct several representations of ROC curves. In the function ROC.graphic(), these curves are used to create boxplot around the values of the sensitivity for several values of the specificity. If, for a given value of the specificity, several sensitivities were available, the median was used to represent this point and the different values of the sensitivity were used to draw a boxplot.
}
\value{
A ROC.graphic object with the following elements
  \item{at}{Either a vector or a 3 rows matrix (if mult.cutoff=FALSE or TRUE resp.) representing the x-axis of the curve for median of the sensitivity (the red line).}
  \item{median}{Either a vector or a 3 rows matrix (if mult.cutoff=FALSE or TRUE resp.) representing the y-axis of the curve for median of the sensitivity (the red line).}
  \item{auc}{Either a single value or a vector of length 3 (if mult.cutoff=FALSE or TRUE resp.) standing for the AUC of the model(s).}
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
  roc1 <- ROC(dat, grp, n=200, healthy="healthy", mult.cutoff=FALSE)
  roc2 <- ROC(dat, grp, n=200, healthy="healthy", mult.cutoff=TRUE)
  ROC.graphic.ktsp(roc1)
  ROC.graphic.ktsp(roc2)
 }
}
\keyword{misc}
