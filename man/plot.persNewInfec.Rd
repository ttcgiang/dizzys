\name{plot.persNewInfec}

\alias{plot}

\title{Plotting Kaplan Meier Survival Curve}

\description{Plotting Kaplan–Meier survival curve based on persistence value or adding one or more straight lines through the current
     plot based on the slot 'pop' of seir object.}

\usage{plot(object, x, y, type, col, xlim, ylim, curvetype, vilabline,...)}

\arguments{
  \item{object}{a seir object with the value of the slot "persistence".}

  \item{x, y}{the names or the positions of the columns in the slot 'pop' of the seir object. This arguments are corresponding to the x and y coordinates of the plot. In this function, we only plot lines 2D, 'x, y' are one of the names or the positions of the columns ("time","S","E","P","R","N") or (1,2,3,4,5,6)}
  
  \item{type}{what type of plot should be drawn. Possible types are "p", "l", "b", "c", ect.}

  \item{col}{the colors for lines. Multiple colors can be specified so that each line can be given its own color. If there are fewer colors than number of lines they are recycled in the standard fashion. Lines will all be plotted in the first colour specified.}

  \item{xlim}{the x limits of the plot.}

  \item{ylim}{the y limits of the plot.}

  \item{curvetype}{type of curve for plot. There are two value for 'curvetype'. If curvetype=="KM" (Kaplan–Meier), it means that we will plot Kaplan–Meier survival curve based on the slot 'persistence' of the object. In the other side, if curvetype=="population", it means that we will plot all lines by using the slot 'pop' of the object with  straight lines.}

  \item{vilabline}{numeric. What subpopulations are chosen to add straight lines to plot, from the subpopulation set. The subpopulation set is a series of natural number from one to n, where n is the number of subpopulation in a metapopualtion of n subpopulation. Normally, 'vilabline' is NULL, we don't add any straight lines to the plot. In the other side, 'vilabline' is a good numeric, it means that 'vilabline' is in the series (1,...,n), we only can add the lines of the subpopulation in 'vilabline' to the plot.}

  \item{\dots}{other arguments.}

}
\details{This function allows us to plot two types of curve to the plot. Plotting Kaplan–Meier survival curve based on persistence value or adding one or more straight lines through the current plot based on the slot 'pop' of seir object.}

\references{ David G. Kleinbaum and Mitchel Klein, Survival Analysis.}

\author{TRAN Thi Cam Giang}

\seealso{Kaplan–Meier survival curve}

\examples{
	p<-persNewInfec(globSEIRSimulNewInfec(nbVilles=5,N=1e5,nbCONTACT0=100,duration=365*30))
	plot.persNewInfec(p)
	x11()
	plot.persNewInfec(p,curvetype="pop",col=c("green","blue"),vilabline=c(1,3))
}

\keyword{ Kaplan–Meier curve }
\keyword{ Kaplan–Meier estimator }% __ONLY ONE__ keyword per line
