\name{plot.seirNewInfec}

\alias{plot}

\title{Plotting 2D/3D a seir Object }

\description{Generic function of 'seir' class for plotting 2D/3D of seir objects.}

\usage{
	plot(object,x="time",y=4,z=NULL,pop=c(),col="black",type="l",unitTIME=1,
				proj=list(),add=F,xlim=NULL,ylim=NULL,zlim=NULL,xlab=x,ylab="#infected",zlab=z)
}

\arguments{
  \item{object}{ seir object.}
  \item{x, y, z}{the names or the positions of the columns in the slot 'pop' of the seir object. This arguments are corresponding to the x, y and z coordinates of the plot. Normally, if 'z' is NULL, we only can add a line 2D to a plot 2D, on the other side, if 'z' is one of the names or the positions of the columns ("time","S","E","P","R","N") or (1,2,3,4,5,6), we can add a line 3D to a plot 3D.}
  \item{pop}{numeric. What subpopulations are chosen to add from the subpopulation set. The subpopulation set is a series of natural number from one to n, where n is the number of subpopulation in a metapopualtion of n subpopulation. Normally, 'pop' is NULL, it allows us to add all lines of the n subpopulations to the plot. In the other side, 'pop' is a good numeric, it means that 'pop' is in the series (1,...,n), we only can add the lines of the subpopulation in 'pop' to the plot.}
  \item{col}{the colors for lines. Multiple colors can be specified so that each line can be given its own color. If there are fewer colors than number of lines they are recycled in the standard fashion. Lines will all be plotted in the first colour specified.}
  \item{type}{what type of plot should be drawn. Possible types are "p", "l", "b", "c", ect.}
  \item{unitTIME}{the unit of time. It is 'day' corresponding to unitTIME=1, is 'week' corresponding to unitTIME=7, is 'month' corresponding to unitTIME=30 and is 'year' corresponding to unitTIME=365. Normally, unitTIME is equal to 1, however we can change it.}
  \item{proj}{list of the plane names. This argument is only used in the plot 3D. It allows us to add the projection of the lines 3D on the planes. It is necessary to give the value of 'proj' corresponding to the x, y and z axes.}
  \item{add}{adding lines of the object to a current plot, if add=T. In contrast, if add=F, we plot lines of the object to a new plot. By default, add=F.}
  \item{xlim, ylim, zlim}{limits to use for the coordinates.}
  \item{xlab, ylab, zlab}{labels for the coordinates.}
  \item{\dots}{further graphical parameters such as 'lwd', 'lty'}
}

\details{Moreover, for the x, y, and z coordinates, the slot 'pop' of the seir object is a list of the data frames corresponding to each subpopulation in a metapopulation of n subpopulations. The each data frame often has six columns that are 'time' being a time column, at time t, 'S' being the number of susceptibles, 'E' being the number of exposed individuals, 'P' being the number of infectives, 'R' being the number of removed individuals with immunity and 'N' being the number of population, it is the sum of S, E,P,R.

In the set of arguments, there are some arguments we can give each subpopulation each value's.

  \bold{S, E, I, R, N}: the initial values of state variables for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own state variables. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.

  \bold{mu, nbCONTACT0, nbCONTACT1, probVISITER, probINFECTER, sigma, gamma, phi}: the parameters of simulation for 'nbVilles' subpopulations. Multiple values can be specified so that each subpopulation can be given its own parameters. If there are fewer values than subpopulations they are recycled in the standard fashion. Subpopulations will all be simulated in the first value specified.
}

\references{Generic X-Y Plotting in R.
	     plot.default in R.}

\author{TRAN Thi Cam Giang}

\seealso{ 'plot.seirNewInfec' in 'dizzys' package.}

\examples{
	obj<-globSEIRSimulNewInfec(nbVilles=3, N=1e6, nbCONTACT0=100)
	plot(obj,col=c("red","blue"),lwd=2,xlab="time (day)", ylab="number of infectives")
	pause()
	plot(obj,z="S",col=c("red","blue"),lwd=2,xlab="time (day)", ylab="number of infectives",zlab="number of susceptible")
	pause()
	#plot(obj,z="S",col=c("red","blue"),lwd=2,proj=list(c("time","P"),c("time","S")),box=F,xlab="time (day)", ylab="number of infectives",zlab="number of susceptible")
}

\keyword{ seir model}
