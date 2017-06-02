\name{lines.seirNewInfec}

\alias{lines}

\title{Add Connected Line Segments of an seir object to a Plot 2D/3D}

\description{A generic function taking coordinates given of the slot 'pop' of a seir object in various ways and joining the corresponding points with line segments.}

\usage{lines(object, x, y, z, pop, col, type, unitTIME, proj,...)}

\arguments{

  \item{object}{a seir object, this is the important parameter.}

  \item{x, y, z}{the names or the positions of the columns in the slot 'pop' of the seir object. This arguments are corresponding to the x, y and z coordinates of the plot. Normally, if 'z' is NULL, we only can add a line 2D to a plot 2D, on the other side, if 'z' is one of the names or the positions of the columns ("time","S","E","P","R","N") or (1,2,3,4,5,6), we can add a line 3D to a plot 3D.}

  \item{pop}{numeric. What subpopulations are chosen to add from the subpopulation set. The subpopulation set is a series of natural number from one to n, where n is the number of subpopulation in a metapopualtion of n subpopulation. Normally, 'pop' is NULL, it allows us to add all lines of the n subpopulations to the plot. In the other side, 'pop' is a good numeric, it means that 'pop' is in the series (1,...,n), we only can add the lines of the subpopulation in 'pop' to the plot.}

  \item{col}{the colors for lines. Multiple colors can be specified so that each line can be given its own color. If there are fewer colors than number of lines they are recycled in the standard fashion. Lines will all be plotted in the first colour specified.}

  \item{unitTIME}{the unit of time. It is 'day' corresponding to unitTIME=1, is 'week' corresponding to unitTIME=7, is 'month' corresponding to unitTIME=30 and is 'year' corresponding to unitTIME=365. Normally, unitTIME is equal to 1, however we can change it.}
  \item{type}{what type of plot should be drawn. Possible types are: "l", "b", "s",...}

  \item{proj}{list of the plane names. This argument is only used in the plot 3D. It allows us to add the projection of the lines 3D on the planes. It is necessary to give the value of 'proj' corresponding to the x, y and z axes.}

  \item{...}{further graphical parameters such as 'lwd', 'lty'.}

}

\details{The function 'lines' is the partial analogue of plot.seir with the argument add=T.
	Moreover, for the x, y, and z coordinates, the slot 'pop' of the seir object is a list of the data frames corresponding to each subpopulation in a metapopulation of n subpopulations. The each data frame often has six columns that are 'time' being a time column, at time t, 'S' being the number of susceptibles, 'E' being the number of exposed individuals, 'P' being the number of infectives, 'R' being the number of removed individuals with immunity and 'N' being the number of population, it is the sum of S, E, P, R.}

\references{Murrell, P. (2005) R Graphics. Chapman & Hall/CRC Press.}

\author{TRAN Thi Cam Giang}

\seealso{}

\examples{
#creating a plot
#adding a line to the plot
	globSEIRSimulNewInfec(nbVilles=2)->obj
	globSEIRSimulNewInfec(nbVilles=1)->obj1
	#2D
	plot(obj,col="red")
	lines(obj1,col="blue",lwd=2)
	#3D
	plot(obj,z="S",col="red",proj=list(c("time","P")))
	lines(obj1,z="S",col="blue",proj=list(c("time","P")))
}
\keyword{ lines }
\keyword{ projection on plane }% __ONLY ONE__ keyword per line
