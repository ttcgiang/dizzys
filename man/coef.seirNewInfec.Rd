\name{coef.seirNewInfec}

\alias{coef}

\title{Coefficient of seirNewInfec Object}

\description{coef is a generic function which extracts coefficients of parameters from objects returned by modeling functions.}

\usage{
coef(object, ...)
}

\arguments{
  \item{object}{an object for which the extraction of model coefficients is meaningful.}
  \item{\dots}{other arguments.}
}

\details{The function uses an object of class 'seirNewInfec' for extracting coefficients of parameters. There are two parameter types to show. One is the initial values of state variables. The other is the values of parameters used in simulation.}

\value{'coef.seirNewInfec' showing the initial values of state variables and the values of parameters used in simulation for all cities.}

\references{Extract Model Coefficients in R.}

\author{TRAN Thi Cam Giang}

\seealso{The generic function 'coef' in R}

\examples{
	seirobj1<-globSEIRSimulNewInfec(N=1e7)
	coef(seirobj1)

	seirobj2<-globSEIRSimulNewInfec(nbVilles=3, N=c(1e7,1e6))
	coef(seirobj2)
}

\keyword{coefficients }
\keyword{seirNewInfec object }
