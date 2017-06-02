\name{pop.seirNewInfec}

\alias{pop.seirNewInfec}

\title{Extract Values of State Variables of each City according to Time.}

\description{Generic function allowing to extract number of S, I, R individuals of each population in a metapopulation of n subpopulations according  to time.}

\usage{pop.seirNewInfec(object, subset, fct,...)}

\arguments{
  \item{object}{a seir object.}
  \item{subset}{numeric. What subpopulations are chosen to extract from the subpopulation set. The subpopulation set is a series of natural number from one to n, where n is the number of subpopulation in a metapopualtion of n subpopulation. Normally, 'subset' is NULL, it allows us to extract all n subpopulations. In the other side, 'subset' is a good numeric, it means that 'subset' is in the series (1,...,n), we only can extract the subpopulations in 'subset' argument.}
  \item{fct}{this is a function such as "sum", "mean",etc. With this function, we can do it with the extracted data from subpopulations.}
  \item{\dots}{other parameters.}
}
\details{Result returned is dependent on the values of the parameter 'subset' and 'fct'. 

	if subset=NULL and fct=NULL, the result is the list of the 'pop' slot in the object given.

	if subset=NULL and fct is a function, the result is a matrix after doing the 'fct' function on the list of the 'pop' slot in the object given.

	if subset!=NULL and fct=NULL, the result is a list of the subpopulations extracted from the list of the 'pop' slot in the object given.

	if subset!=NULL and fct is a function, the result is a matrix after doing the 'fct' function on the list of the subpopulations extracted in the object given.
}
\value{Viewing in detail in the 'detail' part.}

\references{Norman Matloff (2009), The Art of R Programming.}

\author{TRAN Thi Cam Giang}

\seealso{
	'pop.seirNewInfec' fucntion in 'dizzys' package.	
}

\examples{
	obj<-globSEIRSimulNewInfec(nbVilles=3, N=1e6)
	tpobj<-pop(obj) 
	class(tpobj)
	tpobj<-pop(obj,fct="sum")
	class(tpobj)
	tpobj<-pop(obj,subset=c(1,2),fct="sum")	
	class(tpobj)
}
