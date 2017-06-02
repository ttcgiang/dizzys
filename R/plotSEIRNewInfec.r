# plot 2D/3D of an seir object
require(graphics)
# this displays the correlation between two in variables (S, E, I, R, P and t)
# S : number of susceptible over time
# E : number of exposed over time
# I : number of infected during one unit of time
# R : number of recovered over time
# P : number of prevalent people over time
# t : time
# x, y,z  are 'number' or 'name' 
#
plot.seirNewInfec <- function(object,x="time",y=4,z=NULL,pop=c(),col="black",type="l",unitTIME=1,
				proj=list(),add=F,xlim=NULL,ylim=NULL,zlim=NULL,xlab=x,ylab="#infected",zlab=z,...){
	# checking coordinates
	if(is.null(x)||is.null(y)) stop("invalid x,y coordinates!")	
	#population
	colPop <- col
	pops.stoch <- object@pop
	#number of population
	nbVilles <- object@nbVilles	
	for(i in 1:nbVilles) pops.stoch[[i]]$time <-pops.stoch[[i]]$time/unitTIME	
	namecol <-  names(pops.stoch[[1]])
	nbcol <- ncol(pops.stoch[[1]])
	indcol <- c(1:nbcol)
	# finding coordinates due to the initial value of x, y, z
	# x
	if(is.character(x)) indx <- which(namecol==x)
	else
		indx <- which(indcol==x)	
	# y
	if(is.character(y)) indy <- which(namecol==y)
	else
		indy <- which(indcol==y)
	# z
	if(is.character(z)) indz <- which(namecol==z)
	else
		indz <- which(indcol==z)
	#
	nameallcol <-c(namecol[indx],namecol[indy],namecol[indz])
	indallcol <-c(indcol[indx],indcol[indy],indcol[indz])
	#plotting
	if(is.integer0(indx)||is.integer0(indy)) stop("Not exist the x or y columne!")
	else
	{
		#calculating xlim and ylim
		pop1 <- pops.stoch[[1]]
		ragX <- range(c(pop1[,indx]))
		minX <- min(pop1[,indx])
		minY <- min(pop1[,indy])
		ragY <- range(c(pop1[,indy]))
		for(i in 1:nbVilles){
			popi <- pops.stoch[[i]]
			ragX<-range(cbind(ragX,range(popi[,indx])))
			minX <- min(minX,popi[,indx])
			ragY<-range(cbind(ragY,range(popi[,indy])))
			minY <- min(minY,popi[,indy])
		}
		
		taugr<-0.04
		if(is.null(xlim)) xlim<-ragX+c(-1,1)*taugr*diff(ragX) else xlim<-xlim
		if(is.null(ylim)) ylim<-ragY+c(-1,1)*taugr*diff(ragY) else ylim<-ylim

		# plotting in 2D
		faire <- FALSE
		if(is.null(z)){# 2D plot
			if(!add) {# creating a new plot
				with(pop1,plot(pop1[,indx],pop1[,indy],type="n",xlim=xlim,ylim=ylim,
								xlab=xlab,ylab=ylab,zlab=zlab,col=col[1],...))
				faire <- TRUE
			}
			else  if(length(dev.list()>=1)) faire <- TRUE	
			if(faire){#adding lines in a plot given
				if(is.null(pop)){#plotting all
					#all curves
					#make coulor for the populations choosen	
					col<-resizeVectorNewInfec(col,nbVilles)				
					for(i in 1:nbVilles){
						popi <- pops.stoch[[i]]
						with(popi,lines(popi[,indx],popi[,indy],type=type,col=col[i],...))
					}					
				}
				else{#plotting with pop
					#make coulor for the populations choosen					
					for(i in 1:length(pop)){
						if(!is.integer0(which(c(1:nbVilles)==pop[i]))){
							popi <- pops.stoch[[pop[i]]]						
							if(!is.null(colPop)){
								#print(is.na(colPop[i]))
								if(is.na(colPop[i])) coli <- col[1]
								else coli <- colPop[i]
								#print(coli)
							with(popi,lines(popi[,indx],popi[,indy],type=type,col=coli,...))
							}
							else
							with(popi,lines(popi[,indx],popi[,indy],type=type,col=col[1],...))
						}
						else
							print(paste("Invalid city ",pop[i],"in the parameter pop"))
					}	
				}
			}
			else stop("Not exist any graphic window!")
		}
		else # plotting in 3D	
		{	
			# zlim for the coordinate 'z'
			ragZ <- range(c(pop1[,indz]))
			minZ <- min(pop1[,indz])
			for(i in 1:nbVilles){popi <- pops.stoch[[i]];	ragZ<-range(cbind(ragZ,range(popi[,indz])));
						minZ <- min(minZ,popi[,indz])}
			if(is.null(zlim)) zlim=ragZ+c(-1,1)*taugr*diff(ragZ) else zlim<-zlim

			if(!add) {# creating a new plot
				with(pop1,plot3d(pop1[,indx],pop1[,indy],pop1[,indz],type="n",xlim=xlim,
							ylim=ylim,zlim=zlim,xlab=xlab,ylab=ylab,zlab=zlab,...))
				faire <- TRUE
			}
			else faire <- TRUE	
			if(faire){#adding new lines in a plot given
				if(is.null(pop)){#plotting all
					#all curves
					#make coulor for the populations choosen
					col<-resizeVectorNewInfec(col,nbVilles)					
					for(i in 1:nbVilles){
						popi <- pops.stoch[[i]]					
						with(popi,lines3d(popi[,indx],popi[,indy],popi[,indz],
											col=col[i],...))
					}						
				}
				else{#plotting with pop
					#make coulor for the populations choosen					
					for(i in 1:length(pop)){
						if(!is.integer0(which(c(1:nbVilles)==pop[i]))){
							popi <- pops.stoch[[pop[i]]]
							#print(head(popi))
							if(!is.null(colPop)){
								#print(is.na(colPop[i]))
								if(is.na(colPop[i])) coli <- col[1]
								else coli <- colPop[i]
								with(popi,lines3d(popi[,indx],popi[,indy],popi[,indz],
				                  col=coli,...))						
							}
							else
								with(popi,lines3d(popi[,indx],popi[,indy],popi[,indz],
											col=col[1],...))
						}
						else
							print(paste("Invalid city ",pop[i],"in the parameter pop"))
						
					}	
				}
			}
			#finding projection of trajectory on the plane  		
			tauxPrj <-(1-0.01)
			nbProj <- length(proj)		
			if(nbProj>=1){	
				for(i in 1:nbProj){
					proji <- proj[[i]]
					nproji <- length(proji)
					if(nproji==2){
						if(is.character(proji[1])) {indProj1 <- which(nameallcol==proji[1])}
						else indProj1 <- which(indallcol==proji[1])

						if(is.character(proji[2])) indProj2 <- which(nameallcol==proji[2])
						else indProj2 <- which(indallcol==proji[2])

						if(is.integer0(indProj1)||is.integer0(indProj2)) 
							stop("invalid projection parametre!")
						else
						if(indProj1==indProj2)
							stop("invalid projection parametre!")
						else{
							if(is.null(pop)){#plotting all	
								col<-resizeVectorNewInfec(col,nbVilles)			
								for(i in 1:nbVilles){						
								popi <- pops.stoch[[i]]
								X<-popi[,indx];Y<-popi[,indy];Z<-popi[,indz]
								X1<-X*0+tauxPrj*minX
								Y1<-Y*0+tauxPrj*minY
								Z1<-Z*0+tauxPrj*minZ
								if((indProj1==1)||(indProj2==1)) X1<-X
								if((indProj1==2)||(indProj2==2)) Y1<-Y
								if((indProj1==3)||(indProj2==3)) Z1<-Z
								lines3d(X1,Y1,Z1,col=col[i],...)
								}
							}
							else{
							#make coulor for the populations choosen					
								for(i in 1:length(pop)){
								popi <- pops.stoch[[pop[i]]]
								X<-popi[,indx];Y<-popi[,indy];Z<-popi[,indz]
								X1<-X*0+tauxPrj*minX; Y1<-Y*0+tauxPrj*minY;Z1<-Z*0+tauxPrj*minZ
								if((indProj1==1)||(indProj2==1)) X1<-X
								if((indProj1==2)||(indProj2==2)) Y1<-Y
								if((indProj1==3)||(indProj2==3)) Z1<-Z
								#print(head(popi))
								if(!is.null(colPop)){
								#print(is.na(colPop[i]))
								if(is.na(colPop[i])) coli <- col[1]
								else coli <- colPop[i]
								lines3d(X1,Y1,Z1,col=coli,...)
								}
								else
								lines3d(X1,Y1,Z1,col=col[1],...)
							}					
						}	
					}			
				}
			#else
			#	print("invalid projection parametre!")
					
				}
			}				
		}
	}
}
#####################################################################################
#
# adding a new lines on a plot given
# this function is used in the 2D/3D plot
#
lines.seirNewInfec <- function(object,x="time",y=4,z=NULL,pop=c(),col="black",type="l",unitTIME=1,
				proj=list(),...){
	plot.seirNewInfec(object=object,x=x,y=y,z=z,pop=pop,col=col,type=type,unitTIME=unitTIME,
				proj=proj,add=T,...)
}
#############################################################################
