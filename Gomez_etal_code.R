### The following is a toolkit used to simulate data and to perform ML estimation
### under the poisson mixture model.
###	The functions perform data simulation and ML estimation using the optim function
require(MASS)
#######################################################################
##################							 ##########################
################## Single species estimation ##########################
##################							 ##########################
#######################################################################

## nsites = number of point counts sampled
## nvisits = number of visits to each point count
## lambda = conditional mean abundance (number of individuals) given N>0
# Poisson mixture data simulator.
# This simulator is restricted to find at least one individual per sampling event.
# Sampling event is defined by visitin n sites t times. 

rpois.mix	<-	function(nsites,nvisits,p,lambda){
	
	sum.N	<-	0

	while(sum.N==0){

		N	<-	rpois(nsites,lambda)
		sum.N	<-	sum(N)
		
	}

	y	<-	matrix(NA,nrow=nsites,ncol=nvisits)
	
	sum.y	<-	0
	while(sum.y==0){

		for(i in 1:nvisits){
			y[,i]	<-	rbinom(nsites,N,p)
		}
	
		sum.y	<-	sum(y)
	}

	dimnames(y)	<-	list(nsites=1:nsites,nVisits=1:nvisits)
	return(y)
}


# Simulation of the poison mixture model for K species. 
# the probability of detection for each species comes from a beta distribtuion
# with mean p and dispersion tau, where the parameters of the beta distribution are defined as
# alpha = p*tau ; beta = tau*(1-p) 

rpoisbeta.mix	<-	function(nsites,nvisits,alpha,beta,lambdas){
	
	nspp		<-	length(lambdas)
	spp.counts	<-	array(dim=c(nsites,nvisits,nspp),dimnames=list(sites=1:nsites,visits=1:nvisits,species=lambdas))

	p.vec	<-	rbeta(nspp,alpha,beta)

	for (s in 1:nspp){
	
		sum.abund	<-	0
	
		while(sum.abund==0){
	
			lat.abund	<-	rpois(nsites,lambdas[s])
			sum.abund	<-	sum(lat.abund)
	
		}
	
		sum.spp	<-	0
	
		while(sum.spp==0){
			for(n in 1:nsites){

				spp.counts[n,,s]	<-	rbinom(nvisits,lat.abund[n],p.vec[s])
		
			}
		sum.spp	<-	sum(spp.counts[,,s])
		}

	}
		return(spp.counts)	

}


####### The following functions are used to fit the poison mixture and the beta binomial models

#### Negative log likelihood for the poisson mixture model
#### guess = c(p,lambda)
#### counts= observed number of counts in n sites and r visits. A matrix of n x r
#### K = a value large enough to averge over all the possible values that the lambda can take on

pois.mix.nll	<-	function(guess,counts,K){
	
	p		<-	guess[1]#1/(1+exp(-guess[1]))
	lambda	<-	guess[2]#exp(guess[2])
	
	if(p<0|p>1|lambda<0){nll<-9999999}else{
	
	nsites		<-	nrow(counts)
	nvisits		<-	ncol(counts)
	
	# log likelihood of the detection process	
	det.proc	<-	rep(NA,nsites)
	
	for(i in 1:nsites){
		
		ith.site	<-	counts[i,]
		max.count	<-	max(ith.site)
		k.vec		<-	max.count:K
		ith.k		<-	rep(NA,length(k.vec))
		
		for (j in 1:length(k.vec)){
			
			N.j			<-	k.vec[j]
			ith.k[j]	<-	prod(dbinom(ith.site,N.j,p))*dpois(N.j,lambda)						
			
		}
	
	det.proc[i]	<-	sum(ith.k)
	
	}
	
	like	<-	prod(det.proc)
	nll		<-	-log(like)
	}
	return(nll)
}

##### Maximum Likelihood estimation for the poisson mixture model
##### p = detection probability, value to initialize the optimization
##### lambda = mean number of individuals, value to initialize the optimization 
##### counts = observed number of counts in n sites and r visits. A matrix of n x r
##### K = a value large enough to averge over all the possible values that the lambda can take on
##### method = one of the options defined in the optim function

pois.mix.est	<-	function(p,lambda,counts,K,method="Nelder-Mead"){
		
	guess	<-	c(p,lambda)
	
	ml.ests	<-	optim(guess,pois.mix.nll,method=method,counts=counts,K=K)
	
	p.hat		<-	ml.ests$par[1]
	lambda.hat	<-	ml.ests$par[2]
	negLogLike	<-	ml.ests$value
	AIC			<-	2*length(ml.ests$par)+2*negLogLike
	
	results		<-	c(p=p.hat,lambda=lambda.hat,negLogLike=negLogLike,AIC=AIC)
	return(results)	
}

#### Other useful functions 


#### Checking for the if a number falls in a range 

range.check	<-	function(x,interval){
	interval[1]<=x&x<=interval[2]
}

#### The following fucntions transform the parameters for Ml estimation. 
expit	<-	function(xx) {a<-1/(1 + exp(-xx));return(a)}
logit	<-	function(pp) {a<-log(pp) - log(1-pp);return(a)}

#######################################################################
##################							 ##########################
################## Multi-species estimation  ##########################
##################							 ##########################
#######################################################################

#### The Beta n-mixture model for jags
sink("beta_nmixture.txt")
cat("model{

		# Priors
		# Priors for detectability
		a	~	dunif(0,10)
		b	~	dunif(0,10)
		
		# Priors for abundance
		for(i in 1:nspp){
			lambda[i]	~	dgamma(0.01,0.01)
		}
		
		# likelihood
		
		# Loop along clones		
		
		for(k in 1:K){
		
			# Estimate abundance for each speices
		
			for (i in 1:nspp){
		
				# Sample detection probability from distribution
		
				p[i,k] ~ dbeta(a,b)
				
				# Loop along sites
				for (j in 1:nsites){
					
					# Abundance porcess
					
					N[i,j,k] ~ dpois(lambda[i])
					
					for(t in 1:nvisits){
					
						# Detection process
						counts[j,t,i,k]~dbin(p[i,k],N[i,j,k])
				
					}
				}
			}
		}	
	}
")
sink()
##### Beta N-mixture model with extra level of Hierarchy
sink("beta_nmixture2.txt")
cat("
	model{

		# Priors
		# Priors for detectability
		a	~	dunif(0,10)
		b	~	dunif(0,10)
		# Priors for abundance
		mu.a0 ~ dnorm(0,0.001)	# intercept 
		sigma.a0 ~ dunif(0,10)	# intercept
		tau.a0	<-	pow(sigma.a0,-2)

		# likelihood
		
		# Loop along clones		
		
		for(k in 1:K){
		
			# Estimate abundance for each speices
		
			for (i in 1:nspp){
		
				# Sample abundance from distribution
				a0[i,k] ~ dnorm(mu.a0,tau.a0)
				log(lambda[i,k]) <- a0[i,k]			

				# Sample detection probability from distribution
						
				p[i,k] ~ dbeta(a,b)
				
				# Loop along sites
				for (j in 1:nsites){
					
					# Abundance porcess
					
					N[i,j,k] ~ dpois(lambda[i,k])
					
					for(t in 1:nvisits){
					
						# Detection process
						counts[j,t,i,k]~dbin(p[i,k],N[i,j,k])
				
					}
				}
			}
		}	
	}")
sink()

##### Kalman estimates for the beta N-mixture model
sink("beta_nmixture_kalman.txt")
cat("
	model{

		# Priors
		# Priors for detectability
		# For sampling there are no priors but the actual MLEs are used		
		# Priors for abundance
		
		# likelihood
		
		# Loop along clones		
		
		for(k in 1:K){
		
			# Estimate abundance for each speices
		
			for (i in 1:nspp){
		
				# Sample detection probability from distribution
		
				p[i,k] ~ dbeta(a,b)
				
				# Loop along sites
				for (j in 1:nsites){
					
					# Abundance porcess
					
					N[i,j,k] ~ dpois(lambda[i])
					
					for(t in 1:nvisits){
					
						# Detection process
						counts[j,t,i,k]~dbin(p[i,k],N[i,j,k])
				
					}
				}
			}
		}	
	}
",fill=TRUE)
sink()

##### Yamaura et al 2016 model for jags
sink("yamaura.txt")
cat("
	model{
	# prior distributions on community level estimates
	# mean value
	# parameter related to abundance
		mu.a0 ~ dnorm(0,0.001)	# intercept

	# parameter related to detectability
		mu.p0	~ dnorm(0,0.001)
	
	# standard deviation
		# parameter related to abundance
        sigma.a0 ~ dunif(0,10)	# intercept
        
		# parameter related to detectability
		sigma.p0 ~ dunif(0,10)

		# create precision
		tau.a0 <- pow(sigma.a0,-2)
		tau.p0 <- pow(sigma.p0,-2)

		
		# Loop along clons		
		for(k in 1:K){

			# Estimate parameters fo each species

			for(i in 1:nspp){
		
				# generating parameters of each species related to abundance
				a0[i,k] ~ dnorm(mu.a0,tau.a0)#I(-10,10)

				# Estimating lambda from N which is a latent variable
				log(lambda[i,k]) <- a0[i,k]	# equation (4)	
				
				# generating parameters of each species related to detectability
				p0[i,k] ~ dnorm(mu.p0,tau.p0)#I(-10,10)
				
				#values of p for each species (i) and each site (j)
				p[i,k] <- 1/(1+exp(-(p0[i,k])))	# equation (5)

				

				# Loop along sites (point counts)


				for(j in 1:nsites){
						
					N[i,j,k] ~ dpois(lambda[i,k])	# latent abundance of each species in each year			
					
					# Loop along replicate counts
					for(t in 1:nvisits){
		
						# detection process model
						counts[j,t,i,k] ~ dbin(p[i,k],N[i,j,k])	# detection process in each site			
					
					}
				
				}			
			
			}
		
		}
}",fill=TRUE)
sink()

# For Kalman estimates of the Yamaura model
sink()
cat("
	model{
		# Values of mu.a0, sigma.a0, mu.p0 and sigma.p0 are given by the parameter estimates
		# using data clonning
		# create precision
		tau.a0 <- pow(sigma.a0,-2)
		tau.p0 <- pow(sigma.p0,-2)

		
		# Loop along clons		
		for(k in 1:K){

			# Estimate parameters fo each species

			for(i in 1:nspp){
		
				# generating parameters of each species related to abundance
				a0[i,k] ~ dnorm(mu.a0,tau.a0)#I(-10,10)
				
				# Estimating lambda from N which is a latent variable
				log(lambda[i,k]) <- a0[i,k]	# equation (4)
				
				# generating parameters of each species related to detectability
				p0[i,k] ~ dnorm(mu.p0,tau.p0)#I(-10,10)
				
				#values of p for each species (i) and each site (j)
				p[i,k] <- 1/(1+exp(-(p0[i,k])))	# equation (5)


				# Loop along sites (point counts)


				for(j in 1:nsites){
					
					N[i,j,k] ~ dpois(lambda[i,k])	# latent abundance of each species in each year			
					
					# Loop along replicate counts
					for(t in 1:nvisits){
		
						# detection process model
						counts[j,t,i,k] ~ dbin(p[i,k],N[i,j,k])	# detection process in each site			
					
					}
				
				}			
			
			}
		
		}
	}
",fill=TRUE)
sink()

########################################################################
#########													   #########
######### Assesing the performance of the Beta N-mixture model ######### 
#########													   #########
########################################################################

######### 			Bias Benchmark assesment				   #########

# Load the required packages
library(rjags)
library(snow)


# Define parameters required. These are also predefined in the .RData file with the counts
# The following uses parallel computation to run the simulation and assumes that the computer has
# 20 cores to run parallel. To run the code appropriately find the number of cores in the computer and 
# adjust nreps and star.index accordingly. 

load('bias.RData') # counts with low mean detection probability

nsites	<-	25
nvisits	<-	3
area.50m	<-	pi*0.50^2
alpha		<-	0.25*4.5
beta		<-	0.75*4.5
abundances	<-	c(1,2,3,4,5,7,10,15,25,40,55,65,75,85,100)
nreps		<-	20
K			<-	20
full.res	<-	list()
index	<-	1	# index is taken as an outargument from the batch/shell file used to run the program

start.index	<-	seq(0,500,by=20)
tmp.spp.counts.ls	<-	loP.counts[start.index[index]:(start.index[index]+nreps)]

### model

simulations<-function(r){
	
	# Generate data for each simulation

	spp.counts	<-	tmp.spp.counts.ls[[r]]
	nspp		<-	length(abundances)		

	# Arrange in Clones

	spp.clons	<-	array(NA,dim=c(dim(spp.counts),K))

	for (j in 1:K){
	
		spp.clons[,,,j]<-spp.counts
	}

	max.count	<-	apply(spp.clons,c(3,1,4),max)	
	
	# Gomez et al in review model	

	data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

	inits	<-	function(){}

	params	<-	c("lambda","a","b")
	
	jm	<-	jags.model(file="beta_nmixture.txt"
					,data=data_list,n.chains=2,n.adapt=4000,inits=list(lambda=runif(nspp,0.1,10),a=1,b=1,N=(max.count)))
					
	out	<-	coda.samples(model=jm,variable.names=params,n.iter=15000,thin=20)		

	hats.vec<-summary(out)$statistics[,1]
	list(hats.vec,out)
	
}
	
cl<-makeCluster(nreps,"SOCK")
clusterEvalQ(cl,library(rjags))
clusterExport(cl,list("alpha","beta","nvisits","nsites","abundances","area.50m","K","tmp.spp.counts.ls"))
par.samples<-clusterApply(cl,1:nreps,simulations)
stopCluster(cl)

# For mid values of P run the same previous function but change

tm.spp.counts.ls	<-	midP.counts[start.index[index]:(start.index[index]+nreps)]

alpha	<-	0.5*4.5
beta	<-	0.5*4.5

# For High values of P the change

tmp.spp.counts.ls	<-	highP.counts[start.index[index]:(start.index[index]+nreps)]

alpha	<-	0.75*4.5
beta	<-	0.25*4.5

######### Comparison to other community abundance models  #########
# Load Packages

library(colorRamps)
library(plotrix)

#Functions required to simulate the Negative Binomial (Poisson-Gamma) process

plot.circle	<-	function(center,radius,lty=1,lwd=1,col=NA,border="black"){
	
	xc	<-	center[1]
	yc	<-	center[2]
	angle<-	seq(0,2*pi,length=1000)
	X	<-	(cos(angle)*radius)+xc
	Y	<-	(sin(angle)*radius)+yc
	
	polygon(X,Y,col=col,lwd=lwd,lty=lty,border=border)
	
		
}

euclidean.dist	<-	function(x1,y1,x2,y2){
	
	sqrt((x2-x1)^2+(y2-y1)^2)
	
}

spatial.rnd	<-	function(n,distance,A,radius,n.inits){

	sqrt.A		<-	sqrt(A)
	coords.rnd	<-	matrix(runif(2*n.inits,min=0+radius,max=sqrt.A-radius),ncol=2)
	coords		<-	coords.rnd[1,]
	is.far		<-	TRUE
	pos			<-	1
	for(i in 1:(n-1)){
		ith.pos	<-	pos[i]
		for(j in ith.pos:n.inits){
			tmp.coords	<-	coords.rnd[j,]
			all.dist	<-	dist(rbind(coords,coords.rnd[j,]))
			min.dist	<-	min(all.dist)
			if(min.dist>distance){pos<-c(pos,j);break}
		}
		if(length(pos)<i+1){stop("increase n.inits")}
		else{coords<-coords.rnd[pos,]}
	}
		dimnames(coords)	<-	list(1:nrow(coords),c("x","y"))
		return(coords)		

}
		
Obs.negbin <- function(A,a,k,qu,visits,det.prob,alpha,beta,hetero=TRUE,plot.it,scale.lab,type.qu=c("rect","circ"),distance=NULL,n.inits=NULL){
	
	# to find the scale of the gamma parameter
	inv.beta <- beta^(-1);
	
	# Draw lambdas
  if(hetero==FALSE){lambdas <- rep(k,alpha/beta)}else{
	  lambdas <- rgamma(n=k,shape=alpha,scale=inv.beta);lambdas	<-	sort(lambdas)}
	  
	# Draw the latent variable N from poisson distribution with mean lambda[k]
	# There has to be at least 1 individual of each species in the plot
	
	pois.draw <- rep(0,k);
	for(i in 1:k){pois.draw[i]<- 0;while(pois.draw[i]==0){pois.draw[i]<-rpois(n=1,lambda=lambdas[i]*A)}}
	
	nis <- pois.draw; # totals per species drawn randomly
 	Ntot <- sum(nis);     # total number of individuals in A


 	sqrt.A <- sqrt(A); 
 	sqrt.a <- sqrt(a);
 	
 	#randomly locating individuals in A
 	xy.coords <- matrix(runif(n=Ntot*2,min=0,max=sqrt.A),nrow=Ntot,ncol=2); 
 	xs <- xy.coords[,1];
 	ys <- xy.coords[,2];
 	
 	# select colors by species
 	cls <- blue2red(k);
 	marks <- rep(cls[1:k],nis);
 	spp.numtag <- rep(1:k,nis);
 
 	# Plot the quadrant with the distribution of the individuals

	# Create empty matrix to store number of indviduals of each species in each quadrant
	n.sampled <- matrix(rep(0,qu*k),ncol=k,nrow=qu);
	# Empty vector to store in which quadrat is each individual present
 	which.quadrat <- rep(0,Ntot); 
	
	if(is.null(distance)){distance	<-	0}
	if(is.null(n.inits)){n.inits	<-	1000}

	if(type.qu=="rect"){	
		# locates at random the lower left corners of the quadrats
		 	
		 centroids	<-	spatial.rnd(n=qu,distance=distance,A=A,radius=(sqrt.a/2),n.inits=n.inits)
		 low.left.corners <- centroids-sqrt.a/2;
		 quadrats.xs <- matrix(rep(0,4*qu),nrow=qu,ncol=4);
		 quadrats.ys <- quadrats.xs;
		 	
			
		for(i in 1:qu){
		
			## Draws the corners and counts the number of individuals of each species found in each one
			llc <- low.left.corners[i,];
			x.left<- llc[1];
			y.bott<- llc[2];
			x.right<- llc[1]+sqrt.a;
			y.top  <- llc[2]+sqrt.a;
			      		
			for(j in 1:Ntot){
			
				is.inbox <- ((xs[j]>=x.left)&&(xs[j]<=x.right))&&((ys[j]>=y.bott)&&(ys[j]<=y.top));
	   			if(is.inbox==T){	
	       			which.quadrat[j]<- i;
					which.spp  <- spp.numtag[j];
					n.sampled[i,which.spp] <- n.sampled[i,which.spp] + 1;
    	   		}	
			}
	
  		}
  	}else if(type.qu=="circ"){
		
		# Find the radius and randomly locate the centroids of the plots
		# Centroids can follow a distance rule in which case, distance has to be specified
		# If not specified it will default to 0
		
  		radius		<-	sqrt(a/pi)
 		centroid	<-	spatial.rnd(n=qu,distance=distance,A=1,n.inits=n.inits,radius=radius)
 		
  		for(i in 1:qu){				 
			for(j in 1:Ntot){
				
				# Find if an individual is in a point
				
				is.inpt	<-	(sqrt((xs[j]-centroid[i,1])^2+(ys[j]-centroid[i,2])^2)<=radius)

				if(is.inpt==TRUE){
					
					which.quadrat[j]	<-	i
					which.spp	<-	spp.numtag[j]
					n.sampled[i,which.spp]	<-	n.sampled[i,which.spp] + 1
					
				}
			}
		}
			
	}
		
	counts	<-	array(0,dim=c(qu,visits,k))
	
	for(i in 1:k){
		for(j in 1:qu)
			counts[j,,i]	<-	rbinom(3,n.sampled[j,i],prob=det.prob[i])
	}
	
	
	# For Plot
	
	if (plot.it == TRUE)
 	{	
 		if(is.null(scale.lab)){scale.lab==TRUE}
 		if(scale.lab==TRUE){
 			my.layout	<-	matrix(1:2,ncol=2)
 			layout(my.layout,widths=c(10,1))
 			par(mai=c(0.85,0.85,0.35,0.1),oma=c(1.45,1.25,0.5,0))
 			}
 			else{
 				par(mai=c(0.85,0.85,0.35,0.70),oma=c(1.45,1.25,0.5,0.5));}

   		plot(xs,ys,type="n", main="Study Area",xlab="West-East distance (Kms)", ylab="North-South distance (Kms)"
   				,cex=1.5,cex.lab=1.5, cex.axis=1.25,cex.main=1.5
   				,xlim=c(0,sqrt.A),ylim=c(0,sqrt.A));

   		for(i in 1:Ntot)
   		{
			points(xs[i],ys[i],pch=20,col=marks[i]);
   		}
 		
 		for(i in 1:qu){
 			if(type.qu=="rect"){
 				
 				llc <- low.left.corners[i,];
				x.left<- llc[1];
				y.bott<- llc[2];
				x.right<- llc[1]+sqrt.a;
				y.top  <- llc[2]+sqrt.a;
 				
 				rect(xleft=x.left,ybottom=y.bott,xright =x.right,ytop=y.top,density=NA, col=rgb(0.3,0.2,0.1,alpha=0.1), border="blue");
 				
 			}else if(type.qu=="circ"){
 				
 				 plot.circle(centroid[i,],radius,col=rgb(0.3,0.2,0.1,alpha=0),border="blue")
 				
 			}	
 		}
 		
 		if(scale.lab==TRUE){
 			par(mai=c(0.85,0,0.35,0.5))
 			mat	<-	matrix(1:k,nrow=1)
			image(mat,col=cls,axes=F)
			axis(4,at=seq(from=0,to=1,length=3),labels=c(round(min(lambdas)),round(median(lambdas)),round(max(lambdas))))
 		}
 	
 	}

	
	
	
	#Calculate Nq's = number of species after counting the individuals in the q cuadrats

	Nq.vec <- apply(n.sampled, 1,FUN=function(x){sum(x>0)});
	Ntots.persp <- apply(n.sampled,2,sum);
	result <- list(N=n.sampled,counts=counts,Nq.vec=Nq.vec,Ntot=Ntots.persp,lambdas=lambdas);
	return (result);

} #end of the function

# Simulation of the 500 count data sets for 27 species. Sampling was simulated assiming that the species abundance distrbution of the 
# community follows a gamma distribution with parameters 0.65 and 0.33. This parameters are the best fit parameters for the abundance
# distribution of a 100 ha plot in Panama. Robinson et al 2000. Detection probabilites are assumed to follow a uniform distribution
# between 0 and 1. Sampling of the 100 ha plot was performed with 3 visits to 25, 50 meter radius plots randomly located but with a minimum 
# distance among centroids of the plots o 150 meters. 

spp.counts	<-	list()
lambda.mat	<-	matrix(NA,ncol=30,nrow=500)
p.mat		<-	matrix(NA,ncol=30,nrow=500)
# Simulating an observation
for(i in 1:500){
	p.s		<-	runif(30)
	sim.sa <- Obs.negbin(A=1,a=.0078,k=30,qu=25,visits=3,det.prob=p.s,alpha=0.65,beta=0.033,plot.it=FALSE,scale.lab=FALSE,type.qu="circ",distance=0.15,n.inits=10000);	
	spp.counts.ls[[i]]	<-	sim.sa$counts
	lambda.mat[i,]	<-	sim.sa$lambdas
	p.mat[i,]		<-	p.s
}


## Estimation of the paramters for the 500 simulations under the gamma, uniform model
## The following uses parallel computation to run the simulation and assumes that the computer has
## 20 cores to run parallel. To run the code appropriately find the number of cores in the computer and 
## adjust nreps and star.index accordingly. 

load("comparison.RData")
library(rjags)
library(snow)

nreps		<-	20
K			<-	20
full.res	<-	list()
index		<-	1 # We used a batch file to loop this function such that the index was an outside argument.

start.index	<-	seq(0,500,by=20)
tmp.spp.counts.ls	<-	spp.counts.ls[start.index[index]:(start.index[index]+nreps)]

### model

simulations<-function(r){
	
	# Generate data for each simulation

	spp.counts	<-	tmp.spp.counts.ls[[r]]
	nspp		<-	dim(spp.counts)[3]		

	# Arrange in Clones

	spp.clons	<-	array(NA,dim=c(dim(spp.counts),K))

	for (j in 1:K){
	
		spp.clons[,,,j]<-spp.counts
	}

	max.count	<-	apply(spp.clons,c(3,1,4),max)	
	
	# Gomez et al in review model	

	data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

	params	<-	c("lambda","a","b")
	
	jm	<-	jags.model(file="beta_nmixture.txt"
					,data=data_list,n.chains=2,n.adapt=1000,inits=list(lambda=runif(nspp,0.1,10),a=1,b=1,N=(max.count)))
			 		
	if(is.null(jm)){out	<-	NULL}
	else{
		out	<-	tryCatch(coda.samples(model=jm,variable.names=params,n.iter=5000,thin=1),error=function(e){NULL})
		
		}		
	if(is.null(out)){hats.vec	<-	rep(NA,32)}
	else{hats.vec<-summary(out)$statistics[,1]}

	
	list(hats.vec,out)
	
}
	
cl<-makeCluster(nreps,"SOCK")
clusterEvalQ(cl,library(rjags))
clusterExport(cl,list("K","tmp.spp.counts.ls"))
par.samples<-clusterApply(cl,1:nreps,simulations)
stopCluster(cl)

# To obtain the values of p given the data, replace the prior values in the model with the values of the parameter
# estimates and re-run the process using the original dataset without clonning. In this way, you will 
# be sampling from the posterior distribution given the data. See below for examples with the real data.

# Following is the function used to estimate the parameters under the Normal N-mixture model 

nreps		<-	20
K			<-	1
full.res	<-	list()
index		<-	1 # We used a batch file to loop this function such that the index was an outside argument.

start.index	<-	seq(0,500,by=20)
tmp.spp.counts.ls	<-	spp.counts.ls[start.index[index]:(start.index[index]+nreps)]

simulations<-function(r){
	
	# Generate data for each simulation

	spp.counts	<-	tmp.spp.counts.ls[[r]]
	nspp		<-	dim(spp.counts)[3]		

	# Arrange in Clones

	spp.clons	<-	array(NA,dim=c(dim(spp.counts),K))

	for (j in 1:K){
	
		spp.clons[,,,j]<-spp.counts
	}

	max.count	<-	apply(spp.clons,c(3,1,4),max)	
	
	# Gomez et al in review model	

	data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

	params	<-	c("lambda","mu.a0","sigma.a0","mu.p0","sigma.p0","p")
	
	jm	<-	jags.model(file="yamaura.txt"
					,data=data_list,n.chains=2,n.adapt=1000,inits=list(lambda=runif(nspp,0.1,10),a=1,b=1,N=(max.count)))
			 		
	if(is.null(jm)){out	<-	NULL}
	else{
		out	<-	tryCatch(coda.samples(model=jm,variable.names=params,n.iter=5000,thin=20),error=function(e){NULL})
		
		}		
	if(is.null(out)){hats.vec	<-	rep(NA,64)}
	else{hats.vec<-summary(out)$statistics[,1]}

	
	list(hats.vec,out)
	
}
	
cl<-makeCluster(nreps,"SOCK")
clusterEvalQ(cl,library(rjags))
clusterExport(cl,list("K","tmp.spp.counts.ls"))
par.samples<-clusterApply(cl,1:nreps,simulations)
stopCluster(cl)

######### Example using real data  #########

#### Estimation of abundance of 26 species using both normal-nmixture and beta-nmixture.  
#### The data are saved on a separate .RData file called at the beggining of the code

load("real.RData")

# Estimating the abundance of the species using Beta mixture model
# This uses data cloning with 20 clons.

K		<-	20

spp.counts	<-	UIF.counts
spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
	
	spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list		<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

inits			<-	function(){list(lambda=runif(nspp,0.1,10),a=2,b=2,N=max.count)}

params		<-	c("lambda","a","b")

jm			<-	jags.model(file="beta_nmixture.txt"
					,data=data_list,n.chains=2,n.adapt=1000,inits=inits)
out1			<-	coda.samples(model=jm,variable.names=params,n.iter=20000,thin=20)
	
hats			<-	summary(out1)$statistics[,1]

prop_sigma	<-	1.96*(sqrt(summary(out1)$statistics[,2]^2*20))
upper		<-	(hats+prop_sigma)*100/0.78
lower		<-	(hats-prop_sigma)*100/0.78

out1.sum		<-	summary(out1)$statistics

####### Kalman estimates of p for each species.
K	<-	1

spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
	
	spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K
					,lambda=hats[3:28],a=hats[1],b=hats[2])	

params	<-	c("p")
	
jm.kal.beta	<-	jags.model(file="beta_nmixture_kalman.txt"
						,data=data_list,n.chains=1,n.adapt=1
						,inits=	list(N=max.count))
					
out.kal.beta	<-	coda.samples(model=jm.kal.beta,variable.names=params,n.iter=50000,thin=1)

kal.sum.beta	<-	summary(out.kal.beta)$statistics

# Estimation using the multispecies abundance estimation of Yamaura et al 2016.
# This code uses Bayesian estimation as suggested in the orginial study.
# see below for ML estimation

K	<-	1 

spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
	
	spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

params	<-	c("mu.p0","sigma.p0","mu.a0","sigma.a0","lambda")
	
jm	<-	jags.model(file="yamaura.txt"
					,data=data_list,n.chains=2,n.adapt=1000
					,inits=	list(mu.a0=rnorm(1),mu.p0=rnorm(1),
							sigma.a0=runif(n=1,min=0.5,max=1),
							sigma.p0=runif(n=1,min=0.5,max=1),
							N=max.count))
					
out	<-	coda.samples(model=jm,variable.names=params,n.iter=50000,thin=20)

out.sum	<-	summary(out)$statistics

Cred.int	<-	summary(out)$quantiles[,c(1,5)] # When done using Bayesian Estimation

hats.mods	<-	round(cbind(Cred.int[1:26,1]*100/0.78,out.sum[1:26,1]*100/0.78,Cred.int[1:26,1]*100/0.78
							,lower[3:28],hats[3:28]*100/0.78,upper[3:28]),1)
rownames(hats.mods)	<-	dimnames(spp.counts)[[3]]
colnames(hats.mods)	<-	c("2.5","Mean","97.5","2.5","MLE","97.5")

########## ML estimation using Data Cloning of the Normal N-mixture model

K	<-	20

spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
	
	spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K)	

params	<-	c("mu.p0","sigma.p0","mu.a0","sigma.a0","lambda")
	
jm	<-	jags.model(file="yamaura.txt"
					,data=data_list,n.chains=2,n.adapt=1000
					,inits=	list(mu.a0=rnorm(1),mu.p0=rnorm(1),
							sigma.a0=runif(n=1,min=0.5,max=1),
							sigma.p0=runif(n=1,min=0.5,max=1),
							N=max.count))
					
YDC.out	<-	coda.samples(model=jm,variable.names=params,n.iter=20000,thin=20)

YDC.out.sum	<-	summary(out)$statistics

# Obtaining Kalman estimates for lambda and p for real data under the Normal N-Mixture Model
# This needs to be performed using the estimation of mu and sigma based on MLE (i.e. Data Cloning)

K	<-	1

spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
	
	spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list	<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K
					,mu.a0=YDC.out.sum["mu.a0",1],mu.p0=YDC.out.sum["mu.p0",1]
					,sigma.a0=YDC.out.sum["sigma.a0",1],sigma.p0=YDC.out.sum["sigma.p0",1])	

params	<-	c("lambda","p")
	
jm.kal	<-	jags.model(file="yamaura_kalman.txt"
					,data=data_list,n.chains=2,n.adapt=1000
					,inits=	list(N=max.count))

out.kal	<-	coda.samples(model=jm.kal,variable.names=params,n.iter=20000,thin=20)
kal.sum	<-	summary(out.kal)$statistics
	

# Results of lambda estimates from the two multi-species models and different estimation approaches.

real.data.hats	<-	cbind(YBayesian=out.sum[1:26,1]*100/0.78,YDC=kal.sum[,1]*100/0.78,BDC=out1.sum[3:28,1]*100/0.78)

# Model selection between Normal and Beta N mixture model using Data Cloning Likelihood Ratio
# Sample from the beta N-mixture model

K	<-	1

spp.clons		<-	array(NA,dim=c(dim(spp.counts),K))
nspp			<-	dim(spp.counts)[3]

for (j in 1:K){
	
	spp.clons[,,,j]	<-	spp.counts
}

max.count		<-	apply(spp.clons,c(3,1,4),max)	

data_list		<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K
						,a=hats[1],b=hats[2],lambda=hats[3:28])

params	<-	c("N","p")

beta.sample.LR	<-	jags.model(file="beta_nmixture_kalman.txt"
								,data=data_list,n.chains=1,n.adapt=1000
								,inits=	list(N=max.count))
out.beta.LR		<-	coda.samples(model=beta.sample.LR,variable.names=params,n.iter=20000,thin=20)

#### Organize the samples
nsamps	<-	dim(out.beta.LR[[1]])[1]
npcounts	<-	dim(spp.counts)[1]
nreps	<-	dim(spp.counts)[2]
p.lat		<-	out.beta.LR[[1]][,(2470-25):2470]
n.lat		<-	array(NA,dim=c(nsamps,npcounts,nspp))
start		<-	seq(1,(2470-26),by=nspp)
end		<-	seq(26,(2470-26),by=nspp)
for(i in 1:npcounts){
	ith.start	<-	start[i]
	ith.end	<-	end[i]
	n.lat[,i,]	<-	out.beta.LR[[1]][,ith.start:ith.end]
}

##### likelihood for the beta model

alpha	<-	hats[1]
beta	<-	hats[2]
lambdas	<-	hats[3:28]
logpxy.spp	<-	rep(0,nspp)
logpxy.hat2	<-	rep(0,nsamps)

for(i in 1:nsamps){

	for(ii in 1:nspp){
		
		Yloglikes	<-	matrix(0,ncol=nreps,nrow=npcounts)

		for(j in 1:nreps){
			Yloglikes[,j]	<-	dbinom(UIF.counts[,j,ii],n.lat[i,,ii],p.lat[i,ii],log=T)
		}

		Y.logL	<-	sum(Yloglikes)
		X.logL	<-	sum(dpois(n.lat[i,,ii],lambdas[ii],log=T))+dbeta(p.lat[i,ii],alpha,beta,log=T)
		logpxy.spp[ii]	<-	Y.logL+X.logL
		
	}
	
	logpxy.hat2[i]	<-	sum(logpxy.spp)	
	
}

###### likelihood for the normal model

mu.a0	<-	YDC.out.sum["mu.a0",1]
mu.p0	<-	YDC.out.sum["mu.p0",1]
sigma.a0<-	YDC.out.sum["sigma.a0",1]
sigma.p0<-	YDC.out.sum["sigma.p0",1]

data_list		<-	list(counts=spp.clons,nspp=dim(spp.counts)[3],nsites=dim(spp.counts)[1],nvisits=dim(spp.counts)[2],K=K
						,mu.a0=mu.a0,mu.p0=mu.p0,sigma.a0=sigma.a0,sigma.p0=sigma.p0)

params	<-	c("lambda","N","p")

yamaura.sample.LR	<-	jags.model(file="yamaura_kalman.txt",data=data_list,n.chains=1,n.adapt=1000,inits=list(N=max.count))
out.yamaura.LR		<-	coda.samples(model=yamaura.sample.LR,variable.names=params,n.iter=20000,thin=20)

# Organize samples
p.lat		<-	out.yamaura.LR[[1]][,(2496-25):2496]
lam.lat		<-	out.yamaura.LR[[1]][,2445:2470]
n.lat		<-	array(NA,dim=c(nsamps,npcounts,nspp))
start		<-	seq(1,2444,by=nspp)
end		<-	seq(26,2444,by=nspp)
for(i in 1:npcounts){
	ith.start	<-	start[i]
	ith.end	<-	end[i]
	n.lat[,i,]	<-	out.yamaura.LR[[1]][,ith.start:ith.end]
}

##### likelihood for the normal model

logpxy.spp	<-	rep(0,nspp)
logpxy.hat1	<-	rep(0,nsamps)

for(i in 1:nsamps){

	for(ii in 1:nspp){
		
		Yloglikes	<-	matrix(0,ncol=nreps,nrow=npcounts)

		for(j in 1:nreps){
			Yloglikes[,j]	<-	dbinom(UIF.counts[,j,ii],n.lat[i,,ii],p.lat[i,ii],log=T)
		}

		Y.logL	<-	sum(Yloglikes)
		detect.Lik	<-	dnorm(qlogis(p.lat[i,ii]),mu.p0,sigma.p0,log=TRUE)
		abund.Lik	<-	dnorm(log(lam.lat[i,ii]),mu.a0,sigma.a0,log=TRUE)
		X.logL	<-	sum(dpois(n.lat[i,,ii],lambdas[ii],log=T))+detect.Lik+abund.Lik
		logpxy.spp[ii]	<-	Y.logL+X.logL
		
	}
	
	logpxy.hat1[i]	<-	sum(logpxy.spp)	
	
}

# Likelihood Ratio
llratio	<-	(logpxy.hat2-logpxy.hat1)

# Delta AIC
delta.AIC	<-	-2*(mean(llratio))+2*(28-4)
