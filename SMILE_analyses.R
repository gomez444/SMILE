# Testing the SMILE functions
# Load SMILE source code before running this script.

####################################################################################
# Case 1: No age structuring, infection probability is not seasonal
tau=10
b=0.001

# Case 1.0: Low Dispersion. Mean dispersion effort is 0.1
theta=100
sim1.0	<-	smile1(b=b,tau=tau,theta=theta,years=10)

# Case 1.1: Medium Dispersion. Mean disperion effort is 1
theta=10
sim1.1	<-	smile1(b=b,tau=tau,theta=theta,years=10)

# Case 1.2: High Dispersion. Mean dispersion effort is 10
theta = 1
sim1.2	<-	smile1(b=b,tau=tau,theta=theta,years=10)

####################################################################################
# Case 2: No age structuring, infection probability has seasonal forcing
tau=10
b0=-30
b1=0.85
period=3*52

# Case 2.0: Low dispersion
theta=100
sim2.0	<-	smile2(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

# Case 2.1: Intermediate dispersion
theta=10
sim2.1	<-	smile2(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

# Case 2.2: High dispersion
theta = 1
sim2.2	<-	smile2(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

#######################################################################################
# Case 3: With births and no seasonal infection probability
tau=10
b=0.001

# Case 3.0: Low Dispersion. Mean dispersion effort is 0.1
theta=100
sim3.0	<-	smile3(b=b,tau=tau,theta=theta,years=10)

# Case 3.1: Medium Dispersion. Mean dispersion effort is 1
theta = 10
sim3.1	<-	smile3(b=b,tau=tau,theta=theta,years=10)

# Case 3.2: Medium Dispersion. Mean dispersion effort is 10
theta = 1
sim3.2	<-	smile3(b=b,tau=tau,theta=theta,years=10)

#######################################################################################
# Case 4: With births and seasonal forcing of infection
tau=10
b0=-30
b1=0.85
period=3*52

# Case 4.0: Low Dispersion. Mean dispersion effort is 0.1
theta=100
sim4.0	<-	smile4(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

# Case 4.1: Medium Dispersion. Mean dispersion effort is 1
theta = 10
sim4.1	<-	smile4(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

# Case 4.2: Medium Dispersion. Mean dispersion effort is 10
theta = 1
sim4.2	<-	smile4(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

#######################################################################################
# Case 5: With births seasonal forcing, natural deaths and stochasticity to the disease
# dynamics.
tau=10
b0=-30
b1=0.85
period=3*52

# Case 5.0: Low Dispersion. Mean dispersion effort is 0.1
theta=10
sim5.0		<-	smile5(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)

stochastic.sims	<-	list()
LIZ.stoch		<-	matrix(0,nrow=100,ncol=520)
for(i in 1:100){

	stochastic.sims[[i]]	<-	smile5.stoch(b0=b0,b1=b1,period=period,tau=tau,theta=theta,years=10)
	LIZ.stoch[i,]			<-	stochastic.sims[[i]]$LIZ
	
}

SMILE.obs	<-	stochastic.sims
hat.mat			<-	array(0,dim=c(nreps,5,5),dimnames=list(1:nreps
															 ,c("tau","theta","b0","b1","loglik")
															 ,c("SMILE","SMIL","SML","SL","L")))

#######################################################################################
# Using all information
for(i in 1:nreps){
	
	hat.mat[i,,1]	<-	SMILE.param.estim(b0=b0.start,b1=b1.start
										,theta=theta.start,tau=tau.start,SMILE.obs=SMILE.obs[[i]],period=period)
								}

#######################################################################################								
# Using SMIL
SMIL.obs	<-	lapply(SMILE.obs,function(x)x[1:4])
for(i in 1:nreps){
	
	hat.mat[i,,2]	<-	SMILE.param.estim(b0=b0.start,b1=b1.start
										,theta=theta.start,tau=tau.start,SMILE.obs=SMILE.obs[[i]],period=period)
								}

#######################################################################################
# Using SML
SML.obs		<-	lapply(SMILE.obs,function(x)x[c(1,2,4)])
for(i in 1:nreps){
	
	hat.mat[i,,3]	<-	SMILE.param.estim(b0=b0.start,b1=b1.start
										,theta=theta.start,tau=tau.start,SMILE.obs=SML.obs[[i]],period=period)
								}

#######################################################################################
# Using SL
SL.obs		<-	lapply(SMILE.obs,function(x)x[c(1,4)])
for(i in 1:nreps){
	
	hat.mat[i,,4]	<-	SMILE.param.estim(b0=b0.start,b1=b1.start
										,theta=theta.start,tau=tau.start,SMILE.obs=SL.obs[[i]],period=period)
								}

#######################################################################################
# Using L only
L.obs		<-	lapply(SMILE.obs,function(x)x[4])
for(i in 1:nreps){
	
	hat.mat[i,,5]	<-	SMILE.param.estim(b0=b0.start,b1=b1.start
										,theta=theta.start,tau=tau.start,SMILE.obs=L.obs[[i]],period=period)
								}

rel.bias	<-	array(0,dim=c(100,4,5),dimnames=list(1:nreps
													,c("tau","theta","b0","b1")
													,c("SMILE","SMIL","SML","SL","L")))
for(i in 1:5){
	for(j in 1:100){
		
		rel.bias[j,,i]	<-	(hat.mat[j,1:4,i]-c(10,10,-30,0.85))/c(10,10,-30,0.85)
		
	}
}



############################################################################################################
######## Now I am going to take out one time series at a time
hat.mat1	<-	array(0,dim=c(nreps,5,5),dimnames=list(1:nreps
										,c("tau","theta","b0","b1","loglik")
										,c("MILE","SILE","SMLE","SMIE","SMIL")))
rel.bias1	<-	array(0,dim=c(100,4,5),dimnames=list(1:nreps
						,c("tau","theta","b0","b1")
						,c("MILE","SILE","SMLE","SMIE","SMIL")))

for(j in 1:5){
	jth.obs	<-	lapply(SMILE.obs,function(x)x[-j])
	for(i in 1:nreps){
		hat.mat1[i,,j]	<-	SMILE.param.estim(b0=b0.start,b1=b1.start
											,theta=theta.start,tau=tau.start
											,SMILE.obs=jth.obs[[i]],period=period)
		rel.bias1[i,,j]	<-	(hat.mat1[i,1:4,j]-c(10,10,-30,0.85))/c(10,10,-30,0.85)
		
	}	
}

############################################################################################################
###### Now taking two time series at a time
which.sel	<-	expand.grid(1:5,1:5)
which.sel	<-	which.sel[-which(which.sel[,1]==which.sel[,2]),]

hat.mat2	<-	array(0,dim=c(nreps,5,nrow(which.sel)))
rel.bias2	<-	array(0,dim=c(nreps,4,nrow(which.sel)))


pb	<-	txtProgressBar(min=0,max=nrow(which.sel),style=3)
for(j in 1:nrow(which.sel)){
	to.del	<-	as.numeric(which.sel[j,])
	jth.obs	<-	lapply(SMILE.obs,function(x)x[-to.del])


	for(i in 1:nreps){
		hat.mat2[i,,j]	<-	tryCatch(SMILE.param.estim(b0=b0.start,b1=b1.start
											,theta=theta.start,tau=tau.start
											,SMILE.obs=jth.obs[[i]],period=period)
											,error=function(e) rep(NA,5))
		rel.bias2[i,,j]	<-	(hat.mat2[i,1:4,j]-c(10,10,-30,0.85))/c(10,10,-30,0.85)
	}
	setTxtProgressBar(pb,j)
}
close(pb)

all			<-	unlist(strsplit("SMILE",split=""))
info.names	<-	NULL
for(i in 1:nrow(which.sel)){
	to.del	<-	as.numeric(which.sel[i,])
	jth.nam	<-	all[-to.del]
	info.names	<-	c(info.names,paste0(jth.nam[1],jth.nam[2],jth.nam[3]))
}

dimnames(hat.mat2)	<-	list(1:nreps
							,c("tau","theta","b0","b1","loglik")
							,info.names)
dimnames(rel.bias2)	<-	list(1:nreps
							,c("tau","theta","b0","b1")
							,info.names)
							
hat.mat2.1	<-	hat.mat2[,,which(!duplicated(info.names))]
rel.bias2.1<-	rel.bias2[,,which(!duplicated(info.names))]
info.names	<-	info.names[which(!duplicated(info.names))]

############################################################################################################

# Estimating the outbreaks for Montana
raw.dat				<-	read.delim("Montana_TS.txt")
montana.ts			<-	rep(0,52)
montana.ts[30:36]	<-	tapply(raw.dat$Deaths,raw.dat$Week,sum)

n.years			<-	seq(1,20,length=100)
montana.hats	<-	matrix(NA,nrow=100,ncol=5,dimnames=list(n.years,c("tau","theta","b0","b1","loglik")))

for(j in 1:100){
	montana.hats[j,]	<-	SMILE.param.estim(b0=1,b1=1,theta=log(1+0.1),tau=log(1+0.1),period=n.years[j]*52
											,SMILE.obs=list(LIZ=montana.ts))
}

rel.lik	<-	exp(montana.hats[,"loglik"])/max(exp(montana.hats[,"loglik"]))
best.mod<-	which.max(rel.lik)
best.per<-	n.years[best.mod]

montana.stoch	<-	array(NA,dim=c(100,1*52))

montana.pred	<-	smile5(b0=montana.hats[best.mod,"b0"],b1=montana.hats[best.mod,"b1"],period=best.per*52
								,theta=montana.hats[best.mod,"theta"],tau=montana.hats[best.mod,"tau"],years=1)$LIZ

for(i in 1:100){
	montana.stoch[i,]	<-	smile5.stoch(b0=montana.hats[best.mod,"b0"],b1=montana.hats[best.mod,"b1"],period=best.per*52
										,theta=montana.hats[best.mod,"theta"],tau=montana.hats[best.mod,"tau"],years=1)$LIZ
}

# Calculation of R0
E.montana	<-	smile5(b0=montana.hats[best.mod,"b0"],b1=montana.hats[best.mod,"b1"],period=best.per*52
								,theta=montana.hats[best.mod,"theta"],tau=montana.hats[best.mod,"tau"],years=1)$Environment
b.montana	<-	b.season(b0=montana.hats[best.mod,"b0"],b1=montana.hats[best.mod,"b1"],period=best.per*52,t=1:52)
R0.montana	<-	r0(b=b.montana,E=E.montana,theta=montana.hats[best.mod,"theta"],tau=montana.hats[best.mod,"tau"])

# Exploring the parameter space for R0 assuming high number of b and E and varying theta and tau

b.vec		<-	b.season(b0,b1,3*52,1:520)
E.vec		<-	sim5.0$Environment

theta.vec	<-	seq(0.01,100,length=1000)
tau.vec		<-	seq(0.01,100,length=1000)
b.high		<-	0.01
E.high		<-	800


r0.mat1		<-	matrix(NA,nrow=1000,ncol=1000)

for(i in 1:1000){
	
	for(j in 1:1000){
		
		r0.mat1[i,j]	<-	r0(b.high,E.high,theta.vec[i],tau.vec[j])
		
	}
	
}

# Assuming low b and low E

b.low	<-	0.0001
E.low	<-	800

r0.mat3	<-	matrix(NA,nrow=1000,ncol=1000)

for(i in 1:1000){
	
	for(j in 1:1000){
		
		r0.mat3[i,j]	<-	r0(b.low,E.low,theta.vec[i],tau.vec[j])
		
	}
	
}

# Assuming intermiediate b and E

b.int	<-	0.005
E.int	<-	800

r0.mat2	<-	matrix(NA,nrow=1000,ncol=1000)

for(i in 1:1000){
	
	for(j in 1:1000){
		
		r0.mat2[i,j]	<-	r0(b.int,E.int,theta.vec[i],tau.vec[j])
		
	}
	
}

r0.mat1[r0.mat1>1]	<-	1
r0.mat1[r0.mat1<1]	<-	0

r0.mat2[r0.mat2>1]	<-	1
r0.mat2[r0.mat2<1]	<-	0

r0.mat3[r0.mat3>1]	<-	1
r0.mat3[r0.mat3<1]	<-	0