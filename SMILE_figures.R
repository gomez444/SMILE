library(plotrix)
library(diagram)

## Conceptual diagram of the SMILE model
png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/Figure1.png",family="Times"
			,width=8,height=8,units="in",res=300,type="cairo")
par(mar=c(0,0,0,0),oma=c(0.5,0.5,0.5,0.5))
plot(0,0,xlim=c(0,9.5),ylim=c(0,5.5),axes=FALSE,type="n")
# Circles and text of compartments
# Susceptible
selfarrow(c(3.85,4.9),curve=c(0.35,0.35*0.6),lwd=1,code=1,path="L",arr.pos=0.58)
text(3.5,4.9,expression(1-lambda(t)),cex=0.9)
draw.circle(3.5,4.5,0.5,col="grey",border="grey")
text(3.5,4.5,"S",cex=2,col="black")

# Infected
draw.circle(6.5,4.5,0.5,col="grey",border="grey")
text(6.5,4.5,"I",cex=2,col="black")

# Immune
selfarrow(c(4.6,2.5),curve=c(0.35,0.35*0.6),lwd=1,code=1,path="L",arr.pos=0.85)
text(4.5,2.5,expression(1-alpha),pos=2)
draw.circle(5,2.5,0.5,col="grey",border="grey")
text(5,2.5,"M",cex=2,col="black")

# Carcasses
draw.circle(8,2.5,0.5,col="grey",border="grey")
text(8,2.5,"L",cex=2,col="black")

# Environment
selfarrow(c(7.6,0.5),curve=c(0.35,0.35*0.6),lwd=1,code=1,path="L",arr.pos=0.85)
text(7.4,0.5,expression(gamma),pos=2)
draw.circle(8,0.5,0.5,col="grey",border="grey")
text(8,0.5,"E",cex=2,col="black")

# Arrows
# S to I
Arrows(4.1,4.5,5.9,4.5)
text(5,4.5,expression(lambda(t)),pos=3)
# I to M
Arrows(6.5,4.15,5,2.85)
text(5.8,3.5,expression(zeta),pos=2)
# I to L
Arrows(6.5,4.15,8,2.85)
text(7.3,3.5,expression(1-zeta),pos=4)
# L to E
Arrows(8,2.15,8,0.85)
text(8,1.5,expression(psi),pos=2)
# M to S
Arrows(5,2.85,3.5,4.15)
text(4.2,3.5,expression(alpha),pos=4)
# S to death
Arrows(3.1,4.2,2,3.5,lty=2)
text(2,3.5,expression(sigma),adj=c(2,2))
# M to death
Arrows(5,2.15,5,1.15,lty=2)
text(5,1.1,expression(sigma),pos=1)

# Reprod to S
Arrows(1.1,4.5,2.9,4.5,lty=2)
text(1.1,4.5,expression(rho(N[t])),pos=2)

# E to lambda
segments(8.5,0.5,9,0.5)
segments(9,0.5,9,5.1)
segments(9,5.1,4.6,5.1)
Arrows(4.6,5.1,4.6,4.7)

dev.off()

## Plotting the case of intermediate dispersal comparing LIZ and suceptibles

png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/Figure2.png",family="Times"
			,width=7.5,height=7,units="in",res=300,type="cairo")
#tiff("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/SMILE_LIZ_suc.tiff",width=9,height=7,units="in",res=300
#	,compression="lzw",type="cairo",family="times")
par(mfrow=c(2,2),mar=c(2,2,1,2),oma=c(2,2,2,2),mgp=c(1.75,0.5,0),tcl=-0.2,cex.axis=1.25)
plot(sim1.1$LIZ,type="l",xaxt="n",xlab="",ylab="",bty="u")
axis(1,seq(0,520,by=52*2),labels=seq(0,10,by=2))
par(new=TRUE)
plot(sim1.1$Suceptibles,type="l",axes=FALSE,xlab="",ylab="",bty="u",col="red",lty=2)
axis(4,ylim=range(sim1.0$Suceptibles))
mtext("A)",adj=0.05)

plot(sim2.1$LIZ,type="l",xaxt="n",xlab="",ylab=""
	,ylim=range(c(sim2.0$LIZ,sim2.1$LIZ,sim2.2$LIZ)),bty="u")
axis(1,seq(0,520,by=52*2),labels=seq(0,10,by=2))
par(new=TRUE)
plot(sim2.1$Suceptibles,type="l",axes=FALSE,xlab="",ylab="",bty="u",col="red",lty=2)
axis(4,ylim=range(sim2.0$Suceptibles))
mtext("B)",adj=0.05)

plot(sim3.1$LIZ,type="l",xaxt="n",xlab="",ylab=""
	,ylim=range(c(sim3.0$LIZ,sim3.1$LIZ,sim3.2$LIZ)),bty="u")
axis(1,seq(0,520,by=52*2),labels=seq(0,10,by=2))
par(new=TRUE)
plot(sim3.1$Suceptibles,type="l",axes=FALSE,xlab="",ylab="",bty="u",col="red",lty=2)
axis(4,ylim=range(sim3.0$Suceptibles))
mtext("C)",adj=0.05)

plot(sim4.1$LIZ,type="l",xaxt="n",xlab="",ylab=""
	,ylim=range(c(sim3.0$LIZ,sim3.1$LIZ,sim3.2$LIZ)),bty="u")
axis(1,seq(0,520,by=52*2),labels=seq(0,10,by=2))
par(new=TRUE)
plot(sim4.1$Suceptibles,type="l",axes=FALSE,xlab="",ylab="",bty="u",col="red",lty=2)
axis(4,ylim=range(sim4.1$Suceptibles))
mtext("D)",adj=0.05)

mtext("Time (Years)",1,outer=TRUE,cex=1.5)
mtext("Number of Deaths",2,outer=TRUE,cex=1.5)
mtext("Suceptibles",4,outer=TRUE,cex=1.5)

par(mar=c(0,0,0,0),fig=c(0,1,0,1),oma=c(0,0,0,0),new=TRUE)
plot(0,0,type="n",axes=FALSE)
legend("top",legend=c("Deaths","Suceptibles"),lty=c(1,2),col=c("black","red"),bty="n",horiz=TRUE,cex=1.25)
dev.off()

# Boxplots with estimates of b0, b1, thau and theta using different sources of information
png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/Figure3.png",family="Times"
			,width=7,height=7,units="in",type="cairo",res=300)
#tiff("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/SMILE_rel_bias.tiff",width=9,height=9,units="in",res=300
#	,compression="lzw",type="cairo",family="times")
	
fig.lab	<-	c("A)","B)","C)","D)")
#c(expression(tau),expression(theta),expression(b[0]),expression(b[1]))

par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,0.5,0),tcl=-0.2,oma=c(2,2,1,1),cex.axis=1.25)
for(i in 1:4){
	y.lim	<-	range(boxplot(rel.bias[,i,],plot=FALSE)$stats)
	plot(0,0,type="n",bty="l",ylim=y.lim,xlim=c(0.5,5.5),xaxt="n",xlab="",ylab="")
	boxplot(rel.bias[,i,],outline=FALSE,xaxt="n",col="grey",axes=FALSE,add=TRUE)
	axis(1,at=1:5,labels=c("SMILE","SMIL","SML","SL","L"),cex.axis=0.9)
	mtext(fig.lab[i],side=3,adj=0.1)
}
mtext("Relative Bias",side=2,outer=TRUE,cex=1.5)
mtext("Type of Information",side=1,outer=TRUE,cex=1.5)
dev.off()

# Montana outbreak
png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/Figure4.png",family="Times",width=7,height=7,units="in",res=300,type="cairo")
par(mar=c(3,3,1,1),oma=c(1,1,1.5,1),mgp=c(1.75,0.5,0),tcl=-0.2)
plot(0,0,type="n",bty="l",xaxt="n",xlim=c(0,53),ylim=c(0,max(montana.stoch))
	,xlab="Time (Weeks)",ylab="Number of Deaths",cex.lab=1.5)

for(i in 1:100){
	
	points(montana.stoch[i,],type="l",col="grey")
	
}
axis(1,at=seq(0,52,by=13))
points(montana.ts,type="l",lwd=2)
points(montana.pred,type="l",col="red",lwd=2)
par(fig=c(0,1,0,1),mar=c(0,3,0,1),oma=c(0,0,0,0),new=TRUE)
plot(0,0,type="n",axes=FALSE,ylab="")
legend("top",c("Observed","Deterministic\nPrediction","Stochastic\nSimulation"),lty=1,col=c("black","red","grey"),horiz=TRUE,bty="n")
dev.off()

#Figure S1
png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/FigureS1.png",width=5,height=5,units="in",res=299,type="cairo")
mat	<-	matrix(c(1,1,2,2
				,1,4,2,5
				,1,1,2,2
				,0,3,3,7
				,0,3,6,7
				,0,3,3,7)
				,ncol=4,nrow=6,byrow=TRUE)
layout(mat,heights=c(1,1,0.25,1,1,0.25))
par(mar=c(3,3,1,1),oma=c(2,2,1,1))

plot(0,0,type="n",xlim=c(0,100),ylim=c(0,100),xaxs="i",yaxs="i")
polygon(x=c(0,100/8,0),y=c(0,100,100),col="black",border="black")
polygon(x=c(0,100,100,100/8),y=c(0,0,100,100),col="grey",border="grey")
abline(a=0,b=1,col="white")
mtext("A)",side=3,adj=0.05,cex=0.7)
box()

plot(0,0,type="n",xlim=c(0,100),ylim=c(0,100),xaxs="i",yaxs="i")
polygon(x=c(0,100/4,0),y=c(0,100,100),col="black",border="black")
polygon(x=c(0,100,100,100/4),y=c(0,0,100,100),col="grey",border="grey")
abline(a=0,b=1,col="white")
mtext("B)",side=3,adj=0.05,cex=0.7)
box()

plot(0,0,type="n",xlim=c(0,100),ylim=c(0,100),xaxs="i",yaxs="i")
polygon(x=c(0,0,100,100),y=c(0,100,100,100/12.5),col="black",border="black")
polygon(x=c(0,100,100),y=c(0,100/12.5,0),col="grey",border="grey")
abline(a=0,b=1,col="white")
mtext("C)",side=3,adj=0.05,cex=0.7)
box()

# Insets

par(mar=c(2.5,2,0,2),mgp=c(0.5,0.1,0),tcl=-0.1,cex.axis=0.5,cex.lab=0.5)
plot(b.vec,type="l",bty="l",xlab="Years",ylab="b",xaxt="n",yaxt="n")
axis(1,at=seq(0,520,length=3),labels=c(0,5,10),mgp=c(0,0,0))
axis(2,at=c(0,0.005,0.01),labels=c("0","0.005","0.01"))
abline(h=b.high,lty=2)

plot(b.vec,type="l",bty="l",xlab="Years",ylab="b",xaxt="n",yaxt="n")
abline(h=b.int,lty=2)
axis(1,at=seq(0,520,length=3),labels=c(0,5,10),mgp=c(1,0,0))
axis(2,at=c(0,0.005,0.01),labels=c("0","0.005","0.01"))

plot(b.vec,type="l",bty="l",xlab="Years",ylab="b",xaxt="n",yaxt="n",col="white",col.lab="white")
abline(h=b.low,lty=2,col="white")
axis(1,at=seq(0,520,length=3),labels=c(0,5,10),mgp=c(1,0,0),col="white",col.axis="white")
axis(2,at=c(0,0.005,0.01),labels=c("0","0.005","0.01"),col="white",col.axis="white")
box(bty="l",col="white")

par(mar=c(0,0,0,0))
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
legend("center",legend=c(expression(R[0]<1),expression(R[0]>1)),fill=c("black","grey"),border=c("black","grey"),bty="n",cex=1.5)
mtext(expression(tau),1,cex=1.5,outer=TRUE)
mtext(expression(theta),2,cex=1.5,outer=TRUE)
dev.off()

# Figure S2
png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/FigureS2.png",family="Times"
			,width=7,height=7,units="in",res=300,type="cairo")
	
#fig.lab	<-	c(expression(tau),expression(theta),expression(b[0]),expression(b[1]))
fig.lab	<-	c("A)","B)","C)","D)")	
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,0.5,0),tcl=-0.2,oma=c(2,2,1,1),cex.axis=1.25)
for(i in 1:4){
	y.lim	<-	range(boxplot(rel.bias1[,i,],plot=FALSE)$stats)
	plot(0,0,type="n",bty="l",ylim=y.lim,xlim=c(0.5,5.5),xaxt="n",xlab="",ylab="")
	boxplot(rel.bias1[,i,],outline=FALSE,xaxt="n",col="grey",axes=FALSE,add=TRUE)
	axis(1,at=1:5,labels=dimnames(rel.bias1)[[3]])
	mtext(fig.lab[i],side=3,adj=0.1)
}
mtext("Relative Bias",side=2,outer=TRUE,cex=1.5)
mtext("Type of Information",side=1,outer=TRUE,cex=1.5)
dev.off()


# Figure S3
png("~/Dropbox/JuanPabloThesis/PostDoc/SMILE/Figures/FigureS3.png",family="Times"
			,width=7,height=7,units="in",res=300,type="cairo")

fig.lab	<-	c("A)","B)","C)","D)")	
#fig.lab	<-	c(expression(tau),expression(theta),expression(b[0]),expression(b[1]))

par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,0.5,0),tcl=-0.2,oma=c(2.5,2,1,1),cex.axis=1.25)
for(i in 1:4){
	y.lim	<-	range(boxplot(rel.bias2.1[,i,],plot=FALSE)$stats)
	plot(0,0,type="n",bty="l",ylim=y.lim,xlim=c(0.5,11.5),xaxt="n",xlab="",ylab="")
	boxplot(rel.bias2.1[,i,],outline=FALSE,xaxt="n",col="grey",axes=FALSE,add=TRUE)
	axis(1,at=1:10,labels=dimnames(rel.bias2.1)[[3]],las=2)
	mtext(fig.lab[i],side=3,adj=0.1)
}
mtext("Relative Bias",side=2,outer=TRUE,cex=1.5)
mtext("Type of Information",side=1,outer=TRUE,cex=1.5,line=0.5)
dev.off()



