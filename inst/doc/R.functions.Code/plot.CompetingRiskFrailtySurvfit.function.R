"plot.CompetingRiskFrailtySurvfit" <-
function(x,...)
{

#further plotting parameters
cex.main<-0.9
tcl<- -0.1
cex.axis<-0.7
las<-1


#object values
p<-x$p
factor.names<-x$factor.names
grid.frame<-x$grid.frame
varying.list<-x$varying.list
deviation.list<-x$deviation.list
risk.names<-x$risk.names


#number of Competing.Risks
L<-length(risk.names)



#plot baselines
for (l in 1:L)
{
assign(paste("alpha.Baseline.",l,sep=""),varying.list[[l]][,1])
assign(paste("deviation.Baseline.",l,sep=""),deviation.list[[l]][,1])
}
y.range.Baseline<-range(c(unlist(mget(paste("alpha.Baseline.",1:L,sep=""),envir=as.environment(-1)))-unlist(mget(paste("deviation.Baseline.",1:L,sep=""),envir=as.environment(-1))),unlist(mget(paste("alpha.Baseline.",1:L,sep=""),envir=as.environment(-1)))+unlist(mget(paste("deviation.Baseline.",1:L,sep=""),envir=as.environment(-1)))))

for (l in 1:L)
{
plot(grid.frame[[l]],varying.list[[l]][,1],xlab="Duration time (t)",ylab="",main=paste("Baseline.",risk.names[l],sep=""),ylim=y.range.Baseline,cex.main=cex.main,tcl=tcl,cex.axis=cex.axis,las=las,type="n", ...)
vector.minus<-get(paste("alpha.Baseline.",l,sep=""))-get(paste("deviation.Baseline.",l,sep=""))
vector.plus<-get(paste("alpha.Baseline.",l,sep=""))+get(paste("deviation.Baseline.",l,sep=""))
polygon(cbind(c(grid.frame[[l]],grid.frame[[l]][length(grid.frame[[l]]):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.frame[[l]],get(paste("alpha.Baseline.",l,sep="")))
lines(grid.frame[[l]],vector.minus,cex=0.08,col=3)
lines(grid.frame[[l]],vector.plus,cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.frame[[l]],get(paste("alpha.Baseline.",l,sep="")),type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)
}



#plot covariates
if (p > 0)
{
for (k in 1:p)
 {
  for (l in 1:L)    
   {
alpha<-as.matrix(varying.list[[l]][,-1])        #varying coefficients of covariates
deviation<-as.matrix(deviation.list[[l]][,-1])  #deviations
for (i in 1:p) {assign(paste("alpha.",factor.names[i],sep=""),alpha[,i]);assign(paste("deviation.",factor.names[i],sep=""),deviation[,i])}
assign(paste("y.range.",l,sep=""),range(c(unlist(mget(paste("alpha.",factor.names,sep=""),envir=as.environment(-1)))-unlist(mget(paste("deviation.",factor.names,sep=""),envir=as.environment(-1))),unlist(mget(paste("alpha.",factor.names,sep=""),envir=as.environment(-1)))+unlist(mget(paste("deviation.",factor.names,sep=""),envir=as.environment(-1))))))
  }
 }
y.range<-range(mget(paste("y.range.",1:L,sep=""),envir=as.environment(-1)))
}

if (p > 0)
for (k in 1:p)
{
for (l in 1:L)    
  {
alpha<-as.matrix(varying.list[[l]][,-1])        #varying coefficients of covariates
deviation<-as.matrix(deviation.list[[l]][,-1])  #deviations
for (i in 1:p) {assign(paste("alpha.",factor.names[i],sep=""),alpha[,i]);assign(paste("deviation.",factor.names[i],sep=""),deviation[,i])}
plot(grid.frame[[l]],get(paste("alpha.",factor.names[k],sep="")),xlab="Duration time (t)",ylab="",main=paste(factor.names[k],".",risk.names[l],sep=""),ylim=y.range,cex.main=cex.main,tcl=tcl,cex.axis=cex.axis,las=las,type="n", ...)
vector.minus<-get(paste("alpha.",factor.names[k],sep=""))-get(paste("deviation.",factor.names[k],sep=""))
vector.plus<-get(paste("alpha.",factor.names[k],sep=""))+get(paste("deviation.",factor.names[k],sep=""))
polygon(cbind(c(grid.frame[[l]],grid.frame[[l]][length(grid.frame[[l]]):1]),c(vector.minus,vector.plus[length(vector.plus):1])),col="grey")
lines(grid.frame[[l]],get(paste("alpha.",factor.names[k],sep="")))
lines(grid.frame[[l]],get(paste("alpha.",factor.names[k],sep=""))-get(paste("deviation.",factor.names[k],sep="")),cex=0.08,col=3)
lines(grid.frame[[l]],get(paste("alpha.",factor.names[k],sep=""))+get(paste("deviation.",factor.names[k],sep="")),cex=0.08,col=3)
abline(h=0,lty=3,cex=0.05)
par(new=TRUE)
plot(grid.frame[[l]],get(paste("alpha.Baseline.",l,sep="")),type="n",xlab="",ylab="",bty="o",xaxt="n",yaxt="n",...)
  }
}


#dev.off()
#invisible()
}

