"print.CompetingRiskFrailtySurvfit" <-
function(x,...)
{
  
L<-x$L
M.list<-x$M.list
p<-x$p
risk.names<-x$risk.names


#fixed components of varying coefficients
matrix.fixed<-matrix(0,nrow=(p+1)*2,ncol=L)
colnames(matrix.fixed)<-risk.names
if (p > 0) rownames(matrix.fixed)<-c("Baseline.intercept","Baseline.slope",paste(rep(x$factor.names,rep(2,p)),c(".intercept",".slope"),sep="")) else rownames(matrix.fixed)<-c("Baseline.intercept","Baseline.slope")
for (k in 1:(p+1)) for (l in 1:L) {matrix.fixed[2*(k-1)+1,l]<-x$fixed.coef.optim[(k+(l-1)*(p+1))];matrix.fixed[2*k,l]<-x$fixed.coef.optim[(k+L*(p+1)+(l-1)*(p+1))]}

cat("\n\n optimal Fixed parameters of varying coefficients: \n\n")
print(matrix.fixed)
cat("\n\n")

#penalties of random components
matrix.penalty<-matrix(0,nrow=p+1,ncol=L)
colnames(matrix.penalty)<-risk.names
if (p > 0) rownames(matrix.penalty)<-c("penalty.Baseline",paste("penalty",x$factor.names,sep=".")) else rownames(matrix.penalty)<-c("penalty.Baseline")
for (k in 1:(p+1)) for (l in 1:L) matrix.penalty[k,l]<-x$penalty.varying.optim[k+(l-1)*(p+1)]

cat("\n\n optimal Penalty parameters of varying coefficients: \n\n")
print(matrix.penalty)
cat("\n\n")


#mxixture weights of mixture of gamma densites
index.matrix<-expand.grid(eval(parse(text=paste("list(",paste("Column.",1:L,"=1:M.list[[",1:L,"]]",sep="",collapse=","),")"))))
multi.index<-sapply(1:dim(index.matrix)[1],FUN=function(i) paste(index.matrix[i,],sep="",collapse=""))
mixtures<-cbind(as.numeric(multi.index),x$mixture.weights)
dimnames(mixtures)<-list(1:prod(unlist(M.list)),c("index","mixture.weights"))

cat("\n\n optimal mixture weights: \n\n")
print(mixtures)
cat("\n\n")

#optimal values of aic, df of mixtures and marginal log.lik
aic.df.log.lik<-c(x$aic.optim,x$df.weights.optim,x$log.lik.margin.optim)
names(aic.df.log.lik)<-c("aic","df.weights","log.lik.margin")

cat("\n\n AIC, degrees of freedom, marginal log likelihood : \n\n")
print(aic.df.log.lik)
cat("\n\n")


}  #end print.CompetingRiskFrailtySurvfit()

