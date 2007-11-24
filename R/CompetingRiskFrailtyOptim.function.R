CompetingRiskFrailtyOptim<-function(data.set=data.set,control=control,form=formula,risk.names=risk.names)  
{

  
######################################################################
######################################################################
#Description of  functions, which to be called later##################
######################################################################
######################################################################  
  
####################################################################################
#generalized inverse of a matrix (more stable with LINPACK, as with LAPACK)#########
####################################################################################
ginverse<-function(X,tol=1e-100)
{
Xsvd<-svd(X,LINPACK=TRUE)
if (is.complex(X)) Xsvd$u<-Conj(Xsvd$u)
Positive<-Xsvd$d > max(tol * Xsvd$d[1], 0)
if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
else if (!any(Positive)) array(0, dim(X)[2:1])
else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}
#################################################################
#value of the Gamma function,  Striling-formula##################
#################################################################
gamma.stirling.ratio<-function(x,y)
{
exp(y-x)*exp((x-0.5)*log(x)-(y-0.5)*log(y))*(1+1/(12*x)+1/(288*x^2)-139/(51840*x^3)-571/(2488320*x^4))/(1+1/(12*y)+1/(288*y^2)-139/(51840*y^3)-571/(2488320*y^4))
}  
#######################################################################################################################
#specify knots for truncated polynomial slpine basis (s. Ngo/Wand: Smoothing with mixed model software, 2003)##########
#######################################################################################################################
default.knots<-function(x,num.knots)
{
if (missing(num.knots)) num.knots <- max(5,min(floor(length(unique(x))/4),35))
return(quantile(unique(x),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))]))
} 
#################################################
#create artificial poisson data##################
#################################################
survival.to.poisson.frailty<-function(time=time, status.risk=status.risk, x=NULL)
{  
event.time<-time[status.risk==1]
#specify number of integration points (not to be too large beacuse of computation)
if (length(unique(event.time)) < 30) grid<-c(unique(sort(event.time)),max(time+1)) else grid<-c(sort(sample(unique(time[status.risk==1]),30)),max(time+1))
m<-length(grid)
grid.minus1<-c(0,grid[1:(m-1)])
grid.plus1<-c(grid[2:m],max(time)+2)
#artificial {0,1}-poisson.data
yt.list<-lapply(seq(time),FUN=function(ij) c(rep(0,sum(grid<time[ij])),status.risk[ij]))
Mt.list<-lapply(seq(time),FUN=function(ij) sum(grid<time[ij])+1)                
Mt<-unlist(Mt.list)
#inflated design.matrix
Xt<-lapply(seq(time),FUN=function(ij) kronecker(matrix(1,Mt[ij],1),t(matrix(as.matrix(x[ij,])))))
Xt<-eval(parse(text=paste("rbind(",paste("Xt[[",1:length(Xt),"]]",sep="",collapse=","),")",sep="")))
#offset
Ot.list<-lapply(seq(time),FUN=function(ij) log(0.5*(apply(cbind(grid.plus1[1:Mt[ij]],rep(time[ij],Mt[ij])),MARGIN=1,FUN=min) - apply(cbind(grid.minus1[1:Mt[ij]],rep(time[ij],Mt[ij])),MARGIN=1,FUN=min))))
#beachte: diese berechnung von offset weicht geringfgig am rechten rand von der angegebenen formel ab                
return(list(grid=grid,y.list=yt.list,m.list=Mt.list,o.list=Ot.list,x=Xt))   
}


############################################################
#Convex.program(via primal-dual optimization)###############
#optimize mixture weights c.m w.r.t. to constraints######### 
############################################################
primal.dual<-function(penalty.cc,cc,lambda,nu)    
{
###############################################
#ziel.funktion definieren, (= -log.lik)########
###############################################
funktion.wert<-function(cc) {-sum(apply(matrix(cc,nrow=nrow(Delta.i.m),ncol=length(cc),byrow=TRUE)*Delta.i.m,1,FUN=function(zeile) log(sum(zeile)))) + 0.5*penalty.cc*t(cc)%*%Nabla.penalty%*%(cc)}
  
################################
#residuen.funktion definieren###
################################
residuen<-function(cc,lambda,nu)  #residuen.vektor bestimmen
  {
#Delta.i.m ist matrix mit dim=c(N,M)
Delta.i.m.quer<-(1/apply(matrix(cc,nrow=nrow(Delta.i.m),ncol=length(cc),byrow=TRUE)*Delta.i.m,1,FUN=function(x) sum(x)))*Delta.i.m
Nabla.f.0<-apply(Delta.i.m.quer,2,FUN=function(x) -sum(x)) + penalty.cc*Nabla.penalty%*%cc
Nabla.2.f.0<-matrix(0,M,M);for (i in seq(ID.unique)) Nabla.2.f.0<-Nabla.2.f.0 + outer(Delta.i.m.quer[i,],Delta.i.m.quer[i,])
Nabla.2.f.0<-Nabla.2.f.0 + penalty.cc*Nabla.penalty         
#ableitungen der i.constraints in der form ...<0
Nabla.f.i<-diag(-1,length(cc)) #D.f=Nabla.f.i  
Nabla.2.f.i.summe<-matrix(0,ncol=ncol(Nabla.2.f.0),nrow=nrow(Nabla.2.f.0))

r.dual<-Nabla.f.0+t(Nabla.f.i)%*%lambda+t(A)%*%nu
r.cent<- -diag(lambda)%*%f.i-rep(1/tt,length(lambda))
r.pri<-A%*%cc-b.constr
c(r.dual,r.cent,r.pri)
  }
#########################
#start.werte#############
#########################
mu<-1000
epsilon.feas<-1e-07;epsilon<-1e-06  #abbruch.kriterien
m<-length(cc)
eta.dach<- t(cc)%*%lambda  + 1e-08
##############################################
#parameter fr line.search-back.tracking######
##############################################
alpha.backtracking<-0.01
beta.backtracking<-0.5
###########################
#constraints###############
###########################
#i.constraints in der form ...<0
f.i<- -cc  #i.contsriants: cc >=0
#e.constraints
A<-rbind(rep(1,M))
b.constr<-c(1)
repeat
  {
f.i<- -cc
tt<-2*mu/eta.dach
#cat("t= ",tt,"\n\n")

#####################################
#Ableitungen bestimmen###############
#####################################
#Delta.i.m als matrix darstellen mit dim=c(N,M)
Delta.i.m.quer<-(1/apply(matrix(cc,nrow=nrow(Delta.i.m),ncol=length(cc),byrow=TRUE)*Delta.i.m,1,FUN=function(x) sum(x)))*Delta.i.m
Nabla.f.0<-apply(Delta.i.m.quer,2,FUN=function(x) -sum(x)) + penalty.cc*Nabla.penalty%*%cc
Nabla.2.f.0<-matrix(0,M,M);for (i in seq(ID.unique)) Nabla.2.f.0<- Nabla.2.f.0+outer(Delta.i.m.quer[i,],Delta.i.m.quer[i,])
Nabla.2.f.0<-Nabla.2.f.0+penalty.cc*Nabla.penalty
Nabla.f.i<-diag(-1,length(cc)) #D.f=Nabla.f.i
Nabla.2.f.i.summe<-matrix(0,ncol=ncol(Nabla.2.f.0),nrow=nrow(Nabla.2.f.0))

############################
#residuen.vektoren##########
############################
r.dual<-Nabla.f.0+t(Nabla.f.i)%*%lambda+t(A)%*%nu
r.cent<- -diag(lambda)%*%f.i-rep(1/tt,length(lambda))
r.pri<-A%*%cc-b.constr
#cat("r.dual= ",sqrt(sum(r.dual^2)),"\n")
#cat("r.cent= ",sqrt(sum(r.cent^2)),"\n")
#cat("r.pri= ",sqrt(sum(r.pri^2)),"\n")
############################  
#wert der ziel.funktion f.0#
############################
f.0<-funktion.wert(cc)
#############################################
#l�ung durch reduzierung####################
#############################################
#berechnung fr delta.cc und delta.nu (s. Harville,p.468)
H<-Nabla.2.f.0+Nabla.2.f.i.summe+diag(lambda/cc)
#H<-Nabla.2.f.0+Nabla.2.f.i.summe-t(Nabla.f.i)%*%qr.solve(diag(f.i))%*%diag(lambda)%*%Nabla.f.i
b<- -(r.dual+t(Nabla.f.i)%*%ginverse(diag(f.i),tol=1e-100)%*%r.cent)
d<- -r.pri
#k<-1000
#W<-ginverse(H+k*t(A)%*%A,tol=1e-100)  #hier U=Id, aber zur stabilit� nehme z.B U=k*Id mit k=1000
W<-qr.solve(H+t(A)%*%A,tol=1e-100)
TT<-ginverse(A%*%W%*%t(A),tol=1e-100)%*%(A%*%W%*%b-d)
#delta.nu<-as.vector(TT+k*d)
delta.nu<-as.vector(TT+d)
delta.cc<-as.vector(W%*%b-W%*%t(A)%*%TT)
delta.lambda<-as.vector(-ginverse(diag(f.i),tol=1e-100)%*%diag(lambda)%*%Nabla.f.i%*%delta.cc+ginverse(diag(f.i),tol=1e-100)%*%r.cent)
#cat("delta.lambda= ",delta.lambda,"\n")
#cat("delta.nu= ",delta.nu,"\n")
#cat("delta.cc= ",delta.cc,"\n\n")
if  (sqrt(sum(r.pri^2)) < epsilon.feas & sqrt(sum(r.dual^2)) < epsilon.feas & sqrt(sum(r.dual^2)) < epsilon.feas & eta.dach < epsilon) break
##################
#line.search######
##################
s<-1
while (any(cc+s*delta.cc < 0)) s<-s*beta.backtracking
while (sqrt(sum(residuen(cc+s*delta.cc,lambda+s*delta.lambda,nu+s*delta.nu)^2)) > (1-alpha.backtracking*s)*sqrt(sum(residuen(cc,lambda,nu)^2))) s<-s*beta.backtracking
#cat("s= ",s, "\n\n")
#######################
#updates###############
#######################
lambda<-abs(lambda+s*delta.lambda)  #technische hilfe, da theoretisch lambda >0 gilt
nu<-nu+s*delta.nu
cc<-cc+s*delta.cc
#cat("lambda= ",lambda,"\n")
#cat("nu= ",nu,"\n")
#cat("cc= ",cc,"\n")   
eta.dach<- t(cc)%*%lambda  + 1e-12
#cat("eta.dach= ",eta.dach,"\n\n")
 }  #ende repeat
list(cc=as.vector(cc),lambda=as.vector(lambda),nu=as.vector(nu))

} #end of function primal.dual()



#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


names(data.set)[1]<-"ID"



N<-length(data.set[,1])                 #number of observations
N.Cluster<-length(unique(data.set$ID))  #number of clusters
ni<-table(data.set$ID)                  #number of observations in each cluster
ID.unique<-unique(data.set$ID)
names(ID.unique)<-1:length(ID.unique)



#################################################################
#define for each risk type the spline knots for spline bases#####
#################################################################
L<-length(grep("status",names(data.set),fix=TRUE))
Status.names<-names(data.set[,grep("status",names(data.set),fix=TRUE),drop=FALSE])
knots.t<-list()
for (l in 1:L)
{if (length(grep("num.knots",names(control)))==0)  knots.t[[l]]<-default.knots(data.set$time[eval(parse(text=paste("data.set$",Status.names[l],sep="")))==1])
else knots.t[[l]]<-default.knots(data.set$time[eval(parse(text=paste("data.set$",Status.names[l],sep="")))==1],control$num.knots[l])
}

K.t<-lapply(knots.t,FUN=function(x) length(x))



############################################################
#indicator matrix with elementes d.ijl as values############
############################################################
D.ijl<-data.set[,grep("status",names(data.set),fix=TRUE),drop=FALSE] 




##############################################################################################################
#create a model.matrix for (co)variables######################################################################
#N.B.!!! the  reference categories for factors should be defined in the analizied data set bevore applied#####
##############################################################################################################
if (length(attributes(terms(form))$term.labels) == 0)
model.matrix.x<-model.matrix(~numeric(nrow(data.set)))[,1,drop=FALSE] else
model.matrix.x<-model.matrix(formula(paste("~",paste(attributes(terms(form))$term.labels,collapse="+"))),data=data.set)



                
#########################################################################
#artificial poisson data(for each competing risk)########################
#########################################################################                 
grid.points.list<-list()   #t-points, at which the integral is being approximated
y.poisson.list<-list()     #poisson.data
m.list<-list()             #number of poisson data
Design.variables<-list()   #inflated covariates matrix "model.matrix.x"
offset.list<-list()        #offset.parameter

for (l in 1:L)
{
help.surv<-survival.to.poisson.frailty(time=data.set$time,status.risk=eval(parse(text=paste("data.set$",Status.names[l],sep=""))),x=model.matrix.x)
grid.points.list[[l]]<-help.surv$grid 
y.poisson.list[[l]]<-help.surv$y.list
m.list[[l]]<-help.surv$m.list
offset.list[[l]]<-help.surv$o.list

Design.variables[[l]]<-help.surv$x
colnames(Design.variables[[l]])<-colnames(model.matrix.x);colnames(Design.variables[[l]])[1]<-"baseline"
}


p<-ncol(Design.variables[[1]])-1     #number of varaibles (including factor levels by factor variables) in the defined design matrix


               
#define survival time for expanded data set                           
Z.time.list<-lapply(1:L,FUN=function(l) sapply(data.set$time,FUN=function(x) c(grid.points.list[[l]][x>grid.points.list[[l]]],x)))
Z.time<-lapply(1:L,FUN=function(l) unlist(Z.time.list[[l]]))                                          


#############################################################################################
#############################################################################################
#define the start values for beta parameter by means of glm regression for y.poisson#########
#############################################################################################
#############################################################################################

if (p > 0)
  {
#for covariates
variables.time<-lapply(1:L,FUN=function(l) Design.variables[[l]][,2:ncol(Design.variables[[l]]),drop=FALSE]*Z.time[[l]])
for (l in 1:L) colnames(variables.time[[l]])<-paste("variables.time.",colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],sep="")

beta.start<-lapply(1:L,FUN=function(l) coef(glm(as.formula(paste("unlist(y.poisson.list[[l]])~",paste(c("Z.time[[l]]",paste(colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],sep="",collapse="+"),paste(colnames(variables.time[[l]]),sep="",collapse="+"),paste("offset(unlist(offset.list[[",l,"]]))",sep="")),collapse="+"))),data=data.frame(Design.variables[[l]],variables.time[[l]]),family=poisson)))
for (l in 1:L) names(beta.start[[l]])<-c(paste(c("beta.t.intercept.Baseline.","beta.t.slope.Baseline."),"Risk.",l,sep=""),paste("beta.t.intercept.",colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],".Risk.",l,sep=""),paste("beta.t.slope.",colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],".Risk.",l,sep=""))
  } else
{
beta.start<-lapply(1:L,FUN=function(l) coef(glm(as.formula(paste("unlist(y.poisson.list[[l]])~",paste(c("Z.time[[l]]","offset(unlist(offset.list[[l]]))"),collapse="+"),sep="")),data=data.frame(Design.variables[[l]]),family=poisson)))  
for (l in 1:L) names(beta.start[[l]])[1:2]<-c(paste("beta.t.intercept.Baseline.Risk.",l,sep=""),paste("beta.t.slope.Baseline.Risk.",l,sep=""))  
}


#coefficients at Baseline                      
beta.t.intercept.Baseline<-lapply(1:L,FUN=function(l) beta.start[[l]][1])    #for Intercept
beta.t.slope.Baseline<-lapply(1:L,FUN=function(l) beta.start[[l]][2])        #for t-Trend
#coefficients at covariables                  
beta.t.intercept.list<-lapply(1:L,FUN=function(l) list())
beta.t.slope.list<-lapply(1:L,FUN=function(l) list())
     

for (l in 1:L)
   {
beta.t.intercept.list[[l]][[1]]<-beta.t.intercept.Baseline[[l]]
beta.t.slope.list[[l]][[1]]<-beta.t.slope.Baseline[[l]]
if (p >0)
  {
for (i in 3:(2+p)) beta.t.intercept.list[[l]][[i-1]]<-beta.start[[l]][i]
for (i in (2+p+1):(2+2*p)) beta.t.slope.list[[l]][[i-(1+p)]]<-beta.start[[l]][i]                     
   }
}
beta.t.intercept<-list()
beta.t.slope<-list()                         
for (l in 1:L)
  {
beta.t.intercept[[l]]<-unlist(beta.t.intercept.list[[l]])
beta.t.slope[[l]]<-unlist(beta.t.slope.list[[l]])
}                     



############################################################################
#initialize penalty parameters through variances of random effects##########
############################################################################
variance.penalty.t<-list()                 
for (l in 1:L)
  {
if (p >0)
  {
variance.penalty.t[[l]]<-rep(10e-1,p+1)
names(variance.penalty.t[[l]])<-c(paste("variance.t.Baseline.Risk.",l,sep=""),paste(paste("variance.t.",colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],sep=""),".Risk.",l,sep=""))
  }  else
{
variance.penalty.t[[l]]<-10e-1   
names(variance.penalty.t[[l]])<-paste("variance.t.Baseline.Risk.",l,sep="")
}
}



###############################################
#create penalty matrices#######################
###############################################
Lambda.t<-lapply(1:L,FUN=function(l) eval(parse(text=paste("diag(c(",paste(paste("0,0,rep(1/variance.penalty.t[[l]][",1:length(variance.penalty.t[[l]]),sep=""),"],K.t[[l]])",sep="",collapse=","),"))",sep=""))))
D.t <-lapply(1:L,FUN=function(l) diag(rep(c(0,0,rep(1,K.t[[l]])),p+1)))



                               
############################################
#initialize random effects##################
############################################
u.t<-list()
for (l in 1:L)
{
u.t.list<-list()
if (p > 0)
 {
u.t.list[[1]]<-rep(0,K.t[[l]]); names(u.t.list[[1]])<-paste("u.t.Baseline.",1:K.t[[l]],".Risk.",l,sep="") 
for (i in 1:p)
  {u.t.list[[i+1]]<-rep(0,K.t[[l]])
   names(u.t.list[[i+1]])<-paste("u.t.",rep(colnames(Design.variables[[l]])[i+1],K.t[[l]]),".",1:K.t[[l]],".Risk.",l,sep="")
  }
 }  else
 { 
u.t.list[[1]]<-rep(0,K.t[[l]])
names(u.t.list[[1]])<-paste("u.t.Baseline.",1:K.t[[l]],".Risk.",l,sep="")
 }
u.t[[l]]<-u.t.list
}




####################################
#combine in theta vector############
####################################
teta.t<-list()
teta<-list()                 
for (l in 1:L)
   {              
teta.t.Baseline<-c(beta.t.intercept.list[[l]][[1]],beta.t.slope.list[[l]][[1]],u.t[[l]][[1]])
teta.t.list<-list()
teta.t.list[[1]]<-teta.t.Baseline                     
if (p >0) for (i in 1:p) teta.t.list[[i+1]]<-c(beta.t.intercept.list[[l]][[i+1]],beta.t.slope.list[[l]][[i+1]],u.t[[l]][[i+1]])
teta.t[[l]]<-unlist(teta.t.list)
teta[[l]]<-c(teta.t[[l]])
}                      





########################################################################
########################################################################
#construct design matrices for the model################################
######################################################################## 
########################################################################
Design.matrix.t.list<-list()
Basis.t.list<-list()
variables.t<-list()
for (l in 1:L)
{
Basis.t.list[[l]]<-outer(Z.time[[l]],knots.t[[l]],FUN="-")
Basis.t.list[[l]]<-Basis.t.list[[l]]*(Basis.t.list[[l]]>0)
if ( p > 0)
 {
variables.t[[l]]<-Design.variables[[l]][,2:ncol(Design.variables[[l]]),drop=FALSE]
Design.matrix.t.list[[l]]<-cbind(1,Z.time[[l]],Basis.t.list[[l]],eval(parse(text=paste("cbind(",paste("variables.t[[l]][,",1:p,"]*cbind(1,Z.time[[l]],Basis.t.list[[l]])",sep="",collapse=","),")"))))
 } else
Design.matrix.t.list[[l]]<-cbind(1,Z.time[[l]],Basis.t.list[[l]])
}





########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#################   EM-Algorithm  ######################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################




########################################################################
#Initialize parameters of mixture distributions for frailties###########
########################################################################
M.list<-list()                #number M.1,...,M.L of grid functions on the grid M.1.x...x.M.L
mu.list<-list()               #means of Gamma mixture densities =1
sigma.quadrat.l.ml<-list()    #variance of Gamma mixture densities set to
alpha.l.ml<-list()            #shape parameter of gamma mixture densities
beta.l.ml<-list()             #scale parameter of gamma mixture densities

for (l in 1:L)
{
M.list[[l]]<-7
mu.list[[l]]<-rep(1,M.list[[l]])
sigma.quadrat.l.ml[[l]]<-c(0.001,0.05,0.15,0.4,0.6,0.8,0.99)^2
alpha.l.ml[[l]]<-1/sigma.quadrat.l.ml[[l]]
beta.l.ml[[l]]<-alpha.l.ml[[l]]
}

M<-prod(unlist(M.list))


#set the vector of multiindices (expand to array) to order the mixture coefficients c.m
index.matrix<-expand.grid(eval(parse(text=paste("list(",paste("Column.",1:L,"=1:M.list[[",1:L,"]]",sep="",collapse=","),")"))))
#expnad shape and scale parameters accordingly
alpha.l.ml.matrix<-as.matrix(expand.grid(eval(parse(text=paste("list(",paste("Column.",1:L,"=alpha.l.ml[[",1:L,"]]",sep="",collapse=","),")")))))
beta.l.ml.matrix<-as.matrix(expand.grid(eval(parse(text=paste("list(",paste("Column.",1:L,"=beta.l.ml[[",1:L,"]]",sep="",collapse=","),")")))))


##############################################################################################################
#define Penalty Matrix on mixture weights c.m (marginal weights are being penalized)##########################
##############################################################################################################
sigma.quadrat.l.ml.matrix<-as.matrix(expand.grid(eval(parse(text=paste("list(",paste("Column.",1:L,"=sigma.quadrat.l.ml[[",1:L,"]]",sep="",collapse=","),")")))))
sigma.penalty<-apply(sigma.quadrat.l.ml.matrix,1,sum)
Nabla.penalty<-diag(sigma.penalty^2)



################################################################################
################################################################################
#define penalty sequence for penalty parameter on mixture weights c.m###########
################################################################################
################################################################################
penalty.sequence<-c(0.001,0.1,1,10,100,1000,5000,10000,100000)
#penalty.sequence<-c(0.001,0.1)



#####################################################################################################################
#initialize list of frailty values for stopping of the EM-algorithm for each given value of penalty.sequence#########
#####################################################################################################################
w.list<-list()

#########################################################################################
#initialize vector of AIC, df.weigths and log.lik.margin values,#########################
#to define an optimum from them w.r.t. to the  penalty sequence##########################
#########################################################################################
aic.vector<-c()
degree.freedom.weights.vector<-c()
log.lik.margin.vector<-c()
aic.and.log.lik.and.df.EM<-list()
##################################################################################################################
#initialize lists for the estimating parameters of the model,#####################################################
#each element of a list corresponds to a penalty value from the penalty sequence##################################
##################################################################################################################
cc.list<-list()                  #mixture weigths
for (l in 1:L) assign(paste("alpha.t.risk.",l,".list",sep=""),vector(mode="list",length=length(penalty.sequence)))  #varying coefficients
deviation.t.list<-list()         #deviation
frailty.estimate.list<-list()    #frailties
##################################################################################################################
#initialize lists for resulting estimates of fixed and random components of varying coefficients of the model#####
#as well as penalty values for their random parts,################################################################
#each element of a list corresponds to a penalty value from the penalty sequence##################################
##################################################################################################################
estimates.fix<-list()
estimates.random<-list()
penalty.t.list<-list()
##################################################################################################################
#initialize lists for resulting estimates of degrees.of.freedom for varying coefficients##########################
##################################################################################################################
#df.varying<-list()



##############################################################
##############################################################
#loop for penalty parameter from  penalty sequence############
##############################################################
##############################################################
for (penalty.cc in penalty.sequence)
  {

if (control$print.penalty.mixture) cat("\n\n","penalty mixture =",penalty.cc,"\n\n")
    
if (control$print.estimates)
  {
   cat("\n\n","start values: fixed parameters of the varying coefficients","\n\n")
   print(unlist(beta.t.intercept.list))
   print(unlist(beta.t.slope.list))
   cat("","\n\n")
   #print(penalty)
  }
   
#################################
#initialize following values#####
#################################
help.aic.and.log.lik.and.df.EM<-NULL #writes out (in matrix form) the values of aic, marginal log.lik and df's from each EM iteration
w.list<-list()  #writes out frailties from each EM iteration for comparision for stopping criterion
w.list[[1]]<-lapply(1:L,FUN=function(l) lapply(seq(unique(data.set$ID)),FUN=function(i) 1 ))
log.lik.margin<-c()   #marginal log likelihood of the model (given penalty values from the penalty sequence)


variance.em<-list()  #writes out the variances(=1/penalty) of random components of varying coefficients (from each EM iteration)
variance.em[[1]]<-unlist(variance.penalty.t)   




##########################################################
#start values for the convex optimization routine######### (can be defined out of the penalty.cc loop)
##########################################################
cc<-rep(1,M)               #has to be "strong" feasible, i.e. fulfills i.constraints x >0 (but not necessaraly e.constraints Ax=b)
lambda<-rep(1,length(cc))  #fulfills lambda >0  (= number of i.constraints)
nu<-1                      #can be of any initial value (= number of e.constraints)



###############################################################
###############################################################
######################EM.loop##################################
###############################################################
###############################################################
for (em in 1:control$niter.EM)  
  {

if (control$print.EM)  cat("\n\n","EM.iteration =",em,"\n\n") 


#####################################
############E-step###################
#####################################


##########################################################################
#calculate ingredients of the marginal density of poisson data############
##########################################################################

#lambda.ijl (updated from each M-step of EM-algorithm)
lambda.ijl.list<-list()
lambda.ijl.K.ij.list<-list()
lambda.sum<-matrix(0,nrow=length(seq(ID.unique)),ncol=L)
for (l in 1:L)
{
#build vector lambda.ijlk, k=1...K.ij for all ij; for each given risik l
lambda.ijl<-exp(Design.matrix.t.list[[l]]%*%teta[[l]]+unlist(offset.list[[l]]))
#for each cluster i define aggregated number of poisson data (through all spells j=1...ni)
help.lapply<-lapply(ID.unique,FUN=function(i) {help.number<-cumsum(unlist(m.list[[l]][data.set$ID==i]));help.number[length(help.number)]})
#indiziere geeignet fr jedes cluster
help.cut<-cumsum(unlist(help.lapply))   
help.matrix<-cbind(c(0,help.cut[-length(help.cut)])+1,help.cut)  #bildet die indizes.matrix fr das jeweilige individuum
indizes<-lapply(seq(ID.unique),FUN=function(i) help.matrix[i,1]:help.matrix[i,2]) #jeweilige index.vektoren
#zum i-ten Cluster zugeh�iger vektor der lambda.ijl
lambda.ijl.list[[l]]<-lapply(indizes,FUN=function(ind) lambda.ijl[ind])
#summierte lambda ber j und K fr jedes t-te Cluster
lambda.sum[,l]<-sapply(lambda.ijl.list[[l]],FUN=function(x) sum(x))
#die k=K.ij bei jeweiligen spells des i-ten individuums
lambda.ijl.K.ij<-lambda.ijl[cumsum(unlist(m.list[[l]]))]
#die zum i-ten Cluster geh�en,  zusamen.gefasst
lambda.ijl.K.ij.list[[l]]<-lapply(ID.unique,FUN=function(i) lambda.ijl.K.ij[data.set$ID==i])
}
lambda.prod<-unlist(lapply(seq(ID.unique),FUN=function(i) {product.l<-c(); for (l in 1:L) product.l[l]<-prod(lambda.ijl.K.ij.list[[l]][[i]]^D.ijl[data.set$ID==as.numeric(ID.unique[names(ID.unique)==i]),l]);prod(product.l)}))


#updated shape and scale parameters of gamma distributions
alpha.l.ml.given.delta.i<-lapply(ID.unique,FUN=function(i) {help.matrix<-c();for (l in 1:L) {help.matrix<-cbind(help.matrix,sum(D.ijl[data.set$ID==i,l])+alpha.l.ml.matrix[,l])};help.matrix})
beta.l.ml.given.delta.i<-lapply(seq(ID.unique),FUN=function(i) {help.matrix<-c();for (l in 1:L) {help.matrix<-cbind(help.matrix,lambda.sum[i,l]+beta.l.ml.matrix[,l])};help.matrix})
#expand to the 3.dim array
alpha.l.ml.expand.to.array<-array(alpha.l.ml.matrix,c(M,L,length(seq(ID.unique))))
alpha.l.ml.given.delta.i.expand.to.array<-array(unlist(lapply(seq(ID.unique),FUN=function(i) alpha.l.ml.given.delta.i[[i]])),c(M,L,length(seq(ID.unique))))
beta.l.ml.expand.to.array<-array(beta.l.ml.matrix,c(M,L,length(seq(ID.unique))))
beta.l.ml.given.delta.i.expand.to.array<-array(unlist(lapply(seq(ID.unique),FUN=function(i) beta.l.ml.given.delta.i[[i]])),c(M,L,length(seq(ID.unique))))


#calculate ratio of gamma functions durch stirling.approximation
#(some neglegible differences between stirling approximation and direct calculation)
#through stirling approximation
#gamma.ratio<-ifelse(alpha.l.ml.given.delta.i.expand.to.array > 10 & alpha.l.ml.expand.to.array > 10,gamma.stirling.ratio(alpha.l.ml.given.delta.i.expand.to.array,alpha.l.ml.expand.to.array),gamma(alpha.l.ml.given.delta.i.expand.to.array)/gamma(alpha.l.ml.expand.to.array))
#direct through Gamma(z+1)=z*Gamma(z) 
d.array<-alpha.l.ml.given.delta.i.expand.to.array-alpha.l.ml.expand.to.array
gamma.ratio<-array(apply(cbind(seq(alpha.l.ml.expand.to.array)),1,FUN=function(i) prod(alpha.l.ml.expand.to.array[i]:(alpha.l.ml.expand.to.array[i]+d.array[i]-1))),dim=dim(alpha.l.ml.expand.to.array))
gamma.ratio<-ifelse(d.array==0,1,gamma.ratio)
oben<-alpha.l.ml.expand.to.array*log(beta.l.ml.expand.to.array)
unten<-alpha.l.ml.given.delta.i.expand.to.array*log(beta.l.ml.given.delta.i.expand.to.array)
exp.ratio<-exp(oben-unten)
exp.gamma<-exp.ratio*gamma.ratio
exp.gamma.prod<-t(apply(exp.gamma,c(1,3),"prod"))


##########################################################
#combine to the matrix Delta (i-->rows,m-->columns)#######
##########################################################
Delta.i.m<-exp.gamma.prod*matrix(lambda.prod,ncol=ncol(exp.gamma.prod),nrow=nrow(exp.gamma.prod),byrow=FALSE)


############################
#optimal weights c.m########
############################
#(in the next EM iteration the first value will be set nearby that one already calculated from the preceeding iteration, to reduce computation)
cc<-primal.dual(penalty.cc,as.vector(cc+1e-05),lambda,nu)$cc


######################################################
#marginal weights c.(l)ml#############################
######################################################
Delta.i.m.quer<-(1/apply(matrix(cc,nrow=nrow(Delta.i.m),ncol=length(cc),byrow=TRUE)*Delta.i.m,1,FUN=function(x) sum(x)))*Delta.i.m
cc.tilde.i<-matrix(cc,nrow=nrow(Delta.i.m),ncol=length(cc),byrow=TRUE)*Delta.i.m.quer             
c.l.ml.i<-lapply(seq(ID.unique),FUN=function(i) lapply(1:L,FUN=function(l) lapply(1:M.list[[l]],FUN=function(m.l) sum(cc.tilde.i[i,index.matrix[,l]==m.l]))))


##########################################################
#define the degrees of freedom from proposed relation#####
##########################################################
Nabla.2.f.0<-matrix(0,M,M);for (i in seq(ID.unique)) Nabla.2.f.0<- Nabla.2.f.0+outer(Delta.i.m.quer[i,],Delta.i.m.quer[i,])
I<- Nabla.2.f.0  #not penalized fisher information (here implicit multiplied with "-1", because Nabla.2.f.0 calculated for -log.lik)
H<- Nabla.2.f.0+penalty.cc*Nabla.penalty  #penalized fisher information 
degree.freedom.weights<-sum(diag(qr.solve(H,tol=1e-100)%*%I))
#marginal log.lik will be calculated later; hence aic=-log.lik.margin+degree.freedom.weights


############################
#frailties estimates########
############################
w<-lapply(1:L,FUN=function(l) lapply(seq(unique(data.set$ID)),FUN=function(i) sum(unlist(c.l.ml.i[[i]][[l]])*unique(alpha.l.ml.given.delta.i[[i]][,l])/unique(beta.l.ml.given.delta.i[[i]][,l]))))
w.list[[em+1]]<-w
difference.frailty<-c(); for (l in 1:L) difference.frailty[l]<-sum((unlist(w.list[[em+1]][[l]])-unlist(w.list[[em]][[l]]))^2)
#cat("sum(w.alt-w.neu)^2= ",sum(abs.differenz.frailty),"\n\n")



#########################################################
#####################End of E-step of EM#################
#########################################################



#update the offset parameter of the model
#(depending on the value of frailty parameter w)
#this builds the link between frailty and theta
help.lapply<-list()  #for each cluster i define aggregated number of poisson data(through all spells j=1...ni)
for (l in 1:L) help.lapply[[l]]<-unlist(lapply(ID.unique,FUN=function(i) {help.number<-cumsum(unlist(m.list[[l]][data.set$ID==i]));help.number[length(help.number)]}))
#for each risk l a vector of frailty values
offset.frailty<-lapply(1:L,FUN=function(l) unlist(offset.list[[l]]) + rep(log(unlist(w[[l]])),help.lapply[[l]]))





###########################################
############# M-step of EM ################
###########################################


##################################################################################################################
#approximate penalized fisher.type matrix (for frailty=1), to calculate approx. marginal log likilehood###########
##################################################################################################################
I.t.frailty.1<-list()

index.vector<-list()
for (l in 1:L)
{index.vector[[l]]<-(2+1):(2+K.t[[l]])
if (p > 0) for (k in 1:p) index.vector[[l]]<-c(index.vector[[l]],(2+K.t[[l]]+2*k+1+(k-1)*K.t[[l]]):(2+K.t[[l]]+2*k+k*K.t[[l]]))
}



variance.epoch<-list()   #writes out the variances(=1/penalty) of random components of varying coefficients (from each epoch iteration)

#############################################################
#fr jeden Ausfall.typ l=1...L###############################
#############################################################
for (l in 1:L)
{
variance.epoch.l<-list()   #initialize list of variances(=1/lambda); each element of this list is a value from current epoch
variance.epoch.l[[1]]<-unlist(variance.penalty.t[[l]])


for (epoch in 1:control$niter.epoch)
 {

#penalized score for theta.t
S.t<-crossprod(Design.matrix.t.list[[l]],as.vector(unlist(y.poisson.list[[l]])-exp(Design.matrix.t.list[[l]]%*%teta.t[[l]]+unlist(offset.frailty[[l]]))))
S.t.pen<-S.t-Lambda.t[[l]]%*%D.t[[l]]%*%teta.t[[l]]
#penalized observed fisher.type matrix(multipliced  with "-1")
I.t<-crossprod(Design.matrix.t.list[[l]],Design.matrix.t.list[[l]]*as.vector(exp(Design.matrix.t.list[[l]]%*%teta.t[[l]]+unlist(offset.frailty[[l]]))))
I.t.pen<-I.t+Lambda.t[[l]]%*%D.t[[l]]


##########################################################################################
#approximation of penalised fisher.type matrix(with frailties=1)##########################
##########################################################################################
I.t.frailty.1[[l]]<-crossprod(Design.matrix.t.list[[l]],Design.matrix.t.list[[l]]*as.vector(exp(Design.matrix.t.list[[l]]%*%teta.t[[l]]+unlist(offset.list[[l]]))))
I.t.frailty.1[[l]]<-I.t.frailty.1[[l]][index.vector[[l]],index.vector[[l]]] + Lambda.t[[l]][index.vector[[l]],index.vector[[l]]]


#update theta.t
inverse.t<-qr.solve(I.t.pen,tol=1e-100)
teta.t.old<-teta.t[[l]]
teta.t.new<-teta.t.old+inverse.t%*%S.t.pen
names(teta.t.new)<-names(teta.t[[l]])    
teta.t[[l]]<-teta.t.new
names(teta.t[[l]])<-names(teta.t.new)    
teta[[l]]<-c(teta.t[[l]])
#update fixed components beta
beta.t.intercept[[l]]<-teta.t[[l]][grep("beta.t.intercept",names(teta.t[[l]]),fixed=TRUE)]
beta.t.slope[[l]]<-teta.t[[l]][grep("beta.t.slope",names(teta.t[[l]]),fixed=TRUE)]
#update random components u
u.t[[l]]<-teta.t[[l]][grep("u.t",names(teta.t[[l]]),fixed=TRUE)]




#########################################################################
#fixpoint iteration for update of lambda penalties#######################
#just one update w.r.t lambda=F(lambda)##################################
#########################################################################

#old values of variances(=1/lambda) 
variance.penalty.t.old<-variance.penalty.t[[l]]
variance.penalty.old<-c(variance.penalty.t.old)

#define first the submatrices of the (penalized) fisher.information, correspondig to the random coefficients
index.vector.t<-(2+1):(2+K.t[[l]])  #for u.t.Baseline
if (p > 0) for (k in 1:p) index.vector.t<-c(index.vector.t,(2+K.t[[l]]+2*k+1+(k-1)*K.t[[l]]):(2+K.t[[l]]+2*k+k*K.t[[l]]))
I.t.pen.u.t<-I.t.pen[index.vector.t,index.vector.t]    #penalizied fisher.type matrix for t
inverse.I.t.pen.u.t<-qr.solve(I.t.pen.u.t,tol=1e-100)

#update of penalties
variance.penalty.t.new<-c()
#baseline
variance.penalty.t.new[1]<-(sum(diag(inverse.I.t.pen.u.t[1:K.t[[l]],1:K.t[[l]]]))+sum(u.t[[l]][grep("Baseline",names(u.t[[l]]),fixed=TRUE)]^2))/K.t[[l]]
#covariates
if (p > 0)  for (k in 1:p)
variance.penalty.t.new[k+1]<-(sum(diag(inverse.I.t.pen.u.t[(k*K.t[[l]]+1):((k+1)*K.t[[l]]),(k*K.t[[l]]+1):((k+1)*K.t[[l]])]))+sum(u.t[[l]][(k*K.t[[l]]+1):((k+1)*K.t[[l]])]^2))/K.t[[l]]
Lambda.t[[l]]<-eval(parse(text=paste("diag(c(",paste(paste("0,0,rep(1/variance.penalty.t.new[",1:length(variance.penalty.t.new),sep=""),"],K.t[[l]])",sep="",collapse=","),"))",sep="")))


variance.penalty.t[[l]]<-variance.penalty.t.new
variance.epoch.l[[epoch+1]]<-variance.penalty.t[[l]]




#stoping rule
if (sum((teta.t.new-teta.t.old)^2) < control$tol.epoch & sum((variance.epoch.l[[epoch+1]]-variance.epoch.l[[epoch]])^2) < control$tol.variance) break    


}  #end for (epoch in 1:niter.epoch)



variance.epoch[[l]]<-variance.epoch.l[[epoch+1]]



#degrees.of.freedom for varying coefficients
#assign(paste("df.",l,".Baseline",sep=""),sum(diag((qr.solve(I.t.pen,tol=1e-100)%*%I.t)[1:(2+K.t[[l]]),1:(2+K.t[[l]])])))
#if (p > 0) for (k in 1:p)
#{assign(paste("df.",l,".",colnames(Design.variables[[l]])[k+1],sep=""),sum(diag((qr.solve(I.t.pen,tol=1e-100)%*%I.t)[(2+K.t[[l]]+2*(k-1)+(k-1)*K.t[[l]]+1):(2+K.t[[l]]+2*k+k*K.t[[l]]),(2+K.t[[l]]+2*(k-1)+(k-1)*K.t[[l]]+1):(2+K.t[[l]]+2*k+k*K.t[[l]])])))
#}





} #end for (l in 1:L)




#################################################################################################
#marginal log likelihood (all random components integrated out from the model)###################
#(s. supplemented paper)#########################################################################
#################################################################################################
I.part<-sum(apply(matrix(cc,nrow=nrow(Delta.i.m),ncol=length(cc),byrow=T)*Delta.i.m,1,FUN=function(zeile) log(sum(zeile))))
II.part<-0.5*sum(sapply(1:L,FUN=function(l) as.vector(-t(u.t[[l]])%*%Lambda.t[[l]][index.vector[[l]],index.vector[[l]]]%*%u.t[[l]]) + sum(log(eigen(Lambda.t[[l]][index.vector[[l]],index.vector[[l]]],only.values=TRUE)$values)) - sum(log(eigen(I.t.frailty.1[[l]],only.values=TRUE)$values))))

log.lik.margin<-c(log.lik.margin,I.part+II.part)
if (control$print.log.lik) cat("\n\n","log.lik.margin =",I.part+II.part,"\n\n")

#AIC of the model
aic<- -(I.part+II.part)+degree.freedom.weights

#save the values of AIC, df and marginal log.lik from each EM iteration
help.aic.and.log.lik.and.df.EM<-rbind(help.aic.and.log.lik.and.df.EM,c(aic,degree.freedom.weights,I.part+II.part))
aic.and.log.lik.and.df.EM[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-help.aic.and.log.lik.and.df.EM




#stopping rule for EM iteration
variance.em[[em+1]]<-unlist(variance.epoch)
if (sum((teta.t.new-teta.t.old)^2) < control$tol.epoch & sum((variance.em[[em+1]]-variance.em[[em]])^2) < control$tol.variance & sum(difference.frailty) < control$tol.frailty) break



} #end for (em in 1:control$niter.EM)



###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#######################################End of EM Algorithm#################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################



#save the last value of AIC, df and log.lik.margin from the EM algorithm as an optimal value, given the mixture penalty from penalty sequence
aic.vector<-c(aic.vector,aic.and.log.lik.and.df.EM[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]][em,1])
degree.freedom.weights.vector<-c(degree.freedom.weights.vector,aic.and.log.lik.and.df.EM[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]][em,2])
log.lik.margin.vector<-c(log.lik.margin.vector,aic.and.log.lik.and.df.EM[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]][em,3])


#penalty(=1/lambda)
penalty.t<-list()
for (l in 1:L)
{
penalty.t[[l]]<-1/variance.penalty.t[[l]]
names(penalty.t[[l]])[1]<-paste("penalty.t.Baseline.Risk.",l,sep="")
if (p > 0) names(penalty.t[[l]])[-1]<-paste(paste("penalty.t.",colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],sep=""),".Risk.",l,sep="")
}


#print of resulting estimates: theta and penalty
if (control$print.estimates)
{
cat("\n\n","resulting estimates: beta and penalties","\n\n")
print(unlist(beta.t.intercept))
print(unlist(beta.t.slope))
cat("","\n\n")
print(unlist(penalty.t))
}




###############################################################################
###############################################################################
#calculate variances of the resulting estimates for theta (cp. Louis,1982)#####
###############################################################################
###############################################################################
co.variance.teta<-list()
variance.teta<-list()
offset.frailty<-lapply(1:L,FUN=function(l) unlist(offset.list[[l]]) + rep(log(unlist(w[[l]])),unlist(help.lapply[[l]])))
for (l in 1:L)
{
#complete information, already calculated
fisher.complete<-crossprod(Design.matrix.t.list[[l]],Design.matrix.t.list[[l]]*as.vector(exp(Design.matrix.t.list[[l]]%*%teta.t[[l]]+unlist(offset.frailty[[l]]))))
#(conditional) variance of frailty w, given poisson data delta.l
variance.w<-unlist(lapply(seq(ID.unique),FUN=function(i) sum(unlist(c.l.ml.i[[i]][[l]])*unique(alpha.l.ml.given.delta.i[[i]][,l])/unique(beta.l.ml.given.delta.i[[i]][,l]^2))))
variance.w<-rep(variance.w,unlist(help.lapply[[l]]))  #inflated, according to dimension of the Design.matrix
#missing information
vector.diag<-sqrt(variance.w)*as.vector(exp(Design.matrix.t.list[[l]]%*%teta.t[[l]]+unlist(offset.list[[l]])))
fisher.missing<-crossprod(vector.diag*Design.matrix.t.list[[l]])

fisher.observed<-fisher.complete-fisher.missing  
fisher.observed.pen<-fisher.observed+diag(c(diag(Lambda.t[[l]])))  #penalized observed fisher type matrix

#sandwich estimator of covariance matrix of theta
co.variance.teta[[l]]<-qr.solve(fisher.observed.pen,tol=1e-100)%*%fisher.observed%*%qr.solve(fisher.observed.pen,tol=1e-100)
variance.teta[[l]]<-diag(co.variance.teta[[l]])
names(variance.teta[[l]])<-c(names(teta[[l]]))
}



###################################
###Grid for time###################
###################################
length.grid<-1000
grid.t.list<-list()
B.grid.t.list<-list()
for (l in 1:L)
{
grid.t.list[[l]]<-seq(min(Z.time[[l]]),max(Z.time[[l]]),le=length.grid)
B.grid.t.list[[l]]<-outer(grid.t.list[[l]],knots.t[[l]],FUN="-")
B.grid.t.list[[l]]<-B.grid.t.list[[l]]*(B.grid.t.list[[l]]>0)
}



#############################
#Confidence Bands############
#############################
deviation.t<-list()  #list of L elements, each element is a list of p+1 deviations for: baseline, covariate.1,...,covariate.p
C.t.grid<-list()   
variance.t.Baseline<-list()
deviation.t.Baseline<-list()
for (l in 1:L)
{
deviation.t.matrix<-matrix(0,nrow=length(grid.t.list[[l]]),ncol=p+1)
colnames(deviation.t.matrix)<-paste("deviation.",colnames(Design.variables[[l]]),sep="")
C.t.grid[[l]]<-cbind(1,grid.t.list[[l]],B.grid.t.list[[l]])
#Baseline
co.variance.teta.Baseline<-co.variance.teta[[l]][1:(2+K.t[[l]]),1:(2+K.t[[l]])]
variance.t.Baseline[[l]]<-apply(C.t.grid[[l]],1,FUN=function(help.row) t(help.row)%*%co.variance.teta.Baseline%*%help.row)  
deviation.t.Baseline[[l]]<-qnorm(0.975)*sqrt(variance.t.Baseline[[l]])   #conf.level alpha=0.05
deviation.t.matrix[,1]<-deviation.t.Baseline[[l]]
#covariates
if (p > 0) for (k in 1:p)
  {
assign(paste("co.variance.teta.",colnames(Design.variables[[l]])[k+1],sep=""),co.variance.teta[[l]][(2+K.t[[l]]+2*k-1+(k-1)*K.t[[l]]):(2+K.t[[l]]+2*k+k*K.t[[l]]),(2+K.t[[l]]+2*k-1+(k-1)*K.t[[l]]):(2+K.t[[l]]+2*k+k*K.t[[l]])])
assign(paste("variance.t.",colnames(Design.variables[[l]])[k+1],sep=""),apply(C.t.grid[[l]],1,FUN=function(help.row) t(help.row)%*%get(paste("co.variance.teta.",colnames(Design.variables[[l]])[k+1],sep=""))%*%help.row))
assign(paste("deviation.t.",colnames(Design.variables[[l]])[k+1],sep=""),qnorm(0.975)*sqrt(get(paste("variance.t.",colnames(Design.variables[[l]])[k+1],sep=""))))
deviation.t.matrix[,k+1]<-get(paste("deviation.t.",colnames(Design.variables[[l]])[k+1],sep=""))
  }
deviation.t[[l]]<-deviation.t.matrix
}  #end for(l in 1:L)





######################################################################
#varying coefficients for Baseline and covariates#####################
######################################################################
for (l in 1:L) assign(paste("estimate.beta.t.risk.",l,sep=""),eval(parse(text=paste("c(beta.t.intercept[[",l,"]],beta.t.slope[[",l,"]])",sep=""))))
#varying coefficients with values on the grid
for (l in 1:L) assign(paste("alpha.t.risk.",l,sep=""),vector(mode="list",length=p+1))
#Baseline
for (l in 1:L)
{
help.vector<-as.vector(cbind(1,grid.t.list[[l]])%*%get(paste("estimate.beta.t.risk.",l,sep=""))[grep("Baseline",names(get(paste("estimate.beta.t.risk.",l,sep=""))),fixed=TRUE)]+B.grid.t.list[[l]]%*%u.t[[l]][grep("Baseline",names(u.t[[l]]),fixed=TRUE)]) 
eval(parse(text=paste("alpha.t.risk.",l,"[[1]]","<-","help.vector",sep="")))   
}
#covariates
if (p > 0)
for (l in 1:L)
  {
  for (k in 1:p)
  {
help.vector<-as.vector(cbind(1,grid.t.list[[l]])%*%c(beta.t.intercept[[l]][grep(paste("beta.t.intercept.",colnames(Design.variables[[l]])[k+1],sep=""),names(beta.t.intercept[[l]]),fixed=TRUE)],beta.t.slope[[l]][grep(paste("beta.t.slope.",colnames(Design.variables[[l]])[k+1],sep=""),names(beta.t.slope[[l]]),fixed=TRUE)])+B.grid.t.list[[l]]%*%u.t[[l]][grep(paste("u.t.",colnames(Design.variables[[l]])[k+1],sep=""),names(u.t[[l]]),fixed=TRUE)])
eval(parse(text=paste("alpha.t.risk.",l,"[[",k+1,"]]","<-","help.vector",sep="")))   
  }  
}


########################################
#resulting estimates of frailties#######
########################################
frailty.estimate<-w.list[[em]]

##########################################################################
#resulting estimates of fixed and random components of varying.coef#######
#as well as penalty values for their random parts#########################
##########################################################################
estimates.fix[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-c(unlist(beta.t.intercept),unlist(beta.t.slope))
estimates.random[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-u.t
penalty.t.list[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-unlist(penalty.t)

                 

##############################################################################################
#fill out the list components according to the penalty value of mixture weights###############
##############################################################################################
cc.list[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-cc
for (l in 1:L) eval(parse(text=paste("alpha.t.risk.",l,".list","[[",seq(penalty.sequence)[penalty.sequence==penalty.cc],"]]","<-","alpha.t.risk.",l,sep="")))
deviation.t.list[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-deviation.t
frailty.estimate.list[[seq(penalty.sequence)[penalty.sequence==penalty.cc]]]<-frailty.estimate

                      
}
#####################################################
#####################################################
#end for loop (penalty.cc in penalty.sequence)#######
#####################################################
#####################################################





####################################################################
#define the optimal value of penalty of mixture weights#############
####################################################################
penalty.cc.optim<-penalty.sequence[aic.vector==min(aic.vector)]
position.optim<-seq(penalty.sequence)[penalty.sequence==penalty.cc.optim]
if (control$print.estimates) cat("\n\n\n","optimal penalty of mixture weights  =",penalty.cc.optim,"\n\n\n")





################################################################
#define according optimial values of model components###########
################################################################
#optimal mixture weights
cc.optim<-cc.list[[position.optim]]
#deviation as a list with data frames (w.r.t. to each risk) as elements of it
deviation.optim<-deviation.t.list[[position.optim]]
names(deviation.optim)<-paste("deviation.Risk.",1:L,sep="")
#optimal frailties
frailty.estimate.optim<-frailty.estimate.list[[position.optim]]
frailty.list<-list()
for (l in 1:L) frailty.list[[l]]<-unlist(frailty.estimate.optim[[l]])
names(frailty.list)<-paste("frailty.Risk.",1:L,sep="")
#optimal AIC,df.weights andlog.lik.margin values
aic.optim<-aic.vector[aic.vector==min(aic.vector)]
degree.freedom.weights.optim<-degree.freedom.weights.vector[aic.vector==min(aic.vector)]
log.lik.margin.optim<-log.lik.margin.vector[aic.vector==min(aic.vector)]
#varying coefficients as a list with data frames (w.r.t. to each risk) as elements of it
varying.list<-list()
for (l in 1:L)
{
assign(paste("alpha.t.risk.",l,".optim",sep=""),data.frame(get(paste("alpha.t.risk.",l,".list",sep=""))[[position.optim]]))
if (p > 0) var.coef.names<-c(paste("Baseline.Risk.",l,sep=""),paste(colnames(Design.variables[[l]])[2:ncol(Design.variables[[l]])],".Risk.",l,sep=""))
else var.coef.names<-c(paste("Baseline.Risk.",l,sep=""))
eval(parse(text=paste("names(alpha.t.risk.",l,".optim",")","<-var.coef.names",sep="")))
varying.list[[l]]<-get(paste("alpha.t.risk.",l,".optim",sep=""))
}
#grid for plotting
grid.frame<-data.frame(grid.t.list)
names(grid.frame)<-paste("grid.Risk.",1:L,sep="")

#optimal values of fixed and random components of varying coefficients as well as penalties for their random parts
fixed.coef.optim<-estimates.fix[[position.optim]]
random.coef.optim<-estimates.random[[position.optim]]
penalty.varying.optim<-penalty.t.list[[position.optim]]


#optimal values of degrees.of.freedom for varying coefficients
#df.varying.optim<-df.varying[[position.optim]]


if (p > 0) factor.names<-colnames(Design.variables[[1]])[2:length(colnames(Design.variables[[1]]))] else factor.names<-character(0)




#################################################################
#give the supplemented names to the Risks########################
#################################################################
#fixed.coef (beta)
names(fixed.coef.optim)<-paste(sapply(1:length(fixed.coef.optim),FUN=function(i) strsplit(names(fixed.coef.optim),".Risk.")[[i]][1]),rep(risk.names,rep(p+1,L),sep="."))
#random.coef (u)
for (l in 1:L) names(random.coef.optim[[l]])<-paste(sapply(1:length(random.coef.optim[[l]]),FUN=function(i) strsplit(names(random.coef.optim[[l]]),".Risk.")[[i]][1]),risk.names[l],sep=".")
#penalty.varying.optim (lambda)
for (l in 1:L) names(penalty.varying.optim)[(1+(l-1)*(p+1)):(l*(p+1))]<-paste(sapply(1:(p+1),FUN=function(i) strsplit(names(penalty.varying.optim)[(1+(l-1)*(p+1)):(l*(p+1))],".Risk.")[[i]][1]),risk.names[l],sep=".")
#grid.frame
names(grid.frame)<-paste(sapply(1:length(grid.frame),FUN=function(i) strsplit(names(grid.frame),".Risk.")[[i]][1]),risk.names,sep=".")
#varying.list
for (l in 1:L) colnames(varying.list[[l]])<-paste(sapply(1:length(varying.list[[l]][1,]),FUN=function(i) strsplit(colnames(varying.list[[l]]),".Risk.")[[i]][1]),risk.names[l],sep=".")
#deviation.list
names(deviation.optim)<-paste(sapply(1:length(deviation.optim),FUN=function(i) strsplit(names(deviation.optim),".Risk.")[[i]][1]),risk.names,sep=".")
for (l in 1:L) colnames(deviation.optim[[l]])<-paste(colnames(deviation.optim[[l]]),risk.names[l],sep=".")
#frailty.list
names(frailty.list)<-paste(sapply(1:length(frailty.list),FUN=function(i) strsplit(names(frailty.list),".Risk.")[[i]][1]),risk.names,sep=".")



list(L=L,M.list=M.list,fixed.coef.optim=fixed.coef.optim, random.coef.optim=random.coef.optim, penalty.varying.optim=penalty.varying.optim, penalty.weights.optim=penalty.cc.optim, grid.frame=grid.frame, varying.list=varying.list, deviation.list=deviation.optim, frailty.list=frailty.list, mixture.weights=cc.optim, aic.optim=aic.optim, df.weights.optim=degree.freedom.weights.optim, log.lik.margin.optim=log.lik.margin.optim, p=p, factor.names=factor.names, risk.names=risk.names)
#aic.vector=aic.vector,df.weights.vector=degree.freedom.weights.vector,log.lik.margin.vector=log.lik.margin.vector


}  #end of function CompetingRiskFrailtyOptim(...)
