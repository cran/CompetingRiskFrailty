"CompetingRiskFrailtySurvfitCreate" <-
function(formula=formula, data=data, na.action=na.fail, control=control, risk.names=risk.names)  
{    
call <- match.call()


#check if formula missing
if (missing(formula)) stop("formula must be specified")

#matching of control parameters, if some absent, the defaults will be setup
if (missing(control)) control<-list()  
names.missing<-setdiff(names(CompetingRiskFrailtySurvfitControl()),names(control))
if (length(names.missing) >0) for (i in 1:length(names.missing)) control[[names.missing[i]]]<-CompetingRiskFrailtySurvfitControl()[[names.missing[i]]]     


#check for user specified formula parameter which can be: a direct call formula=CompetingRiskFrailtySurv(...), a name formula="name_of_CompetingRiskFrailtySurv_object" or formula="name_of_CompetingRiskFrailtySurv_object"~1 or ~x1+...+xp
if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("CompetingRiskFrailtySurv")) || inherits(formula, "CompetingRiskFrailtySurv"))
      {
        formula <- eval(parse(text = paste(deparse(call[[2]],width.cutoff=500),1, sep = "~")))  
        environment(formula) <- parent.frame()  #if formula=CompetingRiskFrailtySurv(...)~1, means there is only baseline in the model, without any covariates
    }

#check if (otherway) formula not correctly specified
if (!inherits(formula, "formula")) stop("formula is not an object of class formula")



#check for missing values NA, which are not allowed in the model
if (missing(data)) {if (any(is.na(eval(formula[[2]])))) stop("NA's must be removed or imputed from data bevore applied")} else
if ( any(is.na(data)) || any(is.na(eval(formula[[2]]))) ) stop("NA's must be removed or imputed from data bevore applied")



#construct a frame to be evaluated 
m <- match.call()
Terms <- terms(formula, "variables")
m$formula <- Terms
m$control<-NULL
m$risk.names<-NULL
m[[1]] <- as.name("model.frame")
m <- eval(m, parent.frame()) 
n <- nrow(m)
Y <- model.extract(m, "response")   #Y has named columns: ID, time, Risk.indicators as status.1,...,status.k
if (!inherits(Y,"CompetingRiskFrailtySurv")) stop("Response must be a survival object with specified survival time, birth time and censor indicator")
ll <- attr(Terms, "term.labels")
if (length(ll) == 0)  {X <- as.data.frame(rep(1, n));names(X)<-"Intercept"}  else X <- m[ll]  #X consists of covariates


data.set<-data.frame(Y,X)



##########################################################
#call to internal routine for optimization################
##########################################################
temp <- CompetingRiskFrailtyOptim(data.set=data.set,control=control,form=formula,risk.names=risk.names)   #X=covariates, Y=(ID,time,status), control=list with control parameters
        class(temp) <- "CompetingRiskFrailtySurvfit"
#    temp$call <- call


    temp
}

