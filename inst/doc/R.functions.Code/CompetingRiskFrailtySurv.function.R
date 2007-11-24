"CompetingRiskFrailtySurv" <-
function(ID=ID,surv.time=surv.time,status=status)
{
status<-as.matrix(status)

if (!is.numeric(surv.time)) stop("not numeric times arguments")
if (!is.numeric(status)) stop("censoring indicators must be binary variables: 1 for death, 0 for censored")
  
    Surv.object <- cbind(ID,surv.time,status)
    dimnames(Surv.object) <- list(NULL, c("ID","time", paste("status.",1:ncol(status),sep="")))
    #dimnames(Surv.object) <- list(NULL, c("ID","time", colnames(status)))
    class(Surv.object) <- c("CompetingRiskFrailtySurv","matrix")
    Surv.object

}

