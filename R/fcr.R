#' Fuzzy Censored Regression
#'
#' @description  It gives the estimates of Fuzzy Censored Regression by using the iterative based method of Newton Raphson.
#'
#' @usage fcr(Y,X,lower=c(L1,L2,L3),upper=c(U1,U2,U3),type)
#'
#' @param x is an input fuzzy number
#' @param y is an output fuzzy number
#' @param lower is a set contains lower end censored observations of left, centre and right of triangular fuzzy number.Default lower limit is zero "0".
#' @param upper is a set contains Upper end censored observations of left, centre and right of triangular fuzzy number.Default upper limit is infinity "Inf".
#' @param type is model to be used (1 or 2)
#' @param 1 "Fuzzy output, Fuzzy input and Fuzzy Parameters"
#' @param 2 "Fuzzy output, Crisp input and Fuzzy Parameters"

#' @author Mohsin Shahzad
#'
#' @examples If given Triangular Fuzzy NUmber
#'
#' library("VGAM")
#'
#' x_left<-c(1.5,3.0,4.5,6.5,8.0,9.5,10.5,12.0)
#' x_centre<-c(2.0,3.5,5.5,7.0,8.5,10.5,11.0,12.5)
#' x_right<-c(2.5,4.0,6.5,7.5,9.0,11.5,11.5,13.0)
#'
#' y_left<-c(3.5,5.0,6.5,6.0,8.0,7.0,10.0,9.0)
#' y_centre<-c(4.0,5.5,7.5,6.5,8.5,8.0,10.5,9.5)
#' y_right<-c(4.5,6.0,8.5,7.0,9.0,9.0,11.0,10.0)
#'
#' X<-cbind(x_left,x_centre,x_right)
#' Y<-cbind(y_left,y_centre,y_right)
#'
#' fcr(Y,X,lower=c(2.0,3.5,5.0),upper=c(8.5,9.0,10.0),type=1)
#'
#'
fcr<-function(y,x,lower=c(L1=0,L2=0,L3=0),upper=c(U1=Inf,U2=Inf,U3=Inf),type){

if (type==1){
    Il<-rep(0,1)		#Initial Guess
  al=bl=sigl=Il
  model_l<-vglm(y[,1]~x[,1],tobit(Lower=lower[1],Upper=upper[1]))	#model estimation
  coefl<-coef(model_l,matrix=TRUE)[,1]		# coefficients "a" and "b"
  logsigmal<-coef(model_l,matrix=TRUE)[1,2]		#logsigma

  Im<-rep(0,1)		#Initial Guess
  am=bm=sigm=Im
  model_m<-vglm(y[,2]~x[,2],tobit(Lower=lower[2],Upper=upper[2]))	#model estimation
  coefm<-coef(model_m,matrix=TRUE)[,1]		# coefficients "a" and "b"
  logsigmam<-coef(model_m,matrix=TRUE)[1,2]

  Ir<-rep(0,1)		#Initial Guess
  ar=br=sigr=Ir
  model_r<-vglm(y[,3]~x[,3],tobit(Lower=lower[3],Upper=upper[3]))	#model estimation
  coefr<-coef(model_r,matrix=TRUE)[,1]		# coefficients "a" and "b"
  logsigmar<-coef(model_r,matrix=TRUE)[1,2]
}
  if (type==2){
    Il<-rep(0,1)		#Initial Guess
    al=bl=sigl=Il
    model_l<-vglm(y[,1]~x,tobit(Lower=lower[1],Upper=upper[1]))	#model estimation
    coefl<-coef(model_l,matrix=TRUE)[,1]		# coefficients "a" and "b"
    logsigmal<-coef(model_l,matrix=TRUE)[1,2]		#logsigma

    Im<-rep(0,1)		#Initial Guess
    am=bm=sigm=Im
    model_m<-vglm(y[,2]~x,tobit(Lower=lower[2],Upper=upper[2]))	#model estimation
    coefm<-coef(model_m,matrix=TRUE)[,1]		# coefficients "a" and "b"
    logsigmam<-coef(model_m,matrix=TRUE)[1,2]

    Ir<-rep(0,1)		#Initial Guess
    ar=br=sigr=Ir
    model_r<-vglm(y[,3]~x,tobit(Lower=lower[3],Upper=upper[3]))	#model estimation
    coefr<-coef(model_r,matrix=TRUE)[,1]		# coefficients "a" and "b"
    logsigmar<-coef(model_r,matrix=TRUE)[1,2]
  }

  lim<-cbind(lower,upper)
  ans<-cbind(coefl,coefm,coefr)
  sigmaa<-rbind(logsigmal,logsigmam,logsigmar)
  dimnames(lim)=list(c("Left","Centre","Right"),c("Lower","Upper"))
  dimnames(ans)=list(c("Intercept","Slope"),c("Left","Centre","Right"))
  dimnames(sigmaa)=list(c("Left","Centre","Right"))
  return(list(c("Fuzzy Censored Regression"),limits=lim,coefficient=ans,logsigma=sigmaa))
}
