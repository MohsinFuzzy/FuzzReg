#' Fuzzy Linear Regression
#'
#' @description  It gives the estimates of fuzzy regression by using the method of least squares.
#'
#' @usage flr(Y,X,type)
#'
#' @param x is an input fuzzy number
#' @param y is an output fuzzy number
#' @param type is model to be used (1 or 2)
#' @param 1 "Fuzzy output, Fuzzy input and Fuzzy Parameters"
#' @param 2 "Fuzzy output, Crisp input and Fuzzy Parameters"

#' @author Mohsin Shahzad
#'
#' @examples If given Triangular Fuzzy NUmber
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
#' flr(Y,X,type=1)
#'
#'

flr<-function (Y,X,type){
  if (type==1){
    lm1<-lm(Y[,1]~X[,1])
    lm2<-lm(Y[,2]~X[,2])
    lm3<-lm(Y[,3]~X[,3])
  }
  if (type==2){
    lm1<-lm(Y[,1]~X)
    lm2<-lm(Y[,2]~X)
    lm3<-lm(Y[,3]~X)
  }
  coef<-cbind(lm1$coefficients,lm2$coefficients,lm3$coefficients)
  dimnames(coef)=list(c("Intercept","Slope"),c("Left","Centre","Right"))
  return(list(c("Fuzzy Linear Regression"),coefficients=coef))
}
