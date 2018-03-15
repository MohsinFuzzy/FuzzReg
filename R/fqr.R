#' Fuzzy Quantile Regression
#'
#' @description  It gives the estimates of fuzzy quantile regression using the method of Weighted Least Absolute Deviation (WLAD). It converts the input variables into Linear Programming Problem (LPP) and uses the Simplex Algorithm to solve the LPP.
#'
#' @usage fqr(X,y_left,y_centre,y_right,t,type)
#'
#' @param X is an input fuzzy number
#' @param y_left is an left output fuzzy number
#' @param y_centre is an centre of output fuzzy number
#' @param y_right is an right output fuzzy number
#' @param t is a spesified quantile ranges from 0 to 1 i.e t=[0,1]
#' @param type spesifies the model (1 or 2)
#' @param 1 "Fuzzy output, Fuzzy input and Fuzzy Parameters"
#' @param 2 "Fuzzy output, Crisp input and Fuzzy Parameters"
#'
#' @author Mohsin Shahzad
#' @examples If given Triangular Fuzzy NUmber
#'
#' library ("lpSolve")
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
#'
#' t<-0.5
#' fqr(X,y_left,y_centre,y_right,t,type=1)
#'
fqr  <-  function(x,Yl,Ym,Yr, t, type, intercept=TRUE) {
  if (type==1){
    if (intercept) Xl  <-  cbind(1, x[,1]) else X <-  cbind(x)
    Nl   <-  length(Yl)
    nl  <-  nrow(Xl)
    stopifnot(nl == Nl)
    pl  <-  ncol(Xl)
    cl  <-  c(rep(t, nl), rep(1-t, nl), rep(0, 2*pl))
    Al  <- cbind(diag(nl), -diag(nl), Xl, -Xl)
    resl  <-  lp("min", cl, Al, "=", Yl, compute.sens=1)
    sol_l <- resl$solution
    coef1l  <-  sol_l[(2*nl+1):(2*nl+2*pl)]
    coefl <- numeric(length=pl)
    for (i in seq(along=coefl)) {
      coefl[i] <- (if(coef1l[i]<=0)-1 else +1) *  max(coef1l[i], coef1l[i+pl])
    }
    if (intercept) Xm  <-  cbind(1, x[,2]) else X <-  cbind(x)
    Nm   <-  length(Ym)
    nm  <-  nrow(Xm)
    stopifnot(nm == Nm)
    pm  <-  ncol(Xm)
    cm  <-  c(rep(t, nm), rep(1-t, nm), rep(0, 2*pm))
    Am  <- cbind(diag(nm), -diag(nm), Xm, -Xm)
    resm  <-  lp("min", cm, Am, "=", Ym, compute.sens=1)
    sol_m <- resm$solution
    coef1m  <-  sol_m[(2*nm+1):(2*nm+2*pm)]
    coefm <- numeric(length=pm)
    for (j in seq(along=coefm)) {
      coefm[j] <- (if(coef1m[j]<=0)-1 else +1) *  max(coef1m[j], coef1m[j+pm])
    }
    if (intercept) Xr <-  cbind(1, x[,3]) else X <-  cbind(x)
    Nr   <-  length(Yr)
    nr  <-  nrow(Xr)
    stopifnot(nr == Nr)
    pr  <-  ncol(Xr)
    cr  <-  c(rep(t, nr), rep(1-t, nr), rep(0, 2*pr))
    Ar  <- cbind(diag(nr), -diag(nr), Xr, -Xr)
    resr  <-  lp("min", cr, Ar, "=", Yr, compute.sens=1)
    sol_r <- resr$solution
    coef1r  <-  sol_r[(2*nr+1):(2*nr+2*pr)]
    coefr <- numeric(length=pr)
    for (k in seq(along=coefr)) {
      coefr[k] <- (if(coef1r[k]<=0)-1 else +1) *  max(coef1r[k], coef1r[k+pr])
    }
  }
  if (type==2){
    if (intercept) Xl  <-  cbind(1, x) else X <-  cbind(x)
    Nl   <-  length(Yl)
    nl  <-  nrow(Xl)
    stopifnot(nl == Nl)
    pl  <-  ncol(Xl)
    cl  <-  c(rep(t, nl), rep(1-t, nl), rep(0, 2*pl))
    Al  <- cbind(diag(nl), -diag(nl), Xl, -Xl)
    resl  <-  lp("min", cl, Al, "=", Yl, compute.sens=1)
    sol_l <- resl$solution
    coef1l  <-  sol_l[(2*nl+1):(2*nl+2*pl)]
    coefl <- numeric(length=pl)
    for (i in seq(along=coefl)) {
      coefl[i] <- (if(coef1l[i]<=0)-1 else +1) *  max(coef1l[i], coef1l[i+pl])
    }
    if (intercept) Xm  <-  cbind(1, x) else X <-  cbind(x)
    Nm   <-  length(Ym)
    nm  <-  nrow(Xm)
    stopifnot(nm == Nm)
    pm  <-  ncol(Xm)
    cm  <-  c(rep(t, nm), rep(1-t, nm), rep(0, 2*pm))
    Am  <- cbind(diag(nm), -diag(nm), Xm, -Xm)
    resm  <-  lp("min", cm, Am, "=", Ym, compute.sens=1)
    sol_m <- resm$solution
    coef1m  <-  sol_m[(2*nm+1):(2*nm+2*pm)]
    coefm <- numeric(length=pm)
    for (j in seq(along=coefm)) {
      coefm[j] <- (if(coef1m[j]<=0)-1 else +1) *  max(coef1m[j], coef1m[j+pm])
    }
    if (intercept) Xr <-  cbind(1, x) else X <-  cbind(x)
    Nr   <-  length(Yr)
    nr  <-  nrow(Xr)
    stopifnot(nr == Nr)
    pr  <-  ncol(Xr)
    cr  <-  c(rep(t, nr), rep(1-t, nr), rep(0, 2*pr))
    Ar  <- cbind(diag(nr), -diag(nr), Xr, -Xr)
    resr  <-  lp("min", cr, Ar, "=", Yr, compute.sens=1)
    sol_r <- resr$solution
    coef1r  <-  sol_r[(2*nr+1):(2*nr+2*pr)]
    coefr <- numeric(length=pr)
    for (k in seq(along=coefr)) {
      coefr[k] <- (if(coef1r[k]<=0)-1 else +1) *  max(coef1r[k], coef1r[k+pr])
    }
  }
  result<-cbind(coefl,coefm,coefr)
  dimnames(result)=list(c("Intercept","Slope"),c("Left","Centre","Right"))
  return(list(c("Fuzzy Quantile Regression"),Quantile=t,coefficients=result))
}
