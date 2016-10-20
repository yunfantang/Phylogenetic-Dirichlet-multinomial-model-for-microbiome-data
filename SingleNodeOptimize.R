
library(HMP)
library(nloptr)

loglik <- function(x,t){
  l <- 0
  ts <- sum(t)
  nc <- ncol(x)
  for(j in 1:nrow(x)){
    if(sum(x[j,])!=0)
      l <- l - sum(unlist(lapply(list(1:sum(x[j,])),ll,t=ts) ))
    # maybe not    lij <- 0 ## NEW LINE: This initializes lij to zero for each row. Otherwise lij is a cummulant. ##
    for(i in 1:nc){
      if(x[j,i]==0) lij <- 0
      else lij <- sum(unlist(lapply(list(1:(x[j,i])),ll,t=t[i])))
      l <- l + lij
    }
  }
  l
}
ulk <- function(x,lambda,kappa){
  S <- numeric(2)
  x1 = x[,1]; x2 = x[,2]
  for(j in 1:nrow(x)){
    Sk = numeric(2); Sl = 0
    if(x1[j]!=0) Sk[1] = sum(unlist(lapply(list(1:x1[j]),uu,t=kappa*lambda)))
    if(x2[j]!=0) Sk[2] = sum(unlist(lapply(list(1:x2[j]),uu,t=(1-kappa)*lambda)))
    S[2] = S[2] + lambda*(Sk[1]-Sk[2])
    if(x1[j]+x2[j]!=0) Sl = sum(unlist(lapply(list(1:(x1[j]+x2[j])),uu,t=lambda)))
    S[1] = S[1] + kappa*Sk[1]+(1-kappa)*Sk[2]- Sl
  }
  S
}

ll <- function(x,t) log(t+x-1)
uu <- function(x,t) 1/(t+x-1)
ff <- function(x,t) 1/(t+x-1)^2
xloglik = function(z,x){
  pi = z[1]; nux = z[2]; 
  -loglik(x,nux*c(pi,1-pi)) 
}
xu = function(z,x){
  Sx = ulk(x,z[2],z[1])
  -c(Sx[2],Sx[1])
}

OptimizeNode = function(x, algo="NLOPT_LD_SLSQP"){
  DM_node = DM.MoM(x)
  opts = list("algorithm"=algo,"print_level"=1,"xtol_rel"=1e-6)
  EstNode_SLSQP = nloptr(x0=c(DM_node$pi[1], sum(DM_node$gamma)),eval_f=xloglik,eval_grad=xu,lb=c(0,0),ub=c(1,Inf),opts=opts,x=x)
  gradnorm = sum(xu(EstNode_SLSQP$solution,x)^2)
  c(-EstNode_SLSQP$objective, gradnorm)
}



  
 