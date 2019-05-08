#' @title Aggreagated conditional socre test. 
#' 
#' @description 
#' Test for rare variant effects.
#' 
#' 
#' @details  
#' \code{ACST} provides p-value of a rare variant test.
#' 
#' @param y: A numeric vector of 1 (case) and 0 (control).
#' @param X: A numeric matrix. Each low corresponds to a vector of numbers of variants (0,1 or 2) in each indivisual.
#' @param B: A number of permutations for computing p-value. 
#' @param C: A value indicating independent (C=0) or correlated (C=1) structures among score staitstics. 
#' @param rho.set: A numeric vector of location of grids of a correlation parameter.
#' @param mc: A number of Monte Carlo samples for evaluating p-value in each correlation parameter. 
#' @return p-value (C=0), p-value and optimal correlation parameter (C=1).
#' 


ACST=function(y,X,B=200,C=0,rho.set=seq(-0.1,0.1,by=0.01),mc=2000){
  library(MASS)
  ind=ifelse(apply(X,2,sum)>0,1,0)
  X=X[,ind==1]
  M=dim(X)[2]; n=length(y)
  
  score=function(y,X){
    num1=(1:n)[y==1]; num0=(1:n)[y==0]
    n1=length(num1); n2=length(num0)
    res=c()
    for(k in 1:M){
      a=length((1:n1)[X[num1,k]==2])
      b=length((1:n1)[X[num1,k]==1])
      m1=length((1:n)[X[,k]==2])
      m2=length((1:n)[X[,k]==1])
      m3=length((1:n)[X[,k]==0])
      
      num=2*a+b-n1/n*(2*m1+m2)
      SD=sqrt(4*m1*m3+m1*m2+m2*m3)
      res[k]=num/SD*n*sqrt(n-1)/sqrt(n1*n2)
    }
    return(res)
  }
  
  if(C==0){
    origin.score=score(y,X)
    stat=sum(origin.score^2) 
    perm.stat=c()
    for(b in 1:B){
      yb=y[sample(1:n,n)];  perm.score=score(yb,X)
      perm.stat[b]=sum(perm.score^2)
    }
    pvalue=mean(ifelse(perm.stat>stat,1,0))
    return(pvalue)
  }
  
  if(C==1){ 
    set.seed(1)
    slen=length(rho.set)
    C=list()
    for(s in 1:slen){
      C[[s]]=solve((1-rho.set[s])*diag(M)+rho.set[s]*matrix(1,M,M))
    }
    
    calc.stat=function(y,X){
      cmat=cor(X); SS=score(y,X)
      pval=c()
      for(s in 1:slen){
        st=t(SS)%*%C[[s]]%*%SS; st=as.vector(st)
        gam=eigen(C[[s]]%*%cmat)$values
        rn=apply(matrix(rchisq(mc*M,1),M,mc)*gam,2,sum)
        pval[s]=mean(ifelse(rn>st,1,0))
      }
      return(list(min(pval),pval))
    }
    
    res=calc.stat(y,X)
    stat=res[[1]]; opt.rho=rho.set[which.min(res[[2]])]
    perm=c()
    for(b in 1:B){
      yb=y[sample(1:n,n)]
      perm[b]=calc.stat(yb,X)[[1]]
    }
    pvalue=mean(ifelse(perm<stat,1,0))
    res=c(pvalue,opt.rho)
    names(res)=c("p-value","optimal-rho")
    return(res)
  }
}








