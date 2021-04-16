library(phytools)
library(ape)
library(expm)

# function uses numerical optimization to solve for the stationary distribution
# written by Liam J. Revell 2013

#   ____________________________________________________________________________
#   Tool functions                                                          ####

statdist<-function(Q){
  foo<-function(theta,Q){
    Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
    sum((Pi%*%Q)^2)
  }
  k<-nrow(Q)
  if(nrow(Q)>2){
    fit<-optim(rep(1/k,k-1),foo,Q=Q,control=list(reltol=1e-16))
    return(setNames(c(fit$par[1:(k-1)],1-sum(fit$par[1:(k-1)])),rownames(Q)))
  } else {
    fit<-optimize(foo,interval=c(0,1),Q=Q)
    return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
  }
}

#' Title
#'
#' @param m
#' @param q
#' @param index.matrix
#'
#' @return
#' @export
#'
#' @examples
makeQ<-function(m,q,index.matrix){
  Q<-matrix(0,m,m)
  Q[]<-c(0,q)[index.matrix+1]
  diag(Q)<-0
  diag(Q)<--rowSums(Q)
  Q
}

#' Title
#'
#' @param m number of character states
#' @param q number of param
#' @param index.matrix
#'
#' @return Transition matrix in the correct form for Felsensteins pruning algorithm.
#' @export
#'
#' @examples
makeQ10X<-function(m,q,index.matrix){ #For binary optimization for COVID
  Q<-matrix(0,m,m)
  Q[]<-c(0,q)[index.matrix+1]
  diag(Q)<-0
  Q[2,1] = 10*Q[1,2]
  diag(Q)<--rowSums(Q)
  Q
}

## wraps around expm
## written by Liam Revell 2011, 2017
#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#' @importFrom expm expm
#' @examples
EXPM<-function(x,...){
  e_x<-if(isSymmetric(x)) matexpo(x) else expm(x,...)
  dimnames(e_x)<-dimnames(x)
  e_x
}


#' Title
#'
#' @param tree Phylogenetic tree
#' @param x Tip character states that are used to inform the transition matrix
#' @param model Model matrix ordered where unique numbers indicate individual parameters.
#'
#' @return
#' @export
#'
#' @examples
FitTransition10x<-function(tree,x,model=NA){

  output.liks<-TRUE
  q.init<-length(unique(x))/sum(tree$edge.length)
  opt.method<-"nlminb"
  min.q<-1e-12
  N<-Ntip(tree)
  M<-tree$Nnode

  if(is.matrix(x)){
    x<-x[tree$tip.label,]
    m<-ncol(x)
    states<-colnames(x)
  } else {
    x<-to.matrix(x,sort(unique(x)))
    x<-x[tree$tip.label,]
    m<-ncol(x)
    states<-colnames(x)
  }
  pi<-setNames(rep(1/m,m),states)

  rate<-model
  k<-max(rate)
  Q<-matrix(0,m,m)

  #Set transition matrix to the one provided by the user.
  index.matrix<-rate
  tmp<-cbind(1:m,1:m)
  rate[tmp]<-0
  rate[rate==0]<-k+1
  liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))

  pw<-reorder(tree,"pruningwise")

  #Felsensteins pruning algorithm:
  lik<-function(Q,output.liks=FALSE,pi){
    if(any(is.nan(Q))||any(is.infinite(Q))) return(1e50)
    comp<-vector(length=N+M,mode="numeric")
    parents<-unique(pw$edge[,1])
    root<-min(parents)
    for(i in 1:length(parents)){
      anc<-parents[i]
      ii<-which(pw$edge[,1]==parents[i])
      desc<-pw$edge[ii,2]
      el<-pw$edge.length[ii]
      v<-vector(length=length(desc),mode="list")
      for(j in 1:length(v)){
        v[[j]]<-EXPM(Q*el[j])%*%liks[desc[j],]
      }
      vv<-if(anc==root) Reduce('*',v)[,1]*pi else Reduce('*',v)[,1]
      comp[anc]<-sum(vv)
      liks[anc,]<-vv/comp[anc]
    }
    if(output.liks)return(liks[1:M+N,,drop=FALSE])
    logL<--sum(log(comp[1:M+N]))
    return(if(is.na(logL)) Inf else logL)
  }

  #Optimize the fit to find the transition matrix:
  if(length(q.init)!=k) q.init<-rep(q.init[1],k) #Remember to reset q.init if starting from here.
  q.init = c(q.init)


  #fit<-nlminb(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),lower=c(rep(0,k),-1e50),upper=rep(1e50,k+1), control=list(iter.max=400, eval.max=400, trace=5))
  #Unconstrained fit
  fit_unconstrained<-optim(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),method="Nelder-Mead",control=list(trace=5, maxit=5000), hessian =T)
  #Constrained fit.
  fit_constrained10X<-optim(q.init,function(p) lik(makeQ10X(m,p,index.matrix),pi=pi),method="Nelder-Mead",control=list(trace=5, maxit=5000), hessian =T)
  fit_constrained10X_2<-optim(q.init,function(p) lik(makeQ10X(m,p,index.matrix),pi=pi),method="BFGS", control=list(trace=5, maxit=5000), hessian =T)
  fit<-nlminb(q.init,function(p) lik(makeQ10X(m,p,index.matrix),pi=pi),lower=rep(0,k),upper=rep(1e50,k),hessian = T)

  #Print the estimated matrix
  print("Unconstrained we get: ")
  print(makeQ(m, fit_unconstrained$par, index.matrix))
  print("Constrained to 10X using Nelder Mead")
  print(makeQ10X(m,fit_constrained10X$par,index.matrix))
  print("Constrained to 10X using BFGS")
  print(makeQ10X(m,fit_constrained10X_2$par,index.matrix))
  print("Constrained to 10X using nlminb")
  print(makeQ10X(m,fit$par,index.matrix))


#   ____________________________________________________________________________
#   The sections below are not needed as we only want to print the matrix   ####
#   for our purposes.                                                       ####

  #To calculate the variance of the parameter estimate we take the inverse of the Hessian matrix
  #E.G.
  # fisher_info = solve(fit_constrained10X_2$hessian)
  # prop_sigma<-sqrt(diag(fisher_info))
  # prop_sigma
  # upper<-fit_constrained10X_2$par+1.96*prop_sigma
  # lower<-fit_constrained10X_2$par-1.96*prop_sigma
  # c(lower,fit_constrained10X_2$par,upper)
  #And the interpretation:
  #1/c(lower,fit_constrained10X_2$par,upper)

  # fit = fit_constrained10X
  #
  # opt.method="optim"
  # obj<-list(logLik=
  #             if(opt.method=="optim") -fit$value else -fit$objective,
  #           rates=fit$par,
  #           index.matrix=index.matrix,
  #           states=states,
  #           pi=pi,
  #           method=opt.method)
  #
  # if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
  #                                  pi=pi)
  #
  # lik.f<-function(q) -lik(q,output.liks=FALSE,pi=pi)
  # obj$lik<-lik.f
  # class(obj)<-"fitMk"
  # return(obj)
}

# mod = rbind(c(1,2), c(3,4))
# fit10X = fitTransition_10X(Covid19_tree_norway, Locations,model=mod, pi=c(0,1))

