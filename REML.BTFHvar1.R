###############################################################################
###############################################################################
###
###       Fecha de creación: 02/10/2023
###                                                  
###       Ultima modificación: 02/10/2023
###       
###       Autor: Esteban Cabello García
###
###       Trabajo: BTFH-var1 - REML

library(expm)
library(diagonals)
library(dplyr)

V2d_f <- function(V,tp,Phi){
  V2d_elem<- list() 
  cindx <- merge(1:tp, 1:tp)
  for(k in 1:dim(cindx)[1]){ ## We calculete each cov(u_t1,u_t2) and cov(u_t2,u_t1)
    if(cindx[k,1] < cindx[k,2])
      V2d_elem[[k]] <- V%*%(Phi%^%(cindx[k,2]-cindx[k,1]))
    else 
      V2d_elem[[k]] <- (Phi%^%(cindx[k,1]-cindx[k,2]))%*%V
  }
  V2d_row <- list()
  for(k in 1:tp){
    j <- split_vector(1:(tp^2), steps = tp)
    V2d_row[[k]] <-  suppressMessages(bind_cols(V2d_elem[j[[k]]], col_names = F))
  }
  
  V2d <- as.matrix(Reduce(rbind,V2d_row))
  V2d <- V2d[,-dim(V2d)[2]]
  return(V2d)
}

V2da_78f <- function(V,V2dt,tp,Phi,Phi78){
  V2d_elem<- list() 
  cindx <- merge(1:tp, 1:tp)
  for(k in 1:dim(cindx)[1]){ ## We calculete each cov(u_t1,u_t2) and cov(u_t2,u_t1)
    if(cindx[k,1] < cindx[k,2])
      V2d_elem[[k]] <- V%*%(Phi%^%(cindx[k,2]-cindx[k,1])) + (cindx[k,2]-cindx[k,1])*V2dt%*%(Phi%^%(cindx[k,2]-cindx[k,1]-1))%*%Phi78
    else if(cindx[k,1] > cindx[k,2]) 
      V2d_elem[[k]] <- (Phi%^%(cindx[k,1]-cindx[k,2]))%*%V + (cindx[k,1]-cindx[k,2])*(Phi%^%(cindx[k,1]-cindx[k,2]-1))%*%Phi78%*%V2dt
    else
      V2d_elem[[k]] <- V
  }
  V2d_row<- list()
  for(k in 1:tp){
    j <- split_vector(1:(tp^2), steps = tp)
    V2d_row[[k]] <-  Reduce(rbind,V2d_elem[j[[k]]])
  }
  
  V2d <- as.matrix(Reduce(cbind,V2d_row))
  return(V2d)
}

## X = agregated data X (matrix form)
## ydi = target variable (matrix)
## sigma.ini = initial values for sigma (theta parameters)
## Ve = sampling variance diag(1<d<D,diag(1<t<T)) LIST


REML.BTFHvar1 <- function(X,y, sigma.ini, Ve, PRECISION = 0.0001, MAXITER = 40){
  
  kit <- 0
  FLAG <- 0
  diff <- rep(PRECISION + 1, 8)
  sigmau1 <- sigma.ini
  
  col1 <- matrix(1,nrow = tp)
  tcol1col1 <- tcrossprod(col1,col1)
  tX <- t(X)
  ty <- t(y)
  
  Xd <- list() ### Aggregated X by province (Xd)
  j <- 1
  for(i in seq(1,2*D*tp,by = 2*tp)){
    Xd[[j]] <- X[i:(i+2*tp-1),]
    j <- j+1
  }
  
  tp <- dim(Ve[[1]])[1]/2
  D <- length(y)/(2*tp)
  
  while (any(diff > rep(PRECISION, 8)) & (kit < MAXITER)){

    kit <- kit + 1
    
    #print(kit)
    #print(sigmau1)
    
    V1d <- matrix(c(sigmau1[1,],sigmau1[3,]*sqrt(sigmau1[1,]*sigmau1[2,]),sigmau1[3,]*sqrt(sigmau1[1,]*sigmau1[2,]),sigmau1[2,]), byrow = T, nrow = 2)
    
    Phi <- diag(c(sigmau1[7,],sigmau1[8,]))
    V2dt <- matrix(c(sigmau1[4,]/(1-sigmau1[7,]^2), sigmau1[6,]*sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,]), sigmau1[6,]*sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,]),sigmau1[5,]/(1-sigmau1[8,]^2) ), nrow = 2, byrow = T)
    V2d <- V2d_f(V = V2dt,tp = tp, Phi = Phi)
    
    
    Vd <- lapply(Ve, function(x) kronecker(tcol1col1, V1d) + V2d + x)
    Vd.inv <- lapply(Vd,solve)
    Qd <- Map("%*%",lapply(Xd,t),Map("%*%", Vd.inv, Xd))
    Q  <- Reduce("+",Qd)
    
    #### Cálculo de la matriz P mediante Pd1d2_(d1,d2 =1,...,D)
    Pd1 <- lapply(Map("%*%", Vd.inv,Xd), function(x) x%*%solve(Q))
    Pd2 <- Map("%*%", lapply(Xd,t),Vd.inv)
    Pdd_row <- list()
    Pdd_ij <- list()
    for(j in 1:D){
      h <- 0
      for(k in 1:D){
        h <- h + 1
        if(j == k){
          Pdd_ij[[h]] <- Vd.inv[[j]]-Pd1[[j]]%*%Pd2[[j]]
        }
        else{
          Pdd_ij[[h]] <- -1*Pd1[[j]]%*%Pd2[[k]]
        }
      }
      Pdd_row[[j]] <- Reduce(cbind2,Pdd_ij)
    }
    P <- Reduce(rbind2,Pdd_row)
    
    
    V2dt1 <- matrix(c(1, sigmau1[3,]*sqrt(sigmau1[2,])/(2*sqrt(sigmau1[1,])), sigmau1[3,]*sqrt(sigmau1[2,])/(2*sqrt(sigmau1[1,])), 0), ncol = 2, byrow = T)
    V2dt2 <- matrix(c(0, sigmau1[3,]*sqrt(sigmau1[1,])/(2*sqrt(sigmau1[2,])), sigmau1[3,]*sqrt(sigmau1[1,])/(2*sqrt(sigmau1[2,])), 1), ncol = 2, byrow = T)
    V2dt3 <- matrix(c(0, sqrt(sigmau1[1,]*sigmau1[2,]), sqrt(sigmau1[1,]*sigmau1[2,]), 0), ncol = 2, byrow = T)
    
    Vda1 <- kronecker(tcol1col1, V2dt1)
    Vda2 <- kronecker(tcol1col1, V2dt2)
    Vda3 <- kronecker(tcol1col1, V2dt3)
    
    V2dt4 <- matrix(c(1/(1-sigmau1[7,]^2), sigmau1[6,]*sqrt(sigmau1[5,])/(2*sqrt(sigmau1[4,])*(1-sigmau1[7,]*sigmau1[8,])), sigmau1[6,]*sqrt(sigmau1[5,])/(2*sqrt(sigmau1[4,])*(1-sigmau1[7,]*sigmau1[8,])), 0), ncol = 2, byrow = T)
    V2dt5 <- matrix(c(0, sigmau1[6,]*sqrt(sigmau1[4,])/(2*sqrt(sigmau1[5,])*(1-sigmau1[7,]*sigmau1[8,])), sigmau1[6,]*sqrt(sigmau1[4,])/(2*sqrt(sigmau1[5,])*(1-sigmau1[7,]*sigmau1[8,])),1/(1-sigmau1[8,]^2)), ncol = 2, byrow = T)
    V2dt6 <- matrix(c(0, sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,]), sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,]), 0), ncol = 2, byrow = T)
    Vda4 <- V2d_f(V = V2dt4,tp = tp, Phi = Phi)
    Vda5 <- V2d_f(V = V2dt5,tp = tp, Phi = Phi)
    Vda6 <- V2d_f(V = V2dt6,tp = tp, Phi = Phi)
    
    
    V2dt7 <- matrix(c(2*sigmau1[7,]*sigmau1[4,]/(1-sigmau1[7,]^2)^2, sigmau1[8,]*sigmau1[6,]*sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,])^2, sigmau1[8,]*sigmau1[6,]*sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,])^2,0), ncol = 2, byrow = T)
    V2dt8 <- matrix(c(0, sigmau1[7,]*sigmau1[6,]*sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,])^2, sigmau1[7,]*sigmau1[6,]*sqrt(sigmau1[4,]*sigmau1[5,])/(1-sigmau1[7,]*sigmau1[8,])^2, 2*sigmau1[8,]*sigmau1[5,]/(1-sigmau1[8,]^2)^2), ncol = 2, byrow = T)
    
    Phi7 <- diag(c(1,0)) 
    Phi8 <- diag(c(0,1)) 
    
    Vda7 <- V2da_78f(V = V2dt7, V2dt = V2dt, tp = tp, Phi = Phi, Phi78 = Phi7)
    Vda8 <- V2da_78f(V = V2dt8, V2dt = V2dt, tp = tp, Phi = Phi, Phi78 = Phi8)
    
    
    Vda <- list(Vda1, Vda2, Vda3, Vda4, Vda5, Vda6, Vda7, Vda8)
    Va  <- lapply(Vda, function(x) as(Reduce(bdiag,rep(list(x),D)),"dgCMatrix"))
    
    yP  <- ty %*% P
    Py  <- P %*% y
    PVa <- lapply(Va,function(x) P%*%x)
    S   <- matrix(sapply(Va,function(x) -0.5*sum(diag(P%*%x)) + as.vector(0.5*yP%*%x%*%Py)),ncol = 1)
    
    
    u <- c()
    for(j in 1:8){ 
      u <- c(u,1:j)
    }
    
    PVa.a <- PVa[rep(1:8,times = 1:8)]
    PVa.b <- PVa[u]
    
    preF0 <- sapply(Map("%*%", PVa.a, PVa.b), function(x) 0.5*sum(diag(x)))
    F0 <- matrix(0, ncol = 8, nrow = 8, byrow = T) 
    F0[upper.tri(F0,diag = T)] <- preF0
    F0[lower.tri(F0,diag = T)] <- t(F0)[lower.tri(t(F0),diag = T)]
    
    F.inv <- try(solve(F0))

    
    if(any(class(F.inv) == "try-error") | kit == MAXITER){
      stop("Not inversible matrix or not convergence")
      
    }
    else {
      
      sigmau <- sigmau1 + F.inv%*%S
      
      if(any(sigmau[c(1,2,4,5),] < 0) | any(abs(sigmau[c(3,6,7,8),]) > 1)){
        warning("Out of parametric space")
        # stop("Out of parametric space")
        break
      }
      
      
      diff <- abs(sigmau-sigmau1)
      
      V1d <- matrix(c(sigmau[1,],sigmau[3,]*sqrt(sigmau[1,]*sigmau[2,]),sigmau[3,]*sqrt(sigmau[1,]*sigmau[2,]),sigmau[2,]), byrow = T, nrow = 2)
      Phi <- diag(c(sigmau[7,],sigmau[8,]))
      V2dt <- matrix(c(sigmau[4,]/(1-sigmau[7,]^2), sigmau[6,]*sqrt(sigmau[4,]*sigmau[5,])/(1-sigmau[7,]*sigmau[8,]), sigmau[6,]*sqrt(sigmau[4,]*sigmau[5,])/(1-sigmau[7,]*sigmau[8,]),sigmau[5,]/(1-sigmau[8,]^2) ), nrow = 2, byrow = T)
      V2d <- V2d_f(V = V2dt,tp = tp, Phi = Phi)
      Vd <- lapply(Ve, function(x) kronecker(tcol1col1, V1d) + V2d + x)
      Vd.inv <- lapply(Vd,solve)
      V.inv <- as(bdiag(Vd.inv),"dgCMatrix")
      
      sigmau1 <- sigmau
      betaREML <- as.matrix(solve(tX%*%V.inv%*%X)%*%tX%*%V.inv%*%y)
      
      
      
      
    }
    
    
  }
  
  
  
  return(list(Sigma = sigmau1, Beta = betaREML, V.inv = V.inv, Fisher.inv = F.inv, FLAG = FLAG, ITER = kit))
}


