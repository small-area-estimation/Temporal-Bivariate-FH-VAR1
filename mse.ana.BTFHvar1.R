###############################################################################
###############################################################################                                               
###       Analytical approximation of the MSE matrix 
###       
###       Autor: Esteban Cabello García
###
###       Work: BTFH-var1 

library(pracma)
library(mvtnorm)
library(expm)
library(ramify)
library(diagonals)
library(sae)
library(dplyr)
library(vctrs)

mse.ana.BTFHvar1 <- function(X, Sigma, F.inv, Ved){

  # X: Design matrix. Matrix.
  # Sigma: Variance components. Vector.
  # F.inv: Inverse of the Fisher information matrix.
  # Ved: Sampling error covariance matrices. List.


  tp <- dim(Ved[[1]])[1]/2
  D <- length(Ved)
  
  #### Creating V_{1,d} V_{2,d} and V_{u1} V_{u2} matrices

  V1d <- matrix(c(Sigma[1,],Sigma[3,]*sqrt(Sigma[1,]*Sigma[2,]),Sigma[3,]*sqrt(Sigma[1,]*Sigma[2,]),Sigma[2,]), byrow = T, nrow = 2)
  Phi <- diag(c(Sigma[7,],Sigma[8,]))
  V2dt <- matrix(c(Sigma[4,]/(1-Sigma[7,]^2), Sigma[6,]*sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,]), Sigma[6,]*sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,]),Sigma[5,]/(1-Sigma[8,]^2) ), nrow = 2, byrow = T)
  V2d <- V2d_f(V = V2dt,tp = tp, Phi = Phi)
  Vu1 <- Reduce(bdiag,rep(list(V1d),D))
  Vu2 <- Reduce(bdiag,rep(list(V2d),D))
  Vu <- bdiag(Vu1,Vu2)
  
  
  ###### First derivatives for V1,d and V2,d

  V2dt1 <- matrix(c(1, Sigma[3,]*sqrt(Sigma[2,])/(2*sqrt(Sigma[1,])), Sigma[3,]*sqrt(Sigma[2,])/(2*sqrt(Sigma[1,])), 0), ncol = 2, byrow = T)
  V2dt2 <- matrix(c(0, Sigma[3,]*sqrt(Sigma[1,])/(2*sqrt(Sigma[2,])), Sigma[3,]*sqrt(Sigma[1,])/(2*sqrt(Sigma[2,])), 1), ncol = 2, byrow = T)
  V2dt3 <- matrix(c(0, sqrt(Sigma[1,]*Sigma[2,]), sqrt(Sigma[1,]*Sigma[2,]), 0), ncol = 2, byrow = T)
  
  col1 <- matrix(1,nrow = tp)
  tcol1col1 <- tcrossprod(col1,col1)
  
  Vda1 <- kronecker(tcol1col1, V2dt1)
  Vda2 <- kronecker(tcol1col1, V2dt2)
  Vda3 <- kronecker(tcol1col1, V2dt3)
  
  V2dt4 <- matrix(c(1/(1-Sigma[7,]^2), Sigma[6,]*sqrt(Sigma[5,])/(2*sqrt(Sigma[4,])*(1-Sigma[7,]*Sigma[8,])), Sigma[6,]*sqrt(Sigma[5,])/(2*sqrt(Sigma[4,])*(1-Sigma[7,]*Sigma[8,])), 0), ncol = 2, byrow = T)
  V2dt5 <- matrix(c(0, Sigma[6,]*sqrt(Sigma[4,])/(2*sqrt(Sigma[5,])*(1-Sigma[7,]*Sigma[8,])), Sigma[6,]*sqrt(Sigma[4,])/(2*sqrt(Sigma[5,])*(1-Sigma[7,]*Sigma[8,])),1/(1-Sigma[8,]^2)), ncol = 2, byrow = T)
  V2dt6 <- matrix(c(0, sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,]), sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,]), 0), ncol = 2, byrow = T)
  Vda4 <- V2d_f(V = V2dt4,tp = tp, Phi = Phi)
  Vda5 <- V2d_f(V = V2dt5,tp = tp, Phi = Phi)
  Vda6 <- V2d_f(V = V2dt6,tp = tp, Phi = Phi)
  
  V2dt7 <- matrix(c(2*Sigma[7,]*Sigma[4,]/(1-Sigma[7,]^2)^2, Sigma[8,]*Sigma[6,]*sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,])^2, Sigma[8,]*Sigma[6,]*sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,])^2,0), ncol = 2, byrow = T)
  V2dt8 <- matrix(c(0, Sigma[7,]*Sigma[6,]*sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,])^2, Sigma[7,]*Sigma[6,]*sqrt(Sigma[4,]*Sigma[5,])/(1-Sigma[7,]*Sigma[8,])^2, 2*Sigma[8,]*Sigma[5,]/(1-Sigma[8,]^2)^2), ncol = 2, byrow = T)
  
  Phi7 <- diag(c(1,0)) 
  Phi8 <- diag(c(0,1)) 
  
  Vda7 <- V2da_78f(V = V2dt7, V2dt = V2dt, tp = tp, Phi = Phi, Phi78 = Phi7)
  Vda8 <- V2da_78f(V = V2dt8, V2dt = V2dt, tp = tp, Phi = Phi, Phi78 = Phi8)
  
  ### W_id, W_ud, V, W_i, W_u, R and L(i)
  
  Wid <- list(Vda1,Vda2,Vda3,Vda4,Vda5,Vda6,Vda7,Vda8)
  Wud <- kronecker(tcrossprod(col1,col1), V1d) + V2d
  Wi  <- lapply(Wid, function(x) as(Reduce(bdiag,rep(list(x),D)),"dgCMatrix"))
  Wu  <- Reduce(bdiag,rep(list(Wud),D))

  I2TD <- diag(1,2*tp*D)

  Vd <- lapply(Ved, function(x) x + Wud)
  V <- Reduce(bdiag, Vd)
  V.inv <- solve(V)
  R <- Wu%*%V.inv
  Li <- lapply(Wi, function(x) (I2TD-R)%*%x%*%V.inv)
  
  
  ##### g1 component
  
  Z1 <- Reduce(bdiag,rep(list(Reduce(rbind,rep(list(diag(2)),tp))),D))
  Z2d<- Reduce(bdiag,rep(list(diag(c(1,1))),tp))
  Z2 <- diag(rep(1,2*D*tp))
  Z  <- cbind(Z1,Z2)
  tZ <- t(Z)
  
  Tm <- Vu - Vu%*%tZ%*%V.inv%*%Z%*%Vu
  
  G1 <- Z%*%Tm%*%tZ
  
  
  ###### g2 component
  
  tX <- t(X)
  Q <- solve(tX%*%V.inv%*%X)
  Ve <- bdiag(Ved)
  Ve.inv <- solve(Ve)
  G2 <- (X-Z%*%Tm%*%tZ%*%Ve.inv%*%X)%*%Q%*%(tX-tX%*%Ve.inv%*%Z%*%Tm%*%tZ)
  
  
  ###### g3 component 
  
  G3 <- list()
  
  h <- 0
  for(i in 1:8){
    for(j in 1:8){
      G3[[h+1]] <- F.inv[i,j]*Li[[i]]%*%V%*%t(Li[[j]])
      h <- h+1
    }
  }
  
  G3 <- Reduce("+", G3)
  
  
  #### MSE
  
  mse <- diag(G1 + G2 + 2*G3)
  mse <- matrix(mse, ncol = 2, byrow = T)
  
  return(list(G1 = G1, G2 = G2, G3 = 2*G3, mse.ana = mse))
  
}



