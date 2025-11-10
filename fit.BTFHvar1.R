###############################################################################
###############################################################################
###
###       EBLUPs and Tables of the Regression and variance parameters  
###
###       Autor: Esteban Cabello García
###       Work: TBFH-var1


fit.BTFHvar1 <- function(X, y, D, tp, sigma, beta, V.inv, F.inv){
  
  ### X: Design matrix
  ### y: Target variables. Column Vector.
  ### D: domains. Integer.
  ### tp: time periods. Integer.
  ### sigma: Variance components estimates. Vector.
  ### beta: Regression parameters. Vector
  ### V.inv: inverse of V = V_u + V_e. Matrix.
  ### F.inv: inverse of the information Fisher. Matrix.
  
  tX <- t(X)
  col1 <- matrix(1,nrow = tp)
  tcol1col1 <- tcrossprod(col1,col1)
  
  V1d <- matrix(c(sigma[1,],sigma[3,]*sqrt(sigma[1,]*sigma[2,]),sigma[3,]*sqrt(sigma[1,]*sigma[2,]),sigma[2,]), byrow = T, nrow = 2)
  Phi <- diag(c(sigma[7,],sigma[8,]))
  V2dt <- matrix(c(sigma[4,]/(1-sigma[7,]^2), sigma[6,]*sqrt(sigma[4,]*sigma[5,])/(1-sigma[7,]*sigma[8,]), sigma[6,]*sqrt(sigma[4,]*sigma[5,])/(1-sigma[7,]*sigma[8,]),sigma[5,]/(1-sigma[8,]^2) ), nrow = 2, byrow = T)
  V2d <- V2d_f(V = V2dt,tp = tp, Phi = Phi)
  
  
  Z1 <-  Reduce(bdiag,rep(list(Reduce(rbind,rep(list(diag(2)),tp))),D))
  tZ1 <- t(Z1)
  Z2 <-  Reduce(bdiag,rep(list(Reduce(bdiag,rep(list(diag(2)),tp))),D))
  tZ2 <- t(Z2)
  Vu1 <- as.matrix(Reduce(bdiag,rep(list(V1d),D)))
  Vu2 <- as.matrix(Reduce(bdiag,rep(list(V2d),D)))
  
  u1E <- Vu1 %*% tZ1 %*% V.inv %*% (y-X%*%beta)
  u2E <- Vu2 %*% tZ2 %*% V.inv %*% (y-X%*%beta)
  
  
  # Eblup's
  muE <- X %*% beta + Z1 %*% u1E + Z2 %*% u2E
  
  ########## Table for sigma ########
  alpha <- 0.05
  se.theta <- sqrt(diag(F.inv))
  t.val <- sigma/se.theta
  pv <- 2 * pnorm(as.vector(abs(t.val)), lower.tail = FALSE)
  lim.inf <- sigma-qnorm(1-0.05/2)*se.theta 
  lim.sup <- sigma+qnorm(1-0.05/2)*se.theta 
  coefsigma <- cbind(sigma, se.theta, t.val, lim.inf,lim.sup,pv)
  
  colnames(coefsigma) = c("Variances", "std.error", "t.statistics", "lim.inf","lim.sup",
                          "p.value")
  
  ## Table for beta
  
  Q.inv <- solve(tX%*%V.inv%*%X)
  se.b <- sqrt(diag(Q.inv))
  t.val <- beta/se.b
  pv <- 2 * pnorm(as.vector(abs(t.val)), lower.tail = FALSE)
  lim.inf <- beta - qnorm(1-0.05/2)*se.b
  lim.sup <- beta + qnorm(1-0.05/2)*se.b
  coefbeta <- cbind(beta, se.b, t.val, lim.inf, lim.sup, pv)
  colnames(coefbeta) = c("beta", "std.error", "t.statistics", "lim.inf", "lim.sup",
                         "p.value")
  
  return(list(sigma = coefsigma, beta = coefbeta, Eblups = muE, Q.inv = Q.inv))
}







