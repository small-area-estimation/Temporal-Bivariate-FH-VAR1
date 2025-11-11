###############################################################################
###############################################################################                                               
###       Parametric Bootstrap estimator of the MSE matrix 
###       
###       Author: Esteban Cabello García
###
###       Work: TBFH-var1 


mse.boot.BTFHvar1 <- function(X, D, tp, beta, sigma, Ve, Bsize){
  
  ### X: design matrix
  ### D: number of domains. Integer.
  ### tp: number of time periods. Integer.
  ### beta: Regression parameters. Vector.
  ### sigma: Sigma estimates. Vector
  ### Ve: Variances-covariances of the sampling errors. List.
  ### Bsize: bootstrap iterations. Vector.
  
  B.last <- tail(Bsize, 1)
  
  sigma.hat <- sigma
  beta.hat <- beta
  
  sigma1.hat <- sigma.hat[1]
  sigma2.hat <- sigma.hat[2]
  rho1.hat <- sigma.hat[3]
  sigma3.hat <- sigma.hat[4]
  sigma4.hat <- sigma.hat[5]
  rho2.hat <- sigma.hat[6]
  phi1.hat <- sigma.hat[7]
  phi2.hat <- sigma.hat[8]
  
  
  ##### Creating V_{1,d} matrix
  V1d.hat <- matrix(c(sigma1.hat,
                  rho1.hat*sqrt(sigma1.hat*sigma2.hat),
                  rho1.hat*sqrt(sigma1.hat*sigma2.hat),
                  sigma2.hat), byrow = T, nrow = 2)
  
  ##### Creating V_{2,d} matrix
  Phi.hat <- diag(c(phi1.hat,phi2.hat))
  
  V2dt.hat <- matrix(c(sigma3.hat/(1-phi1.hat^2), 
                   rho2.hat*sqrt(sigma3.hat*sigma4.hat)/(1-phi1.hat*phi2.hat), 
                   rho2.hat*sqrt(sigma3.hat*sigma4.hat)/(1-phi1.hat*phi2.hat), 
                   sigma4.hat/(1-phi2.hat^2) ), nrow = 2, byrow = T)
  
  V2d.hat <- V2d_f(V = V2dt.hat,tp = tp, Phi = Phi.hat)
  
  
  Z1 <-  Reduce(bdiag,rep(list(Reduce(rbind,rep(list(diag(2)),tp))),D))
  Z1t <- t(Z1)
  Z2d <- Reduce(bdiag,rep(list(diag(c(1,1))),tp))
  Z2 <-  diag(rep(1,2*D*tp))
  Z2t <- t(Z2)
  
  
  col1 <- rep(1,tp)
  tcol1col1 <- tcrossprod(col1,col1)
  tX <- t(X)
  
  Vedt <- Ve[[1]]
  
  mu.ast <- mu.ast.hat <- eblup.ast.hat <- matrix(0, ncol = B.last, nrow = 2*D*tp)
  mse.eblup.boot <- mse.dir.boot <- as.data.frame(matrix(0, ncol = length(Bsize), nrow = 2*D*tp))
  dif.mu.eblup <- dif.mu.dir <- 0
  
  b <- 0
  BadTot_b <- BadTot_2 <- 0
  

  while(b < B.last){
    b <- b+1
    
    u1d.ast <- mvrnorm(D, mu= rep(0,2), Sigma = V1d.hat) ### Domain random-effects
    u2d.ast <- mvrnorm(D, mu = rep(0,2*tp), Sigma = V2d.hat) ### Subdomain random-effects
    edt.ast <- mvrnorm(D, mu = rep(0,2*tp), Sigma = Vedt) ### Sampling errors
    
    
    #################
    #### Generation of target variables y_{dt}
    ################
    
    u1.ast <- matrix(t(u1d.ast), ncol = 1, byrow = T)
    u2.ast <- matrix(t(u2d.ast), ncol = 1, byrow = T)
    e.ast <- matrix(t(edt.ast), ncol = 1, byrow = T)
    
    mu.ast[,b] <- as.matrix(X%*%beta.hat+ Z1%*%u1.ast + Z2%*%u2.ast)
    y.ast <- mu.ast[,b] + e.ast

    
    
    fit.boot <- try(REML.BTFHvar1(X, y.ast, sigma.ini = sigma.hat, Ve = Ve, PRECISION = 10^-3, MAXITER = 60), TRUE)
    
    
    if(inherits(fit.boot,"try-error")){
      BadTot_2 <- BadTot_2 + 1
      cat("\t Bootstrap_sample", b, " rejected by try-error\n")
      b <- b-1
    }
    else {
      cat("Fitting Bootstrap sample number", b, "\n")
      

      sigma.ast.hat <- fit.boot$Sigma
      beta.ast.hat <- fit.boot$Beta
      V.inv.ast.hat <- fit.boot$V.inv
      F.inv.ast.hat <- fit.boot$Fisher.inv
      
      eblup.ast.hat[,b] <- as.numeric(fit.BTFHvar1(X, y.ast, D, tp, sigma = sigma.ast.hat, beta = beta.ast.hat, V.inv = V.inv.ast.hat, F.inv = F.inv.ast.hat)$Eblups)
        

      dif.mu.eblup <- dif.mu.eblup + (eblup.ast.hat[,b] - mu.ast[,b])^2
      dif.mu.dir <- dif.mu.dir + (y.ast - mu.ast[,b])^2
      
    }
    
    if(any(Bsize == b)){
      
      indx.B <- which(Bsize == b)
      mse.eblup.boot[,indx.B] <- dif.mu.eblup/b
      mse.dir.boot[,indx.B] <- dif.mu.dir/b
    
    }
    
  }
  
  colnames(mse1.boot) <- colnames(mean.deltaEB.square) <- colnames(mean.g1g2.boot) <- paste0("B=",as.character(Bsize))
  

  return(list(mse.eblup.boot = mse.eblup.boot, mse.dir.boot = mse.dir.boot, BadTot_2))
  
}
