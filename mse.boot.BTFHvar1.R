###############################################################################
###############################################################################
###
###       Fecha de creación: 04/03/2025
###                                                  
###       Ultima modificación: 05/03/2025
###       
###       Autor: Esteban Cabello García, Agustín Pérez, Lola Esteban
###
###       Trabajo: BTFH-var1 - Bootstrap


mse.BTFHvar1.boot <- function(X, D, tp, beta, sigma, Ve, Bsize, k){
  
  ### Bsize: vector con nº iteraciones bootstrap
  
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
  
  V2d.hat <- V2d_f_opt(V = V2dt.hat,tp = tp, Phi = Phi.hat)
  
  
  Z1 <-  Reduce(bdiag,rep(list(Reduce(rbind,rep(list(diag(2)),tp))),D))
  Z1t <- t(Z1)
  Z2d <- Reduce(bdiag,rep(list(diag(c(1,1))),tp))
  Z2 <-  diag(rep(1,2*D*tp))
  Z2t <- t(Z2)
  
  
  col1 <- rep(1,tp)
  tcol1col1 <- tcrossprod(col1,col1)
  tX <- t(X)
  
  Vedt <- Ve[[1]]
  
  mu.ast <- mu.ast.hat <- eblup.ast.hat <- blup.ast.hat <- matrix(0, ncol = B.last, nrow = 2*D*tp)
  mse1.boot <- mean.deltaEB.square <- mean.g1g2.boot <- as.data.frame(matrix(0, ncol = length(Bsize), nrow = 2*D*tp))
  dif.mu.blup <- dif.mu.eblup <- dif.mu.eblup.blup <- sum.g1g2.boot <- 0
  
  b <- 0
  BadTot_b <- BadTot_2 <- 0
  

  while(b < B.last){
    b <- b+1
    
    u1d.ast <- mvrnorm(D, mu= rep(0,2), Sigma = V1d.hat) ### Domain random-effects
    u2d.ast <- mvrnorm(D, mu = rep(0,2*tp), Sigma = V2d.hat) ### Subdomain random-effects
    edt.ast <- mvrnorm(D, mu = rep(0,2*tp), Sigma = Vedt) ### Sampling errors
    
    
    #################
    #################
    #### Generation of target variables y_{dt}
    ################
    ################
    
    u1.ast <- matrix(t(u1d.ast), ncol = 1, byrow = T)
    u2.ast <- matrix(t(u2d.ast), ncol = 1, byrow = T)
    e.ast <- matrix(t(edt.ast), ncol = 1, byrow = T)
    
    mu.ast[,b] <- as.matrix(X%*%beta.hat+ Z1%*%u1.ast + Z2%*%u2.ast)
    y.ast <- mu.ast[,b] + e.ast
    ty.ast <- t(y.ast)
    
    
    fit.boot <- try(REML.BTFHvar1_opt(X, y.ast, sigma.ini = sigma.hat, Ve = Ve, PRECISION = 10^-3, MAXITER = 60), TRUE)
    
    
    if(inherits(fit.boot,"try-error")){
      BadTot_2 <- BadTot_2 + 1
      # excepcion <- c(excepcion, b)
      write.table(data.frame(class(fit.boot), D, k, b), file=paste0("Results/WARNING_BOOTSTRAP B=", B.last, ".txt"), append=TRUE, col.names=FALSE, row.names=FALSE)
      cat("\t Muestra_bootstrap", b, " rechazada por try-error\n")
      b <- b-1
    }
    else {
      cat("Iteration", k, ".Fitting Bootstrap sample number", b, "\n")
      

      sigma.ast.hat <- fit.boot$Sigma
      beta.ast.hat <- fit.boot$Beta
      V.inv.ast.hat <- fit.boot$V.inv
      F.inv.ast.hat <- fit.boot$Fisher.inv
      
      eblup.ast.hat[,b] <- fit.BTFHvar1(X, y.ast, D, tp, sigma = sigma.ast.hat, beta = beta.ast.hat, V.inv = V.inv.ast.hat, F.inv = F.inv.ast.hat)$Eblups
      blup.ast.hat[,b] <- fit.blup.BTFHvar1(X, y.ast, D, tp, sigma = sigma.hat, beta = beta.hat, Ve = Ve)$Blups
        
      mse.ana.boot <- mse.ana.BTFHvar1(X, sigma.ast.hat, F.inv.ast.hat, Ved)
      
      sum.g1g2.boot <- sum.g1g2.boot + (mse.ana.boot$G1 + mse.ana.boot$G2)

      dif.mu.blup <- dif.mu.blup + (blup.ast.hat[,b] - mu.ast[,b])^2
      dif.mu.eblup <- dif.mu.eblup + (eblup.ast.hat[,b] - mu.ast[,b])^2
      dif.mu.eblup.blup <- dif.mu.eblup.blup + (eblup.ast.hat[,b] - blup.ast.hat[,b])^2
      
    }
    
    if(any(Bsize == b)){
      
      indx.B <- which(Bsize == b)
      mse1.boot[,indx.B] <- dif.mu.eblup/b
      mean.deltaEB.square[,indx.B] <- dif.mu.eblup.blup/b
      mean.g1g2.boot <- sum.g1g2.boot/b
      

    }
    
  }
  
  colnames(mse1.boot) <- colnames(mean.deltaEB.square) <- colnames(mean.g1g2.boot) <- paste0("B=",as.character(Bsize))
  
  #cat("Fitting of Bootstrap sample", b, "in iteration i =",i, " finished \n")
  return(list(mse1.boot = mse1.boot, deltaEB.2 = mean.deltaEB.square, g1g2.boot = mean.g1g2.boot, BadTot_2))
  
}
