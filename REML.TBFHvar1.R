###############################################################################
###############################################################################
###
###       REML Fisher-scoring fitting algorithm
###
###       Author: Esteban Cabello García
###
###       Work: TBFH-VAR1

library(Matrix)
  
  V2d_f <- function(V, tp, Phi) {
  phi_diag <- diag(Phi)
  n_phi <- length(phi_diag)
  cindx <- expand.grid(1:tp, 1:tp)
  V2d_elem <- vector("list", nrow(cindx))
  
  for (k in 1:nrow(cindx)) {
    diff_val <- cindx[k, 2] - cindx[k, 1]
    
    
    if (diff_val > 0) {
      V2d_elem[[k]] <- V %*% diag(phi_diag^diff_val, nrow = n_phi)
    } else if (diff_val < 0) {
      V2d_elem[[k]] <- diag(phi_diag^(-diff_val), nrow = n_phi) %*% V
    } else {
      V2d_elem[[k]] <- V
    }
  }
  
  V2d_row <- lapply(1:tp, function(i) do.call(cbind, V2d_elem[cindx[, 1] == i]))
  V2d <- do.call(rbind, V2d_row)
  return(V2d)
}


V2da_78f <- function(V, V2dt, tp, Phi, Phi78) {
  phi_diag <- diag(Phi)
  n_phi <- length(phi_diag)
  cindx <- expand.grid(1:tp, 1:tp)
  V2d_elem <- vector("list", nrow(cindx))
  
  for (k in 1:nrow(cindx)) {
    diff_val <- cindx[k, 2] - cindx[k, 1]
    
    if (diff_val > 0) {
      V2d_elem[[k]] <- V %*% diag(phi_diag^diff_val, nrow = n_phi) + 
        diff_val * V2dt %*% diag(phi_diag^(diff_val - 1), nrow = n_phi) %*% Phi78
    } else if (diff_val < 0) {
      ndiff <- -diff_val
      V2d_elem[[k]] <- diag(phi_diag^ndiff, nrow = n_phi) %*% V + 
        ndiff * diag(phi_diag^(ndiff - 1), nrow = n_phi) %*% Phi78 %*% V2dt
    } else {
      V2d_elem[[k]] <- V
    }
  }
  
  V2d_row <- lapply(1:tp, function(i) do.call(cbind, V2d_elem[cindx[, 1] == i]))
  V2d <- do.call(rbind, V2d_row)
  return(V2d)
}


REML.BTFHvar1 <- function(X, y, sigma.ini, Ve, PRECISION = 0.0001, MAXITER = 40) {
  
  kit <- 0
  FLAG <- 0
  diff <- rep(PRECISION + 1, 8)
  sigmau1 <- sigma.ini
  
  tp <- dim(Ve[[1]])[1] / 2
  D <- length(y) / (2 * tp)
  
  col1 <- matrix(1, nrow = tp)
  tcol1col1 <- tcrossprod(col1, col1)
  tX <- t(X)
  ty <- t(y)
  
  Phi7 <- diag(c(1, 0))
  Phi8 <- diag(c(0, 1))
  
  V.inv_sparse <- NULL
  betaREML <- NULL
  
  while (any(diff > rep(PRECISION, 8)) & (kit < MAXITER)) {
    kit <- kit + 1
    
    # Construction of the variance-covariance matrix V
    V1d <- matrix(c(sigmau1[1, ], sigmau1[3, ] * sqrt(sigmau1[1, ] * sigmau1[2, ]), 
                    sigmau1[3, ] * sqrt(sigmau1[1, ] * sigmau1[2, ]), sigmau1[2, ]), byrow = T, nrow = 2)
    
    Phi <- diag(c(sigmau1[7, ], sigmau1[8, ]))
    
    V2dt <- matrix(c(sigmau1[4, ] / (1 - sigmau1[7, ]^2), 
                     sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ]), 
                     sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ]), 
                     sigmau1[5, ] / (1 - sigmau1[8, ]^2)), nrow = 2, byrow = T)
    
    V2d <- V2d_f(V = V2dt, tp = tp, Phi = Phi)
    
    Vd <- lapply(Ve, function(x) kronecker(tcol1col1, V1d) + V2d + x)
    Vd.inv <- lapply(Vd, solve)
    
    V.inv_sparse <- as(bdiag(Vd.inv), "dgCMatrix")
    Q_inv <- solve(tX %*% V.inv_sparse %*% X)
    
    P_sparse <- V.inv_sparse - (V.inv_sparse %*% X %*% Q_inv %*% tX %*% V.inv_sparse)
    P <- as.matrix(P_sparse)
    
    # Construction of the derivatives of the variance components
    V2dt1 <- matrix(c(1, sigmau1[3, ] * sqrt(sigmau1[2, ]) / (2 * sqrt(sigmau1[1, ])), 
                      sigmau1[3, ] * sqrt(sigmau1[2, ]) / (2 * sqrt(sigmau1[1, ])), 0), ncol = 2, byrow = T)
    V2dt2 <- matrix(c(0, sigmau1[3, ] * sqrt(sigmau1[1, ]) / (2 * sqrt(sigmau1[2, ])), 
                      sigmau1[3, ] * sqrt(sigmau1[1, ]) / (2 * sqrt(sigmau1[2, ])), 1), ncol = 2, byrow = T)
    V2dt3 <- matrix(c(0, sqrt(sigmau1[1, ] * sigmau1[2, ]), sqrt(sigmau1[1, ] * sigmau1[2, ]), 0), ncol = 2, byrow = T)
    
    Vda1 <- kronecker(tcol1col1, V2dt1)
    Vda2 <- kronecker(tcol1col1, V2dt2)
    Vda3 <- kronecker(tcol1col1, V2dt3)
    
    V2dt4 <- matrix(c(1 / (1 - sigmau1[7, ]^2), 
                      sigmau1[6, ] * sqrt(sigmau1[5, ]) / (2 * sqrt(sigmau1[4, ]) * (1 - sigmau1[7, ] * sigmau1[8, ])), 
                      sigmau1[6, ] * sqrt(sigmau1[5, ]) / (2 * sqrt(sigmau1[4, ]) * (1 - sigmau1[7, ] * sigmau1[8, ])), 0), ncol = 2, byrow = T)
    V2dt5 <- matrix(c(0, 
                      sigmau1[6, ] * sqrt(sigmau1[4, ]) / (2 * sqrt(sigmau1[5, ]) * (1 - sigmau1[7, ] * sigmau1[8, ])), 
                      sigmau1[6, ] * sqrt(sigmau1[4, ]) / (2 * sqrt(sigmau1[5, ]) * (1 - sigmau1[7, ] * sigmau1[8, ])), 
                      1 / (1 - sigmau1[8, ]^2)), ncol = 2, byrow = T)
    V2dt6 <- matrix(c(0, sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ]), 
                      sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ]), 0), ncol = 2, byrow = T)
    
    Vda4 <- V2d_f(V = V2dt4, tp = tp, Phi = Phi)
    Vda5 <- V2d_f(V = V2dt5, tp = tp, Phi = Phi)
    Vda6 <- V2d_f(V = V2dt6, tp = tp, Phi = Phi)
    
    V2dt7 <- matrix(c(2 * sigmau1[7, ] * sigmau1[4, ] / (1 - sigmau1[7, ]^2)^2, 
                      sigmau1[8, ] * sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ])^2, 
                      sigmau1[8, ] * sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ])^2, 0), ncol = 2, byrow = T)
    V2dt8 <- matrix(c(0, 
                      sigmau1[7, ] * sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ])^2, 
                      sigmau1[7, ] * sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ])^2, 
                      2 * sigmau1[8, ] * sigmau1[5, ] / (1 - sigmau1[8, ]^2)^2), ncol = 2, byrow = T)
    
    Vda7 <- V2da_78f(V = V2dt7, V2dt = V2dt, tp = tp, Phi = Phi, Phi78 = Phi7)
    Vda8 <- V2da_78f(V = V2dt8, V2dt = V2dt, tp = tp, Phi = Phi, Phi78 = Phi8)
    
    Vda <- list(Vda1, Vda2, Vda3, Vda4, Vda5, Vda6, Vda7, Vda8)
    Va <- lapply(Vda, function(x) as(Matrix::bdiag(rep(list(x), D)), "dgCMatrix"))
    
    yP <- ty %*% P
    Py <- P %*% y
    PVa <- lapply(Va, function(x) as.matrix(P %*% x))
    
    S <- matrix(sapply(1:8, function(idx) {
      -0.5 * sum(diag(PVa[[idx]])) + as.vector(0.5 * yP %*% Va[[idx]] %*% Py)
    }), ncol = 1)
    
    u <- unlist(lapply(1:8, function(j) 1:j))
    PVa.a <- PVa[rep(1:8, times = 1:8)]
    PVa.b <- PVa[u]
    
    #Property: tr(AB) = sum(A * t(B)), with * the Hadamard product.
    preF0 <- sapply(1:length(PVa.a), function(idx) {
      0.5 * sum(PVa.a[[idx]] * t(PVa.b[[idx]])) 
    })
    
    F0 <- matrix(0, ncol = 8, nrow = 8, byrow = T) 
    F0[upper.tri(F0, diag = T)] <- preF0
    F0[lower.tri(F0, diag = T)] <- t(F0)[lower.tri(t(F0), diag = T)]
    
    F.inv <- try(solve(F0), silent = TRUE)
    
    if (inherits(F.inv, "try-error") | kit == MAXITER) {
      stop("Not inversible matrix or not convergence")
    } else {
      
      sigmau <- sigmau1 + F.inv %*% S
      
      if (any(sigmau[c(1, 2, 4, 5), ] < 0) | any(abs(sigmau[c(3, 6, 7, 8), ]) > 1)) {
        warning("Out of parametric space")
        break
      }
      
      diff <- abs(sigmau - sigmau1)
      sigmau1 <- sigmau
      
      V1d_new <- matrix(c(sigmau1[1, ], sigmau1[3, ] * sqrt(sigmau1[1, ] * sigmau1[2, ]), 
                          sigmau1[3, ] * sqrt(sigmau1[1, ] * sigmau1[2, ]), sigmau1[2, ]), byrow = T, nrow = 2)
      Phi_new <- diag(c(sigmau1[7, ], sigmau1[8, ]))
      V2dt_new <- matrix(c(sigmau1[4, ] / (1 - sigmau1[7, ]^2), 
                           sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ]), 
                           sigmau1[6, ] * sqrt(sigmau1[4, ] * sigmau1[5, ]) / (1 - sigmau1[7, ] * sigmau1[8, ]), 
                           sigmau1[5, ] / (1 - sigmau1[8, ]^2)), nrow = 2, byrow = T)
      
      V2d_new <- V2d_f(V = V2dt_new, tp = tp, Phi = Phi_new)
      Vd_new <- lapply(Ve, function(x) kronecker(tcol1col1, V1d_new) + V2d_new + x)
      
      Vd.inv_new <- lapply(Vd_new, solve)
      V.inv_sparse <- as(bdiag(Vd.inv_new), "dgCMatrix")
      
      betaREML <- as.matrix(solve(tX %*% V.inv_sparse %*% X) %*% tX %*% V.inv_sparse %*% y)
    }
  }
  
  return(list(Sigma = sigmau1, Beta = betaREML, V.inv = V.inv_sparse, Fisher.inv = F.inv, FLAG = FLAG, ITER = kit))
}
  



