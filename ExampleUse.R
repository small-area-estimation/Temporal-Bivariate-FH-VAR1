###############################################################################
###############################################################################
###
###       Example of use of the TBFH-var1 model with data
###       
###       Author: Esteban Cabello García
###
###       Work: TBFH-var1 - Aplicacion a datos reales


library(saery)

source("REML.TBFHvar1.R")
source("fit.TBFHvar1.R")
source("mse.boot.TBFHvar1.R")
source("mse.ana.TBFHvar1.R")

### Reading the data
df <- readRDS("example_data.rds")

### Initial values for the REML algorithm. Univariate fits with AR(1) structure.
D <- length(unique(df$domain)) * length(unique(df$group))
tp <- length(unique(df$time))
md <- rep(length(unique(df$time)),D)

X1 <- as.matrix(df$X1)
X2 <- as.matrix(df$X2)

marginalY1 <- fit.saery.AR1(X = X1, ydi = df$Y1, D = D, md = md, sigma2edi = df$varY1, conf.level = 0.95)
marginalY2 <- fit.saery.AR1(X = X2, ydi = df$Y2, D = D, md = md, sigma2edi = df$varY2, conf.level = 0.95)

theta_1 <- marginalY1$SIGMA[1:2,1]
theta_2 <- marginalY2$SIGMA[1:2,1]
phi_1 <- marginalY1$SIGMA[3,1]
phi_2 <- marginalY2$SIGMA[3,1]

sigma.initial <- matrix(c(theta_1[1], theta_2[1], 0, 
                          theta_1[2], theta_2[2], 0, 
                          phi_1, phi_2), ncol = 1)


### Construction of the matrices Ve and Ved
varcov_row <- matrix(apply(df,1, function(x) x[c(6,8,8,7)]),ncol = 4,byrow = T)
Vedt <- rep(list(list()),D*tp)
for(i in 1:(D*tp)){
  Vedt[[i]] <- matrix(varcov_row[i,],ncol = 2,byrow = T)
}
indx <- split_vector(1:(D*tp),steps =D)

Ved <- list()
for(i in 1:length(indx)){
  Ved[[i]] <- bdiag(Vedt[indx[[i]]])
}

Ve <- bdiag(Vedt)

#### Preparation for REML algorithm. Construct X (explanatory variables) and y (target vector) matrices.

formula <- list(f1 = Y1 ~ X1,
                f2 = Y2 ~ X2)


X1 <- model.matrix(formula[[1]], df)
X2 <- model.matrix(formula[[2]], df)

x.matrix <- do.call(rbind, lapply(1:nrow(X1), function(i) {
  as.matrix(Matrix::bdiag(X1[i, , drop = FALSE], X2[i, , drop = FALSE]))
}))

y <- matrix(unlist(t(df[,c("Y1","Y2")])), ncol = 1, nrow = 2*D*tp)


############ REML algorithm ############ 
model <- REML.BTFHvar1(X = x.matrix, y, sigma.ini = sigma.initial, 
                            Ve = Ved, PRECISION = 0.001, MAXITER = 20)



beta.hat <- model$Beta
sigma.hat <- model$Sigma
Fish.inv <- model$Fisher.inv
V.inv <- model$V.inv

fitted <- fit.BTFHvar1(x.matrix, y, D, tp, sigma.hat, beta.hat, V.inv, Fish.inv, conf.level = 0.95)

fitted$sigma
fitted$beta
head(fitted$Eblups)

### MSE estimation
mse.ana <- mse.ana.BTFHvar1(X = x.matrix, Sigma = sigma.hat, F.inv = Fish.inv, Ved)
mse.ana.y1 <- mse.ana$mse[,1]
mse.ana.y2 <- mse.ana$mse[,2]

### The procedure can be slow. Be patient :)
mse.boot <- mse.boot.BTFHvar1(X = x.matrix, D, tp, beta.hat, sigma.hat, Ved, Bsize = c(10,50))
mse.boot.eblup <- mse.boot$mse.eblup.boot
mse.boot.dir <- mse.boot$mse.dir.boot


