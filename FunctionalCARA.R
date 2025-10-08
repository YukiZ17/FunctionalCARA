library(rARPACK);library(MFPCA);library(tidyverse);library(MASS);require(zipfR);require(fOptions);require(Matrix)
library(pec);library(fda);library(survAUC);library(Triangulation);library(BPST);library(survival)

source("ufun.R")


############# Imaging Data ###########
V.est = readRDS('V.est') # vertices for triangles for simulation setting
Tr.est = readRDS('Tr.est') # triangles
Z = readRDS('Z') # expanded grid points

############# Hyperparameter Setting ###############
d.est = 2; r = 1
theta = 0.01 # an example for a pre-specified theta (needs to be tuned)
npc = 2 # number of basis 
ncr = 40 # dimension of the simulation images
lambda = c(0, 1, 10, 10^2, 10^3, 10^6) # tunning parameter for bernstein
#ind.inside=inVT(V.est, Tr.est, Z[,1], Z[,2])$ind.inside


############### Optimal Design ####################

## Simulation parameter

N <- 500
gamma.cara.1.N <- c()
alpha.cara.1.N <- c()
T.cara.1.N <- c()
Y.cara.1.N <- c()
P.cara.1.N <- c()

gamma.cara.2.N <- c()
alpha.cara.2.N <- c()
T.cara.2.N <- c()
Y.cara.2.N <- c()
P.cara.2.N <- c()

gamma.cara.3.N <- c()
alpha.cara.3.N <- c()
T.cara.3.N <- c()
Y.cara.3.N <- c()
P.cara.3.N <- c()

gamma.eq.N <- c()
alpha.eq.N <- c()
T.eq.N <- c()
Y.eq.N <- c()
P.eq.N <- c()

## Parameter inputs
n.sample <- 500
n.burnin <- 100
n.batch <- 40
batch <- (n.sample-n.burnin)/n.batch
alpha = 0.25
gamma = 0


for (j in 1:N){
  
  ## Result Storage
  T.all.cara.1 <- c()
  Y.all.cara.1 <- c()
  Alpha.batch.cara.1 <- c()
  b.1.batch.cara.1 <- c()
  b.2.batch.cara.1 <- c()
  Gamma.batch.cara.1 <- c()
  
  T.all.cara.2 <- c()
  Y.all.cara.2 <- c()
  Alpha.batch.cara.2 <- c()
  b.1.batch.cara.2 <- c()
  b.2.batch.cara.2 <- c()
  Gamma.batch.cara.2 <- c()
  
  T.all.cara.3 <- c()
  Y.all.cara.3 <- c()
  Alpha.batch.cara.3 <- c()
  b.1.batch.cara.3 <- c()
  b.2.batch.cara.3 <- c()
  Gamma.batch.cara.3 <- c()
  
  
  ##### Generating data
  a <- train_dat.id$Z1[1,];na <- which(is.na(a))
  
  basis <- eFun(argvals = seq(0, 1, length.out=40), M = 2,  type = "Poly")@X
  
  phi.1 <- basis[1,] %*% t(basis[2,]); dim(phi.1) <- c(1,1600)
  phi.2 <- basis[2,] %*% t(basis[1,]); dim(phi.2) <- c(1,1600)
  phi.3 <- basis[2,] %*% t(basis[2,]); dim(phi.3) <- c(1,1600)
  
  w1 <- matrix(rnorm(n.sample,0,sqrt(10)),n.sample,1)
  w2 <- matrix(rnorm(n.sample,0,sqrt(8)),n.sample,1)
  w3 <- matrix(rnorm(n.sample,0,sqrt(4)),n.sample,1)
  
  X<- w1%*%phi.1+w2%*%phi.2+w3%*%phi.3 
  X[,na] <- NA
  X.inside <- X
  X.inside[,na] <- 0
  
  
  ## Burnin -- Equal:
  
  t.burnin <- rbinom(n.burnin,size=1,prob = 0.5)
  
  T.all.cara.1 <- c(T.all.cara.1,t.burnin)
  T.all.cara.2 <- c(T.all.cara.2,t.burnin)
  T.all.cara.3 <- c(T.all.cara.3,t.burnin)
  
  
  pi.burnin <- rep(alpha,n.burnin)+ gamma*t.burnin +(X.inside[1:n.burnin,]%*%t(phi.3)/921)+(X.inside[1:n.burnin,]%*%t(phi.2)/921)+(X.inside[1:n.burnin,]%*%t(phi.1)/921) #divided by 1600? 921?
  Y.outcome.burnin<- rbinom(n.burnin,1,(1/(1+exp(-pi.burnin))))
  
  Y.all.cara.1 <- c(Y.all.cara.1,Y.outcome.burnin)
  Y.all.cara.2 <- c(Y.all.cara.2,Y.outcome.burnin)
  Y.all.cara.3 <- c(Y.all.cara.3,Y.outcome.burnin)
  
  
  sup_fun.binary=fpca_img.binary(type = "bernstein", Y = X.inside[1:n.sample,], n_length = n.sample,
                                   theta = NA , lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                   r = r, Z = Z, ncr = ncr)
  
  
  o_tfs.bin = as.matrix(sup_fun.binary[[1]])
  
  
  b.basis.bin = as.matrix(sup_fun.binary[[3]]); b.scores.bin = as.matrix(sup_fun.binary[[4]])
  
  score_sup.bin = t(o_tfs.bin) %*% (b.basis.bin%*%t(b.scores.bin))
  
  
  if(npc == 2){
    score_sup.bin = t(score_sup.bin)
  }else{
    score_sup.bin = cbind(score_sup.bin[1,])
  }
  

  id <- 1:n.burnin
  
  train <- data.frame(id)
  train$y <- Y.outcome.burnin
  train$score1 <- score_sup.bin[1:n.burnin,1]
  train$score2 <- score_sup.bin[1:n.burnin,2]
  train$t <- t.burnin
  
  
  glm_fit_burnin = glm(y~score1+score2+t,data=train,family="binomial")
  Alpha.batch.cara.1 <- c(Alpha.batch.cara.1,glm_fit_burnin$coefficients[1])
  b.1.batch.cara.1 <- c(b.1.batch.cara.1,glm_fit_burnin$coefficients[2])
  b.2.batch.cara.1 <- c(b.2.batch.cara.1,glm_fit_burnin$coefficients[3])
  Gamma.batch.cara.1 <- c(Gamma.batch.cara.1,glm_fit_burnin$coefficients[4])
  
  Alpha.batch.cara.2 <- c(Alpha.batch.cara.2,glm_fit_burnin$coefficients[1])
  b.1.batch.cara.2 <- c(b.1.batch.cara.2,glm_fit_burnin$coefficients[2])
  b.2.batch.cara.2 <- c(b.2.batch.cara.2,glm_fit_burnin$coefficients[3])
  Gamma.batch.cara.2 <- c(Gamma.batch.cara.2,glm_fit_burnin$coefficients[4])
  
  Alpha.batch.cara.3 <- c(Alpha.batch.cara.3,glm_fit_burnin$coefficients[1])
  b.1.batch.cara.3 <- c(b.1.batch.cara.3,glm_fit_burnin$coefficients[2])
  b.2.batch.cara.3 <- c(b.2.batch.cara.3,glm_fit_burnin$coefficients[3])
  Gamma.batch.cara.3 <- c(Gamma.batch.cara.3,glm_fit_burnin$coefficients[4])
  
  
  ################ Equal ############### 
  
  t.equal <- rbinom(n.sample,size=1,prob = 0.5)
  
  pi.equal <- rep(alpha,n.sample)+ gamma*t.equal +(X.inside%*%t(phi.3)/921)+(X.inside%*%t(phi.2)/921) +(X.inside%*%t(phi.1)/921)
  
  Y.outcome.equal<- rbinom(n.sample,1,(1/(1+exp(-pi.equal))))
  
  id <- 1:n.sample
  
  train <- data.frame(id)
  train$y <- Y.outcome.equal
  train$score1 <- score_sup.bin[,1]
  train$score2 <- score_sup.bin[,2]
  train$t <- t.equal
  
  
  glm_fit_equal = glm(y~score1+score2+t,data=train,family="binomial")
  
  gamma.eq.N <- c(gamma.eq.N,glm_fit_equal$coefficients[4])
  alpha.eq.N <- c(alpha.eq.N,glm_fit_equal$coefficients[1])
  Y.eq.N <- c(Y.eq.N,sum(Y.outcome.equal))
  T.eq.N <- c(T.eq.N,sum(t.equal))
  wald.coef <- summary(glm_fit_equal)$coefficients 
  p.eq <- wald.coef["t","Pr(>|z|)"]
  P.eq.N <- c(P.eq.N,p.eq)
  
  
  ################## Adaptive CARA 3 Types ####################
  
  ###### CARA 1: Classical #########
  
  for (k in 1:batch){
    
    
    l.k <- n.burnin+1+(k-1)*n.batch
    u.k <- n.burnin+k*n.batch
    X.k <- X.inside[l.k:u.k,]
    
    t.k <- rbinom(n.batch,size=1,prob =1/(1+exp(-Gamma.batch.cara.1[k])))
    
    T.all.cara.1 <- c(T.all.cara.1,t.k)
    
    pi.k <- rep(alpha,n.batch)+ gamma*t.k +(X.k%*%t(phi.3)/921)+(X.k%*%t(phi.2)/921)+(X.k%*%t(phi.1)/921) 
    Y.outcome.k<- rbinom(n.batch,1,(1/(1+exp(-pi.k))))
    Y.all.cara.1 <- c(Y.all.cara.1,Y.outcome.k)
    
    
    id <- 1:u.k
    
    train <- data.frame(id)
    train$y <- Y.all.cara.1
    train$score1 <- score_sup.bin[id,1]
    train$score2 <- score_sup.bin[id,2]
    train$t <- T.all.cara.1
    
    
    
    glm_fit = glm(y~score1+score2+t,data=train,family="binomial")
    Alpha.batch.cara.1 <- c(Alpha.batch.cara.1,glm_fit$coefficients[1])
    b.1.batch.cara.1 <- c(b.1.batch.cara.1,glm_fit$coefficients[2])
    b.2.batch.cara.1 <- c(b.2.batch.cara.1,glm_fit$coefficients[3])
    Gamma.batch.cara.1 <- c(Gamma.batch.cara.1,glm_fit$coefficients[4])
    
 
    
  }
  
  gamma.cara.1.N <- c(gamma.cara.1.N,Gamma.batch.cara.1[length(Gamma.batch.cara.1)])
  alpha.cara.1.N <- c(alpha.cara.1.N,Alpha.batch.cara.1[length(Alpha.batch.cara.1)])
  T.cara.1.N <- c(T.cara.1.N,sum(T.all.cara.1))
  Y.cara.1.N <- c(Y.cara.1.N,sum(Y.all.cara.1))
  wald.coef <- summary(glm_fit)$coefficients 
  p.cara.1 <- wald.coef["t","Pr(>|z|)"]
  P.cara.1.N <- c(P.cara.1.N,p.cara.1)
  
  
  
  ###### CARA 2: Optimal #########  
  
  for (k in 1:batch){
    
    
    l.k <- n.burnin+1+(k-1)*n.batch
    u.k <- n.burnin+k*n.batch
    X.k <- X.inside[l.k:u.k,]
    
    logit.p.a.k.hat <- Alpha.batch.cara.2[k]+b.1.batch.cara.2[k]*score_sup.bin[l.k:u.k,1]+b.2.batch.cara.2[k]*score_sup.bin[l.k:u.k,2]+Gamma.batch.cara.2[k]
    logit.p.b.k.hat <- Alpha.batch.cara.2[k]+b.1.batch.cara.2[k]*score_sup.bin[l.k:u.k,1]+b.2.batch.cara.2[k]*score_sup.bin[l.k:u.k,2]
    p.a.k.hat <- 1/(1+exp(-logit.p.a.k.hat))
    p.b.k.hat <- 1/(1+exp(-logit.p.b.k.hat))
    
    t.k <- rbinom(n.batch,size=1,prob = sqrt(p.a.k.hat)/(sqrt(p.a.k.hat)+sqrt(p.b.k.hat)))
    T.all.cara.2 <- c(T.all.cara.2,t.k)
    
    pi.k <- rep(alpha,n.batch)+ gamma*t.k +(X.k%*%t(phi.3)/921)+(X.k%*%t(phi.2)/921)+(X.k%*%t(phi.1)/921) 
    Y.outcome.k<- rbinom(n.batch,1,(1/(1+exp(-pi.k))))
    Y.all.cara.2 <- c(Y.all.cara.2,Y.outcome.k)
    
    
    id <- 1:u.k
    
    train <- data.frame(id)
    train$y <- Y.all.cara.2
    train$score1 <- score_sup.bin[id,1]
    train$score2 <- score_sup.bin[id,2]
    train$t <- T.all.cara.2
    
    
    
    glm_fit = glm(y~score1+score2+t,data=train,family="binomial")
    Alpha.batch.cara.2 <- c(Alpha.batch.cara.2,glm_fit$coefficients[1])
    b.1.batch.cara.2 <- c(b.1.batch.cara.2,glm_fit$coefficients[2])
    b.2.batch.cara.2 <- c(b.2.batch.cara.2,glm_fit$coefficients[3])
    Gamma.batch.cara.2 <- c(Gamma.batch.cara.2,glm_fit$coefficients[4])
    
 
    
  }
  
  gamma.cara.2.N <- c(gamma.cara.2.N,Gamma.batch.cara.2[length(Gamma.batch.cara.2)])
  alpha.cara.2.N <- c(alpha.cara.2.N,Alpha.batch.cara.2[length(Alpha.batch.cara.2)])
  T.cara.2.N <- c(T.cara.2.N,sum(T.all.cara.2))
  Y.cara.2.N <- c(Y.cara.2.N,sum(Y.all.cara.2))
  wald.coef <- summary(glm_fit)$coefficients 
  p.cara.2 <- wald.coef["t","Pr(>|z|)"]
  P.cara.2.N <- c(P.cara.2.N,p.cara.2)
  
  
  
  
  ###### CARA 3: Neyman #########  
  
  for (k in 1:batch){
    
    
    l.k <- n.burnin+1+(k-1)*n.batch
    u.k <- n.burnin+k*n.batch
    X.k <- X.inside[l.k:u.k,]
    
    logit.p.a.k.hat <- Alpha.batch.cara.3[k]+b.1.batch.cara.3[k]*score_sup.bin[l.k:u.k,1]+b.2.batch.cara.3[k]*score_sup.bin[l.k:u.k,2]+Gamma.batch.cara.3[k]
    logit.p.b.k.hat <- Alpha.batch.cara.3[k]+b.1.batch.cara.3[k]*score_sup.bin[l.k:u.k,1]+b.2.batch.cara.3[k]*score_sup.bin[l.k:u.k,2]
    p.a.k.hat <- 1/(1+exp(-logit.p.a.k.hat))
    p.b.k.hat <- 1/(1+exp(-logit.p.b.k.hat))
    q.a.k.hat <- 1-p.a.k.hat
    q.b.k.hat <- 1-p.b.k.hat
    
    
    t.k <- rbinom(n.batch,size=1,prob = sqrt(p.b.k.hat*q.b.k.hat)/(sqrt(p.a.k.hat*q.a.k.hat)+sqrt(p.b.k.hat*q.b.k.hat)))
    T.all.cara.3 <- c(T.all.cara.3,t.k)
    
    pi.k <- rep(alpha,n.batch)+ gamma*t.k +(X.k%*%t(phi.3)/921)+(X.k%*%t(phi.2)/921)+(X.k%*%t(phi.1)/921) #divided by 1600? 921?
    Y.outcome.k<- rbinom(n.batch,1,(1/(1+exp(-pi.k))))
    Y.all.cara.3 <- c(Y.all.cara.3,Y.outcome.k)
    
    
    id <- 1:u.k
    
    train <- data.frame(id)
    train$y <- Y.all.cara.3
    train$score1 <- score_sup.bin[id,1]
    train$score2 <- score_sup.bin[id,2]
    train$t <- T.all.cara.3
    
    
    
    glm_fit = glm(y~score1+score2+t,data=train,family="binomial")
    Alpha.batch.cara.3 <- c(Alpha.batch.cara.3,glm_fit$coefficients[1])
    b.1.batch.cara.3 <- c(b.1.batch.cara.3,glm_fit$coefficients[2])
    b.2.batch.cara.3 <- c(b.2.batch.cara.3,glm_fit$coefficients[3])
    Gamma.batch.cara.3 <- c(Gamma.batch.cara.3,glm_fit$coefficients[4])
    
    #print(k)
    
  }
  
  gamma.cara.3.N <- c(gamma.cara.3.N,Gamma.batch.cara.3[length(Gamma.batch.cara.3)])
  alpha.cara.3.N <- c(alpha.cara.3.N,Alpha.batch.cara.3[length(Alpha.batch.cara.3)])
  T.cara.3.N <- c(T.cara.3.N,sum(T.all.cara.3))
  Y.cara.3.N <- c(Y.cara.3.N,sum(Y.all.cara.3))
  wald.coef <- summary(glm_fit)$coefficients 
  p.cara.3 <- wald.coef["t","Pr(>|z|)"]
  P.cara.3.N <- c(P.cara.3.N,p.cara.3)
  
  
  
  print(j)  
  
}



### Summary
mean(gamma.cara.1.N)
var(gamma.cara.1.N)
mean(alpha.cara.1.N)
var(alpha.cara.1.N)
sum(P.cara.1.N<0.05)
sum(P.cara.1.N<0.10)

mean(T.cara.1.N)
sd(T.cara.1.N)
mean(Y.cara.1.N)
sd(Y.cara.1.N)



mean(gamma.cara.2.N)
var(gamma.cara.2.N)
mean(alpha.cara.2.N)
var(alpha.cara.2.N)
sum(P.cara.2.N<0.05)
sum(P.cara.2.N<0.10)

mean(T.cara.2.N)
sd(T.cara.2.N)
mean(Y.cara.2.N)
sd(Y.cara.2.N)



mean(gamma.cara.3.N)
var(gamma.cara.3.N)
mean(alpha.cara.3.N)
var(alpha.cara.3.N)
sum(P.cara.3.N<0.05)
sum(P.cara.3.N<0.10)

mean(T.cara.3.N)
sd(T.cara.3.N)
mean(Y.cara.3.N)
sd(Y.cara.3.N)




mean(gamma.eq.N)
var(gamma.eq.N)
mean(alpha.eq.N)
var(alpha.eq.N)
sum(P.eq.N<0.05)
sum(P.eq.N<0.10)

mean(T.eq.N)
sd(T.eq.N)
mean(Y.eq.N)
sd(Y.eq.N)

