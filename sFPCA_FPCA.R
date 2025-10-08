########################## sFPCA versus FPCA #######################

## Simulation parameter

N <- 600

theta.sfpca <- c(0.00001,0.0001,0.001,0.01,0.1,1)

data.theta <- matrix(data=NA,nrow=N,ncol=5*length(theta.sfpca))



## Parameter inputs
n.sample <- 500
n.burnin <- 100
n.batch <- 40
batch <- (n.sample-n.burnin)/n.batch
alpha = -0.25
gamma = 0.5

alpha.eq.N.1 <- c()
gamma.eq.N.1 <- c()
Y.eq.N.1 <- c()
T.eq.N.1 <- c()
P.eq.N.1 <- c()


for(j in 1:N){
  
  ################ Equal ############### 
  
  t.equal <- rbinom(n.sample,size=1,prob = 0.5)
  
  pi.equal <- rep(alpha,n.sample)+ gamma*t.equal +(2*X.inside%*%t(phi.3)/921)
  
  Y.outcome.equal<- rbinom(n.sample,1,(1/(1+exp(-pi.equal))))
  
  sup_fun.binary=sfpca_img.binary.2(type = "bernstein", Y = X.inside, train_y=Y.outcome.equal,
                                    theta = 1, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                    r = r, Z = Z, ncr = ncr)
  
  
  o_tfs.bin.eq = as.matrix(sup_fun.binary[[1]])
  
  
  b.basis.bin.eq = as.matrix(sup_fun.binary[[3]]); b.scores.bin.eq = as.matrix(sup_fun.binary[[4]])
  
  score_sup.bin.eq = t(t(o_tfs.bin.eq) %*% (b.basis.bin.eq%*%t(b.scores.bin.eq)))
  
  
  
  
  id <- 1:n.sample
  
  train.eq <- data.frame(id)
  train.eq$y <- Y.outcome.equal
  train.eq$score1 <- score_sup.bin.eq[,1]
  train.eq$score2 <- score_sup.bin.eq[,2]
  train.eq$t <- t.equal
  
  
  glm_fit_equal = glm(y~score1+score2+t,data=train.eq,family="binomial")
  
  gamma.eq.N <- c(gamma.eq.N,glm_fit_equal$coefficients[4])
  alpha.eq.N <- c(alpha.eq.N,glm_fit_equal$coefficients[1])
  Y.eq.N <- c(Y.eq.N,sum(Y.outcome.equal))
  T.eq.N <- c(T.eq.N,sum(t.equal))
  
  wald.coef.eq <- summary(glm_fit_equal)$coefficients 
  p.eq <- wald.coef.eq["t","Pr(>|z|)"]
  
  P.eq.N <- c(P.eq.N,p.eq)
  
}



for (j in 1:N){
  
  
  
  
  
  ## Generating data
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
  
  
  
  for (i in 1:length(theta.sfpca)){
    
    ## Result Storage
    
    T.all <- c()
    Y.all <- c()
    Gamma.batch <- c()
    
    
    
    ## Burnin -- Equal:
    
    t.burnin <- rbinom(n.burnin,size=1,prob = 0.5)
    T.all <- c(T.all,t.burnin)
    pi.burnin <- rep(alpha,n.burnin)+ gamma*t.burnin +(2*X.inside[1:n.burnin,]%*%t(phi.3)/921) #+(X.inside[1:n.burnin,]%*%t(phi.2)/921) #+(X.inside[1:n.burnin,]%*%t(phi.1)/921) #divided by 1600? 921?
    Y.outcome.burnin<- rbinom(n.burnin,1,(1/(1+exp(-pi.burnin))))
    Y.all <- c(Y.all,Y.outcome.burnin)
    
    sup_fun.binary=sfpca_img.binary.2(type = "bernstein", Y = X.inside[1:n.burnin,], train_y=Y.outcome.burnin,
                                      theta = theta.sfpca[i], lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                      r = r, Z = Z, ncr = ncr)
    
    
    o_tfs.bin = as.matrix(sup_fun.binary[[1]])
    
    
    b.basis.bin = as.matrix(sup_fun.binary[[3]]); b.scores.bin = as.matrix(sup_fun.binary[[4]])
    
    score_sup.bin = t(o_tfs.bin) %*% (b.basis.bin%*%t(b.scores.bin))
    
    
    if(npc == 2){
      score_sup.bin = t(score_sup.bin)
    }else{
      score_sup.bin = cbind(score_sup.bin[1,])
    }
    
    #score_sup.bin
    
    id <- 1:n.burnin
    
    train <- data.frame(id)
    train$y <- Y.outcome.burnin
    train$score1 <- score_sup.bin[,1]
    train$score2 <- score_sup.bin[,2]
    train$t <- t.burnin
    
    
    glm_fit_burnin = glm(y~score1+score2+t,data=train,family="binomial")
    Gamma.batch <- c(Gamma.batch,glm_fit_burnin$coefficients[4])
    
    
    
    ## Adaptive
    
    
    for (k in 1:batch){
      
      
      l.k <- n.burnin+1+(k-1)*n.batch
      u.k <- n.burnin+k*n.batch
      X.k <- X.inside[l.k:u.k,]
      
      t.k <- rbinom(n.batch,size=1,prob =1/(1+exp(-Gamma.batch[k])))
      T.all <- c(T.all,t.k)
      pi.k <- rep(alpha,n.batch)+ gamma*t.k +(2*X.k%*%t(phi.3)/921) #+(X.k%*%t(phi.2)/921) #+(X.k%*%t(phi.1)/921) #divided by 1600? 921?
      Y.outcome.k<- rbinom(n.batch,1,(1/(1+exp(-pi.k))))
      Y.all <- c(Y.all,Y.outcome.k)
      
      
      
      sup_fun.binary=sfpca_img.binary.2(type = "bernstein", Y = X.inside[1:u.k,], train_y=Y.all,
                                        theta = theta.sfpca[i], lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                        r = r, Z = Z, ncr = ncr)
      
      
      o_tfs.bin = as.matrix(sup_fun.binary[[1]])
      
      
      b.basis.bin = as.matrix(sup_fun.binary[[3]]); b.scores.bin = as.matrix(sup_fun.binary[[4]])
      
      score_sup.bin = t(o_tfs.bin) %*% (b.basis.bin%*%t(b.scores.bin))
      
      
      if(npc == 2){
        score_sup.bin = t(score_sup.bin)
      }else{
        score_sup.bin = cbind(score_sup.bin[1,])
      }
      
      #score_sup.bin
      
      id <- 1:u.k
      
      train <- data.frame(id)
      train$y <- Y.all
      train$score1 <- score_sup.bin[,1]
      train$score2 <- score_sup.bin[,2]
      train$t <- T.all
      
      
      glm_fit = glm(y~score1+score2+t,data=train,family="binomial")
      Gamma.batch <- c(Gamma.batch,glm_fit$coefficients[4])
      
      
      
    }
    wald.coef <- summary(glm_fit)$coefficients 
    p <- wald.coef["t","Pr(>|z|)"]
    data.theta[j,5*i-4] <- glm_fit$coefficients[1] #alpha
    data.theta[j,5*i-3] <- glm_fit$coefficients[4] #gamma
    data.theta[j,5*i-2] <- sum(T.all)
    data.theta[j,5*i-1] <- sum(Y.all)
    data.theta[j,5*i] <- p
    
    print(i)
  }
  
  
  
  
  print(j)
  
  
  
}
