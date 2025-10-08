############### Real Case Study ############

## Simulation parameter

N <- 1000
gamma.cara.1.N <- c()
alpha.cara.1.N <- c()
T.cara.1.N <- c()
Y.cara.1.N <- c()
P.cara.1.N <- c()
age.cara.1.N <- c()
dna.cara.1.N <- c()

gamma.cara.2.N <- c()
alpha.cara.2.N <- c()
T.cara.2.N <- c()
Y.cara.2.N <- c()
P.cara.2.N <- c()
age.cara.2.N <- c()
dna.cara.2.N <- c()

gamma.cara.3.N <- c()
alpha.cara.3.N <- c()
T.cara.3.N <- c()
Y.cara.3.N <- c()
P.cara.3.N <- c()
age.cara.3.N <- c()
dna.cara.3.N <- c()

gamma.eq.N <- c()
alpha.eq.N <- c()
T.eq.N <- c()
Y.eq.N <- c()
P.eq.N <- c()
age.eq.N <- c()
dna.eq.N <- c()

acc1 <- c()
acc2 <- c()
acc3 <- c()
acceq <- c()


## Parameter inputs: from summary data
n.sample <- 1200
n.burnin <- 200
n.batch <- 100
n.test <- 200
batch <- (n.sample-n.burnin)/n.batch
alpha = 1.32
gamma = -0.506
tdna <- -0.759
tage <- -0.0354
theta = 0.1

for (j in 1:N){
  
  ## Result Storage
  T.all.cara.1 <- c()
  Y.all.cara.1 <- c()
  Alpha.batch.cara.1 <- c()
  b.1.batch.cara.1 <- c()
  b.2.batch.cara.1 <- c()
  bint.1.batch.cara.1 <- c()
  bint.2.batch.cara.1 <- c()
  age.batch.cara.1 <- c()
  dna.batch.cara.1 <- c()
  Gamma.batch.cara.1 <- c()
  
  T.all.cara.2 <- c()
  Y.all.cara.2 <- c()
  Alpha.batch.cara.2 <- c()
  b.1.batch.cara.2 <- c()
  b.2.batch.cara.2 <- c()
  bint.1.batch.cara.2 <- c()
  bint.2.batch.cara.2 <- c()
  age.batch.cara.2 <- c()
  dna.batch.cara.2 <- c()
  Gamma.batch.cara.2 <- c()
  
  T.all.cara.3 <- c()
  Y.all.cara.3 <- c()
  Alpha.batch.cara.3 <- c()
  b.1.batch.cara.3 <- c()
  b.2.batch.cara.3 <- c()
  bint.1.batch.cara.3 <- c()
  bint.2.batch.cara.3 <- c()
  age.batch.cara.3 <- c()
  dna.batch.cara.3 <- c()
  Gamma.batch.cara.3 <- c()
  
  
  
  ## Generating data
  a <- train_dat.id$Z1[1,];na <- which(is.na(a))
  
  basis <- eFun(argvals = seq(0, 1, length.out=40), M = 2,  type = "Poly")@X
  
  phi.1 <- basis[1,] %*% t(basis[2,]); dim(phi.1) <- c(1,1600)
  phi.2 <- basis[2,] %*% t(basis[1,]); dim(phi.2) <- c(1,1600)
  phi.3 <- basis[2,] %*% t(basis[2,]); dim(phi.3) <- c(1,1600)
  
  w1 <- matrix(rnorm(n.sample,0,sqrt(4)),n.sample,1)
  w2 <- matrix(rnorm(n.sample,0,sqrt(2)),n.sample,1)
  w3 <- matrix(rnorm(n.sample,0,sqrt(0.5)),n.sample,1)
  
  X<- w1%*%phi.1+w2%*%phi.2+w3%*%phi.3 
  X[,na] <- NA
  X.inside <- X
  X.inside[,na] <- 0
  
  age <- runif(n.sample,10,80)
  dna <- rbinom(n.sample,1,0.5)
  
  w1.test <- matrix(rnorm(n.test,0,sqrt(4)),n.test,1)
  w2.test <- matrix(rnorm(n.test,0,sqrt(2)),n.test,1)
  w3.test <- matrix(rnorm(n.test,0,sqrt(0.5)),n.test,1)
  
  X.test<- w1.test%*%phi.1+w2.test%*%phi.2+w3.test%*%phi.3 
  X.test[,na] <- 0
  
  t.true <- (gamma-4.45*X.test%*%t(phi.2)/921)>0
  
  
  
  ## Burnin -- Equal:
  
  t.burnin <- rbinom(n.burnin,size=1,prob = 0.5)
  
  T.all.cara.1 <- c(T.all.cara.1,t.burnin)
  T.all.cara.2 <- c(T.all.cara.2,t.burnin)
  T.all.cara.3 <- c(T.all.cara.3,t.burnin)
  
  pi.burnin <- rep(alpha,n.burnin)+ gamma*t.burnin+ ((2.55-4.45*t.burnin)*X.inside[1:n.burnin,]%*%t(phi.2)/921)+(-2.15*X.inside[1:n.burnin,]%*%t(phi.1)/921)
  + tdna*dna[1:n.burnin]+tage*age[1:n.burnin]
  
  Y.outcome.burnin<- rbinom(n.burnin,1,(1/(1+exp(-pi.burnin))))
  
  Y.all.cara.1 <- c(Y.all.cara.1,Y.outcome.burnin)
  Y.all.cara.2 <- c(Y.all.cara.2,Y.outcome.burnin)
  Y.all.cara.3 <- c(Y.all.cara.3,Y.outcome.burnin)
  
  
  sup_fun.binary=sfpca_img.binary.2(type = "bernstein", Y = X.inside[1:n.burnin,], train_y=Y.outcome.burnin,
                                    theta = theta, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                    r = r, Z = Z, ncr = ncr)
  
  o_tfs.bin = as.matrix(sup_fun.binary[[1]])
  o_tfs.bin.2 <- o_tfs.bin
  o_tfs.bin.3 <- o_tfs.bin
  
  
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
  train$age <- age[1:n.burnin]
  train$dna <- dna[1:n.burnin]
  
  
  glm_fit_burnin = glm(y~score1+score2+t+age+dna+score1*t+score2*t,data=train,family="binomial")
  Alpha.batch.cara.1 <- c(Alpha.batch.cara.1,glm_fit_burnin$coefficients[1])
  b.1.batch.cara.1 <- c(b.1.batch.cara.1,glm_fit_burnin$coefficients[2])
  b.2.batch.cara.1 <- c(b.2.batch.cara.1,glm_fit_burnin$coefficients[3])
  Gamma.batch.cara.1 <- c(Gamma.batch.cara.1,glm_fit_burnin$coefficients[4])
  age.batch.cara.1 <- c(age.batch.cara.1,glm_fit_burnin$coefficients[5])
  dna.batch.cara.1 <- c(dna.batch.cara.1,glm_fit_burnin$coefficients[6])
  bint.1.batch.cara.1 <- c(bint.1.batch.cara.1,glm_fit_burnin$coefficients[7])
  bint.2.batch.cara.1 <- c(bint.2.batch.cara.1,glm_fit_burnin$coefficients[8])
  
  
  Alpha.batch.cara.2 <- c(Alpha.batch.cara.2,glm_fit_burnin$coefficients[1])
  b.1.batch.cara.2 <- c(b.1.batch.cara.2,glm_fit_burnin$coefficients[2])
  b.2.batch.cara.2 <- c(b.2.batch.cara.2,glm_fit_burnin$coefficients[3])
  Gamma.batch.cara.2 <- c(Gamma.batch.cara.2,glm_fit_burnin$coefficients[4])
  age.batch.cara.2 <- c(age.batch.cara.2,glm_fit_burnin$coefficients[5])
  dna.batch.cara.2 <- c(dna.batch.cara.2,glm_fit_burnin$coefficients[6])
  bint.1.batch.cara.2 <- c(bint.1.batch.cara.2,glm_fit_burnin$coefficients[7])
  bint.2.batch.cara.2 <- c(bint.2.batch.cara.2,glm_fit_burnin$coefficients[8])
  
  Alpha.batch.cara.3 <- c(Alpha.batch.cara.3,glm_fit_burnin$coefficients[1])
  b.1.batch.cara.3 <- c(b.1.batch.cara.3,glm_fit_burnin$coefficients[2])
  b.2.batch.cara.3 <- c(b.2.batch.cara.3,glm_fit_burnin$coefficients[3])
  Gamma.batch.cara.3 <- c(Gamma.batch.cara.3,glm_fit_burnin$coefficients[4])
  age.batch.cara.3 <- c(age.batch.cara.3,glm_fit_burnin$coefficients[5])
  dna.batch.cara.3 <- c(dna.batch.cara.3,glm_fit_burnin$coefficients[6])
  bint.1.batch.cara.3 <- c(bint.1.batch.cara.3,glm_fit_burnin$coefficients[7])
  bint.2.batch.cara.3 <- c(bint.2.batch.cara.3,glm_fit_burnin$coefficients[8])
  
  
  ################## Adaptive CARA 3 Types ####################
  
  ###### CARA 1: Classical #########
  
  for (k in 1:batch){
    
    
    l.k <- n.burnin+1+(k-1)*n.batch
    u.k <- n.burnin+k*n.batch
    X.k <- X.inside[l.k:u.k,]
    
    basis.est = bernstein(Y = X.k, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                          lambda = lambda)
    
    
    b.basis.new.bin = as.matrix(basis.est[[1]]); b.scores.new.bin = as.matrix(t(basis.est[[2]]))
    
    score_sup_predict.bin = t(t(o_tfs.bin) %*% ((b.basis.new.bin%*%t(b.scores.new.bin))))
    
    
    
    logit.p.a.k.hat <- Alpha.batch.cara.1[k]+(b.1.batch.cara.1[k]+bint.1.batch.cara.1[k])*score_sup_predict.bin[,1]+(b.2.batch.cara.1[k]+bint.2.batch.cara.1[k])*score_sup_predict.bin[,2]+Gamma.batch.cara.1[k]
    +age.batch.cara.1[k]*age[l.k:u.k]+dna.batch.cara.1[k]*dna[l.k:u.k]
    logit.p.b.k.hat <- Alpha.batch.cara.1[k]+(b.1.batch.cara.1[k])*score_sup_predict.bin[,1]+(b.2.batch.cara.1[k])*score_sup_predict.bin[,2]
    +age.batch.cara.1[k]*age[l.k:u.k]+dna.batch.cara.1[k]*dna[l.k:u.k]
    
    p.a.k.hat <- 1/(1+exp(-logit.p.a.k.hat))
    p.b.k.hat <- 1/(1+exp(-logit.p.b.k.hat))
    odd.a.hat <- p.a.k.hat/(1-p.a.k.hat)
    odd.b.hat <- p.b.k.hat/(1-p.b.k.hat)
    
    
    t.k <- rbinom(n.batch,size=1,prob = odd.a.hat/(odd.a.hat+odd.b.hat))
    
    T.all.cara.1 <- c(T.all.cara.1,t.k)
    
    
    pi.k <- rep(alpha,n.batch)+ gamma*t.k +((2.55-4.45*t.k)*X.k%*%t(phi.2)/921)+(-2.15*X.k%*%t(phi.1)/921)
    +tdna*dna[l.k:u.k]+tage*age[l.k:u.k]
    
    Y.outcome.k<- rbinom(n.batch,1,(1/(1+exp(-pi.k))))
    
    Y.all.cara.1 <- c(Y.all.cara.1,Y.outcome.k)
    
    sup_fun.binary=sfpca_img.binary.2(type = "bernstein", Y = X.inside[1:u.k,], train_y=Y.all.cara.1,
                                      theta = theta, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                      r = r, Z = Z, ncr = ncr)
    
    
    o_tfs.bin = as.matrix(sup_fun.binary[[1]])
    
    
    b.basis.bin = as.matrix(sup_fun.binary[[3]]); b.scores.bin = as.matrix(sup_fun.binary[[4]])
    
    score_sup.bin = t(o_tfs.bin) %*% (b.basis.bin%*%t(b.scores.bin))
    
    
    if(npc == 2){
      score_sup.bin = t(score_sup.bin)
    }else{
      score_sup.bin = cbind(score_sup.bin[1,])
    }
    
    
    id <- 1:u.k
    
    train <- data.frame(id)
    train$y <- Y.all.cara.1
    train$score1 <- score_sup.bin[id,1]
    train$score2 <- score_sup.bin[id,2]
    train$t <- T.all.cara.1
    train$age <- age[1:u.k]
    train$dna <- dna[1:u.k]
    
    
    
    glm_fit = glm(y~score1+score2+t+age+dna+score1*t+score2*t,data=train,family="binomial")
    Alpha.batch.cara.1 <- c(Alpha.batch.cara.1,glm_fit$coefficients[1])
    b.1.batch.cara.1 <- c(b.1.batch.cara.1,glm_fit$coefficients[2])
    b.2.batch.cara.1 <- c(b.2.batch.cara.1,glm_fit$coefficients[3])
    Gamma.batch.cara.1 <- c(Gamma.batch.cara.1,glm_fit$coefficients[4])
    age.batch.cara.1 <- c(age.batch.cara.1,glm_fit$coefficients[5])
    dna.batch.cara.1 <- c(dna.batch.cara.1,glm_fit$coefficients[6])
    bint.1.batch.cara.1 <- c(bint.1.batch.cara.1,glm_fit$coefficients[7])
    bint.2.batch.cara.1 <- c(bint.2.batch.cara.1,glm_fit$coefficients[8])
    
    
    print(k)
    
  }
  
  gamma.cara.1.N <- c(gamma.cara.1.N,Gamma.batch.cara.1[length(Gamma.batch.cara.1)])
  alpha.cara.1.N <- c(alpha.cara.1.N,Alpha.batch.cara.1[length(Alpha.batch.cara.1)])
  age.cara.1.N <- c(age.cara.1.N,age.batch.cara.1[length(age.batch.cara.1)])
  dna.cara.1.N <- c(dna.cara.1.N,dna.batch.cara.1[length(dna.batch.cara.1)])
  T.cara.1.N <- c(T.cara.1.N,sum(T.all.cara.1))
  Y.cara.1.N <- c(Y.cara.1.N,sum(Y.all.cara.1))
  
  
  wald.coef <- summary(glm_fit)$coefficients 
  p.cara.1 <- wald.coef["t","Pr(>|z|)"]
  
  P.cara.1.N <- c(P.cara.1.N,p.cara.1)
  
  ####### Test #######
  basis.est = bernstein(Y = X.test, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                        lambda = lambda)
  
  
  b.basis.new.bin = as.matrix(basis.est[[1]]); b.scores.new.bin = as.matrix(t(basis.est[[2]]))
  
  score_sup_predict.bin = t(t(o_tfs.bin) %*% ((b.basis.new.bin%*%t(b.scores.new.bin))))
  
  
  t.test.1 <- (glm_fit$coefficients[8]*score_sup_predict.bin[,2]
               +glm_fit$coefficients[7]*score_sup_predict.bin[,1]
               +glm_fit$coefficients[4])>0
  
  acc1 <- c(acc1,sum(t.true == t.test.1))
  
  
  ###### CARA 2: Optimal #########  
  
  for (k in 1:batch){
    
    
    l.k <- n.burnin+1+(k-1)*n.batch
    u.k <- n.burnin+k*n.batch
    X.k <- X.inside[l.k:u.k,]
    
    basis.est = bernstein(Y = X.k, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                          lambda = lambda)
    
    
    b.basis.new.bin = as.matrix(basis.est[[1]]); b.scores.new.bin = as.matrix(t(basis.est[[2]]))
    
    score_sup_predict.bin.2 = t(t(o_tfs.bin.2) %*% ((b.basis.new.bin%*%t(b.scores.new.bin))))
    
    
    logit.p.a.k.hat <- Alpha.batch.cara.2[k]+(b.1.batch.cara.2[k]+bint.1.batch.cara.2[k])*score_sup_predict.bin.2[,1]+(b.2.batch.cara.2[k]+bint.2.batch.cara.2[k])*score_sup_predict.bin.2[,2]+Gamma.batch.cara.2[k]
    +age.batch.cara.2[k]*age[l.k:u.k]+dna.batch.cara.2[k]*dna[l.k:u.k]
    logit.p.b.k.hat <- Alpha.batch.cara.2[k]+(b.1.batch.cara.2[k])*score_sup_predict.bin.2[,1]+(b.2.batch.cara.2[k])*score_sup_predict.bin.2[,2]
    +age.batch.cara.2[k]*age[l.k:u.k]+dna.batch.cara.2[k]*dna[l.k:u.k]
    
    p.a.k.hat <- 1/(1+exp(-logit.p.a.k.hat))
    p.b.k.hat <- 1/(1+exp(-logit.p.b.k.hat))
    
    t.k.2 <- rbinom(n.batch,size=1,prob = sqrt(p.a.k.hat)/(sqrt(p.a.k.hat)+sqrt(p.b.k.hat)))
    T.all.cara.2 <- c(T.all.cara.2,t.k.2)
    
    pi.k <- rep(alpha,n.batch)+ gamma*t.k.2 +((2.55-4.45*t.k.2)*X.k%*%t(phi.2)/921)+(-2.15*X.k%*%t(phi.1)/921)
    +tdna*dna[l.k:u.k]+tage*age[l.k:u.k]
    Y.outcome.k.2<- rbinom(n.batch,1,(1/(1+exp(-pi.k))))
    Y.all.cara.2 <- c(Y.all.cara.2,Y.outcome.k.2)
    
    
    sup_fun.binary.2=sfpca_img.binary.2(type = "bernstein", Y = X.inside[1:u.k,], train_y=Y.all.cara.2,
                                        theta = theta, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                        r = r, Z = Z, ncr = ncr)
    
    
    o_tfs.bin.2 <-  as.matrix(sup_fun.binary.2[[1]])
    
    
    b.basis.bin.2 = as.matrix(sup_fun.binary.2[[3]]); b.scores.bin.2 = as.matrix(sup_fun.binary.2[[4]])
    
    score_sup.bin.2 = t(t(o_tfs.bin.2)%*% (b.basis.bin.2%*%t(b.scores.bin.2)))
    
    
    
    id <- 1:u.k
    
    train.2 <- data.frame(id)
    train.2$y <- Y.all.cara.2
    train.2$score1 <- score_sup.bin.2[id,1]
    train.2$score2 <- score_sup.bin.2[id,2]
    train.2$t <- T.all.cara.2
    train.2$age <- age[1:u.k]
    train.2$dna <- dna[1:u.k]
    
    
    
    glm_fit.2 = glm(y~score1+score2+t+age+dna+score1*t+score2*t,data=train.2,family="binomial")
    Alpha.batch.cara.2 <- c(Alpha.batch.cara.2,glm_fit.2$coefficients[1])
    b.1.batch.cara.2 <- c(b.1.batch.cara.2,glm_fit.2$coefficients[2])
    b.2.batch.cara.2 <- c(b.2.batch.cara.2,glm_fit.2$coefficients[3])
    Gamma.batch.cara.2 <- c(Gamma.batch.cara.2,glm_fit.2$coefficients[4])
    age.batch.cara.2 <- c(age.batch.cara.2,glm_fit.2$coefficients[5])
    dna.batch.cara.2 <- c(dna.batch.cara.2,glm_fit.2$coefficients[6])
    bint.1.batch.cara.2 <- c(bint.1.batch.cara.2,glm_fit.2$coefficients[7])
    bint.2.batch.cara.2 <- c(bint.2.batch.cara.2,glm_fit.2$coefficients[8])
    
    
    print(k)
    
  }
  
  gamma.cara.2.N <- c(gamma.cara.2.N,Gamma.batch.cara.2[length(Gamma.batch.cara.2)])
  alpha.cara.2.N <- c(alpha.cara.2.N,Alpha.batch.cara.2[length(Alpha.batch.cara.2)])
  age.cara.2.N <- c(age.cara.2.N,age.batch.cara.2[length(age.batch.cara.2)])
  dna.cara.2.N <- c(dna.cara.2.N,dna.batch.cara.2[length(dna.batch.cara.2)])
  T.cara.2.N <- c(T.cara.2.N,sum(T.all.cara.2))
  Y.cara.2.N <- c(Y.cara.2.N,sum(Y.all.cara.2))
  
  wald.coef.2 <- summary(glm_fit.2)$coefficients 
  p.cara.2 <- wald.coef.2["t","Pr(>|z|)"]
  
  P.cara.2.N <- c(P.cara.2.N,p.cara.2)
  
  ########### Test ###############
  
  
  basis.est = bernstein(Y = X.test, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                        lambda = lambda)
  
  
  b.basis.new.bin = as.matrix(basis.est[[1]]); b.scores.new.bin = as.matrix(t(basis.est[[2]]))
  
  score_sup_predict.bin.2 = t(t(o_tfs.bin.2) %*% ((b.basis.new.bin%*%t(b.scores.new.bin))))
  
  
  t.test.2 <- (glm_fit.2$coefficients[8]*score_sup_predict.bin.2[,2]
               +glm_fit.2$coefficients[7]*score_sup_predict.bin.2[,1]
               +glm_fit.2$coefficients[4])>0
  
  
  
  acc2 <- c(acc2,sum(t.true == t.test.2))
  
  
  
  ###### CARA 3: Neyman #########  
  
  for (k in 1:batch){
    
    
    l.k <- n.burnin+1+(k-1)*n.batch
    u.k <- n.burnin+k*n.batch
    X.k <- X.inside[l.k:u.k,]
    
    basis.est.3 = bernstein(Y = X.k, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                            lambda = lambda)
    
    
    b.basis.new.bin.3 = as.matrix(basis.est.3[[1]]); b.scores.new.bin.3 = as.matrix(t(basis.est.3[[2]]))
    
    score_sup_predict.bin.3 = t(t(o_tfs.bin.3) %*% ((b.basis.new.bin.3%*%t(b.scores.new.bin.3))))
    
    
    logit.p.a.k.hat.3 <- Alpha.batch.cara.3[k]+(b.1.batch.cara.3[k]+bint.1.batch.cara.3[k])*score_sup_predict.bin.3[,1]+(b.2.batch.cara.3[k]+bint.2.batch.cara.3[k])*score_sup_predict.bin.3[,2]+Gamma.batch.cara.3[k]
    +age.batch.cara.3[k]*age[l.k:u.k]+dna.batch.cara.3[k]*dna[l.k:u.k]
    logit.p.b.k.hat.3 <- Alpha.batch.cara.3[k]+(b.1.batch.cara.3[k])*score_sup_predict.bin.3[,1]+(b.2.batch.cara.3[k])*score_sup_predict.bin.3[,2]
    +age.batch.cara.3[k]*age[l.k:u.k]+dna.batch.cara.3[k]*dna[l.k:u.k]
    
    p.a.k.hat.3 <- 1/(1+exp(-logit.p.a.k.hat.3))
    p.b.k.hat.3 <- 1/(1+exp(-logit.p.b.k.hat.3))
    q.a.k.hat.3 <- 1-p.a.k.hat.3
    q.b.k.hat.3 <- 1-p.b.k.hat.3
    
    
    t.k.3 <- rbinom(n.batch,size=1,prob = sqrt(p.a.k.hat.3*q.a.k.hat.3)/(sqrt(p.a.k.hat.3*q.a.k.hat.3)+sqrt(p.b.k.hat.3*q.b.k.hat.3)))
    T.all.cara.3 <- c(T.all.cara.3,t.k.3)
    
    pi.k.3 <-  rep(alpha,n.batch)+ gamma*t.k.3 +((2.55-4.45*t.k.3)*X.k%*%t(phi.2)/921)+(-2.15*X.k%*%t(phi.1)/921)
    +tdna*dna[l.k:u.k]+tage*age[l.k:u.k]
    Y.outcome.k.3<- rbinom(n.batch,1,(1/(1+exp(-pi.k.3))))
    Y.all.cara.3 <- c(Y.all.cara.3,Y.outcome.k.3)
    
    
    
    sup_fun.binary.3=sfpca_img.binary.2(type = "bernstein", Y = X.inside[1:u.k,], train_y=Y.all.cara.3,
                                        theta = theta, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                                        r = r, Z = Z, ncr = ncr)
    
    
    o_tfs.bin.3 <-  as.matrix(sup_fun.binary.3[[1]])
    
    
    b.basis.bin.3 = as.matrix(sup_fun.binary.3[[3]]); b.scores.bin.3 = as.matrix(sup_fun.binary.3[[4]])
    
    score_sup.bin.3 = t(t(o_tfs.bin.3) %*% (b.basis.bin.3%*%t(b.scores.bin.3)))
    
    
    
    
    id <- 1:u.k
    
    train.3 <- data.frame(id)
    train.3$y <- Y.all.cara.3
    train.3$score1 <- score_sup.bin.3[id,1]
    train.3$score2 <- score_sup.bin.3[id,2]
    train.3$t <- T.all.cara.3
    train.3$age <- age[1:u.k]
    train.3$dna <- dna[1:u.k]
    
    
    
    glm_fit.3 = glm(y~score1+score2+t+age+dna+score1*t+score2*t,data=train.3,family="binomial")
    Alpha.batch.cara.3 <- c(Alpha.batch.cara.3,glm_fit.3$coefficients[1])
    b.1.batch.cara.3 <- c(b.1.batch.cara.3,glm_fit.3$coefficients[2])
    b.2.batch.cara.3 <- c(b.2.batch.cara.3,glm_fit.3$coefficients[3])
    Gamma.batch.cara.3 <- c(Gamma.batch.cara.3,glm_fit.3$coefficients[4])
    
    age.batch.cara.3 <- c(age.batch.cara.3,glm_fit.3$coefficients[5])
    dna.batch.cara.3 <- c(dna.batch.cara.3,glm_fit.3$coefficients[6])
    bint.1.batch.cara.3 <- c(bint.1.batch.cara.3,glm_fit.3$coefficients[7])
    bint.2.batch.cara.3 <- c(bint.2.batch.cara.3,glm_fit.3$coefficients[8])
    
    print(k)
    
  }
  
  gamma.cara.3.N <- c(gamma.cara.3.N,Gamma.batch.cara.3[length(Gamma.batch.cara.3)])
  alpha.cara.3.N <- c(alpha.cara.3.N,Alpha.batch.cara.3[length(Alpha.batch.cara.3)])
  age.cara.3.N <- c(age.cara.3.N,age.batch.cara.3[length(age.batch.cara.3)])
  dna.cara.3.N <- c(dna.cara.3.N,dna.batch.cara.3[length(dna.batch.cara.3)])
  T.cara.3.N <- c(T.cara.3.N,sum(T.all.cara.3))
  Y.cara.3.N <- c(Y.cara.3.N,sum(Y.all.cara.3))
  
  wald.coef.3 <- summary(glm_fit.3)$coefficients 
  p.cara.3 <- wald.coef.3["t","Pr(>|z|)"]
  
  P.cara.3.N <- c(P.cara.3.N,p.cara.3)
  
  ################### Test #######################
  
  
  basis.est = bernstein(Y = X.test, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                        lambda = lambda)
  
  
  b.basis.new.bin = as.matrix(basis.est[[1]]); b.scores.new.bin = as.matrix(t(basis.est[[2]]))
  
  score_sup_predict.bin.3 = t(t(o_tfs.bin.3) %*% ((b.basis.new.bin%*%t(b.scores.new.bin))))
  
  
  t.test.3 <- (glm_fit.3$coefficients[8]*score_sup_predict.bin.3[,2]
               +glm_fit.3$coefficients[7]*score_sup_predict.bin.3[,1]
               +glm_fit.3$coefficients[4])>0
  
  
  
  
  acc3 <- c(acc3,sum(t.true == t.test.3))
  
  
  ################ Equal ############### 
  
  t.equal <- rbinom(n.sample,size=1,prob = 0.5)
  
  pi.equal <- rep(alpha,n.sample)+ gamma*t.equal +((2.55-4.45*t.equal)*X.inside%*%t(phi.2)/921)+(-2.15*X.inside%*%t(phi.1)/921)
  +tdna*dna[l.k:u.k]+tage*age[l.k:u.k]
  
  Y.outcome.equal<- rbinom(n.sample,1,(1/(1+exp(-pi.equal))))
  
  sup_fun.binary=sfpca_img.binary.2(type = "bernstein", Y = X.inside, train_y=Y.outcome.equal,
                                    theta = theta, lambda = lambda, npc = npc, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
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
  
  train.eq$age <- age
  train.eq$dna <- dna
  
  
  glm_fit_equal = glm(y~score1+score2+t+age+dna+score1*t+score2*t,data=train.eq,family="binomial")
  
  gamma.eq.N <- c(gamma.eq.N,glm_fit_equal$coefficients[4])
  alpha.eq.N <- c(alpha.eq.N,glm_fit_equal$coefficients[1])
  age.eq.N <- c(age.eq.N,glm_fit_equal$coefficients[5])
  dna.eq.N <- c(dna.eq.N,glm_fit_equal$coefficients[6])
  Y.eq.N <- c(Y.eq.N,sum(Y.outcome.equal))
  T.eq.N <- c(T.eq.N,sum(t.equal))
  
  wald.coef.eq <- summary(glm_fit_equal)$coefficients 
  p.eq <- wald.coef.eq["t","Pr(>|z|)"]
  
  P.eq.N <- c(P.eq.N,p.eq)
  
  ############ Test ##############
  
  basis.est = bernstein(Y = X.test, V.est = V.est, Tr.est = Tr.est, d.est = d.est, r = r, Z = Z,
                        lambda = lambda)
  
  
  b.basis.new.bin = as.matrix(basis.est[[1]]); b.scores.new.bin = as.matrix(t(basis.est[[2]]))
  
  score_sup_predict.bin.eq = t(t(o_tfs.bin.eq) %*% ((b.basis.new.bin%*%t(b.scores.new.bin))))
  
  
  t.test.eq <- (glm_fit_equal$coefficients[8]*score_sup_predict.bin.eq[,2]
                +glm_fit_equal$coefficients[7]*score_sup_predict.bin.eq[,1]
                +glm_fit_equal$coefficients[4])>0
  
  acceq <- c(acceq,sum(t.true == t.test.eq))
  
  
  print(j)  
  
}





