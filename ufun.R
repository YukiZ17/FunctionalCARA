bernstein = function(Y, V.est, Tr.est, d.est, r, Z, lambda){
  n <- nrow(Y)
  Bfull.est <- basis(V.est,Tr.est,d.est,r,Z)
  B <- Bfull.est$B
  ind.inside <- Bfull.est$Ind.inside
  Q2 <- Bfull.est$Q2
  K <- Bfull.est$K
  Y <- matrix(Y[,ind.inside],nrow=n)
  lambda <- as.matrix(lambda)
  t.area = Bfull.est$tria.all 
  
  this.call <- match.call()
  n <- nrow(Y)
  npix <- ncol(Y)
  J <- ncol(Q2)
  
  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,t(Y))
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  
  gcv_all <- sapply(lambda,FUN=function(Lam){  
    Dlam <- Lam*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    theta <- crossprod(t(lhs.inv),rhs)
    gamma <- crossprod(t(Q2),theta)
    Yhat <- crossprod(t(W),theta)
    res <- t(Y)-Yhat
    sse <- apply(res^2,2,sum)
    if(!flag){
      df <- sum(1/(1+Cval*Lam))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df <- sum(diag(Hmtx))
    }
    gcv <- npix*sse/(npix-df)^2
  })
  gcv_all <- matrix(gcv_all,nrow=n)
  lam.ind <- apply(gcv_all,1,which.min)
  lambdac <- lambda[lam.ind]
  
  theta <- c()
  gamma <- c()
  Yhat <- c()
  df <- c()
  for (i in 1:n){
    lamc.tmp <- lambdac[i]
    Dlam <- lamc.tmp*D
    lhs <- WW+Dlam
    lhs.inv <- chol2inv(chol(lhs));
    rhs.tmp <- as.matrix(rhs[,i],ncol=1)
    theta.tmp <- crossprod(t(lhs.inv),rhs.tmp)
    theta <- cbind(theta,theta.tmp)
    gamma.tmp <- crossprod(t(Q2),theta.tmp) 
    gamma <- cbind(gamma,gamma.tmp)
    Yhat.tmp <- crossprod(t(W),theta.tmp)
    Yhat <- cbind(Yhat,Yhat.tmp)
    if(!flag){
      df.tmp <- sum(1/(1+Cval*lamc.tmp))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df.tmp <- sum(diag(Hmtx))
    }
    df <- c(df,df.tmp)
  }
  
  return(list(B, gamma, lambdac, t.area))
}



sfpca_img.binary = function(type = "bernstein", Y , train_y, train_dat.id=NA, theta, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr){
  

  if(type == 'bernstein'){
    est = bernstein(Y= Y, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                    r = r, Z = Z, lambda = lambda)
    basis = est[[1]];scores = t(est[[2]]);t.area = est[[4]]
  }
  
  org.idx = c(1:ncr*ncr)
  desMat = matrix(0, nrow = ncr*ncr, ncol = ncol(scores))
  
  for(zz in 1:ncol(desMat)){
    
    desMat[ind.inside,zz] = basis[,zz]
  }
  
  B <- aperm(array(desMat, c(ncr, ncr,ncol(scores))), c(3, 1, 2)) 
  mya = array(dim = c(1, ncr, ncr))
  
  for(i in 1){
    mya[i,,] = matrix(Y[i,], ncr, ncr)
  }
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W = MFPCA:::calcBasisIntegrals(B, 2,g@argvals) 
  
  S = t(scores)
  maty = matrix(rep(train_y,each=nrow(S)),nrow=nrow(S))
  M = rowSums(maty*W%*%S)
  MM = as.matrix(M)%*%t(as.matrix(M))
  sqrM = function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  U =theta/length(train_y)*W%*%S%*%t(S)%*%W+(1-theta)*MM/(length(train_y)^2); G = W 
  halfG_inv  = sqrM(G)
  tM = t(halfG_inv)%*%U%*%halfG_inv
  eigen_res = eigen(tM) 
  
  fd_list = lapply(1:npc, function(ipc){
    coef_pc= halfG_inv%*%as.matrix(Real(eigen_res$vectors[,ipc])) 
  })
  
  sup_basis = NULL
  for(k in 1:npc){
    kth = matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis = cbind(sup_basis, t(kth))
  }
  
  return(list(sup_basis, eigen_res$values[1:npc], basis, scores, t.area))
}

sfpca_img.binary.2 = function(type = "bernstein", Y , train_y, train_dat.id=NA, theta, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr){
  

  
  if(type == 'bernstein'){
    est = bernstein(Y= Y, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                    r = r, Z = Z, lambda = lambda)
    basis = est[[1]];scores = t(est[[2]]);t.area = est[[4]]
  }
  
  org.idx = c(1:ncr*ncr)
  desMat = matrix(0, nrow = ncr*ncr, ncol = ncol(scores))
  for(zz in 1:ncol(desMat)){
    desMat[ind.inside,zz] = basis[,zz]
  }
  
  B <- aperm(array(desMat, c(ncr, ncr,ncol(scores))), c(3, 1, 2)) 
  mya = array(dim = c(1, ncr, ncr))
  for(i in 1){
    mya[i,,] = matrix(Y[i,], ncr, ncr)
  }
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W = MFPCA:::calcBasisIntegrals(B, 2,g@argvals) 
  
  S = t(scores)
  maty1 = matrix(rep(train_y,each=nrow(S)),nrow=nrow(S))
  maty0 = matrix(rep(1-train_y,each=nrow(S)),nrow=nrow(S))
  M1 = rowSums(maty1*W%*%S)
  M0 = rowSums(maty0*W%*%S)
  MM1 = as.matrix(M1)%*%t(as.matrix(M1))
  MM0 = as.matrix(M0)%*%t(as.matrix(M0))
  
  sqrM = function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  U =theta/length(train_y)*W%*%S%*%t(S)%*%W+(1-theta)*(MM1/sum(train_y)+MM0/(length(train_y)-sum(train_y))); G = W 
  halfG_inv  = sqrM(G)
  tM = t(halfG_inv)%*%U%*%halfG_inv
  eigen_res = eigen(tM) 
  
  fd_list = lapply(1:npc, function(ipc){
    coef_pc= halfG_inv%*%as.matrix(Real(eigen_res$vectors[,ipc])) 
  })
  
  sup_basis = NULL
  for(k in 1:npc){
    kth = matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis = cbind(sup_basis, t(kth))
  }
  
  return(list(sup_basis, eigen_res$values[1:npc], basis, scores, t.area))
}


fpca_img.binary = function(type = "bernstein", Y ,n_length, train_y=NA, train_dat.id=NA, theta=NA, lambda, npc, V.est, Tr.est, d.est, r, Z, ncr){
  

  
  if(type == 'bernstein'){
    est = bernstein(Y= Y, V.est = V.est, Tr.est = Tr.est, d.est = d.est, 
                    r = r, Z = Z, lambda = lambda)
    basis = est[[1]];scores = t(est[[2]]);t.area = est[[4]]
  }
  
  org.idx = c(1:ncr*ncr)
  desMat = matrix(0, nrow = ncr*ncr, ncol = ncol(scores))
  for(zz in 1:ncol(desMat)){
    desMat[ind.inside,zz] = basis[,zz]
  }
  
  B <- aperm(array(desMat, c(ncr, ncr,ncol(scores))), c(3, 1, 2)) 
  mya = array(dim = c(1, ncr, ncr))
  for(i in 1){
    mya[i,,] = matrix(Y[i,], ncr, ncr)
  }
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W = MFPCA:::calcBasisIntegrals(B, 2,g@argvals) 
  
  S = t(scores)
  
  sqrM = function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  U =W%*%S%*%t(S)%*%W/n_length; G = W 
  halfG_inv  = sqrM(G)
  tM = t(halfG_inv)%*%U%*%halfG_inv
  eigen_res = eigen(tM) 
  
  fd_list = lapply(1:npc, function(ipc){
    coef_pc= halfG_inv%*%as.matrix(Real(eigen_res$vectors[,ipc])) 
  })
  
  sup_basis = NULL
  for(k in 1:npc){
    kth = matrix(fd_list[[k]], nrow = 1) %*% t(basis)
    sup_basis = cbind(sup_basis, t(kth))
  }
  
  return(list(sup_basis, eigen_res$values[1:npc], basis, scores, t.area))
}
