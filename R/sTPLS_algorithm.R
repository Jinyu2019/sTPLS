
sTPLS <- function(Tensor, nfactor=2, c1, c2, beta=2, initNum=5, niter=100, seed0=1, err=0.0001){
  ptm = proc.time()
  DIM = dim(Tensor)
  U = matrix(0, nrow=DIM[1], ncol=nfactor)
  V = matrix(0, nrow=DIM[2], ncol=nfactor)
  W = matrix(0, nrow=DIM[3], ncol=nfactor)
  rownames(U) = dimnames(Tensor)[[1]]
  rownames(V) = dimnames(Tensor)[[2]]
  rownames(W) = dimnames(Tensor)[[3]]
  
  d_iter = list()
  
  n01 <- DIM[1]*c1
  c1 <- round((sqrt(n01)-1)/(sqrt(DIM[1])-1),digits=2)
  
  n02 <- DIM[2]*c2
  c2 <- round((sqrt(n02)-1)/(sqrt(DIM[2])-1),digits=2)
  
  # print(paste0("c1=", c1, ", c2=",c2, ", beta=", beta, ", nfactor=", nfactor))
  
  for(j in 1:nfactor){
    # print(j)
    out <- L1_STPLS_rank1(Tensor, c1, c2, beta, initNum, niter, seed0, err)
    # save results
    d_iter[[j]] = out$objs
    U[,j] = out$u
    V[,j] = out$v
    W[,j] = out$w
    Tensor = update.tensor(Tensor,U[,j],V[,j],W[,j])
  }
  tim = proc.time() - ptm
  return (list(U=U, V=V, W=W, d_iter=d_iter, time = tim[1]))
}

# sTPLS algorithm
L1_STPLS_rank1 = function(Tensor, c1, c2, beta=2, initNum=5, niter=100, seed0=1, err=0.0001){
  ptm = proc.time()
  DIM = dim(Tensor)
  opt_obj = 10^6
  for(init in 1:initNum){
    set.seed(init*seed0)
    u = matrix(rnorm(DIM[1]),ncol=1); u = u/norm(u,'E')
    v = matrix(rnorm(DIM[2]),ncol=1); v = v/norm(v,'E')
    w = matrix(runif(DIM[3],min=0.1,max=0.9),nrow=DIM[3],ncol=1)
    
    objs= c()
    # Iterative algorithm
    for (j in 1:niter){
      XTY = compute.C(Tensor, w, DIM[3])
      
      z.u = XTY%*%v
      u = L1L2_Proj(z.u, c1)
      
      z.v = crossprod(XTY,u)
      v = L1L2_Proj(z.v, c2) 
      
      z.w = compute.zw(Tensor, u, v, DIM[3])
      w = STPLS_proj_w(z.w,beta)
      
      obj  = - t(w)%*%z.w
      objs = c(objs, obj) 
      
      # Algorithm termination condition
      if((j > 20) && (abs(objs[j]- objs[j-1])<= err)) break 
    }
    # print(objs)
    # save optimal solution
    if(obj < opt_obj){
      modularity = mean(abs(Tensor[which(u!=0), which(v!=0), which(w!=0)]))
      opt_obj = obj
      out = list(u=u, v=v, w=w, objs=objs, opt_obj=opt_obj, modularity=modularity)
    }
  }
  # save algorithm running time
  tim = proc.time() - ptm
  out$time = tim[1]
  return(out)
}

# Compute matrix C
compute.C = function(Tensor, w, dim){
  lapply(1:dim,function(i,dat,v1){dat[,,i]*v1[i]},dat=Tensor,v1=w)->res
  Reduce("+", res) -> C
  # C = w[1]*Tensor[,,1]
  # if(dim>2){
  #   for(i in 2:dim){
  #     C = C + w[i]*Tensor[,,i]
  #   }
  # }
  return(C)
}

# Compute vector zw
compute.zw = function(Tensor, u, v, dim){
  sapply(1:dim,function(i,dat,v1,v2){t(v1)%*%dat[,,i]%*%v2},dat=Tensor,v1=u,v2=v)->zw
  zw = matrix(zw,ncol=1)
  # zw = matrix(0, nrow=dim, ncol=1)
  # for(i in 1:dim){
  #   zw[i] = t(u)%*%Tensor[,,i]%*%v
  # }
  return(zw)
}

# Projection function of w
STPLS_proj_w = function(z,beta){
  
  z[z<=0] = 0
  
  z2 = z^(1/(beta-1))
  z2[z<=0] = 0
  
  z3 = z^(beta/(beta-1))
  z3[z<=0] = 0
  
  s0 = sum(z3)^(1/beta)
  w = z2/s0
  return(w)
}

# Projection function of u and v
L1L2_Proj = function(y, c0){
  # min_x ||x-y||_2^2 s.t. ||x||_1<c and ||x||_2=1
  # ICLR 2013 Block Coordinate Descent for Sparse NMF [ncnmf]
  # CJE 2017 A Novel Sparse Penalty for Singular Value Decomposition (Algorithm 4)
  # c0 in [0,1] sparse level
  c = (sqrt(length(y))-1)*c0 + 1; # print(c)
  
  v = abs(y)
  if(c < 1){
    print("Full sparsity. Return zero vector.")
    return(0*y)
  }
  if(c == 1){
    ID = which.max(abs(y))
    y = 0*y
    y[ID] = 1
    return(y)
  }
  if(sum(v)^2/sum(v^2) < c^2){
    print("No sparsity.")
    return(y/sqrt(sum(y^2)))
  }
  
  z = sort(v,decreasing = F) #v[order(v, decreasing = F)]
  p = length(z)
  
  zmat1 = matrix(rep(z,p), nrow=p,ncol=p, byrow = FALSE)
  zmat2 = matrix(rep(z,p), nrow=p,ncol=p, byrow = TRUE)
  zmat = zmat1-zmat2
  zmat[zmat<0]=0
  i = which(colSums(zmat)^2 < c^2*colSums(zmat^2))[1]
  k_star = p - i + 1
  z_top = z[(p-k_star+1):p]
  s1 = sum(z_top)/k_star
  s2 = sum(z_top^2)/k_star
  
  # find opt eta 
  # CJE 2017 A Novel Sparse Penalty for SVD (Algorithm 4)
  eta = s1 - sqrt(c^2*(s2-s1^2)/(k_star-c^2)) # opt eta 4~5
  
  x = v-eta; x[x<0] = 0
  x = sign(y)*(x)
  x = x/sqrt(sum(x^2))
  return(x)
} 

# Update tensor 
update.tensor = function(Tensor,u,v,w){
  zw <- compute.zw(Tensor, u, v, dim(Tensor)[3])
  d <- zw*w
  mat <- u%*%t(v)
  for(k in 1:length(w)){
    Tensor[,,k] = Tensor[,,k] - d[k]*mat
  }
  
  return(Tensor)
}



