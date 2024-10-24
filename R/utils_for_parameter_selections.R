library(parallel)

find_entry <- function(vec, thrd, top=FALSE){
  if(top){
    sort(vec, decreasing = T, index.return = TRUE) -> tmp
    res <- tmp$ix[1:thrd]
  }else{
    i <- which(vec > thrd)
    if(length(i)>0){
      sort(vec[i], decreasing = T, index.return = TRUE) -> tmp
      res <- i[tmp$ix]
    }else{
      res <- which.max(vec)
    }
  }
  return(res)
}

compute_modularity_thresholding_byTop <- function(tensor,U,V,W,uvwthrds){
  modularity_func <- function(ithrd, uvwthrds, uvwlist, iuvwlist, dat){
    
    newvector <- function(x,xval,k){
      y <- matrix(0, nrow = length(xval), ncol = 1)
      y[x[1:k],1] <- xval[x[1:k]]
      return(y)
    }
    
    compute_zw = function(Tensor, u, v){
      zw = apply(Tensor, 3, function(x,a,b){t(a)%*%x%*%b}, a=u, b=v, simplify = T)
      return(matrix(zw,ncol = 1))
    }
    
    mapply(newvector, iuvwlist[[1]], uvwlist[[1]], k=uvwthrds[ithrd,1], SIMPLIFY = F) -> newu.list
    mapply(newvector, iuvwlist[[2]], uvwlist[[2]], k=uvwthrds[ithrd,2], SIMPLIFY = F) -> newv.list
    mapply(newvector, iuvwlist[[3]], uvwlist[[3]], k=uvwthrds[ithrd,3], SIMPLIFY = F) -> neww.list
    
    lapply(iuvwlist[[1]], function(x,k){x[1:k]}, k=uvwthrds[ithrd,1]) -> newiu.list
    lapply(iuvwlist[[2]], function(x,k){x[1:k]}, k=uvwthrds[ithrd,2]) -> newiv.list
    lapply(iuvwlist[[3]], function(x,k){x[1:k]}, k=uvwthrds[ithrd,3]) -> newiw.list
    
    mapply(function(u,v,w,iu,iv,iw){
      modularity = mean(abs(dat[iu, iv, iw]))
      return(list(modularity = modularity))}, 
      newu.list, newv.list, neww.list, 
      newiu.list, newiv.list, newiw.list,
      SIMPLIFY = F) -> res
    
    modvec <- unlist(lapply(res, function(x){x$modularity}))
    result = list(modularity=modvec)
    return(result)
  }
  
  
  nfactor <- ncol(U)
  
  iulist <- apply(abs(U), 2, order, decreasing = TRUE, simplify = F)
  ivlist <- apply(abs(V), 2, order, decreasing = TRUE, simplify = F)
  iwlist <- apply(W, 2, order, decreasing = TRUE, simplify = F)
  
  ulist <- lapply(apply(U,2,list),unlist)
  vlist <- lapply(apply(V,2,list),unlist)
  wlist <- lapply(apply(W,2,list),unlist)
  
  thrd.len = nrow(uvwthrds)
  modularity_mat <- matrix(NA, nrow=thrd.len, ncol=nfactor)
  
  cl <- makeCluster(6) 
  results <- parLapply(cl, X=1:thrd.len, modularity_func, 
                       uvwthrds = uvwthrds, 
                       uvwlist = list(u=ulist,v=vlist,w=wlist),
                       iuvwlist = list(u=iulist,v=ivlist,w=iwlist),
                       dat = tensor)
  
  t(sapply(results, function(x){x$modularity})) -> modularity_mat
  stopCluster(cl) 
  
  return(modularity_mat)
  
}


