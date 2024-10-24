# generate simulation data
rm(list = ls(all=TRUE))

generate_x = function(n,p,q, xclusters,yclusters,signs,rou=0.9){
  x = array(rnorm(n*p), c(n,p))
  y = array(rnorm(n*q), c(n,q))
  for(i in 1:length(xclusters)){
    sigma = matrix(rou,length(xclusters[[i]]),length(xclusters[[i]]))
    diag(sigma) = 1
    ave = rep(2,length(xclusters[[i]]))
    x[,xclusters[[i]]] = mvrnorm(n, ave, sigma)
    beta = matrix(signs[i]*1/length(xclusters[[i]]), length(xclusters[[i]]),
                  length(yclusters[[i]]))
    y[,yclusters[[i]]] = x[,xclusters[[i]]]%*%beta+
      matrix(rnorm(n*length(yclusters[[i]]), mean = 1, sd = 1), nrow=n) 
  }
  res=list(x,y,xclusters,yclusters)
  return(res)
}

isim = 100
data.rou = 0.3
p = 120;q = 150; nTypes = 4; n=50;

set.seed(isim)
dat1 = generate_x(n,p,q,xclusters = list(v1=c(1:40),v2=c(81:120)),
                  yclusters = list(v1=c(1:50),v2=c(101:150)),signs=c(-1,-1),rou=data.rou)
dat2 = generate_x(n,p,q,xclusters = list(v1=c(81:120)),
                  yclusters = list(v1=c(1:50)),signs=c(1,1),rou=data.rou)
dat3 = generate_x(n,p,q,xclusters = list(v1=c(1:40),v2=c(41:80)),
                  yclusters = list(v1=c(1:50),v2=c(51:100)),signs=c(-1,1),rou=data.rou)
dat4 = generate_x(n,p,q,xclusters = list(v1=c(1:40),v2=c(41:80),v1=c(81:120)),
                  yclusters = list(v1=c(1:50),v2=c(51:100),v1=c(101:150)),
                  signs=c(-1,1,-1),rou=data.rou)

Xlist = list(dat1[[1]],dat2[[1]],dat3[[1]],dat4[[1]])
Ylist = list(dat1[[2]],dat2[[2]],dat3[[2]],dat4[[2]])
Xfeatures = paste0("X", 1:p)
Yfeatures = paste0("Y", 1:q)
Conditions = paste0("C", 1:nTypes)

# real comodules
tmodules = list()
tmodules[[1]] = list(v1=dat1[[3]]$v1,v2=dat1[[4]]$v1+p,v3=c(1,3,4))
tmodules[[2]] = list(v1=dat3[[3]]$v2,v2=dat3[[4]]$v2+p,v3=c(3,4))
tmodules[[3]] = list(v1=dat1[[3]]$v2,v2=dat1[[4]]$v2+p,v3=c(1,4))
tmodules[[4]] = list(v1=dat2[[3]]$v1,v2=dat2[[4]]$v1+p,v3=2)

res = list(Xlist = Xlist, Ylist = Ylist, true_comodules = tmodules,
           Xfeatures = Xfeatures, Yfeatures = Yfeatures, 
           Conditions = Conditions)
saveRDS(res, file = "./data/simulated_data_XY.RDS")



