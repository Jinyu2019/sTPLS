

res = readRDS("./data/GDSC_data_XY.RDS")
Xlist = res$Xlist
Ylist = res$Ylist
Xfeatures = res$Xfeatures
Yfeatures = res$Yfeatures
Conditions = res$Conditions

p = dim(Xlist[[1]])[2]
q = dim(Ylist[[1]])[2]
m = length(Xlist)

remove_outlier_cor = function(x,y,thrd=3){
  # x, y are matrix, row is sample
  p = ncol(x); q = ncol(y)
  # remove outliers in each column vector of x and y
  xirow = apply(abs(scale(x)), 2, function(x){which(x>thrd)}, simplify = F)
  yirow = apply(abs(scale(y)), 2, function(x){which(x>thrd)}, simplify = F)
  
  cormat = matrix(0, p,q)
  x1=x; y1=y
  for(i in 1:p){
    x1[xirow[[i]],i] <- NA
  }
  for(i in 1:q){
    y1[yirow[[i]],i] <- NA
  }
  for(id in 1:p){
    for(ig in 1:q){
      cor(x1[,id], y1[,ig], use = "complete") -> cormat[id,ig]
    }
  }
  return(cormat)
}

data.tensor = array(0, c(p,q,m))
dimnames(data.tensor)[[1]] = Xfeatures
dimnames(data.tensor)[[2]] = Yfeatures
dimnames(data.tensor)[[3]] = Conditions

for(i in 1:m){
  data.tensor[,,i] = remove_outlier_cor(Xlist[[i]], Ylist[[i]], thrd=3)
  corxy = data.tensor[,,i]
  id = which.max(apply(corxy,1,max))
  ig = which.max(corxy[id,])
  plot(Xlist[[i]][,id], Ylist[[i]][,ig], 
       xlab = Xfeatures[id], 
       ylab = Yfeatures[ig], 
       main = paste0(Conditions[i],":", round(corxy[id,ig],digits = 2),
                     " vs ", round(cor(Xlist[[i]][,id],Ylist[[i]][,ig],
                                       use = "complete"), digits = 2)))
  
  id = which.min(apply(corxy,1,min))
  ig = which.min(corxy[id,])
  
  plot(Xlist[[i]][,id], Ylist[[i]][,ig], 
       xlab = Xfeatures[id], 
       ylab = Yfeatures[ig], 
       main = paste0(Conditions[i],":", round(corxy[id,ig],digits = 2),
                     " vs ", round(cor(Xlist[[i]][,id],Ylist[[i]][,ig],
                                       use = "complete"), digits = 2)))
}
saveRDS(data.tensor, file = "./data/GDSC_data_tensor.RDS")



