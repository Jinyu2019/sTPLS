# identify_multi_factors: apply sTPLS to input tensor and identify comodules
# parameters
# tensor: input tensor data
# alpha_cu, alpha_cv, beta: respectively for c_u, c_v, beta in the objective function of sTPLS.
# nfactor: the number of identified comodules
# repeatNum: the algorithm is executed in multiple times and finally select the best one
# filepath: the directory where the results of sTPLS are saved
# 
# return: 
# sTPLS_alpha*_*_beta*_nfactor*.RDS; Weights_alpha*_*_beta*_nfactor*.pdf

identify_multi_factors <- function(tensor, beta, alpha_cu, alpha_cv, nfactor, repeatNum=5, filepath="./"){
  foldername <- paste0(filepath, "alpha",alpha_cu,"_",alpha_cv,"_beta",beta,"_nfactor",nfactor,"/")
  if(!dir.exists(foldername)){dir.create(foldername)}
  # start.time <- Sys.time()
  out.sTPLS <- sTPLS(Tensor=tensor, nfactor=nfactor, 
                     c1=alpha_cu, c2=alpha_cv, beta=beta,
                     initNum=repeatNum, niter=100, 
                     seed0=1, err=0.0001)
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # print(time.taken)
  
  ## save results
  params <- list(alpha_cu = alpha_cu, alpha_cv = alpha_cv, beta = beta, nfactor = nfactor)
  filename <- paste0(foldername, "sTPLS_alpha",alpha_cu,"_",alpha_cv,"_beta", beta,"_nfactor", nfactor, ".RDS")
  saveRDS(list(out.sTPLS = out.sTPLS, params = params), file = filename)
  
  W <- out.sTPLS$W
  U <- out.sTPLS$U
  V <- out.sTPLS$V
  objs <- out.sTPLS$d_iter
  
  ## output results into an pdf
  #  objs; U; V; W.
  pdf(paste0(foldername, "Weights_alpha", alpha_cu,"_",alpha_cv,"_beta", beta,"_nfactor", nfactor,".pdf"), 
      onefile = TRUE, width = 16)
  for(i in 1:nfactor){
    q1 <- generate_ggplot_pointline(x = 1:nrow(U), y = U[,i], labs = c("Dim1", "Weights"), isline = FALSE)
    q1_2 <- generate_ggplot_pointline(x = 1:nrow(U), y = sort(U[,i]), labs = c("Dim1", "Weights"), isline = FALSE)
    q2 <- generate_ggplot_pointline(x = 1:nrow(V), y = V[,i], labs = c("Dim2", "Weights"), isline = FALSE)
    q2_2 <- generate_ggplot_pointline(x = 1:nrow(V), y = sort(V[,i]), labs = c("Dim2", "Weights"), isline = FALSE)
    q3 <- generate_ggplot_pointline(x = 1:nrow(W), y = W[,i], labs = c("Dim3", "Weights"),isline = FALSE)
    q4 <- generate_ggplot_pointline(x = 1:length(objs[[i]]), y = objs[[i]], labs = c("Iteration", "Objective value"), 
                                    isline = TRUE, title = paste0("ifactor =",i))
    grid.arrange(q4, q1,q1_2,q2,q2_2,q3, ncol = 3, nrow = 2)
  }
  
  dev.off()
  
}

generate_ggplot_pointline <- function(x, y, labs, isline = TRUE, title = ""){
  res <- data.frame(x = x, y = y)
  if(!isline){
    p <- ggplot(res, aes(x=x,y=y))+geom_point(size = 1)+
      xlab(labs[1])+ylab(labs[2])+ggtitle(title) 
  }else{
    p <- ggplot(res, aes(x=x,y=y))+geom_point(size = 1)+geom_line()+
      xlab(labs[1])+ylab(labs[2])+ggtitle(title) 
  }
  
  q <- p +
    theme(axis.text.x = element_text(size = 10,colour = "black"),
          axis.text.y = element_text(size = 10,colour = "black"),
          plot.title = element_text(size = 12,colour = "black",hjust = 0.5),
          axis.title.x = element_text(size = 11,colour = "black", vjust = 0.1),
          axis.title.y = element_text(size = 11,colour = "black", vjust = 1.5))
  
  return(q)
}