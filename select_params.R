

source("./R/sTPLS_algorithm.R")
source("./R/identify_multi_factors.R")

# select_params: selecting the optimal parameters (beta, alpha_cu, alpha_cv) used in sTPLS.
# parameters
# data_tensor: input data
# parameter_list: data.frame(beta, alpha_cu, alpha_cv, nfactor) recording candidates for parameters (alpha_cu, alpha_cv, beta)
# uvwthrd_list: data.frame(uthrd, vthrd, wthrd) recording candidates for the number of selected members in comodules.
# data_type: folder name where results are saved 
#
# return the optimal parameters

select_params = function(data_tensor, parameter_list, uvwthrd_list, data_type="test"){
  resfolder = paste0("./examples/", data_type, "/select_params/")
  if(!dir.exists(resfolder)){
    dir.create(resfolder)
  }
  
  nparams = nrow(parameter_list)
  for(iparam in 1:nparams){
    beta = parameter_list[iparam, 1]
    alpha_cu = parameter_list[iparam, 2]
    alpha_cv = parameter_list[iparam, 3]
    print(paste0("using parameters: alpha_cu=",alpha_cu,", alpha_cv=",alpha_cv,", beta=",beta,", nfactor=", nfactor))
    rdspath = paste0(resfolder, "alpha",alpha_cu,"_",alpha_cv,"_beta",beta,"_nfactor", nfactor) 
    if(!dir.exists(rdspath)){
      identify_multi_factors(tensor=data_tensor, 
                             beta=beta, nfactor=nfactor,
                             alpha_cu=alpha_cu, alpha_cv=alpha_cv, 
                             filepath=resfolder, repeatNum=5)
    }
  }
  
  print(paste0("All results have been saved in ", resfolder))
  
  # select parameters (alpha_cu, alpha_cv, beta). 
  source("./R/utils_for_parameter_selections.R")
  
  paramresname = paste0(resfolder, "modularity_score_for_parameter_selection.RDS")
  if(!file.exists(paramresname)){
    modularitymat <- list()
    nparams = nrow(parameter_list)
    nthrds = nrow(uvwthrd_list)
    
    for(iparam in 1:nparams){
      beta = parameter_list[iparam, 1]
      alpha_cu = parameter_list[iparam, 2]
      alpha_cv = parameter_list[iparam, 3]
      nfactor = parameter_list[iparam, 4]
      # print(paste0("alpha",alpha_cu,"_",alpha_cv,"_beta",beta,"_nfactor", nfactor))
      foldername = paste0(resfolder, "alpha", alpha_cu, "_", alpha_cv, "_beta", beta, 
                          "_nfactor", nfactor,"/") 
      filename = paste0("sTPLS_alpha", alpha_cu, "_", alpha_cv, "_beta", beta, 
                        "_nfactor", nfactor, ".RDS")
      res = readRDS(paste0(foldername, filename)) 
      modularitymat[[iparam]] = compute_modularity_thresholding_byTop(tensor=data_tensor,
                                                                      U=res$out.sTPLS$U,
                                                                      V=res$out.sTPLS$V,
                                                                      W=res$out.sTPLS$W,
                                                                      uvwthrds=uvwthrd_list)
    }
    
    modmat = sapply(modularitymat, function(x){apply(x,1,median)}) 
    y.values = apply(modmat,2,median)
    iparam = which.max(y.values)
    opt_params = parameter_list[iparam,]
    
    params_select_res = list(parameter_list = parameter_list, 
                             uvwthrd_list = uvwthrd_list,
                             modularity_score = modularitymat, 
                             optimal_parameter = opt_params)
    saveRDS(params_select_res, file = paramresname)
    
  }
  
  params_select_res = readRDS(paramresname)
  opt_params = params_select_res$optimal_parameter
  print(paste0("Select the optimal parameters: alpha_cu=", opt_params[2], 
               ", alpha_cv=", opt_params[3], ", beta=", opt_params[1],
               ", nfactor =", opt_params[4]))
  modularitymat = params_select_res$modularity_score
  modmat = sapply(modularitymat, function(x){apply(x,1,median)}) 
  # row is uvwthrd index; column is parameter index.
  colnames(modmat) = 1:ncol(modmat)
  rownames(modmat) = 1:nrow(modmat)
  y.values = apply(modmat,2,median)
  
  res = data.frame(x = 1:length(y.values), y = y.values)
  p = ggplot(res, aes(x,y))+geom_point(size = 1)+geom_line(linewidth = 0.5)+
    xlab("Parameter index")+ylab("Modularity")+ggtitle("") 
  fontsize = 12
  q = p + scale_x_continuous(breaks = seq(1,length(y.values),by=2)) + theme_bw() +
    theme(axis.text.x = element_text(size = fontsize,colour = "black"),
          axis.text.y = element_text(size = fontsize,colour = "black"),
          plot.title = element_text(size = fontsize+3,colour = "black",hjust = 0.5),
          axis.title.x = element_text(size = fontsize,colour = "black", vjust = 0.1),
          axis.title.y = element_text(size = fontsize,colour = "black", vjust = 1.5))
  print(q)
  ggsave(paste0(resfolder, "/parameter_selection.png"), width = 6, height = 3, dpi = 600, plot = q)
  
  return(opt_params)
}



