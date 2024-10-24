
source("./R/sTPLS_algorithm.R")
source("./R/identify_multi_factors.R")

# run_sTPLS: apply sTPLS to 'data_tensor' when using the optimal parameters (beta, alpha_cu, alpha_cv) and save the results in the directory "./examples/data_type". 
# If parameter selection is performed before, user could copy the folder corresponding to the optimal parameters in the directory "./examples/data_type/select_params/" to the destination path "./examples/data_type/".

# parameters
# data_tensor: input data
# beta, alpha_cu, alpha_cv: parameters in sTPLS
# nfactor: the number of identified comodules
# data_type: folder name where results are saved 

# return the results of sTPLS including weight matrices U,V,W and parameters.

run_sTPLS = function(data_tensor, alpha_cu, alpha_cv, beta, nfactor, data_type="test"){
  resfolder = paste0("./examples/", data_type, "/")
  foldername = paste0("alpha",alpha_cu,"_",alpha_cv,"_beta",beta,"_nfactor",nfactor)
  stpls_rdsname = paste0(resfolder, foldername, "/sTPLS_", foldername,".RDS")
  
  if(file.exists(stpls_rdsname)){
    cat(paste0("The results had been saved in ", resfolder, foldername, ".\nLoad the results when alpha_cu=",alpha_cu,", alpha_cv=",alpha_cv,", beta=", beta, ", nfactor=", nfactor, ".\n"))
  }else{
    print(paste0("Run sTPLS algorithm when alpha_cu=",alpha_cu,", alpha_cv=",alpha_cv,", beta=",beta,", nfactor=", nfactor, "."))
    identify_multi_factors(tensor=data_tensor, 
                           beta=beta, nfactor=nfactor, 
                           alpha_cu=alpha_cu, alpha_cv=alpha_cv, 
                           filepath=resfolder, repeatNum=5)
    print(paste0("The results are saved in ", resfolder, foldername, "."))
  }
  
  dat <- readRDS(stpls_rdsname)
  return(dat$out.sTPLS)
}

# plot_UVW: heatmaps of weight matrices U,V,W.
# parameters
# U,V,W: weight matrices of sTPLS

# return a list of heatmap objects

plot_UVW = function(U,V,W){
  uirow = which(rowSums(abs(U))>0)
  virow = which(rowSums(abs(V))>0)
  wirow = which(rowSums(W)>0)
  nfactor = ncol(U)
  uheatmap = Heatmap(U[uirow,], name = " ",
                     cluster_rows = T, cluster_columns = F,
                     show_row_names=F, show_column_names=T,
                     column_title ="U", column_labels = paste0("M",1:nfactor),
                     width = unit(5, "cm"))
  vheatmap = Heatmap(V[virow,], name = " ",
                     cluster_rows = T, cluster_columns = F,
                     show_row_names=F, show_column_names=T,
                     column_title ="V",column_labels = paste0("M",1:nfactor),
                     width = unit(5, "cm"))
  wheatmap = Heatmap(W[wirow,], name = " ", cluster_rows = T, 
                     cluster_columns = F, show_row_names=T,
                     show_column_names=T, column_title = "W",
                     row_names_side = "left", width = unit(5, "cm"),height=unit(6,"cm"),
                     column_labels = paste0("M",1:nfactor),
                     col = circlize::colorRamp2(c(0.25, 1), c("grey", "red")))
  draw(uheatmap)
  draw(vheatmap)
  draw(wheatmap)
}

determine_comodule = function(weight_list, uvwthrd_list){
  U = weight_list$U
  V = weight_list$V
  W = weight_list$W
  uthrd = uvwthrd_list$uthrd
  vthrd = uvwthrd_list$vthrd
  wthrd = uvwthrd_list$wthrd
  
  print(paste0("Using uvwthrd = ", uthrd, " and ", vthrd, " and ", wthrd, " to select comodule members."))
  zu = scale(abs(U)); zv = scale(abs(V)); zw = scale(W)
  apply(zu,2,function(x){which(x>uthrd)}, simplify = F) -> ulist
  apply(zv,2,function(x){which(x>vthrd)}, simplify = F) -> vlist
  apply(zw,2,function(x){which(x>wthrd)}, simplify = F) -> wlist
  
  cat(paste0("Numbers of selected comodule members:\n for Dim1 - ", 
             paste0(sapply(ulist,length), collapse=", "), ";\n for Dim2 - ",
             paste0(sapply(vlist,length), collapse=", "), ";\n for Dim3 - ",
             paste0(sapply(wlist,length), collapse=", ")), ".\n")
  
  sortdf = function(df, icol, decrease=F){
    order(df[,icol], decreasing=decrease) -> i
    res = df[i,]
    rownames(res) = NULL
    return(res)
  }
  
  nfactor = ncol(U)
  comodules = vector("list", nfactor)
  for(i in 1:nfactor){
    ddf = data.frame(Dim1_index = ulist[[i]],
                     Name = rownames(U)[ulist[[i]]], Weight = U[ulist[[i]],i])
    ddf = sortdf(ddf, icol = 3, decrease=T)
    
    gdf = data.frame(Dim2_index = vlist[[i]], 
                     Name = rownames(V)[vlist[[i]]], Weight = V[vlist[[i]],i])
    gdf = sortdf(gdf, icol = 3, decrease=T)
    
    cdf = data.frame(Dim3_index = wlist[[i]], 
                     Name = rownames(W)[wlist[[i]]], Weight = W[wlist[[i]],i])
    cdf = sortdf(cdf, icol = 3,decrease=T)
    
    comodules[[i]] = list(Dim1 = ddf, Dim2 = gdf, Dim3 = cdf)
  }
  return(comodules)
}

comodule_to_excel = function(comodules, destdir){
  library(openxlsx)
  nfactor = length(comodules)
  wb <- createWorkbook()
  for(i in 1:nfactor){
    sheetname <- paste0("Comodule", i)
    addWorksheet(wb, sheetname)
    writeData(wb, sheetname, x = comodules[[i]]$Dim3[,2:3], startCol = 1, startRow = 1, rowNames = FALSE)
    writeData(wb, sheetname, x = comodules[[i]]$Dim1[,2:3], startCol = 5, startRow = 1, rowNames = FALSE)
    writeData(wb, sheetname, x = comodules[[i]]$Dim2[,2:3], startCol = 9, startRow = 1, rowNames = FALSE)
  }
  
  saveWorkbook(wb, paste0(destdir, "/comodules_members.xlsx"), overwrite = TRUE)
  
}
