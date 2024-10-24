

# plot_comodule_member: show composition of comodules.
# parameters
# tensor: input data
# comodules: lists recording the indexes of comodule members corresponding to three dimension of input tensor

plot_comodule_member = function(tensor, comodules){
  n1 <- dim(tensor)[1]
  n2 <- dim(tensor)[2]
  n3 <- dim(tensor)[3]
  nfactor = length(comodules)
  
  tu = matrix(0, nrow = n1, ncol = nfactor)
  tv = matrix(0, nrow = n2, ncol = nfactor)
  tw = matrix(0, nrow = n3, ncol = nfactor)
  rownames(tu) = dimnames(tensor)[[1]]
  colnames(tu) <- paste0("M",1:nfactor)
  rownames(tv) = dimnames(tensor)[[2]]
  colnames(tv) <- paste0("M",1:nfactor)
  rownames(tw) = dimnames(tensor)[[3]]
  colnames(tw) <- paste0("M",1:nfactor)
  
  all.i1 = c();  all.i2 = c();  all.i3 = c();
  for(ifactor in 1:nfactor){
    d1list = comodules[[ifactor]]$Dim1$Dim1_index
    d2list = comodules[[ifactor]]$Dim2$Dim2_index
    d3list = comodules[[ifactor]]$Dim3$Dim3_index
    
    tu[d1list, ifactor] = 1
    tv[d2list, ifactor] = 1
    tw[d3list, ifactor] = 1
    all.i1 = union(all.i1, d1list)
    all.i2 = union(all.i2, d2list)
    all.i3 = union(all.i3, d3list)
  }
  
  col_fun = colorRamp2(c(0, 1), c("white", "black"))
  font.size = 10
  
  uobj = Heatmap(t(tu[all.i1, ]), cluster_rows = F, cluster_columns = F, 
                 name = "U",col = col_fun, 
                 column_title = "U", 
                 show_heatmap_legend = FALSE, show_column_names = F,
                 column_names_gp = gpar(fontsize = font.size),
                 row_names_gp = gpar(fontsize=font.size), width = 4, height=unit(6,"cm"),
                 border_gp = gpar(col = "black", lty = 1))
  
  vobj = Heatmap(t(tv[all.i2, ]), cluster_rows = F, cluster_columns = F, 
                 name = "V",col = col_fun, 
                 column_title = "V", 
                 show_heatmap_legend = FALSE, show_column_names = F,
                 column_names_gp = gpar(fontsize = font.size),
                 row_names_gp = gpar(fontsize=font.size), width = 4, height=unit(6,"cm"),
                 border_gp = gpar(col = "black", lty = 1))
  
  wobj = Heatmap(t(tw[all.i3,]), cluster_rows = F, cluster_columns = F, 
                 name = "W",col = col_fun, 
                 column_title = "W", 
                 show_heatmap_legend = FALSE, show_column_names = T,
                 column_names_gp = gpar(fontsize = font.size),
                 row_names_gp = gpar(fontsize=font.size), width = 2, height=unit(6,"cm"),
                 border_gp = gpar(col = "black", lty = 1))
  ht_list = uobj+vobj+wobj
  draw(ht_list)
}

# plot_comodule_heatmap: show heatmaps of identified comodules
# parameters
# tensor: input data
# weight_list: list(U,V,W) where U,V,W are weight matrices of sTPLS
# comodules: lists recording the indexes of comodule members corresponding to three dimension of input tensor
# mindex: vector recording the indexes of comodules which are shown 
# includeAllCond: if includeAllCond=T, demonstrate heatmaps across all conditions; if includeAllCond=F, demonstrate heatmaps across selected conditions in comodule

plot_comodule_heatmap = function(tensor, weight_list, comodules, mindex = 1, includeAllCond = T){
  ave_num_dim1 = median(unlist(lapply(comodules,function(x){nrow(x[[1]])})))
  ave_num_dim2 = median(unlist(lapply(comodules,function(x){nrow(x[[2]])})))
  
  if(ave_num_dim2/ave_num_dim1>=5){
    # transpose the first and second dimensions of tensor to show pretty heatmaps
    ttensor = array(NA, dim = c(dim(tensor)[2],dim(tensor)[1],dim(tensor)[3]))
    for(i in 1:dim(tensor)[3]){
      ttensor[,,i] = t(tensor[,,i])
    }
    dimnames(ttensor)[[1]] = dimnames(tensor)[[2]]
    dimnames(ttensor)[[2]] = dimnames(tensor)[[1]]
    dimnames(ttensor)[[3]] = dimnames(tensor)[[3]]
    
    U = weight_list$V[,mindex,drop=F]
    V = weight_list$U[,mindex,drop=F] 
    W = weight_list$W[,mindex,drop=F]
    d1list = d2list = d3list = vector("list", length(mindex))
    for(i in 1:length(mindex)){
      ifactor = mindex[i]
      d1list[[i]] = comodules[[ifactor]]$Dim2$Dim2_index
      d2list[[i]] = comodules[[ifactor]]$Dim1$Dim1_index
      d3list[[i]] = comodules[[ifactor]]$Dim3$Dim3_index
    }
  }else{
    ttensor = tensor
    U = weight_list$U[,mindex,drop=F]
    V = weight_list$V[,mindex,drop=F] 
    W = weight_list$W[,mindex,drop=F]
    d1list = d2list = d3list = vector("list", length(mindex))
    for(i in 1:length(mindex)){
      ifactor = mindex[i]
      d1list[[i]] = comodules[[ifactor]]$Dim1$Dim1_index
      d2list[[i]] = comodules[[ifactor]]$Dim2$Dim2_index
      d3list[[i]] = comodules[[ifactor]]$Dim3$Dim3_index
    }
  }
  
  ntype <- dim(ttensor)[3]
  if(ntype>12){
    CondCol <- scico::scico(ntype, palette = 'lapaz')
  }else{
    CondCol <- brewer.pal(12, "Set3")[1:ntype]
  }
  names(CondCol) <- dimnames(ttensor)[[3]]
  
  sortx <- function(x,ix, decreasing = T){
    order(x[ix], decreasing = decreasing) -> tmp
    return(ix[tmp])
  }
    
  for(i in 1:length(mindex)){
    sortx(x = U[,i], ix = d1list[[i]], decreasing = F) -> iu
    sortx(x = V[,i], ix = d2list[[i]], decreasing = F) -> iv
    
    if(includeAllCond){
      iw = order(W[,i],decreasing = T)
    }else{
      sortx(x = W[,i], ix = d3list[[i]], decreasing = T) -> iw
    }
    newx = tensor2mat(ttensor[iu,iv,iw, drop = FALSE], isscale = F)
    
    CondType <- factor(unlist(lapply(dimnames(ttensor)[[3]][iw], rep, length(iv))), 
                      levels = dimnames(ttensor)[[3]][iw])
    topCol <- HeatmapAnnotation(df = data.frame(Cond = CondType),
                                col = list(Cond = CondCol),
                                annotation_height = unit(0.5,"cm"),
                                annotation_name_rot = 0, 
                                show_legend = TRUE)
    
    comodule_heatmap = Heatmap(newx, name = " ", top_annotation = topCol,
                               show_row_names = FALSE, show_column_names = FALSE, 
                               cluster_rows = F, cluster_columns = F, 
                               row_title = " ", 
                               column_title = paste0("Comodule ",mindex[i]," (",
                                                     paste0(c(length(iu),length(iv),
                                                              length(d3list[[i]])),
                                                            collapse=", "),")"),
                               column_title_gp = gpar(fontsize = 12), 
                               row_title_gp = gpar(fontsize = 12),
                               column_split = CondType, column_gap = unit(0.15,"cm"),
                               show_heatmap_legend = T, 
                               col = colorRamp2(c(-0.5, 0, 0.5), c("blue","white", "red")),
                               width = unit(10,"cm"), height = unit(6,"cm"))
    draw(comodule_heatmap, adjust_annotation_extension = F)
    
  }
  
}

tensor2mat <- function(tendat, isscale = F){
  x <- tendat[,,1]
  if(isscale){
    if(dim(tendat)[3]>1){
      for(i in 2:dim(tendat)[3]){
        cbind(x, scale(tendat[,,i])) -> x
      }
    }
    dat <- x
  }else{
    if(dim(tendat)[3]>1){
      for(i in 2:dim(tendat)[3]){
        cbind(x, tendat[,,i]) -> x
      }
    }
    dat <- x
  }
  
  return(dat)
}
