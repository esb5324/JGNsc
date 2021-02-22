# ---------------------------
# basic functions
# ---------------------------

#' truncate the (precision) matrix
#' 
#' @param z a matrix
#' @param threshold the threshold value, where if a value's absolute value is less than this, it will be truncated to zero.
#' @return a truncated matrix
#' @export
trunc_precision <- function(z, threshold=0.01){
  z[abs(z) <threshold] <-0
  return(z)
}

#' Transform the precision matrix to partial correlation matrix
#'  
#' @param mat the precision matrix
#' @return a partial correlation matrix
#' @export
prec2partialcorr <- function(mat){ 
  denom <- diag(mat)
  n <- ncol(mat)
  out <- mat
  for (i in 1:n){
    for (j in 1:n){
      out[i,j] <- - mat[i,j]/sqrt(denom[i]*denom[j])
    }
  }
  diag(out) <- 0
  return(out)
}

#' create list of matrices according to conditions
#' 
#' @param mtx raw count matrix, genes by samples
#' @param group vector of group identities. For example: c("A","A","B","C")
#' @param min.cells integer, keep genes with at least __ cells expressed. Default is 10.
#' @param geneSet a vector of genes of interest
#' @return a list of matrices by group
#' @export
getObservedList <- function(mtx, group = NULL, min.cells = 10, geneSet = NULL){ 
  if (is.null(geneSet)){
    message("Screening all available genes /n Please make sure group variable matches the samples ID")
    geneSet = rownames(mtx)
  }  
  if (is.null(group)){
    group = rep("OneGroup", ncol(mtx))
  }
  mtx.filter = mtx[rowSums(mtx>0) > min.cells,]
  genes.use = intersect(rownames(mtx.filter), geneSet)
  mtx.filter = mtx.filter[genes.use,]
  gid = unique(group)
  obsList = list()
  for (i in 1:length(gid)){
    obsList[[i]] = mtx.filter[,group == gid[i]]
    names(obsList)[i] <- gid[i]
  }
  message("The generated List dimensions:", lapply(obsList, dim))
  return(obsList)
}

#' get the genes connected to a TARGET GENE under each condition.
#' @param gs character, 'gene selected' or 'TARGET GENE': e.g. "GAPDH", "CD4",...
#' @param pcor0 numeric, the threshold for partial correlations. values greater than pcor0 will be kept.
#' @param partcorr.list list of partial correlation matrices
#' @param conditions a vector of names for each condition. e.g.: if the list has two matrices, here conditions = c("ConditionA","ConditionB")
#' @return a table of partial correlations for genes connected to TARGET GENE under all conditions
#' @export
getCondRelated <- function(gs, pcor0 = 0.0001, partcorr.list, conditions = NULL){
  Grelated <- list()
  for (k in 1:length(partcorr.list)){
    Grelated[[k]] <- names(which(abs(partcorr.list[[k]][gs,])>pcor0))
  }
  gset.temp <- unique(unlist(Grelated))
  xtemp <- NULL
  for (k in 1:length(partcorr.list)){
    xtemp <- cbind(xtemp, partcorr.list[[k]][gset.temp,gs])
  }
  colnames(xtemp) <- gsub(" ", "_", conditions)
  return(xtemp)
}

#' Given a set of genes, get the frequencies of genes in each pathway.
#' @param meta a data.frame including the pathway information of genes
#' @param VARPATH the variable name for pathways
#' @return a table showing the frequencies of genes in each pathway
#' @export
#' 
getPathFreq <- function(meta, VARPATH ="SigmaMiniMap.Term"){
  nn <- nrow(meta)
  paths <- str_trim(unlist(strsplit(meta[,VARPATH], ";")), side = "both")
  paths <- paths[paths !=""]
  paths.count <- count(paths)
  paths.count$pct <- paths.count$freq / nn 
  return(paths.count)
}

#' plot: boxplot for a specific gene by group
#' @param countmat raw count matrix, gene by sample
#' @param group group identities of cells
#' @param gene the name of a TARGET GENE. e.g.: "MYC"
#' @return A plot
#' @export
getGENEbox <- function(countmat, group, gene){
  data = data.frame(Gene = countmat[gene,], 
                    group = group)
  datamelt <- melt(data)
  ggplot(datamelt, aes(x=group, y=value)) + geom_jitter(alpha = 0.3, color = "lightblue")+ geom_boxplot(alpha = 0.5) +
    ggtitle(gene) + xlab("") + ylab("")
}

#' plot: histogram of raw counts by group
#' @param countmat the raw count matrix , gene by sample
#' @param group a vector of group identities for cells/samples
#' @param gene the name of a TARGET GENE. e.g.: "MYC"
#' @param xlabel.pos the position of samplesize/expression rate information
#' @return a plot
#' @export
getGENEhist <- function(countmat, group, gene, xlabel.pos = 500){ 
  data <- data.frame(group = group, 
                     GENEbin = countmat[gene,] > 0,
                     GENEraw = countmat[gene,])
  groupn = table(data$group)
  pct = round(table(data$group, data$GENEbin)[,2] / groupn,2)
  gid = names(groupn)
  for (ii in 1:length(gid)){
    data$pct[data$group == gid[ii]] <- pct[ii]
  }
  p1 = ggplot(data, aes(x = GENEraw)) + geom_histogram(binwidth = 30) +
    facet_wrap(. ~ group , scales = "free_y") + geom_text(aes(x=xlabel.pos, y=min(groupn)/4, label=pct)) +
    xlab("Raw") + ylab("Frequency") + ggtitle(gene)
  print(p1)
}

#' Map the gene set of interest to the known pathways
#' @param partcorr.list a list of estimated partial correlation matrices
#' @param conditions a vector of the names of the corresponding partial correlation matrices
#' @param GeneInterest character. e.g. "MYC"
#' @param pathwayRef a data frame, including the gene and pathway information. e.g. column1 is genes, column2 is the pathways this gene belongs to according to Sigma.Mini.Map
#' @param pathwayRef_geneVariable character, the column name of gene variable
#' @param pathwayRef_pathVariable character, the column name of pathway variable
#' @param threshold the cutoff value for partial correlation, for visualization and GSEA purpose.
#' @return x.connect: the list of genes that connected to the GeneInterest under each condition
#' @return GSEA.table: the final GSEA table
#' @export
Map2Pathways <- function(partcorr.list, conditions = NULL, GeneInterest = NULL,
                         pathwayRef, pathwayRef_geneVariable, pathwayRef_pathVariable,
                         threshold =0){ 
  if (!is.data.frame(pathwayRef)){
    pathwayRef <- as.data.frame(pathwayRef)
  }
  relatedGenes <- getCondRelated(gs=GeneInterest, partcorr.list = partcorr.list, pcor0 = threshold, conditions = conditions)
  x.connect <- pathwayRef[toupper(pathwayRef[,pathwayRef_geneVariable]) %in% rownames(relatedGenes),]
  
  # constant 
  x.connect <- cbind.data.frame(x.connect, relatedGenes[toupper(x.connect[,pathwayRef_geneVariable]),])
  x.connect[,pathwayRef_pathVariable][x.connect[,pathwayRef_pathVariable] ==""] <- "Undefined;"
  x.allgenes <- pathwayRef[toupper(pathwayRef[,pathwayRef_geneVariable]) %in% toupper(rownames(partcorr.list[[1]])),]
  x.allgenes[,pathwayRef_pathVariable][x.allgenes[,pathwayRef_pathVariable] ==""] <- "Undefined;"
  
  cond.connect <- list()
  cond.freqs <- list()
  cond.vars <- gsub(" ", "_", conditions)
  for (m in 1:length(conditions)){
    cond.connect[[m]] <- x.connect[abs(x.connect[,cond.vars[m]]) >0,]
    cond.freqs[[m]] <- getPathFreq(cond.connect[[m]], VARPATH =pathwayRef_pathVariable)
  }
  all.freqs <- getPathFreq(x.allgenes, VARPATH =pathwayRef_pathVariable)
  
  # merge these proportions
  meta.merge <- all.freqs
  for (i in 1:length(cond.freqs)){
    meta.merge <- merge(x = meta.merge, cond.freqs[[i]], by = "x", all = T)
  }
  # TODO: make sure conditions and the data col names match!!!!
  colnames(meta.merge)[2:ncol(meta.merge)] <- c("AllGenes_freq","AllGenes_pct",
                                                paste(rep(cond.vars, each =2), c("freq","pct"), sep = "_"))
  meta.merge[is.na(meta.merge)] <- 0
  nconnect = c(nrow(x.allgenes), unlist(lapply(cond.connect, nrow)))
  names(nconnect) <- c("All", cond.vars)
  
  getFisherPval <- function(subg = 1){
    ft.table <- NULL
    for(ii in 1:nrow(all.freqs)){
      subpath <- all.freqs$x[ii]
      x1 = c(all.freqs[all.freqs$x == subpath,2], nrow(x.allgenes)-all.freqs[all.freqs$x == subpath,2])
      x2 = c(cond.freqs[[subg]][cond.freqs[[subg]]$x == subpath,2], nrow(cond.connect[[subg]]) - cond.freqs[[subg]][cond.freqs[[subg]]$x == subpath,2])
      if (length(x2)>0){
        contigency = rbind(x2, x1-x2)
        ft = fisher.test(contigency)
        ft.table <- rbind(ft.table,  c(round(ft$p.value,4), all.freqs$x[ii]))
      } else {
        x2 <- c(0,nrow(cond.connect[[2]]))
        contigency = rbind(x2, x1-x2)
        ft = fisher.test(contigency)
        ft.table <- rbind(ft.table,  c(round(ft$p.value,4), all.freqs$x[ii]))
      }
    }
    colnames(ft.table) <- c("Fisher_pval","x")
    ft.table <- as.data.frame(ft.table)
    ft.table$Fisher_pval <- as.numeric(ft.table$Fisher_pval)
    return(ft.table)
  }
  
  test.fisher <- lapply(1:length(conditions), getFisherPval)
  # all(test.fisher[[1]]$x ==test.fisher[[2]]$x)
  temp = do.call(cbind, test.fisher)[,- (2*(1:length(conditions)))]
  colnames(temp) <- paste("Fisher_pval",cond.vars, sep = "_")
  
  GSEA.table = cbind(meta.merge, temp)
  return(list(x.connect = x.connect,
              GSEA.table = GSEA.table))
}