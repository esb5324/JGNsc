#' RunJGNsc framework main function
#'
#' @param observed.list the list containing the matrices of K conditions. dim: genes by samples
#' @param mask.rate iterative imputation procedure
#' @param nrep number of iterations in the imputation procedure, default is 50
#' @param warm the number of warm steps in the MCMC procedure
#' @param iter the number of iterations in the MCMC procedure
#' @param min.cell the min number of cells that should a gene should express
#' @param runNetwork logical. Run the joint graphical lasso procedure or not.
#' @param l1.vec if runNetwork=T, the vector of candidate values for the tuning parameter lambda1
#' @param l2.vec if runNetwork=T, the vector of candidate values for the tuning parameter lambda2
#' @param a1 alpha parameter for Beta(alpha,beta) distribution which is the prior for non-dropouts. default is 3. The expected non-dropout rate is a1/(a1+b1) = 0.75 by default.
#' @param b1 beta parameter for Beta(alpha,beta) distribution which is the prior for non-dropouts. default is 1.
#' @param dropThreshold the threshold less than which we define a dropout, default is 0.75. set dropThreshold = 0 for UMI
#' @param AvoidIterImp Defualt if F. To allow the usage for UMI data, we allow users to turn off the iterative imputation by setting AvoidIterImp = T.
#' @return theta.star.npn: the imputed and Gaussian transformed list of matrices
#' @return if runNetwork = T, then it will return the JGL results (the precision matrices can be accessed by result$JGL$theta)
#'         , aic.table (the AIC values with corresponding tuning parameter candidates), partcorr (the list of partial correlation matrices)
#' @export

RunJGNsc <- function(observed.list, warm = 1000, iter = 5000,
                     mask.rate = 0.15, nrep = 50, min.cell = 3, runNetwork = F, l1.vec = NULL,
                     l2.vec = NULL, a1 = 2, b1 = 1, dropThreshold = 0.75, AvoidIterImp = F){
  # observed.list: the list containing the matrices of K conditions. dim: genes by samples
  # mask.rate: iterative imputation procedure
  # nrep: number of iterations in the imputation procedure
  # set dropThreshold = 0 for UMI
  if (dropThreshold == 0){
    AvoidIterImp = T
    print("UMI Mode \n")
  }
  zip.list <- lapply(observed.list, JGNsc_cont_cpp, warm = warm, iter = iter, minCell= min.cell, dropThreshold = dropThreshold, a1 = a1, b1 = b1)
  for(kk in 1:length(zip.list)){
    colnames(zip.list[[kk]]$y.impute) = colnames(observed.list[[kk]])
    gnames = rownames(observed.list[[kk]])[zip.list[[kk]]$keep.gene > min.cell]
    rownames(zip.list[[kk]]$y.impute) = gnames
  }

  g.common <- Reduce(intersect, lapply(zip.list, function(x){rownames(x$y.impute)}))
  theta.star <- lapply(zip.list, function(x){
    y = x$y.impute[g.common,]
    return(y)
  }) #genes by samples
  observed.list2 <- lapply(observed.list, function(x){
    y = x[g.common,]
  })
  if(AvoidIterImp){
    # UMI mode
    # theta.star.npn = theta.star
    theta.star.t <- lapply(theta.star, t)
    theta.star.npn <- lapply(theta.star.t, huge.npn)  # samples by genes
  } else {
    for (k in 1:length(zip.list)){
      nr = nrow(observed.list2[[k]])
      nc = ncol(observed.list2[[k]])
      cat("cond",k,"\n")
      imp.sum <-0
      mask.sum <-0
      x = t(theta.star[[k]])
      for(kk in 1:nrep){
        rmat <- matrix(data = runif(nr*nc) > mask.rate,nrow = nr, ncol = nc) # keep 1-mask.rate,
        mask <- t((observed.list2[[k]]!=0) | rmat)
        x.mask = x*mask
        x.mask.mcimpute = mcImpute_cpp(x.mask, preprocess = T)
        current.imp <- t(x.mask.mcimpute$data)
        x = t(current.imp)
        imp.sum <- imp.sum + current.imp
      }
      imp.mean <- imp.sum / nrep
      rownames(imp.mean) <- rownames(observed.list2[[k]])
      colnames(imp.mean) <- colnames(observed.list2[[k]])
      if (is.null(rownames(observed.list2[[k]]))){
        rownames(imp.mean) <- paste("gene",1:ncol(x), sep = "")
      }
      if (is.null(colnames(observed.list2[[k]]))){
        colnames(imp.mean) <- paste("sample",1:nrow(x), sep = "")
      }
      theta.star[[k]] <- imp.mean
    }
    theta.star.t <- lapply(theta.star, t)
    theta.star.npn <- lapply(theta.star.t, huge.npn)
  }

  # joint lasso:
  if(runNetwork){
    if (is.null(l1.vec)){
      l1.vec <- seq(1,30,by=2)/100
    }
    if (is.null(l2.vec)){
      l2.vec <- seq(1,30,by=2)/100
    }
    jgnsc.aic <- NULL
    for (lam1 in l1.vec){
      for(lam2 in l2.vec){
        cat("tuning parameter:", lam1,", ", lam2,"\n")
        # mat.list.t.gau is: samples by genes
        tps <- AIC_select(mat.list.t.gau = theta.star.npn, lam1 = lam1, lam2 = lam2)
        jgnsc.aic <- rbind(jgnsc.aic, tps)
      }
    }
    JGL.res <- JGL(theta.star.npn, lambda1 = jgnsc.aic[which.min(jgnsc.aic[,3]),1],
                   lambda2 = jgnsc.aic[which.min(jgnsc.aic[,3]),2], return.whole.theta = T)
    partcorr <- lapply(JGL.res$theta, prec2partialcorr)
  } else {
    JGL.res <- NULL
    jgnsc.aic <- NULL
    partcorr <- NULL
  }
  return(list(theta.star.npn = theta.star.npn,
              theta.star.t = theta.star.t,
              JGL = JGL.res,
              aic.table = jgnsc.aic,
              partcorr = partcorr))
}


