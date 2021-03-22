#'---------------------------------------
#' Tuning parameter selection by AIC
#'
#' @param mat.list.t.gau the list of matrices, which contains the Gaussian distributed counts
#' @param lam1 numeric, Tuning parameter for the sparsity
#' @param lam2 numeric, Tuning parameter for the similarities across the conditions
#' @param returnJGL logical, whether return the JGL model results. default is F.
#' @return a list of 1) a vector of tuning parameter values and the AIC value; 2) the JGL result object
#' @export
AIC_select <- function(mat.list.t.gau, lam1 = 0.1, lam2 = 0.1, returnJGL = F){
  res.S <- lapply(mat.list.t.gau, cov)
  res.jgl <- JGL::JGL(mat.list.t.gau, lambda1 = lam1, lambda2 = lam2, return.whole.theta = T)
  aic <- 0; bic <- 0; ebic <- 0
  p <- ncol(res.jgl$theta[[1]])
  for (k in 1:length(mat.list.t.gau)){
    nk <- nrow(mat.list.t.gau[[k]])
    traceST <- sum(diag(res.S[[k]] %*% res.jgl$theta[[k]]))
    Ek <- sum(! res.jgl$theta[[k]] == 0)
    detT <- det(res.jgl$theta[[k]])
    aick <- (traceST*nk - log(detT)*nk + 2*Ek)/1e4
    aic <- aic + aick
    # bick <- (traceST*nk - log(detT)*nk + log(nk)*Ek)/1e4
    # bic <- bic + bick
    # ebick <- (traceST*nk - log(detT)*nk + log(nk)*Ek + 4/2*Ek*log(p))/1e4
    # ebic <- ebic + ebick
  }
  res <- c(lam1, lam2, aic)
  if (returnJGL){
    return(list(res,
                res.jgl))
  }else{
    return(res)
  }
}

#' Derive the Joint Graphical Lasso estimation result by AIC tuning parameter selection
#'
#' @param GauList a list of gaussian distributed matrices
#' @param l1vec a vector of candidate values for lambda1, if NULL, the values will be 1:5/20
#' @param l2vec a vector of candidate values for lambda2, if NULL, the values will be 1:5/50
#' @return aic.vec: a table of AIC values ; jgl.res: the JGL result object
#' @export
getJGLTuningParamResult <- function(GauList, l1vec = NULL, l2vec = NULL){
  # Gaulist: samples by genes
  if (is.null(l1vec)){
    l1vec = 1:5/20
  }
  if (is.null(l2vec)){
    l2vec = 1:5/50
  }
  aic.vec = NULL
  jgl.res = NULL
  for (l1 in l1vec){
    for (l2 in l2vec){
      cat("Tuning parameter l1=", l1,", l2=", l2, "\n")
      t1 <- proc.time()
      tps <- AIC_select(mat.list.t.gau = GauList, lam1 = l1, lam2 = l2, returnJGL = T)
      t2 <- proc.time()
      message(t2 -t1, "\n")

      if (is.null(aic.vec)){
        jgl.res <- tps[[2]]
      } else if (tps[[1]][3] < min(aic.vec)){
        jgl.res <- tps[[2]]
      }
      aic.vec <- rbind(aic.vec, tps[[1]])
    }
  }
  return(list(aic.vec = aic.vec,
              jgl.res = jgl.res))
}

#'---------------------------------------------
#' The StARS tuning parameter selection method from package 'huge'. Adapted to JGL models.
#'
#' @param inputlist a list of matrices with Gaussian distributed counts
#' @param lambda a vector of candidate tuning parameter values. If NULL, the values will be generated according to the data.
#' @param nlambda interger, number of candidate lambda values
#' @param lambda.min.ratio numeric, the ratio of lambda.min / lambda.max, if NULL, it willbe set to 0.1
#' @param stars.subsample.ratio.vec vector, the subsample ratios. If NULL, it will be decided by data
#' @param rep.num integer, the number of subsampling
#' @param stars.thresh the threshold value in StARS
#' @param lam1init the initial value of lambda1 in JGL. This is set when try to optimize over lambda2
#' @param lam2init the initial value of lambda2 in JGL. This is set when try to optimize over lambda1
#' @param verbose logical, whether show processing messages.
#' @return lambda: a vector of candidate lambda used. opt.lambda: the selected lambda value
#' @export
huge_stars_joint <- function(inputlist, lambda = NULL, nlambda = 10, lambda.min.ratio = NULL,
                             stars.subsample.ratio.vec = NULL, rep.num = 20, stars.thresh = 0.1,
                             lam2init = NULL, lam1init = NULL, verbose = T, ...){
  commong <- Reduce(intersect, lapply(inputlist, colnames))
  inputlist <- lapply(inputlist, function(x){
    x[,commong]
  })

  d <- ncol(inputlist[[1]]) #number of genes
  if (is.null(lambda)){
    templam <- lapply(inputlist, function(x){
      x = scale(x)
      S = cor(x)
      if (is.null(lambda)) {
        if (is.null(nlambda))
          nlambda = 10
        if (is.null(lambda.min.ratio))
          lambda.min.ratio = 0.1
        lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
        lambda.min = lambda.min.ratio * lambda.max
        lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
      }
      return(lambda)
    })
    ltemp <- round(unlist(templam),2)
    lambda <- unique(ltemp)[order(unique(ltemp))]
    nlambda <- length(lambda)
  } else {
    nlambda <- length(lambda)
  }

  est.merge = list()
  for(i in 1:nlambda){
    est.merge[[i]] <- list()
    for (k in 1:length(inputlist)){
      est.merge[[i]][[k]] =  Matrix::Matrix(0,d,d) #matrix(0,d,d)#
    }
  }

  if(is.null(stars.subsample.ratio.vec)){
    stars.subsample.ratio.vec <- vector(length = length(inputlist))
    for (k in 1:length(inputlist)){
      nk <- nrow(inputlist[[k]])
      if(nk>144) stars.subsample.ratio = 10*sqrt(nk)/nk
      if(nk<=144) stars.subsample.ratio = 0.8
      stars.subsample.ratio.vec[k] <- stars.subsample.ratio
    }
  }

  for(i in 1:rep.num){
    if(verbose){
      mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/rep.num), "%"), collapse="")
      cat(mes, "\r")
      flush.console()
    }

    templist <- inputlist
    for (k in 1:length(inputlist)){
      ind.samplek = sample(c(1:nk), floor(nk*stars.subsample.ratio.vec[k]), replace=FALSE)
      templist[[k]] <- as.matrix(templist[[k]][ind.samplek,])
    }

    for(j in 1:nlambda){
      cat(j,"... lambda1 =", lambda[j], "\n ")
      if(is.null(lam1init)){
        tmp = try(JGL(templist, lambda1 = lambda[j], lambda2 = lam2init, return.whole.theta = T), silent = T)
      } else if(is.null(lam2init)){
        tmp = try(JGL(templist, lambda1 = lam1init, lambda2 = lambda[j], return.whole.theta = T), silent = T)
      }

      # if (!is.null(tmp)){ # will trigger error when tmp=NULL
      for (k in 1:length(inputlist)){
        if (!is.null(tmp$theta[[k]])){ #just in case sometimes, there's error of JGL...
          temp01 <- tmp$theta[[k]]
          est.merge[[j]][[k]] = est.merge[[j]][[k]] + abs(sign(temp01))
        }
      }
    }
    try(rm(ind.samplek,tmp, templist, temp01), silent = T)
    gc()
  }

  if(verbose){
    mes = "Conducting Subsampling....done.                 "
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }

  est.variability = rep(0,nlambda)
  for(i in 1:nlambda){
    for (k in 1:length(inputlist)){
      est.merge[[i]][[k]] = est.merge[[i]][[k]]/rep.num
      est.variability[i] = est.variability[i] + 4*sum(est.merge[[i]][[k]]*(1-est.merge[[i]][[k]]))/(d*(d-1))
    }

  }
  # -----------
  # est.opt.index = max(which.max(est.variability >= stars.thresh)[1]-1,1) # why >=? shouldn't it be <=?
  # opt.lambda = lambda[est.opt.index]
  est.opt.index2 = max(which.max(est.variability <= stars.thresh)[1]-1,1) # why >=? shouldn't it be <=?
  opt.lambda2 = lambda[est.opt.index2]
  return(list(lambda = lambda,
              opt.lambda = opt.lambda2))
}

#'---------------------------------------------
#' Tuning parameter selection in JGL by StARS
#'
#' @param inputlist a list of matrices with Gaussian distributed counts
#' @param nlamb integer, specify the number of tuning parameter candidate
#' @param rep.num integer, the number of subsampling to perform, default is 15
#' @param lambvec vector, the candidate tuning parameter values for both lambda1 and lambda2
#' @return stars.l1: the stars tuning parameter selection result when fixing lambda2
#' @return stars.l2: the stars tuning parameter selection result when fixing lambda2
#' @return res: the JGL result ;
#' @return runtime: run time
#' @export
stars.select <- function(inputlist, nlamb = 15, rep.num = 15, lambvec = NULL){
  tt0 = proc.time()
  stars.l1 <- huge_stars_joint(inputlist, lambda = lambvec, nlambda = nlamb, lambda.min.ratio = NULL, stars.subsample.ratio.vec = NULL,
                               rep.num = rep.num, stars.thresh = 0.1, lam2init = 0.05, lam1init = NULL, verbose = T)
  stars.l2 <- huge_stars_joint(inputlist, lambda = lambvec, nlambda = nlamb, lambda.min.ratio = NULL, stars.subsample.ratio.vec = NULL,
                               rep.num = rep.num, stars.thresh = 0.1, lam2init = NULL, lam1init = stars.l1$opt.lambda, verbose = T)
  res = try(JGL::JGL(inputlist, lambda1 = stars.l1$opt.lambda , lambda2 = stars.l2$opt.lambda, return.whole.theta = T), silent = T)
  tt1 = proc.time()
  return(list(stars.l1 = stars.l1,
              stars.l2 = stars.l2,
              res = res,
              runtime = tt1 - tt0))
}
