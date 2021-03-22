#' JGNsc framework main function
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
#' @return theta.star.npn: the imputed and Gaussian transformed list of matrices
#' @return if runNetwork = T, then it will return the JGL results (the precision matrices can be accessed by result$JGL$theta)
#'         , aic.table (the AIC values with corresponding tuning parameter candidates), partcorr (the list of partial correlation matrices)
#' @export

JGNsc <- function(observed.list, mask.rate = 0.15, nrep = 20, warm = 100, iter = 200,
                  min.cell = 3, runNetwork = F, l1.vec = NULL,
                  l2.vec = NULL){
  zip.list <- lapply(observed.list, JGNsc.cont, warm = warm, iter = iter, min.cell= min.cell)
  g.common <- Reduce(intersect, lapply(zip.list, function(x){rownames(x$y.impute)}))
  theta.star <- lapply(zip.list, function(x){
    y = x$y.impute[g.common,]
    return(y)
  })
  observed.list2 <- lapply(observed.list, function(x){
    y = x[g.common,]
  })

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
      x.mask.mcimpute = mcImpute.R(x.mask, preprocess = T) #samples by genes input
      current.imp <- t(x.mask.mcimpute)
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
  theta.star.npn <- lapply(theta.star.t, huge::huge.npn)

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
        tps <- AIC_select(mat.list.t.gau = theta.star.npn, lam1 = lam1, lam2 = lam2)
        jgnsc.aic <- rbind(jgnsc.aic, tps)
      }
    }
    JGL.res <- JGL::JGL(theta.star.npn, lambda1 = jgnsc.aic[which.min(jgnsc.aic[,3]),1],
                   lambda2 = jgnsc.aic[which.min(jgnsc.aic[,3]),2], return.whole.theta = T)
    partcorr <- lapply(JGL.res$theta, prec2partialcorr)
  } else {
    JGL.res <- NULL
    jgnsc.aic <- NULL
    partcorr <- NULL
  }
  return(list(theta.star.npn = theta.star.npn,
              JGL = JGL.res,
              aic.table = jgnsc.aic,
              partcorr = partcorr))
}


#' JGNsc framework -- MCMC continuize function
#'
#' @param y matrix. dim: genes by samples
#' @param stepsize the stepsize for the random proposal, default is 0.5
#' @param warm the number of warm steps in the MCMC procedure
#' @param iter the number of iterations in the MCMC procedure
#' @param min.cell the min number of cells that should a gene should express
#' @param a1 the initial value of the parameter from the alpha parameter distribution: alpha ~ Gamma(a1,b1)
#' @param b1 the initial value of the parameter from the alpha parameter distribution: alpha ~ Gamma(a1,b1)
#' @param a2 the initial value of the parameter from the beta parameter distribution: beta ~ Gamma(a2,b2)
#' @param b2 the initial value of the parameter from the beta parameter distribution: beta ~ Gamma(a2,b2)
#' @param a3 the initial value of the parameter from the (non-dropout) p parameter distribution: p ~ Beta(a3,b3)
#' @param b3 the initial value of the parameter from the (non-dropout) p parameter distribution: p ~ Beta(a3,b3)
#' @param c0 a constant number close to zero, for the dropout case, theta value is imputed as c0 first.
#' @return y.cont: imputed and continuized count matrix
#' @return keep.gene: genes kept
#' @export
JGNsc.cont <- function(y, stepsize = 0.5, warm = 50, a2 = 0.001, a3 = 0.001, b2 = 1e8, b3 = 1e8,
                        a1 = 2, b1 = 1, iter = 50, c0 = 0, min.cell = 2, dropThreshold = 0.65){
  # JGNsc.cont2()
  # change this to Rcpp
  # c0: a constant number close to zero, for the dropout case, theta value
  # y is genes by samples
  # remove genes with zero expression (low expression level: at least 20 cells expressed this gene) across samples
  # a1 and b1 are the parameter for the non-dropout event!!!!

  keep.gene <- rowSums(y>0)> min.cell
  y <- y[keep.gene,]
  tau <- colSums(y)/median(colSums(y))
  ng <- nrow(y)
  nsample <- ncol(y)
  z <- matrix(1, ncol = ncol(y), nrow = nrow(y))
  yz = y*z
  # by subject: subj by genes
  yzsum = apply(yz, 1, sum)
  zsum = apply(z,1, sum)

  # set up initial values
  logalpha <- rep(0, ng)
  alpha <- exp(logalpha)
  beta <- seq(1, ng)
  thetaj <- (alpha + yzsum)/(beta + zsum)
  z.total <- matrix(0, ncol = nsample, nrow = ng)

  ave <- matrix(0, ncol = ng, nrow = nsample)
  theta.total <- rep(0, ng)
  npoi.total <- matrix(0, ncol = ng, nrow = nsample )
  pi <- rowMeans(y>0)
  pi.total <- rep(0,ng)

  # write this in parallel mode ???
  # estimation of posterior parameters
  # for a gene, if it expressed in a subject in less than 3 cells, skip impute
  for (ii in 1:(iter+warm)){ #number of iterations, warm: warm up (burn), dump the first warm number of initial results
    cat(ii,"..\n")
    for (gg in 1:ng){ # for each gene
      # ii=1
      # gg=1
      ##### update z #####
      P0 <- ifelse(y[gg,]==0, 1,0)*(1-pi[gg])
      P1 <- exp(-thetaj[gg]*tau)*pi[gg]
      z[gg,] <- sapply(P1/(P1+P0), function(x){
        if (is.na(x)){
          x = 1
        }
        rbinom(1,1, x)
      })
      z[gg,y[gg,]>0] <- 1 #make sure: y>0 -> z=1

      ##### update alpha #####comparing log(f), use the complete function, not propto...
      newlogalpha <- logalpha[gg] + gasdev(1)*stepsize
      alpha[gg] <- exp(logalpha[[gg]])
      newalpha <- exp(newlogalpha)

      newsum <- 0
      newsum <- (a2-1)*newlogalpha - lgamma(newalpha) +
        newalpha*(- b2+log(beta[gg]) + sum(log(thetaj[gg])))

      newsum <- newsum / b2
      sum <- 0
      sum <- (a2-1)*logalpha[gg] - lgamma(alpha[gg]) +
        alpha[gg]*(- b2+ log(beta[gg]) + sum(log(thetaj[gg])))
      sum <- sum / b2

      r <- newsum - sum
      if (is.na(r)){
        accept = 0} else if(r>0){
          accept = 1
        } else {
          un = runif(1)
          if (un < exp(r) & newalpha < alpha[gg]){
            accept = 1
          } else{
            accept = 0
          }}

      if (accept ==1){
        logalpha[gg] = newlogalpha
        alpha[gg] = newalpha}
      # accept_loc[gg] = accept_loc[gg] +1}
      # total_loc[gg] = total_loc[gg] + 1

      ##### update beta #####
      sum <- 0
      sum <- sum(thetaj[gg])
      beta[gg] <- rgamma(1, shape = alpha[gg]+a3, rate = sum + b3)

      ##### update theta*z #####
      temp.yz = y[gg,]*z[gg,]
      # by subject: subj by genes
      temp.yzsum = sum(temp.yz)
      temp.zsum = sum(z[gg,])
      temp.zsumtau = sum(z[gg,]*tau)
      temp.y0sum = sum(y[gg,]>0)

      if (temp.y0sum > 1){ #if only <2 cells in a subject expressed this gene, skip
        thetaj[gg] <- rgamma(1, shape = alpha[gg]+temp.yzsum , rate = beta[gg]+temp.zsumtau )
      } else {
        thetaj[gg] <- 0
      }

      ##### update p #####
      temp.1zy0 = sum((y[gg,]==0)*(1-z[gg,]))
      pi[gg] <- rbeta(1, shape1 = a1 + temp.zsum , shape2 = b1+temp.1zy0 )

    } #end for for each gene gg
    b2 <- b2 + 1/iter
    b3 <- b2
    if (ii > warm){
      theta.total = theta.total + thetaj
      npoi <- thetaj >0
      npoi.total <- npoi.total + npoi
      pi.total <- pi.total + pi
      z.total <- z.total + z
    }
  }

  # estimate the final theta values
  ave <- theta.total/npoi.total
  pi.final <- pi.total/iter
  z.final <- z.total/iter
  ave[is.na(ave)] <- 0
  ave[is.infinite(ave)] <- 0
  avetau <- ave*tau
  # --------
  #imputed theta: y>0 parts are the same, y=0 parts are imputed using Gaussian distribution
  z <- z.final >= dropThreshold
  ave.impute <- y
  for (gg in 1:ng){
    ave.impute[gg,] = y[gg,] + avetau[gg]*(1-z[gg,])*(y[gg,]==0)
  }
  ave.impute[is.na(ave.impute)] <- 0
  ave.impute[is.infinite(ave.impute)] <- 0

  res <- list(thetaij = ave,
              thetatau = avetau,
              pi = pi.final,
              zijc = z ,
              y.impute = ave.impute,
              keep.gene = keep.gene )
  return(res)
}


#' soft threshold function
#' @param x numeric
#' @param thld the threshold value
#' @return a number by soft thresholding
#' @export
softThreshold <- function(x, thld = 0.5){
  z= sign(x)* max(0, abs(x) - thld)
  return(z)
}

#' Matrix completion imputation method
#' @description the original algorithm was implemented in Matlab. Here, we rewrite it in R.
#' @references https://www.frontiersin.org/articles/10.3389/fgene.2019.00009/full
#' @source https://github.com/aanchalMongia/McImpute_scRNAseq
#' @param data a matrix, samples by genes
#' @param preprocess logical, whether to preprocess data: Removing BAD genes, median normalization and log-transformation.
#' @param eps a small positive number that controls convergence, when calculating the Frobenius norm.
#' @param normfac normalize factor, default is 1.
#' @param insweep the max number of iterations when minimizing the rank of matrix given a tuning parameter(lambda) value
#' @param tol a small positive number that controls convergence
#' @param decfac decay factor, controls the rate when searching for the softthresholding parameter.
#' @param min_count for filtering step, minimum number of counts that a gene expressed in a cell
#' @param min_cells for filtering step, minimum number of cells that a gene expressed in.
#' @param verbose logical, show message.
#' @return an imputed matrix
#' @export
mcImpute.R <- function(data, preprocess = T,
                       eps = 1e-12, normfac = 1, insweep = 20,
                       tol = 1e-4, decfac = 0.7, min_count=2,
                       min_cells=3, verbose = F){
  M = data >0
  if (any(is.na(data))){
    message("NA values exist in the data matrix.")
    break
  }
  if (preprocess){
    # Removing BAD genes, median normalization and log-transformation
    gene.filter = colSums(data > min_count) > min_cells
    data.filterg = data[,gene.filter]
    libsize = rowSums(data.filterg)
    data.filterg.norm = apply(data.filterg, 1, function(x){
      x/sum(x)*median(libsize)
    })
    y = log2(data + 1)
  } else {
    y = data
  }
  alpha = 1.1*normfac
  xinit = matrix(0, nrow = nrow(data), ncol = ncol(data))
  x = xinit
  lambda.init = decfac*max(abs(y))
  lambda = lambda.init

  f_current = norm(y - M*x,"f") + lambda*norm(x)
  while (lambda > lambda.init*tol){
    if (verbose){
      cat("Lambda = ", lambda, "\n")
    }
    for (ins in 1:insweep){
      f_previous = f_current
      B = x + (1/alpha)*M*(y - M*x)
      USV = svd(B)
      s = sapply(USV$d, softThreshold, thld = lambda/(2*alpha))
      S = diag(s)
      X = USV$u%*%S%*%t(USV$v)
      X[X<0] = 0
      x =X
      f_current = norm(y - M*x,"f") + lambda*norm(x)
      if (abs(f_current-f_previous)/ abs(f_current+f_previous) < tol){
        break
      }
    }
    if (norm(y - M*x,"f") < eps){
      break
    }
    lambda = decfac*lambda
  }

  if (preprocess){
    resX = round(2^x - 1)
  } else {
    resX = round(x)
  }
  return(resX)
}
