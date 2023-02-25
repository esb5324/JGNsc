library(igraph)

crt_net <- function(nodes.n,net.a,net.m){
  # generate scale free net
  nt <- sample_pa(n = nodes.n,power = net.a,m = net.m,directed = FALSE)
  output <- get.adjacency(nt,names = FALSE,edges = FALSE, sparse = FALSE)
  return(output) # power is the power of preferential attachment, m is the # of edges to add in each time step
}

crt_mat_A1 <- function(A0,U=0){
  # make A pd
  egv1 <- min(eigen(A0)$value)
  egv2 <- NULL
  if (is.matrix(U)) {egv2 <- min(eigen(A0+U)$value)} # if ID module, add U to A
  delta <- abs(min(egv1,egv2))+0.05 # + 0.05 makes sure we get pd
  # min of egv1 and egv2: based on papers
  A1 <- A0+U+delta*diag(dim(A0)[1])
  return(A1)
}

# -------------------------------------
# simulation functions
# -------------------------------------

#' generate the small matrix in block i
#' @param ni the number of rows and columns in block i
#' @param ud the range of uniform distribution that to generate the numbers in the block
#' @export
generateBlki <- function(ni, ud= c(-100:-60, 60:100)/100, runif_threshold, t_net, pd){
  mati <- matrix(0, ni, ni)
  if (t_net=="random"){
  for (i in 1:ni){
    for (j in i:ni){
      mati[i,j] <- ifelse(i!=j & runif(1) > runif_threshold, sample(ud, 1),  #half chance: coexpression of i and j. prob in [-1,-0.6] U [0.6,1]
                          ifelse(i==j,1,0))
    }
  }
  } else if (t_net=="power_law"){
    pl_mat <- crt_new(nodes.n=ni, net.a=1, net.m=1)
  for (i in 1:ni){
    for (j in i:ni){
      mati[i,j] <- ifelse(i!=j & pl_mat[i,j]==1, sample(ud, 1),  #half chance: coexpression of i and j. prob in [-1,-0.6] U [0.6,1]
                          ifelse(i==j,1,0))
    }
  }
    }
  mati0 = mati + t(mati) - diag(1, nrow = ni, ncol = ni)
  mati1 <- mati0 + diag(1, nrow = ni, ncol = ni)
if (pd=="diagplus1"){
  for (i in 1:20){
    # cat(i,".. \n")
    if (matrixcalc::is.positive.definite(mati1)){
      mati1 <- mati1
      break
    } else {
      mati1 <- mati1 + diag(1, nrow = ni, ncol = ni)
   }
    }
    } else if (pd=="eigen"){
      mati1 <- crt_mat_A1(mati1)
      }

  binv = solve(mati1)
  bii <- diag(binv)
  sigmam <- binv
  for (i in 1:ni){
    for (j in 1:ni){
      sigmam[i,j] <- binv[i,j]/sqrt(bii[i]*bii[j])
    }
  }
  res <- list(sigmam = sigmam,
              Bm = mati1)
  return(res)
}

#' check if the list of block dimensions are the same for the different conditions
#' @param xlist a list of vectors
#' @export
check_ident_list <- function(xlist){
  if (length(unique(unlist(lapply(xlist, length)))) == 1){ #length of sublists are the same
    nsame = 0 # identical pairs of sublist
    for (cc in 1:length(xlist)){
      for (dd in 1:length(xlist)){
        if (cc < dd){
          nsame = nsame + all(xlist[[cc]] == xlist[[dd]])
        }
      }
    }
    if (nsame == choose(length(xlist), 2)){
      y = T
    } else {
      message("The elements in sublists are not all identical... \n")
      y = F
    }
  } else {
    message("The length of list elements are different... \n")
    y = F
  }
  return(y)
}

#' generate a list of covariance matrices with designed structures
#' @param nivec.list a list of vectors. Each vector specifies the dimensions of the blocks under a condition.
#' For example, c(20,20) means a covariance matrix is composed of two block matrices, each having 20 genes in the block.
#' @param ud The uniform distribution range. Default is c(-100:-60, 60:100)/100
#' @param structure structure S, Weight W: "Identical S, Identical W", "Identical S, Diff W", "Diff S, Identical W", "Diff S, Diff W".
#' If using the structure "Diff S, Identical W", setting some structure different while keep the rest of structure and weights the same, specify diffblk = list(diffblk1=1, diffblk2=c(1,2),...)
#' @param diffblk a list of vectors indicating the location of blocks that have different structure/weights.
#' For example, nivec.list = list(nivec= c(10,20,30,40, rep(20,5)),nivec2 = rep(20,10)); diffblk = list(1:4,1:5)
#' @return a list of covariance matrices
#' @export
generateSigmaList <- function(nivec.list, ud = c(-100:-60, 60:100)/100,
                              structure = "Identical S, Identical W", diffblk = NULL, blk_runif_threshold,true_net="random",pos_def="diagplus1"){

  sigma.list <- list()
  # check if identical structure
  checkI <- check_ident_list(nivec.list)
  if (structure =="Identical S, Identical W"){
    if (checkI){ #------------------------------------------- if structures are identical
      blklist <- list()
      sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[1]]))
      nblk <- length(nivec.list[[1]])
      for (b in 1:nblk){
        ni <- nivec.list[[1]][b]
        blklist[[b]] <- generateBlki(ni=ni, ud=ud,runif_threshold=blk_runif_threshold,t_net=true_net,pd=pos_def)
        zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[1]][0:(b-1)]))
        zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[1]][(b+1):nblk]),0))
        temp <- blklist[[b]]$sigmam
        sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
      }
      gnames = paste("gene",1:sum(nivec.list[[1]]), sep = "")
      rownames(sigma) <- gnames
      colnames(sigma) <- gnames
      for(ss in 1:length(nivec.list)){
        sigma.list[[ss]] <- sigma
      }
    } else {
      message("nivec.list and the selected network structure do NOT match ...\n")
      break()
    }

  } else if (structure =="Identical S, Diff W"){
    if (checkI){ #------------------------------------- if structures are identical
      for(ss in 1:length(nivec.list)){
        blklist <- list()
        sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[ss]]))
        # sigma <- matrix(0, nrow = ni*nblk, ncol = ni*nblk)
        nblk <- length(nivec.list[[ss]])
        for (b in 1:nblk){
          ni <- nivec.list[[ss]][b]
          blklist[[b]] <- generateBlki(ni=ni, ud=ud,runif_threshold=blk_runif_threshold)
          zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[ss]][0:(b-1)]))
          zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[ss]][(b+1):nblk]),0))
          temp <- blklist[[b]]$sigmam
          sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
        }
        gnames = paste("gene",1:sum(nivec.list[[ss]]), sep = "")
        rownames(sigma) <- gnames
        colnames(sigma) <- gnames
        sigma.list[[ss]] <- sigma
      }
    } else {
      message("nivec.list and the selected network structure do NOT match ...\n")
      break()
    }

  } else if (structure =="Diff S, Identical W"){
    # for the part that structures are the same, weights are the same
    # assume the first ndiff rows have different structure.
    blklist <- list()
    sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[1]]))
    nblk <- length(nivec.list[[1]])
    for (b in 1:nblk){
      ni <- nivec.list[[1]][b]
      blklist[[b]] <- generateBlki(ni=ni, ud=ud,runif_threshold=blk_runif_threshold)
      zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[1]][0:(b-1)]))
      zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[1]][(b+1):nblk]),0))
      temp <- blklist[[b]]$sigmam
      sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
    }
    gnames = paste("gene",1:sum(nivec.list[[1]]), sep = "")
    rownames(sigma) <- gnames
    colnames(sigma) <- gnames
    sigma.list[[1]] <- sigma

    for(ss in 2:length(nivec.list)){ # from the second matrix, only change the diffblk[[ss]] block part
      blklist <- list()
      sigma2 <- sigma
      nblk <- length(nivec.list[[ss]])
      temps <- matrix(0, nrow = 0, ncol = sum(nivec.list[[ss]]))
      for (b in diffblk[[ss]]){
        ni <- nivec.list[[ss]][b]
        blklist[[b]] <- generateBlki(ni=ni, ud=ud,runif_threshold=blk_runif_threshold)
        zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[ss]][0:(b-1)]))
        zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[ss]][(b+1):nblk]),0))
        temp <- blklist[[b]]$sigmam
        temps <- cbind(zeroleft, temp, zeroright)
        diffid <- (sum(nivec.list[[ss]][0:(b-1)])+1) : sum(nivec.list[[ss]][0:b])
        sigma2[diffid,] <- temps
      }
      sigma.list[[ss]] <- sigma2
    }


  } else if (structure == "Diff S, Diff W"){
    for(ss in 1:length(nivec.list)){
      blklist <- list()
      sigma <- matrix(0, nrow = 0, ncol = sum(nivec.list[[ss]]))
      # sigma <- matrix(0, nrow = ni*nblk, ncol = ni*nblk)
      nblk <- length(nivec.list[[ss]])
      for (b in 1:nblk){
        ni <- nivec.list[[ss]][b]
        blklist[[b]] <- generateBlki(ni=ni, ud=ud,runif_threshold=blk_runif_threshold)
        zeroleft <- matrix(0, nrow = ni, ncol = sum(nivec.list[[ss]][0:(b-1)]))
        zeroright <- matrix(0, nrow = ni, ncol = ifelse(b<nblk,sum(nivec.list[[ss]][(b+1):nblk]),0))
        temp <- blklist[[b]]$sigmam
        sigma <- rbind(sigma, cbind(zeroleft, temp, zeroright))
      }
      gnames = paste("gene",1:sum(nivec.list[[ss]]), sep = "")
      rownames(sigma) <- gnames
      colnames(sigma) <- gnames
      sigma.list[[ss]] <- sigma
    }
  } else {
    message("Please choose the joint network structure from: 'Identical S, Identical W', 'Identical S, Diff W', 'Diff S, Identical W', 'Diff S, Diff W'")
  }
  return(sigma.list)
}

#' Map scRNA-seq counts to a known covariance structure
#'
#' @param sigma a covairance matrix
#' @param ngene number of genes
#' @param n number of samples
#' @param a10 the hyperparameter for generating alpha, where alpha~Gamma(a10,b10)
#' @param b10 the hyperparameter for generating alpha, where alpha~Gamma(a10,b10)
#' @param a20 the hyperparameter for generating beta, where beta~Gamma(a20,b20)
#' @param b20 the hyperparameter for generating beta, where beta~Gamma(a20,b20)
#' @param a3 the hyperparameter for generating pij, where pij~beta(a3, b3)
#' @param b3 the hyperparameter for generating pij, where pij~beta(a3, b3)
#' @return a list of i) raw count data matrix; ii) adjacency matrix; iii) the true expression mean level theta matrix; iv) zij: the dropout matrix; v) sigma: the covariance matrix; vi) precision: the precision matrix.
#' @export
CountMap <- function(sigma, ngene, n, a1 = 3,
                      b1 = 1, a20 = 2,  b20 = 3, a30 = 1, b30 = 10){
  # CountMap5()
  precision1 <- solve(sigma)
  mu <- rep(0, ngene)
  adj <- abs(sign(precision1))
  diag(adj) <- 0
  x <- mvtnorm::rmvnorm(n, mu, sigma)
  y <- x

  # -------------------------------------------------
  # generate count data from posterior distribution
  # one subject at a time
  # subject i
  z <- matrix(1, nrow = n, ncol = ngene)
  pij.vec <- rep(1, ngene)
  ycomplete <- x

  alphavec <- rep(0, ngene)
  betavec <- rep(0,ngene)
  scmeanvec <- rep(0,ngene)
  # -----------------------------
  for (j in 1:ngene) {
    dat <- x[, j]
    mu_v <- mean(dat)
    sd_v <- sd(dat)
    p_v <- pnorm(dat, mu_v, sd_v)

    # ---------------------
    # simulate sc from zip -- 1 gene
    # 3. alpha, beta, thetaij
    alphaj <- rgamma(1, shape = a20, rate = b20)
    betaj <- rgamma(1, shape = a30, rate = b30)
    scmeanj <- rgamma(1, shape = alphaj, rate = betaj)
    if (any(scmeanj >=1000) | any(scmeanj <1) ){
      scmeanj <- runif(1, 1,10)
    }
    sctemp <- rpois(n, scmeanj)
    alphavec[j] <- alphaj
    betavec[j] <- betaj
    scmeanvec[j] <- scmeanj

    ytemp <- quantile(sctemp, p_v)
    ycomplete[,j] <- ytemp
    pij <- rbeta(1, shape1 = a1, shape2 = b1)
    zijc <- sapply(1:n, function(x){
      px = ifelse(scmeanvec[j] > 10, 1, 1- ((1-pij) + pij * exp(-scmeanvec[j])))
      rbinom(1, size = 1, prob = px)
    })

    z[,j] <- zijc
    pij.vec[j] <- round(pij,3)
  }

  # -------------------------------
  yout = ycomplete*z

  colnames(yout) <- paste("gene",1:ngene, sep = "")
  rownames(yout) <-  paste("cell",1:n,sep = "")
  colnames(ycomplete) <- paste("gene",1:ngene, sep = "")
  rownames(ycomplete) <- paste("cell",1:n,sep = "")
  colnames(z) <- paste("gene",1:ngene, sep = "")
  rownames(z) <- paste("cell",1:n,sep = "")

  result <- list()
  result$count <- yout
  result$count.nodrop <- ycomplete
  result$zijc <- z
  result$thetaj <- round(scmeanvec,2)
  result$pij <- pij.vec
  result$sigma <- sigma
  result$precision <- precision1
  return(result)
}


#' generate a list of raw count matrices based on the list of covariance matrices
#'
#' @param sigma.list a list of covariance matrices. Please make sure the positive definite property.
#' @param nvec a vector of sample size under each condition.
#' @param ngene integer. the number of genes
#' @param a3 the parameter for the non-dropout rate pij~ beta(a3, b3)
#' @param b3 the parameter for the non-dropout rate pij~ beta(a3, b3)
#' @return a list of count matrices
#' @export
getCountList <- function(sigma.list, nvec = c(500, 500), ngene = NULL, a3 = 2, b3 = 1){
  if (is.null(ngene)){
    ngene = ncol(sigma.list[[1]])
  }
  count.list <- list()
  for (c in 1:length(sigma.list)){
    counti <- CountMap(sigma = sigma.list[[c]], ngene = ngene, n=nvec[c], a3 = a3, b3 = b3)
    count.list[[c]] <- counti
  }
  return(count.list)
}

#' calculate area under ROC curve
#' @param score a score (here we use partial correlation values) for each edge
#' @param bool logical, indicating the existance of an edge in a network.
#' @references https://blog.mbq.me/augh-roc/
#' @return the value of AUROC
#' @export
auroc <- function(score, bool) {
  n1 <- sum(!bool)
  n2 <- sum(bool)
  U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

#' calculate area under precision-recall curve
#' @param z a list of scores (here we use partial correlation values)
#' @param trueadj.list a list of true adjacency matrices
#' @return AUPRC
#' @import PRROC
#' @export
PRROC_pr <- function(z, trueadj.list){
  s1 = c(abs(unlist(z)))[c(unlist(trueadj.list))==1]
  s0 = c(abs(unlist(z)))[c(unlist(trueadj.list))==0]
  pr <- pr.curve(scores.class0 = s1,
                 scores.class1 = s0)
  pr$auc.integral
}
