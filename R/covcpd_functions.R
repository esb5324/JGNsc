#' covcpd change point detection by covariance testing
#' @param X raw count matrix, with genes by samples
#' @param nperm number of permutation
#' @param siglvl significance level, default is 0.05
#' @param k the minimum length of a segment / partition
#' @param maxseg maximum number of segments / partitions
#' @param verbose whether to print out the processing info while runing the funtion.
#' @param search_by the space between candidate change points.
#' @return cps: change points
#' cps_order: change points by the order when they are found
#' candtau: candidate change points that were tested
#' pvals: p-values for the candidate change points
#' pthresh: the learned p-value threshold.
#' @import evd
#' @export 
cpd_cov_dist <- function(X, nperm = 100, siglvl = 0.05, k=30, maxseg = NULL, 
                         verbose = T, search_by = NULL){
  n <- ncol(X)
  cp.1o = cp.1 <- c(1,n) # NULL SET FOR START
  pvec = NULL
  candtau = NULL
  
  # find a cp from the whole interval first
  seg <- 1
  m = length(cp.1)
  if(n/k < 3){
    k = round(n/3)
  }
  if(is.null(maxseg)){
    maxseg = max(3, ceiling(n/k))
  }
  if(is.null(search_by)){
    if(n>1000){
      search_by = 20
    } else if(n>500){
      search_by = 10
    } else if(n <=500 & n>200){
      search_by = 5
    } else {
      search_by = 2
    }
  }
  
  #~~~~~~~~~~~~ 1. learn mu and beta
  if (verbose){
    message("Start learning Gumbel distr. parameters mu and beta... \n")
  }
  if ((n-k)/k <= 5){
    by1 = 0.2
  } else if ((n-k)/k <= 9){
    by1 = 0.4
  } else if((n-k)/k <= 15){
    by1 = 0.6
  } else {
    by1 = 1
  }
  nr =seq(1,(n-k)/k, by1)
  searchx = unique(c(round(n/(1+nr)), n - round(n/(1+nr))))
  searchx = searchx[order(searchx)]
  
  stat.allperm1 = CLXstatPerm(t(X), searchx = searchx, start = 1, nperm = nperm)
  stat.allperm1 =as.data.frame(stat.allperm1)
  stat.allperm1$nratio = pmax(stat.allperm1$V2, n - stat.allperm1$V2)/pmin(stat.allperm1$V2, n - stat.allperm1$V2)
  
  dtnullmed1 = aggregate(stat.allperm1$V1, list(stat.allperm1$nratio), median)
  dtnullvar1 = aggregate(stat.allperm1$V1, list(stat.allperm1$nratio), var)
  dtnullbeta1 = sqrt(6*dtnullvar1$x/pi^2)
  dtnullmu1  = dtnullmed1$x + dtnullbeta1*log(log(2))
  
  dtx1 = data.frame(mu = dtnullmu1, beta=dtnullbeta1, nratio = dtnullmed1$Group.1)
  dtx1$nratio2 = dtx1$nratio^2
  fit1s = lm(mu ~ nratio + nratio2, data = dtx1)
  fit2s = lm(beta ~ nratio + nratio2, data = dtx1)
  
  
  coa = coef(fit1s)
  cob = coef(fit2s)
  func_permmu = function(xx){
    out = coa[1] + xx*coa[2] + xx^2*coa[3]
    return(out)
  }
  func_permbeta = function(xx){
    out = cob[1] + xx*cob[2] + xx^2*cob[3]
    return(out)
  }
  
  # ~~~~~~~~~~~~~~~ 2. learn p-threshold
  if (verbose){
    message("Start learning type-I error threshold by permutation... \n")
  }
  # searchx = seq(k,n-k, by = max(1, round((n-2*k)/10)))
  searchx = seq(k,(n-1 - k), by = search_by) 
  learnp.perm = CLXPermHier(t(X), searchx = searchx, start = 1,nt2 = 4, minband = k, nperm = nperm)
  learnp.perm = as.data.frame(learnp.perm)
  learnp.perm$nratio = pmax(learnp.perm$V2, learnp.perm$V4-learnp.perm$V2)/pmin(learnp.perm$V2, learnp.perm$V4 - learnp.perm$V2)
  learnp.perm$pgum = 1- evd::pgumbel(learnp.perm$V1, loc = func_permmu(learnp.perm$nratio), scale = func_permbeta(learnp.perm$nratio))
  learnp_acrosst2 = aggregate(learnp.perm$pgum, list(t1 = learnp.perm$V2, permindex = learnp.perm$V3), min)
  learnp_acrossperm = aggregate(learnp_acrosst2$x, list(learnp_acrosst2$permindex), min)
  pthresh = quantile(learnp_acrossperm$x, siglvl)  
  
  # start search
  # ~~~~
  if (verbose){
    message("Start searching for change points... \n")
  }
  while (seg < maxseg){ 
    newcps = NULL
    for (jj in 1:(m-1)){
      if (cp.1[jj] > 1){
        start.p = cp.1[jj] +1
      } else {
        start.p = 1
      }
      end.p = cp.1[jj+1] 
      if ((end.p - start.p) <= 2*k){
        message("skip interval (", start.p, ", ", end.p, ")\n")
      } else { # keep finding cp
        
        stat.all <- NULL 
        searchx = seq((start.p + k),(end.p-1 - k), by = search_by) 
        
        for (t1 in searchx){ 
          t2seq = c(seq(t1+1+k, end.p, by = max(round((end.p - t1-1-k)/4), 1)), end.p)
          pt2   = 1
          pgums = NULL
          mns = NULL
          for(tt in 1:length(t2seq)){
            X1 <- as.matrix(t(X[,start.p:t1]))
            X2 <- as.matrix(t(X[,(t1+1):t2seq[tt]]))  
            stat0 = testCov_cpp(X=X1,Y=X2)
            nr = max(nrow(X1),nrow(X2))/min(nrow(X1),nrow(X2))
            pgum = 1- evd::pgumbel(stat0$CLX, loc = func_permmu(nr), scale = func_permbeta(nr))
            pgums = rbind(pgums, c(pgum, stat0$pvalue, t2seq[tt]))
            mns = rbind(mns, c( stat0$CLX,func_permmu(nr)))
            pt2 = min(pgum, pt2)
          }     
          # permutation step for each candidate point
          stat.all <- rbind(stat.all, c(max(abs(mns[,1] - mns[,2])),pt2, t1))
        }
        stat.all <- as.data.frame(stat.all) 
        
        if(length(searchx)>5){
          minp = min(stat.all$V2)
          taus = stat.all$V3[stat.all$V2 == minp] 
          if (length(taus)>1){ 
            stmp = stat.all[stat.all$V3 %in% taus,]
            tau0 <- stmp$V3[which.max(stmp$V1)]
          } else {
            tau0 = taus
          }
          
          if (search_by > 2){
            candidates = seq(tau0 - search_by, tau0 + search_by, by=1)
            stat.all.cand = NULL
            for (t1 in candidates){ 
              t2seq = c(seq(min(t1+1+k, end.p), end.p, by = max(round((end.p - t1-1-k)/5), 2)), end.p)
              pt2   = 1
              pgums = NULL
              mns = NULL
              for(tt in 1:length(t2seq)){
                X1 <- as.matrix(t(X[,start.p:t1]))
                X2 <- as.matrix(t(X[,(t1+1):t2seq[tt]]))  
                stat0 = testCov_cpp(X=X1,Y=X2)
                nr = max(nrow(X1),nrow(X2))/min(nrow(X1),nrow(X2))
                pgum = 1- evd::pgumbel(stat0$CLX, loc = func_permmu(nr), scale = func_permbeta(nr))
                pgums = rbind(pgums, c(pgum, stat0$pvalue, t2seq[tt]))
                mns = rbind(mns, c( stat0$CLX,func_permmu(nr)))
                pt2 = min(pgum, pt2)
              }     
              stat.all.cand <- rbind(stat.all.cand, c(max(abs(mns[,1] - mns[,2])),pt2, t1))
            }
            
            stat.all.cand = as.data.frame(stat.all.cand)
            minp = min(stat.all.cand$V2)
            taus = stat.all.cand$V3[stat.all.cand$V2 == minp] 
            
            if (length(taus)>1){ 
              stmp = stat.all.cand[stat.all.cand$V3 %in% taus,]
              tau0 <- stmp$V3[which.max(stmp$V1)]
            } else {
              tau0 = taus
            }
          }  
          
        } else { 
          minp = min(stat.all$V2)
          tau0 = stat.all$V3[which.min(stat.all$V2)]
        } 
        
        pvec = c(pvec, minp)
        candtau = c(candtau, tau0)
        if (minp < pthresh){
          newcps = c(newcps, tau0)  
        }  
        
      }
      
    } 
    
    cp.1o <- append(cp.1o, newcps)
    cp.1 <- unique(cp.1o)
    cp.1 <- cp.1[order(cp.1)]
    
    if (length(cp.1) > m){
      m = length(cp.1)
      seg = m - 1
    } else {
      break
    }
  }
  return(list(cps = cp.1,
              cps_order = cp.1o,
              pvals = pvec,
              candtau = candtau,
              pthresh = pthresh))
  
}
