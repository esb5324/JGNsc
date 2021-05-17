// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec arma_rowSums(const arma::mat X){
  int       nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  };
  return(out);
};

// [[Rcpp::export]]
arma::vec arma_colSums(const arma::mat X){
  int       nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i));
  };
  return(out);
};

// [[Rcpp::export]]
arma::vec arma_rowMeans(const arma::mat X){
  int       nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = mean(X.row(i));
  };
  return(out);
};

// [[Rcpp::export]]
arma::vec arma_colMeans(const arma::mat X){
  int       nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = mean(X.col(i));
  };
  return(out);
};

// [[Rcpp::export]]
arma::mat BinMatrix(arma::mat x){
  arma::mat out(x.n_rows, x.n_cols, fill::zeros);
  for (unsigned int i=0; i<x.n_rows; i++){
    for (unsigned int j=0; j<x.n_cols; j++){
      if(x(i,j) > 0){
        out(i,j) =1;
      }
    }
  };
  return out;
};

// [[Rcpp::export]]
double MaxAbsMatrix(arma::mat x){
  double maxe = 0;
  for (unsigned int i =0; i<x.n_rows; i++){
    for (unsigned int j =0; j<x.n_cols; j++){
      if (abs(x(i,j)) > maxe){
        maxe = abs(x(i,j));
      }
    }
  };
  return maxe;
};

// [[Rcpp::export]]
double frobeniousNorm(arma::mat x){
  double out1 =0.0, out = 0.0;
  for (unsigned int i = 0; i < x.n_rows; i++){
    for (unsigned int j = 0; j < x.n_cols; j++){
      out += pow(x(i,j), 2);
    }
  }
  out1 = sqrt(out);
  return out1;
};

// [[Rcpp::export]]
double oneNorm(arma::mat x){
  // maximum abs col sum
  return max(arma_colSums(x));
};


// [[Rcpp::export]]
double softThresh_cpp(double x, double threshold = 0.5){
  double z = sign(x) * std::max(0.0, abs(x) - threshold);
  return z;
};


// [[Rcpp::export]]
double gasdev(double r=1){
  double v1 =0.0, v2 =0.0;
  while(r>=1){
    v1 = R::runif(0,1)*2-1;
    v2 = R::runif(0,1)*2-1;
    r = pow(v1,2) + pow(v2,2);
  };
  double fac = sqrt(-2*log(r)/r);
  return v2*fac;
};


// [[Rcpp::export]]
List JGNsc_cont_cpp(arma::mat y, int minCell, int iter = 5000, int warm = 2000, double stepsize = 0.5,
                    double dropThreshold = 0.7,
                    double a1=2, double b1=1,
                    double a2 = 0.001, double b2 = 1e8,
                    double a3 = 0.001, double b3 = 1e8){
  arma::vec keepgene(y.n_rows);
  for (unsigned int i =0; i < y.n_rows; ++i){ //genes
    int ybin = 0;
    for (unsigned int j=0; j < y.n_cols; ++j){ // cells
      if (y(i,j) >0){
        ybin +=1;
      };
      keepgene(i) = ybin;
    };
  };
  arma::mat ysub = y.rows(find(keepgene > minCell));

  int ng      = ysub.n_rows;
  int nsample = ysub.n_cols;
  // calculate library size
  arma::vec tau = arma_colSums(ysub) / median(arma_colSums(ysub));

  // indicator matrix z
  arma::mat Z(ng, nsample, fill::ones);
  arma::mat yz = ysub % Z;
  arma::vec yzsum = arma_rowSums(yz);
  arma::vec zsum = arma_rowSums(Z);

  // SET UP INITIAL VALUES
  arma::vec logalpha = zeros(ng);
  arma::vec alpha = exp(logalpha);
  arma::vec beta = ones(ng);
  arma::vec thetaj = (alpha + yzsum)/(beta + zsum);
  arma::mat zTotal(ng, nsample, fill::zeros);
  arma::vec thetaTotal = zeros(ng);
  arma::vec npoiTotal = zeros(ng);
  arma::vec piTotal = zeros(ng);
  arma::mat ysub_bin(ng, nsample, fill::zeros);
  arma::vec ave(ng);
  arma::mat ave_impute = ysub;

  for (unsigned int i =0; i < ysub.n_rows; ++i){ //genes
    for (unsigned int j=0; j < ysub.n_cols; ++j){ // cells
      if (ysub(i,j) >0){
        ysub_bin(i,j) = 1;
      };
    };
  };
  arma::vec pi_vec = arma_rowMeans(ysub_bin);

  // parallel mode?
  //-------- start mcmc
  arma::vec P0, P1, prob_, thetaj_bin = zeros(ng);
  double    newsum_ = 0.0, sum_ =0.0, newlogalpha, newalpha, r, tempyzsum, tempzsum, tempzsumtau, tempy0sum, temp1zy0;
  int       accept_;
  arma::vec _ones = ones(nsample);

  for (int ii =0; ii < iter+warm; ii++){ // number of iterations, warm: warm up (burn), dump the first warm number of initial results
    for (int gg=0; gg < ng; gg++){
      // ~~~~~ 1. update Z
      P0 = (_ones - ysub_bin.row(gg).t())*(1 - pi_vec(gg));
      P1 = exp(-thetaj(gg)*tau)*pi_vec(gg);
      prob_ = P1/(P1+P0);
      for (int ss=0; ss<nsample; ss++){
        if(R_IsNA(prob_(ss))){
          prob_(ss) = 1;
        };
        Z(gg, ss) = R::rbinom(1,prob_(ss)); //Z(gg,ss)
      };

      // ~~~~~ 2. update alpha
      newlogalpha = logalpha(gg) + stepsize*gasdev(1);
      alpha(gg)   = exp(logalpha(gg));
      newalpha    = exp(newlogalpha);
      newsum_     = (a2-1)*newlogalpha - R::lgamma1p(newalpha) + newalpha*(-b2 + log(beta(gg)) + log(thetaj(gg)));
      newsum_     /= b2;
      sum_        = (a2-1)*logalpha(gg) - R::lgamma1p(alpha(gg)) + alpha(gg)*(-b2 + log(beta(gg)) + log(thetaj(gg)));
      sum_        /= b2;
      r           = newsum_ - sum_;
      if(R_IsNA(r)){
        accept_ = 0;
      } else if (r>0){
        accept_ = 1;
      } else{
        double un_ = R::runif(0,1);
        if(un_ < exp(r) && newalpha < alpha(gg)){
          accept_ = 1;
        } else {
          accept_ = 0;
        }
      };
      if (accept_ == 1){
        logalpha(gg) = newlogalpha;
        alpha(gg) = newalpha;
      };

      // ~~~~~ 3. update beta
      beta(gg) = R::rgamma(alpha(gg)+a3, 1/(thetaj(gg)+b3)); // rgamma is slightly different from R version

      // ~~~~~ 4. update theta*z
      tempyzsum   = sum(ysub.row(gg).t() % Z.row(gg).t());
      tempzsum    = sum(Z.row(gg));
      tempzsumtau = sum(Z.row(gg).t() % tau);
      tempy0sum   = sum(ysub_bin.row(gg));
      if(tempy0sum >1){
        thetaj(gg) = R::rgamma(alpha(gg)+tempyzsum, 1/(beta(gg)+tempzsumtau));
      } else {
        thetaj(gg) = 0;
      };
      if (thetaj(gg)>0){
        thetaj_bin(gg) = 1;
      } else {
        thetaj_bin(gg) = 0;
      };

      // ~~~~~ 5. update p
      temp1zy0   = sum(ysub_bin.row(gg).t() % (ones(nsample) - Z.row(gg).t()));
      pi_vec(gg) = R::rbeta(a1+tempzsum, b1+temp1zy0);
    };
    b2 += 1/iter;
    b3 = b2;
    if(ii > warm){
      thetaTotal += thetaj;
      npoiTotal += thetaj_bin;
      piTotal += pi_vec;
      zTotal += Z;
    }
  };

  // average and estimate final values
  ave = thetaTotal / npoiTotal;
  piTotal /= iter;
  zTotal /= iter;
  for (int gg=0; gg<ng; gg++){
    if (R_IsNA(ave(gg)) || !(arma::is_finite(ave(gg)))  ){
      ave(gg) = 0;
    }
  };

  for (unsigned int i =0; i < Z.n_rows; ++i){ //genes
    for (unsigned int j=0; j < Z.n_cols; ++j){ // cells
      if (zTotal(i,j) > dropThreshold){
        Z(i,j) = 1;
      } else {
        Z(i,j) = 0;
      };
      ave_impute(i,j) = ysub(i,j) + ave(i)*tau(j)*(1-Z(i,j))*(1- ysub_bin(i,j));
    };
  };


  // final output
  return List::create(
    _["ysub"] = ysub,
    _["zij"] = Z,
    _["libsize"] = tau,
    _["thetaj"] = thetaj,
    _["pi_vec"] = piTotal,
    _["y.impute"] = ave_impute,
    _["keep.gene"] = keepgene
  );
};


// [[Rcpp::export]]
List mcImpute_cpp(arma::mat data, bool preprocess = true, double eps=1e-12, double normfac=1, int insweep = 20, double tol=1e-4,
                  double decfac = 0.7, int min_count=1, int min_cells =1, bool verbose = false){
  // data: samples by genes
  arma::mat y, B, U, V, X, datasub, datasubnorm, resX;
  arma::vec geneFilter(data.n_cols), s, libsize;
  double alpha = 1.1*normfac, lambdaInit, lambda, f_current, f_previous;

  if (preprocess){
    //Removing BAD genes, median normalization and log-transformation
    for (unsigned int i =0; i < data.n_cols; ++i){ //genes
      int ybin = 0;
      for (unsigned int j=0; j < data.n_rows; ++j){ // cells
        if (data(j,i) >= min_count){
          ybin +=1;
        };
      };
      geneFilter(i) = ybin;
    };
    datasub = data.cols(find(geneFilter >= min_cells));
    libsize = arma_rowSums(datasub);
    datasubnorm = datasub;
    for (unsigned int i=0; i < datasubnorm.n_rows; i++){
      datasubnorm.row(i) = datasub.row(i) / sum(datasub.row(i)) * median(libsize);
    };
    y = log2(datasubnorm + ones(datasub.n_rows, datasub.n_cols));
  } else {
    y = data;
  };
  arma::mat ybin = BinMatrix(y);
  arma::mat x(y.n_rows, y.n_cols, fill::zeros);
  lambdaInit = decfac * MaxAbsMatrix(y);
  lambda = lambdaInit;

  arma::mat S(y.n_rows, y.n_cols, fill::zeros);

  f_current = frobeniousNorm(y - ybin % x) + lambda * oneNorm(x);
  while (lambda > lambdaInit*tol){
    for (int ins = 0; ins<insweep; ins++){
      f_previous = f_current;
      B = x + (1/alpha)* ybin % (y - ybin % x);
      arma::svd(U, s, V, B);
      for (unsigned int jj = 0; jj < s.size(); jj++){
        s(jj) = softThresh_cpp(s(jj), lambda/(2*alpha));
      };
      S.diag() = s;
      X = U * S * V.t();
      for (unsigned int ii = 0; ii < X.n_elem; ii ++){
        if (X(ii) < 0){
          X(ii) = 0;
        }
      };
      x = X;
      f_current = frobeniousNorm(y - ybin % x) + lambda * oneNorm(x);
      if ((abs(f_current - f_previous)/ abs(f_current + f_previous)) < tol){
        break;
      };
    };
    if (frobeniousNorm(y - ybin % x) < eps){
      break;
    };
    lambda = decfac * lambda;
  };

  if (preprocess){
    resX = round(exp2(x) - 1);
  } else {
    resX = round(x);
  }

  //----------
  return List::create(
    _["data"] = resX,
    _["geneFilter"] = geneFilter >= min_cells
  );
};
