#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h> 
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]


arma::cube impute(const arma::mat& cov, const arma::vec& ftime,const arma::vec& ltime,
                     const arma::vec& fstatus, const arma::vec& F_coef,
                     const arma::vec& S1_coef,const arma::vec& S2_coef,
                     const NumericVector p_inp,const int B){
  
  int p = cov.n_cols;
  int n = ftime.n_elem;
  int nc = F_coef.n_elem;
  int c = p_inp.length()-1;
  
  arma::cube imput_ar(n,p+9,B);

  
  Function rmultinom("rmultinom");
  Function pexp("pexp");
  Function dexp("dexp");

  
  
  for(int i=0;i<B;i++){
    NumericMatrix gen_binom = rmultinom(n,1,p_inp);
    arma::mat temp_binom_re = arma::mat(gen_binom.begin(),gen_binom.nrow(),gen_binom.ncol(),false );
    arma::mat temp_binom = temp_binom_re.t();
    
    
    /* initial setup and parameters to calculate probability */
    imput_ar(span(0,n-1),span(p-c,p-1),span(i,i)) = temp_binom(span(0,n-1),span(1,c));
    imput_ar(span(0,n-1),span(0,p-c-1),span(i,i)) = cov(span(0,n-1),span(0,p-c-1));
    imput_ar(span(0,n-1),span(p),span(i,i)) = ftime;
    imput_ar(span(0,n-1),span(p+1),span(i,i)) = ltime;
    
    imput_ar(span(0,n-1),span(p+2),span(i,i)) = fstatus;
    arma::mat new_cov = imput_ar(span(0,n-1),span(0,p-1),span(i,i));
    imput_ar(span(0,n-1),span(p+3),span(i,i)) =  new_cov * F_coef(span(1,nc-1));
    imput_ar(span(0,n-1),span(p+4),span(i,i)) =  new_cov * S1_coef(span(1,nc-1));
    imput_ar(span(0,n-1),span(p+5),span(i,i)) =  new_cov * S2_coef(span(1,nc-1));
    
    
    arma::mat linkF = exp(imput_ar(span(0,n-1),span(p+3),span(i,i)));
    arma::mat linkS1 = exp(imput_ar(span(0,n-1),span(p+4),span(i,i)));
    arma::mat linkS2 = exp(imput_ar(span(0,n-1),span(p+5),span(i,i)));

    NumericVector F_p_temp = pexp(ftime % linkF,F_coef[0]);
    NumericVector S1_p_temp = pexp(ftime % linkS1,S1_coef[0]);
    NumericVector S2_p_temp = pexp(ftime % linkS2,S2_coef[0]);
    
    NumericVector F_p_temp_l = pexp(ltime % linkF,F_coef[0]);
    NumericVector S1_p_temp_l = pexp(ltime % linkS1,S1_coef[0]);
    NumericVector S2_p_temp_l = pexp(ltime % linkS2,S2_coef[0]);
    
    
    NumericVector F_d_temp = dexp(ftime % linkF,F_coef[0]);
    NumericVector S1_d_temp = dexp(ftime % linkS1,S1_coef[0]);
    NumericVector S2_d_temp = dexp(ftime % linkS2,S2_coef[0]);
    
    arma::mat F_p = (1-F_p_temp)/(1-F_p_temp_l); 
    arma::mat S1_p = (1-S1_p_temp)/(1-S1_p_temp_l); 
    arma::mat S2_p = (1-S2_p_temp)/(1-S2_p_temp_l);
    arma::mat F_d = F_d_temp/(1-F_p_temp_l); 
    arma::mat S1_d = S1_d_temp/(1-S1_p_temp_l); 
    arma::mat S2_d = S2_d_temp/(1-S2_p_temp_l);
    
    
    imput_ar(span(0,n-1),span(p+6),span(i,i)) =  F_p % S2_p% S1_d % linkS1; /* FS2s1  */
    imput_ar(span(0,n-1),span(p+7),span(i,i)) =  F_p % S1_p% S2_d % linkS2; /* FS1s2  */
    imput_ar(span(0,n-1),span(p+8),span(i,i)) =  F_d % S2_p% S1_p % linkF;  /* fS2S1  */
  }    
    
  return imput_ar;
  
  
}

