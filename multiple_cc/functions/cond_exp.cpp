#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h> 
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]

arma::mat cond_expect(const arma::cube& imput_dt, const arma::vec& beta){
  int n = size(imput_dt)[0];
  int p = size(imput_dt)[1]-9;
  int b = size(imput_dt)[2];
  
  arma::mat exp_mt(n,p+1);
  
  arma::vec disease = imput_dt(span(0,n-1),span(p+2,p+2),span(0,0));
  
  arma::vec censor(n);
  censor.elem(find(disease==0)).fill(1);
  
  arma::vec cause1(n);
  cause1.elem(find(disease==1)).fill(1);
  
  arma::vec cause2(n);
  cause2.elem(find(disease==2)).fill(1);
  
  arma::mat expZ(n,p);
  arma::mat expE(n,1);
  

  for(int i=0;i<n;i++){
    arma::mat cov = imput_dt(span(i,i),span(0,p-1),span(0,b-1));
    arma::vec exp_form = exp(cov.t() * beta);
    
    arma::vec FS2s1 = imput_dt(span(i,i),span(p+6,p+6),span(0,b-1));
    arma::vec FS1s2 = imput_dt(span(i,i),span(p+7,p+7),span(0,b-1));
    arma::vec fS2S1 = imput_dt(span(i,i),span(p+8,p+8),span(0,b-1));
    
    double exp_temp = accu(exp_form % FS2s1 * cause1[i])/accu(FS2s1) + 
      accu(exp_form  % FS1s2 * cause2[i])/accu(FS1s2) +
      accu(exp_form  % fS2S1 * censor[i])/accu(fS2S1) ; 
    expE(i,0) = exp_temp;
    
    for(int j=0;j<p;j++){
      arma::vec Z_temp = imput_dt(span(i,i),span(j,j),span(0,b-1));
      double expZ_temp = accu(exp_form % Z_temp % FS2s1 * cause1[i])/accu(FS2s1) + 
        accu(exp_form % Z_temp % FS1s2 * cause2[i])/accu(FS1s2) +
        accu(exp_form % Z_temp % fS2S1 * censor[i])/accu(fS2S1) ; 
      expZ(i,j) = expZ_temp;
      
    }
    
  }
  exp_mt(span(0,n-1),span(0,p-1)) = expZ;
  exp_mt.col(p) = expE;


  
  
  return exp_mt;
  
  
  
  
}

