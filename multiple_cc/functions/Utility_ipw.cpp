#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h> 
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]


Rcpp::List score_info_fun(const arma::mat& cov, const arma::vec& ftime,arma::vec& ltime,
                     const arma::vec& disease, const arma::vec delta,
                     const arma::vec& selvar, const arma::vec& subvar,
                     const arma::vec& kmfit,
                     const arma::vec& beta, const Rcpp::StringVector cal_Meat){
   
  int p = cov.n_cols;
  int n = cov.n_rows; 
  int t = ftime.n_elem;
  
  
  
  double alpha = accu(selvar % (1-delta))/accu(1-delta);
  arma::vec rho = delta + (1-delta)%selvar/alpha;
  
  arma::vec e_time = ftime.elem(find(disease==1));
  
  int etime_len = e_time.n_elem;
  
  arma::mat rho_matrix(n,etime_len);
  rho_matrix.each_col() = rho;
  
  
  /*Setup score matrix and information matrix  */
  arma::vec U_score(p);
  arma::mat Info(p,p);
  arma::mat Meat(p,p);
  
  
  /*Setup counting matrix, risk matrix and weight matrix  */
  
  arma::mat dN_matrix(n,etime_len);
  arma::mat risk_matrix(n,etime_len);
  arma::mat Wi_in_matrix(n,etime_len);
  arma::mat ri_matrix(n,etime_len);
  
  arma::vec bt(t);
  
  for(int i=0;i<t;i++){
    arma::vec bt_vec(n);
    arma::uvec b_loc1 = (ltime <= ftime[i]) % (ftime >= ftime[i]);
    bt_vec.elem(find( b_loc1==1 )).fill(1);
    bt[i] = sum(bt_vec)/n;
  }
  
  //  int i=2;
  
  
  arma::vec ftime_0(t+1);
  arma::vec kmfit_0(t+1);
  ftime_0[0] = 0; ftime_0(span(1,t)) = ftime;
  kmfit_0[0] = 1; kmfit_0(span(1,t)) = kmfit;
  
  
  
  
  
  for(int i=0;i<etime_len;i++){
    
    arma::vec n_risk(n);
    arma::uvec n_loc1 = (ftime >= e_time[i]) % (ltime <= e_time[i]);
    n_risk.elem( find(n_loc1 == 1)).fill(1);
    arma::uvec n_loc2 = (ftime <= e_time[i]) % (disease !=1)  % (delta !=0) %  (ltime <= ftime) ;
    
    n_risk.elem(find(n_loc2==1)).fill(1);
    risk_matrix.col(i) = n_risk;
    
    
    arma::vec n_event(n);
    n_event.elem(find(ftime==e_time[i])).fill(1);
    dN_matrix.col(i) = n_event % disease;
    
    arma::vec ri_elem(n);
    arma::uvec r_loc1 = (ltime <= e_time[i]) % (ftime >= e_time[i]);
    
    ri_elem.elem(find( r_loc1==1  )).fill(1);
    arma::uvec r_loc2 = (ltime <= e_time[i]) % (ftime < e_time[i]) % (delta != 0);
    ri_elem.elem(find(r_loc2==1  )).fill(1);
    
    ri_matrix.col(i) = ri_elem;
    
    arma::vec Wi_in_temp(n);
    
    // weight matrix
   if(ftime[0] == e_time[i] ){
      Wi_in_matrix.col(i) = arma::vec(n,fill::value(1));
      
    }else{
      arma::vec SXi = ftime_0.elem(find(ftime < e_time[i]));
      arma::vec bXi_time= ftime.elem(find(ftime < e_time[i]));
      

      
      
      Wi_in_temp.elem(find(ftime >= e_time[i])).fill(1);

      Rcpp::NumericVector ftime2 = NumericVector(ftime.begin(),ftime.end());
      Rcpp::NumericVector ftime_02 = NumericVector(ftime_0.begin(),ftime_0.end());
      Rcpp::NumericVector bXi_time2 = NumericVector(bXi_time.begin(),bXi_time.end());
      Rcpp::NumericVector SXi2 = NumericVector(SXi.begin(),SXi.end());
      
      
      Rcpp::LogicalVector bXi_loc = in(bXi_time2,ftime2);
      Rcpp::LogicalVector GXi_loc = in(SXi2,ftime_02);
      arma::uvec bXi_loc1 = as<arma::uvec>(wrap(bXi_loc));
      arma::uvec GXi_loc1 = as<arma::uvec>(wrap(GXi_loc));
      
      
      arma::vec bXi = bt.elem(find(bXi_loc1 ==1 ));
      arma::vec GXi = kmfit_0.elem(find(GXi_loc1==1 ));
      
      //int bXi_len = bXi.n_elem;
      //int GXi_len = GXi.n_elem;
      
      //Rprintf("bXi %i", bXi_len);
      //Rprintf("GXi %i",GXi_len);
      

      arma::vec Gt_d = kmfit_0.elem(find(e_time[i] == ftime));
      arma::vec bt_n = bt.elem(find(e_time[i] == ftime));

      
      Wi_in_temp.elem(find(ftime < e_time[i])) = GXi/Gt_d[0] * bt_n[0]/bXi;
      Wi_in_matrix.col(i) = Wi_in_temp;

    }
    
  }
  
  arma::mat Wi_matrix = ri_matrix % Wi_in_matrix;
  
  
  /* Calculate Score function*/
  
  arma::mat S1(etime_len,p);
  arma::mat S0(etime_len,p);
  arma::mat Ehat(etime_len,p);
  
  arma::mat exp_term = exp(beta.t()*cov.t()).t();
  arma::mat rwrho_mt = risk_matrix % Wi_matrix % rho_matrix;
  
  
  
  
  for(int i=0;i<p;i++){
    arma::mat S1_part1 = rwrho_mt.each_col() % cov.col(i);
    arma::mat S1_temp = S1_part1.t() * exp_term;
    S1.col(i) = S1_temp.col(0); /* t*1 */
    arma::mat S0_temp = rwrho_mt.t() * exp_term;
    S0.col(i) = S0_temp.col(0); /* t*1 */
    Ehat.col(i) = S1.col(i)/S0.col(i);/* t*1 */
  
    arma::mat Z_E_hat(n,etime_len);
    arma::mat cov_expand(n,etime_len);
    arma::mat E_expand(n,etime_len);
  
  
    cov_expand.each_col() = cov.col(i);
    E_expand.each_row() = (Ehat.col(i)).t();
    Z_E_hat = cov_expand - E_expand;
  
    arma::mat final_mt = (Z_E_hat % Wi_matrix % dN_matrix);
    final_mt.elem(find_nonfinite(final_mt)).fill(0);
    U_score(i) = accu(final_mt);
  
  }
  
  
  /* Calculate Information matrix */
  
  
  for(int i=0;i<p;i++){
    for(int j=i;j<p;j++){
      arma::mat S2_temp_part1 = rwrho_mt.each_col() % (cov.col(i) % cov.col(j));
      arma::mat S2_temp = S2_temp_part1.t() * exp_term;
      
      arma::mat f_term1 = S2_temp/S0.col(i);
      
      arma::mat f_term2 = S1.col(i) % S1.col(j)/(S0.col(i) % S0.col(i)); /* t*1*/
      arma::mat f_in = f_term1 - f_term2;
      arma::mat f_in_expand(n,etime_len);
      f_in_expand.each_row() = (f_in.col(0)).t();
      arma::mat final_mt = (f_in_expand % Wi_matrix % dN_matrix);
      final_mt.elem(find_nonfinite(final_mt)).fill(0);
      if(i==j){
        Info(i,i) = accu(final_mt);
      }else{
          Info(i,j) = accu(final_mt);
          Info(j,i) = accu(final_mt);
      }
    }  
  }
  
  
  /* Calculate Meat part in variance matrix  */
  if(cal_Meat(0) == "True"){
    
    
    arma::mat dLambda_smt1_expand(n,etime_len);
    dLambda_smt1_expand.each_row() = (1/S0.col(0)).t();
    arma::mat dLambda_sm = sum(dLambda_smt1_expand % Wi_matrix % dN_matrix,0);
    arma::mat dLambda_10(n,etime_len);
    dLambda_10.each_row() = dLambda_sm.row(0);
    
    arma::mat dN_dot(n,t);
    arma::mat risk_dot(n,t);
    
    
    for(int i=0;i<t;i++){
      arma::vec tdot_event(n);
      tdot_event.elem(find(ftime==ftime[i])).fill(1);
      dN_dot.col(i) = tdot_event % delta;
      arma::vec ndot_event(n);
      arma::uvec ndot_loc = (ftime >= ftime[i]) % (ltime <= ftime[i])  ;
      ndot_event.elem(find(ndot_loc == 1)).fill(1);
      risk_dot.col(i) = ndot_event;
      
    }
    
    arma::mat dLambda_dotp = sum(dN_dot,0)/sum(risk_dot,0);
    arma::mat dLambda_dot(n,t);
    dLambda_dot.each_row() = dLambda_dotp;
    
    arma::mat dM_doti = dN_dot - risk_dot % dLambda_dot;
    arma::mat exp_term_expand(n,etime_len);
    exp_term_expand.each_col()= exp_term;
    arma::mat dM_1i = dN_matrix - risk_matrix % exp_term_expand % dLambda_10;
    
    arma::mat pi_u = sum(risk_dot,0);
    
    Rcpp::NumericVector ftime2 = NumericVector(ftime.begin(),ftime.end());
    Rcpp::NumericVector e_time2 = NumericVector(e_time.begin(),e_time.end());
    Rcpp::LogicalVector e_time_loc = in(ftime2,e_time2);
    arma::uvec e_time_loc1 = as<arma::uvec>(wrap(e_time_loc));
    
    arma::mat bt_etime = bt.elem(find(e_time_loc1));
    
    
    /*Construct eta1 eta2 eta3 */
    arma::mat eta1_list(p,n);
    arma::mat eta2_list(p,n);
    arma::mat eta3_list(p,n);
    arma::mat eta4_list(p,n);
    
    for(int i=0;i<p;i++){
      arma::mat Z_E_hat(n,etime_len);
      arma::mat cov_expand(n,etime_len);
      arma::mat E_expand(n,etime_len);
      
      cov_expand.each_col() = cov.col(i);
      E_expand.each_row() = (Ehat.col(i)).t();
      Z_E_hat = cov_expand - E_expand;
      
      arma::mat eta1_temp = Z_E_hat % Wi_matrix % dM_1i;
      eta1_temp.elem(find_nonfinite(eta1_temp)).fill(0);
      
      eta1_list.row(i) = (sum(eta1_temp,1)).t();
      
      arma::mat eta4_first(n,etime_len);
      eta4_first.each_col() = (1- delta);
      arma::mat eta4_mt = eta4_first % Wi_matrix % risk_matrix % Z_E_hat % exp_term_expand % dLambda_10;
      eta4_mt.elem(find_nonfinite(eta4_mt)).fill(0);
      
      eta4_list.row(i) = (sum(eta4_mt ,1)).t();
      
      arma::vec qt(t);
      
      for(int j=1;j<t;j++){
        arma::uvec Xid = find(ftime < ftime[j]);
        arma::uvec tid = find(e_time >= ftime[j]);
        
        arma::mat qt_temp = Z_E_hat(Xid,tid) % Wi_matrix(Xid,tid) % dM_1i(Xid,tid);
        qt_temp.elem(find_nonfinite(qt_temp)).fill(0);
        
        qt(j) = accu(qt_temp);
        
      }
      
      arma::mat qt_expand(n,t);
      qt_expand.each_row() = qt.t();
      arma::mat pi_expand(n,t);
      pi_expand.each_row() = pi_u;
      
      arma::mat eta2_mt = qt_expand / pi_expand % dM_doti;
      eta2_mt.elem(find_nonfinite(eta2_mt)).fill(0);
      eta2_list.row(i) = (sum(eta2_mt,1)).t();
      
      //eta2_list.row(i) = (sum(qt_expand / pi_expand % dM_doti,1)).t();
      
      arma::mat phiI(n,etime_len);
      for(int j=0;j<etime_len;j++){
        arma::vec phiI_vec(n);
        phiI_vec.elem(find(ftime < e_time[j])).fill(1);
        phiI.col(j) = phiI_vec;
        
      }
      
      arma::vec eta_3(n);
      
      for(int j=0; j<n;j++){
        arma::uvec brack1_loc = (ltime[j] <= e_time) % (e_time <= ftime[j]);
        arma::uvec brack2_loc = (ltime[j] <= ftime) % (ftime <= ftime[j]);
        arma::vec brack_in1_vec(etime_len);
        arma::vec brack_in2_vec(n);
        
        brack_in1_vec(find(brack1_loc == 1)) = 1/(n*bt_etime(find(brack1_loc == 1)));
        brack_in2_vec(find(brack2_loc ==1)) = 1/(n*bt(find(brack2_loc==1)));
        
        arma::mat brack_in1 = repelem(brack_in1_vec,1,n);
        brack_in1 = brack_in1.t();
        arma::mat brack_in2 = repelem(brack_in2_vec,1,etime_len);
        
        arma::mat phi1 = Z_E_hat % Wi_matrix % dM_1i % phiI % (brack_in1-brack_in2);
        phi1.elem(find_nonfinite(phi1)).fill(0);
        eta_3[j] = accu(phi1);
        eta3_list.row(i) = eta_3.t();
        
      }
      
      
      
    }
    
    arma::mat eta1_list_sub = eta1_list.each_row() % subvar.t();
    arma::mat eta2_list_sub = eta2_list.each_row() % subvar.t();
    arma::mat eta3_list_sub = eta3_list.each_row() % subvar.t();
    arma::mat eta4_list_sub = eta4_list.each_row() % subvar.t();
    
    
    double alpha_star = accu(subvar)/n;
    
    arma::mat V0 = (eta1_list_sub + eta2_list_sub + eta3_list_sub) * (eta1_list_sub + eta2_list_sub + eta3_list_sub).t() * 1/alpha_star;
    arma::mat V1 =  eta4_list_sub * eta4_list_sub.t() * (1-alpha_star)/alpha_star * 1/alpha_star;
    
    Meat = V0+V1;
    
    

  }else{
    Meat = Meat;
  }
  
  
  

  return Rcpp::List::create(Rcpp::Named("score") = U_score,
                            Rcpp::Named("fisher")=Info,
                            Rcpp::Named("meat")=Meat );
  
  

  
  
}