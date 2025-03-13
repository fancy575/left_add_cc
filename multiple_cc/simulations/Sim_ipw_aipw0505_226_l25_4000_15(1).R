source("Data Generation.R")
source("cp_ipw_function.R")
source("cp_aipw_function.R")
library(cmprsk)

M <- 100

Monte_gen_est <- data.frame(matrix(ncol=4))
Monte_gen_var <- data.frame(matrix(ncol=4))
colnames(Monte_gen_est) <- c("ipw_Z1","ipw_Z2" ,"aipw_Z1","aipw_Z2")
colnames(Monte_gen_var) <- c("ipw_Z1","ipw_Z2","aipw_Z1","aipw_Z2")


for (m  in 1:M) {
 
  set.seed(m*3448+1)
   cp_cc_data <- cp_left_data_gen(size=5400,p=0.55,
                                    beta1=c(0.5,-0.5),beta2 = c(0.5,0.5),
                                    Z_dist=c("categorical","categorical"),phi=c(1,1),
                                    kappa=c(1,2),strata=2, 
                                     gen_parm=list(list(c(0.4,0.6),c(0.6,0.4))
                                               ,list(c(0.3,0.7),c(0.7,0.3) )),
                                    censor_b = 1.0736842, 
                                    trunc_rate = 6.405263,
                                    sub_prob=c(0.15, 0.15)) 

  ipw_fit <- tryCatch(
    {cp_ipw_fit(ite=1000,tol=1e-8,data=cp_cc_data[[3]],time="time",
                            left_time ="left_time",
                            event="disease",covname=c("C1_2","C2_2"),
                            selname = "cc_sel",subname = "sb_ind",
                            failcode = 1,cencode = 0)},
    error = function(e) list(root=c(NA,NA), cov=matrix(rep(NA,4),nrow=2) )
    
  )


  aipw_fit <- tryCatch(
    {cp_aipw_fit(ite=1000,tol=1e-8,data=cp_cc_data[[3]],time="time",
                left_time ="left_time",
                event="disease",covname=c("C1_2","C2_2"),
                selname = "cc_sel",subname = "sb_ind",
                failcode = 1,cencode = 0,mis.cov="C2",mis.cat=c(2),
                B=500)},
    error = function(e) list(root=c(NA,NA), cov=matrix(rep(NA,4),nrow=2) )
    
  )


  if(m==1){
    Monte_gen_est[m,] <- c(ipw_fit$root,aipw_fit$root)
    Monte_gen_var[m,] <- c(ipw_fit$cov[1,1],ipw_fit$cov[2,2],
                           aipw_fit$cov[1,1],aipw_fit$cov[2,2])
    
  }else{
    Monte_gen_est <- rbind(Monte_gen_est,c(ipw_fit$root,aipw_fit$root))
    Monte_gen_var <- rbind(Monte_gen_var,
                           c(ipw_fit$cov[1,1],ipw_fit$cov[2,2],
                             aipw_fit$cov[1,1],aipw_fit$cov[2,2]))
  }


save(Monte_gen_est,Monte_gen_var,file="Monte_1101_iapu0505_226_l25B5_4015strata(1).RData")

}