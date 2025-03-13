library(flexsurv)
library(nnet)
sourceCpp("monte_imput.cpp")
monte_imput <- function(ftime,ltime,fstatus,cov,selvar,B=200,
                        mis.cov="C",cencode=0,failcode=1){
  
  
  alpha <- sum(selvar*( 1-ifelse(fstatus==cencode,0,1)))/ sum(( 1-ifelse(fstatus==cencode,0,1)))
  rho <- ifelse(fstatus==cencode,0,1)  + (1-ifelse(fstatus==cencode,0,1))*selvar/alpha 
  
  
  
  glm_dt <- data.frame(cov[which(selvar==1),],rho=rho[which(selvar==1)],
                       ftime=ftime[which(selvar==1)],ltime=ltime[which(selvar==1)],
                       fstatus=fstatus[which(selvar==1)])
#  glm_dt$diff_time <- glm_dt$ftime - glm_dt$ltime
  glm_form <- paste(mis.cov,"~1",sep="")
  
#  dtime <- ftime - ltime
  
  
  glm_bi <- multinom(as.formula(glm_form),data=glm_dt,weight=rho)
  
  glm_dt$censor_ind <- ifelse(glm_dt$fstatus==cencode,1,0)
  glm_dt$cause1 <- ifelse(glm_dt$fstatus==failcode,1,0)
  glm_dt$cause2 <- ifelse(glm_dt$fstatus!=cencode & glm_dt$fstatus!=failcode,1,0  )
  
  
  
#  aft_form_censor <- as.formula(paste("Surv(ltime,ftime,censor_ind)","~",paste(names(cov)[1:(ncol(cov)-1)],collapse = "+"),sep=""))
#  aft_form_cause1 <- as.formula(paste("Surv(ltime,ftime,cause1)","~",paste(names(cov)[1:(ncol(cov)-1)],collapse = "+"),sep=""))
#  aft_form_cause2 <- as.formula(paste("Surv(ltime,ftime,cause2)","~",paste(names(cov)[1:(ncol(cov)-1)],collapse = "+"),sep=""))
  
  flex_form_censor <- as.formula(paste("Surv(ltime,ftime,censor_ind)","~",paste(names(cov)[1:(ncol(cov)-1)],collapse = "+"),sep=""))
  flex_form_cause1 <- as.formula(paste("Surv(ltime,ftime,cause1)","~",paste(names(cov)[1:(ncol(cov)-1)],collapse = "+"),sep=""))
  flex_form_cause2 <- as.formula(paste("Surv(ltime,ftime,cause2)","~",paste(names(cov)[1:(ncol(cov)-1)],collapse = "+"),sep=""))
  
  
  F_censor <- flexsurvreg(flex_form_censor,weight=rho,data=glm_dt,dist = "exp",control=list(fnscale=1000))
  S_cause1 <- flexsurvreg(flex_form_cause1,weights=rho,data=glm_dt,dist = "exp",control=list(fnscale=1000))
  S_cause2 <- flexsurvreg(flex_form_cause2,weights=rho,data=glm_dt,dist = "exp",control=list(fnscale=1000))
  
  
  F_coef <- c(exp(F_censor$coefficients[1]),
              F_censor$coefficients[2:length(F_censor$coefficient)]) ## rate
  
  S1_coef <- c(exp(S_cause1$coefficients[1]),
               S_cause1$coefficients[2:length(S_cause1$coefficients)]) ## rate
  
  
  S2_coef <- c(exp(S_cause2$coefficients[1]),
               S_cause2$coefficients[2:length(S_cause2$coefficients)]) ## rate
  p_inp <- as.matrix(predict(glm_bi,type="probs"))[1,]
  if(length(p_inp)==1){
    p_inp <- c(1-p_inp,p_inp)
  }
  
  
#  my_imput <- impute(as.matrix(cov)[,1:(ncol(cov)-1)],ftime,fstatus,
#                 F_coef,S1_coef,S2_coef,p_inp,B)
  
  my_imput <- impute(as.matrix(cov)[,1:(ncol(cov)-1)],ftime,ltime,fstatus,
                     F_coef,S1_coef,S2_coef,p_inp,B)
  
  
  return(my_imput)
  
  
  
}


