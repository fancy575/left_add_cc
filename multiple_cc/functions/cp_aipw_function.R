library(survival)
library(Rcpp)
source("monte_imput.R")
sourceCpp("cond_exp.cpp")
sourceCpp("Utility_aug.cpp")


cp_aipw_fit <- function(ite,tol,data,time,left_time,event,covname,
                                selname,subname,failcode=1,cencode=0,
                        mis.cov="C",mis.cat=c(2,3) ,B=100){
  
  ### data preparation 
  ftime_strata <- list()
  ltime_strata <- list()
  fstatus_strata <- list()
  cov_strata <- list()
  selvar_strata <- list()
  subvar_strata <- list()
  G_km_fit_strata <- list()
  disease_strata <- list()
  delta_strata <- list()
  my_imput_strata <- list()
  
  n_total <- 0
  
  
  for(k in 1:length(data)){
    
    
    ftime <- data[[k]][,time]
    ltime <- data[[k]][,left_time]
    fstatus <- data[[k]][,event]
    cov <- data[[k]][,covname]
    selvar <- data[[k]][,selname]
    subvar <- data[[k]][,subname]
    
    
    ### reorder the time
    new_order <- order(ftime)
    ftime <- ftime[new_order]
    ltime <- ltime[new_order]
    cov <- cov[new_order,]
    fstatus <- fstatus[new_order]
    selvar <- selvar[new_order]
    subvar <- subvar[new_order]
    
    
    
    
    ### fit kaplan meier
    
    G_km_data <- as.data.frame(cbind(ltime,ftime,fstatus))
    names(G_km_data) <- c("left_time","time","event")
    G_km_data$event <- ifelse(G_km_data$event!=0,1,0)
    
    G_km <- survfit(Surv(left_time,time,event)~1,data=G_km_data)
    G_km_sum <- summary(G_km,ftime)
    cols <- lapply(c(2,6), function(x) G_km_sum[x])
    G_km_fit <- do.call(data.frame,cols)
    
    disease <- ifelse(fstatus==failcode,1,0)
    delta <- ifelse(fstatus!=cencode,1,0)
    
    
    
    new_cov <- numeric(nrow(cov))
    for(i in 1:(length(mis.cat))){
      new_cov[which(cov[,paste(mis.cov,mis.cat[i],sep="_")]==1)] <- mis.cat[i]
    }
    
    
    cov2 <- cbind(cov[,which(!colnames(cov) %in% paste(mis.cov,mis.cat,sep="_"))],
                  cov[,which(colnames(cov) %in% paste(mis.cov,mis.cat,sep="_") )],
                  new_cov)
    cov2 <- as.data.frame(cov2)
    names(cov2)[1:(ncol(cov2)- length(mis.cat)-1 )] <- names(cov)[which(!names(cov) %in% paste(mis.cov,mis.cat,sep="_") )]
    
    names(cov2)[(ncol(cov2)- length(mis.cat)):(ncol(cov2)-1)] <- names(cov)[which(names(cov) %in% paste(mis.cov,mis.cat,sep="_") )]
    
    names(cov2)[ncol(cov2)] <- mis.cov
    
    
    
    ftime_strata[[k]] <- ftime
    ltime_strata[[k]] <- ltime
    fstatus_strata[[k]] <- fstatus
    cov_strata[[k]] <- cov
    selvar_strata[[k]] <- selvar
    subvar_strata[[k]] <- subvar
    G_km_fit_strata[[k]] <- G_km_fit
    disease_strata[[k]] <- disease
    delta_strata[[k]] <- delta
    
    n_total <- n_total + nrow(data[[k]])
    
    
    my_imput <- monte_imput(ftime =ftime, ltime=ltime,fstatus=fstatus,cov=cov2,
                            selvar=selvar,B=B,
                            mis.cov=mis.cov,cencode=0,failcode=1)
    
    my_imput_strata[[k]] <- my_imput
    
  }

  
  diff <- 1000
  maxite <- 0
  beta <- numeric(ncol(cov))
  
  while(diff > tol){
    
    score <- 0
    fisher <- 0
    cond_strata <- list()
    
    for(k in 1:length(data)){
      
      cond <- cond_expect(my_imput_strata[[k]],beta)
      cond_strata[[k]] <- cond
      
      
      score_info <- score_infoaug_fun(cov=as.matrix(cov_strata[[k]]),ftime=ftime_strata[[k]], ltime=ltime_strata[[k]],
                                      disease=disease_strata[[k]],
                                      delta=delta_strata[[k]],selvar=selvar_strata[[k]],subvar=subvar_strata[[k]],
                                      kmfit=G_km_fit_strata[[k]]$surv,
                                      beta=beta,
                                      cond,cal_Meat=c("False"))
      score <- score + score_info$score-score_info$aug
      fisher <- fisher + score_info$fisher
      
    }
    
    

    dev <- solve(-fisher)
    
    beta1 <- beta - dev %*% matrix(score,nrow=length(beta))
    beta1 <- c(beta1)
    
    diff <- sum(abs(beta1-beta))
    
    if(diff > tol) {
      beta <- beta1
      maxite <- maxite + 1
    }
    
    if(maxite > ite){
      print (paste("root =",beta,"do not converge",sep=" "))
      return (list(root = rep(NA,ncol(cov)),iteration = maxite,eps = diff))
    }
    
    
  }
  
  final_fisher <- 0
  final_meat <- 0
  
  for(k in 1:length(data)){
    final_score_info <- score_infoaug_fun(cov=as.matrix(cov_strata[[k]]),ftime=ftime_strata[[k]],ltime=ltime_strata[[k]],
                                          disease=disease_strata[[k]],
                                          delta=delta_strata[[k]],selvar=selvar_strata[[k]],subvar=subvar_strata[[k]],
                                          kmfit=G_km_fit_strata[[k]]$surv,
                                          beta=beta,
                                          cond_strata[[k]],cal_Meat=c("True"))
    
    ql <- nrow(data[[k]])/n_total
    
    final_fisher <- final_fisher + ql * final_score_info$fisher/nrow(data[[k]])
    final_meat <- final_meat + ql * final_score_info$meat/nrow(data[[k]])
    
  }
  
  
  
  covariance <- solve(final_fisher) %*% final_meat %*% solve(final_fisher)/n_total
  colnames(covariance) <- rownames(covariance) <- covname
  names(beta) <- covname
  
  return (list(root = beta,iteration = maxite,
               eps = diff,cov=covariance,infomation = final_fisher
  ))
  
  
  
}
