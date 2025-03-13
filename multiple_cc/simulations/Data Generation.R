
cp_left_data_gen <- function(size=2600,p=0.22,
                             beta1=c(1,1,0.5),beta2 = c(1,-1,0.5),
                             Z_dist=c("normal","categorical"),phi=c(1,1),
                             kappa=c(1,2),strata=2, gen_parm=list(c(1,0.1),list(c(0.25,0.5,0.25),c(0.2,0.4,0.4))),
                             censor_b = 1.607035, # p=0.12, b=3.535176, 2:4:4 (1:2:0)
                                                  #p=0.22 b=1.607035, 2:2:6 (1:2:0)
                             trunc_rate=4.458458, # 2:4:4 trunc_rate = 1.040040, 50%, trunc_rate=3.162162, 25%
                                                  #2:2:6 trunc_rate = 1.575576, 50%, trunc_rate=4.458458, 25%
                             sub_prob=c(0.15,0.15)){
  
  if(length(phi)!=length(kappa) | length(phi)!=strata){
    return(print("number of strata is not consistent number of phi and kappa setting"))
  }
  
  cp_data_ori <- list()
  cp_data_left <- list()
  cp_data_cc <- list()
  
  for(k in 1:strata){
    
    ### generate covariates
    Z_mt <- matrix(nrow=size,ncol=length(beta1))
    colnames(Z_mt) <- paste("X",1:length(beta1),sep="_")
    
    for(i in 1:length(Z_dist)){
      if(Z_dist[i]=="normal"){
        Z_mt[,i] <- rnorm(size,0,gen_parm[[i]][[k]])
        colnames(Z_mt)[i] <- paste("X",i,sep="_")
      }else if(Z_dist[i] == "categorical"){
        temp <- t(rmultinom(size,1,gen_parm[[i]][[k]] ))
        temp_cov <- matrix(0,nrow = size,ncol=length(gen_parm[[i]][[k]])-1)
        for(j in 2:length(gen_parm[[i]][[k]])){
          temp_cov[which(temp[,j]==1),j-1] <- 1
        }
        colnames(temp_cov) <- paste(paste("C",i,sep=""),2:length(gen_parm[[i]][[k]]),sep="_")
        
        Z_mt[,i:(length(gen_parm[[i]][[k]])+i-2)] <- temp_cov
        colnames(Z_mt)[i:(length(gen_parm[[i]][[k]])+i-2)] <- colnames(temp_cov)
      }
    }
    
    ### generate event time
    
    Z1_multi <- as.vector(beta1 %*% t(Z_mt))
    Z2_multi <- as.vector(beta2 %*% t(Z_mt))
    
    P1 <- 1- (1-p)^(exp(Z1_multi))
    P2 <- (1-p)^(exp(Z1_multi))
    U1 <- runif(size,0,1)
    cause <- ifelse(U1<=P1,1,2)
    
    T <- numeric(size)
    n1 <- length(which(cause==1))
    n2 <- length(which(cause==2))
    T[which(cause==1)] <- (-log(1- 1/p* (1-(1-runif(n1)*P1[which(cause==1)])^(exp(-Z1_multi[which(cause==1)]) ) ))/phi[k])^(1/kappa[k])
    T[which(cause==2)] <- -log(1-runif(n2))/(exp(Z2_multi[which(cause==2)] ))     
    
    ## generate censor 
    C <- runif(size,0,censor_b)
    X <- apply(cbind(C,T), 1, min)
    
    Delta <- ifelse(T>=C,0,1)
    
    cp_data <- data.frame(Z_mt,disease=ifelse(Delta==1,cause,0),time = X,cc_sel=1,sb_ind=1 )
    cp_data$strata <- k
    # prop.table(table(cp_data$disease))
    
    cp_data_ori[[k]] <- cp_data
    
    
    
    ### Generate left truncation time
    
    
    L <- rexp(size,trunc_rate)
    # prop.table(table(L<=X))
    
    cp_data <- data.frame(Z_mt,disease=ifelse(Delta==1,cause,0),left_time = L,time = X,cc_sel=1,sb_ind=1 )
    
    cp_data <- cp_data[which(cp_data$left_time<=cp_data$time),]
    
    cp_data_left[[k]] <- cp_data
    
    
    
    ### generate case-cohort data
    sub_sample <- sample(1:nrow(cp_data),nrow(cp_data)*sub_prob[k],replace = F)
    cp_casecohort <- cp_data
    cp_casecohort[sub_sample,'sb_ind'] <- 1
    cp_casecohort[-sub_sample,'sb_ind'] <- 0
    cp_casecohort[-sub_sample,][which(cp_casecohort[-sub_sample,]$disease == 0),"cc_sel"] <- 0  
    
    cp_data_cc[[k]] <- cp_casecohort
    
    
  }
  
  


  
  return(list(cp_data_ori,cp_data_left,cp_data_cc))
  
  
}
