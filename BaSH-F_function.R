rm(list=ls())
library(MASS)
library(pscl)
library(MCMCpack)
library(GIGrvg)
library(statmod)
library(LCA)
library(LaplacesDemon)
library(Matrix)
library(magic)
library(psych)
library(mvtnorm)

BaSHF <- function(Mod,Sim,mcmc,burn){
  
  dname <- paste("Y_",Mod,".RData",sep="")  ## data
  load(dname)
  Y=Y
  Omega_true <- Omega_true
  
  M <- M                     ## Number of Groups                       
  p <- p                     ## Number of Covariates                      
  K_choice <- 3:8            ## Range of factors
  K_model <- K_choice[Sim]   ## Sim defines the Kth factor       
  n <- n                     ## Number of sample in each group
  N <- sum(n)
  gr_ind <- 1:M
  group <- c(rep(1,n[1]),rep(2,n[2]),rep(3,n[3]))
                      
  prob_a <- 1                      
  prob_b <- 1
  tau_alpha <- 1                      
  tau_beta <- 1
  c <- 0.1                         
  d <- 0.1
  theta <- 1                      
  a <- rep(0.2,M)                      
  tau <- 1
  Nchain <- 3
  thin <- 10
  nthin <- mcmc/thin
  nSwap <- 5
  
  Indicator <- array(NA,dim=c(p,K_model,nthin,Nchain))
  error_variance <- array(NA,dim=c(M,p,nthin,Nchain))
  mixture_prob <- array(NA,dim=c(nthin,Nchain))
  LAMBDA <- array(NA,dim=c(p,K_model,M,nthin,Nchain))
  TAU <- array(NA,dim=c(nthin,Nchain))
  W <- array(NA,dim=c(p,K_model,nthin,Nchain))
  NU <- array(NA,dim=c(nthin,M,Nchain))
  ETA <- array(NA,dim=c(N,K_model,nthin,Nchain))
  count <- array(0,dim=c(mcmc,p,Nchain))
  
  for(Chain in 1:Nchain)
  { 
    set.seed(999*Chain)
    eta <- matrix(c(rnorm(N*K_model,0,1)),ncol=K_model)
    nu <- rlnorm(M,0,1)
    lambda <- array(c(rnorm(p*K_model*M,0,1)),dim=c(p,K_model,M))
    w <- matrix(c(rnorm(p*K_model,0,1)),ncol=K_model)+rnorm(p*K_model,0,2)
    sigma <- matrix(c(rlnorm(p*M,0,1)),ncol=p)
    Z <- matrix(c(rbinom(p*K_model,1,0.5)),ncol=K_model)
    Y_new <- Y[group==Chain,]
    fit <- principal(Y_new,nfactors= K_model)
    lam_temp <- fit$loadings[,1:K_model]+rnorm(p*K_model,0,2)
    lambda_start <- array(NA,dim=c(p,K_model,M))
    for(m in 1:M){lambda_start[,,m] <- fit$loadings[,1:K_model]}
    w_start <- lam_temp
    lambda <- lambda_start
    w <- w_start
    epsilon <- array(NA,dim=c(N,p)) 
    var_term <- array(NA,M)
    temp <- array(NA,M)
    const <- array(NA,M)
    W_varterm <- array(NA,dim=M)
    W_meanterm <- array(NA,M)
    temp1 <- array(NA,M)
    temp2 <- array(NA,M)
    u1 <- array(NA,M)
    u2 <- array(NA,M)
    u3 <- array(NA,M)
    u4 <- array(NA,M)
    u5 <- array(NA,M)
    u6 <- u7 <- u9 <- c()  
    prod1 <- array(NA,M)
    prod2 <- array(NA,M)
    sum_term1 <- array(NA,M)
    Omega_mat <- array(NA,dim=c(p,p,M)) 
    Log_swap <- array(NA,M)
    Log <- array(NA,M)
    
    store_count <- 0
    start_time <- Sys.time()
    for(rep in 1:mcmc)
    {    
     prob <- rbeta(1,prob_a+sum(Z),prob_b+(p*K_model)-sum(Z))
      for(j in 1:p)
      {                   
        for (k in sample(K_model))
        {              
          
          log_mZ <- function(z){
            if(sum(z)!=0){
              u6 <- array(NA,dim=c(sum(z),sum(z),M))
              u7 <- array(NA,dim=c(sum(z),1,M))
              u9 <- array(NA,dim=c(sum(z),sum(z),M))
              for( m in 1:M)
              {
                I_k <- diag(1,sum(z))    
                u1[m] <- (n[m]+sum(z)-sum(z))*log(2*pi)
                u2[m] <- n[m]*log(sigma[m,j])+sum(z)*log(nu[m])
                u3[m] <- (t(Y[group==m,j])%*%Y[group==m,j])/sigma[m,j]
                temp_u4 <- ((t(eta[group==m,as.logical(z)])%*%eta[group==m,as.logical(z)])/sigma[m,j])+(I_k/nu[m])
                u4[m] <- log(det(temp_u4))
                u5[m] <- ((t(Y[group==m,j])%*%eta[group==m,as.logical(z)])/sigma[m,j])%*%solve(temp_u4)%*%((t(eta[group==m,as.logical(z)])%*%Y[group==m,j])/sigma[m,j])
                u6[,,m] <- (nu[m]^(-1))*I_k-((nu[m]^(-1))*solve(temp_u4)*(nu[m]^(-1)))
                u7[,,m] <- ((nu[m]^(-1))*solve(temp_u4)%*%((t(eta[group==m,as.logical(z)])%*%Y[group==m,j])/sigma[m,j]))
              }
              u8 <-  t(apply(u7,c(1,2),sum))%*%solve((I_k/tau)+apply(u6,c(1,2),sum))%*%(apply(u7,c(1,2),sum))  
              I2 <- -0.5*(sum(u1+u2+u3+u4-u5))-0.5*(sum(z)-sum(z))*log(2*pi)-0.5*sum(z)*log(tau)-0.5*log(det((I_k/tau)+
                                                                                                               apply(u6,c(1,2),sum)))+0.5*u8
            }else{
              for( m in 1:M)
              {
                u1[m] <- (n[m]+sum(z)-sum(z))*log(2*pi)
                u2[m] <- n[m]*log(sigma[m,j])+sum(z)*log(nu[m])
                u3[m] <- (t(Y[group==m,j])%*%Y[group==m,j])/sigma[m,j]
              }
              I2 <- -0.5*(sum(u1+u2+u3))
            }
            return(I2)
          }
          
          z0 <- z1 <- Z[j,]
          z0[k] <- 0
          z1[k] <- 1
          
          logA <- log_mZ(z0)+sum(z0)*log(prob)+(K_model-sum(z0))*log(1-prob)
          logB <- log_mZ(z1)+sum(z1)*log(prob)+(K_model-sum(z1))*log(1-prob)
          p_star <- 1/(1+exp((logA-logB)))
          Z[j,k] <- rbinom(1,1,p_star)     
          
        } 
       
           if(sum(which(Z[j,]==1))!=0 & sum(which(Z[j,]==0))!=0)    
        {
          samp1 <- which(Z[j,]==1)       
          samp2 <- which(Z[j,]==0)       
          K1 <- ifelse(length(samp1)>1,sample(samp1,1),samp1)   
          K0 <- ifelse(length(samp2)>1,sample(samp2,1),samp2)   
          Z_old <- Z[j,]     
          Z_swap <- Z[j,]               
          Z_swap[K1] <- 0    
          Z_swap[K0] <- 1    
          
          log_old <- log_mZ(Z_old)+sum(Z_old)*log(prob)+(K_model-sum(Z_old))*log(1-prob)    
          log_new <- log_mZ(Z_swap)+sum(Z_swap)*log(prob)+(K_model-sum(Z_swap))*log(1-prob)  
          p_swap<- 1/(1+exp((log_old-log_new)))  
          
          if(runif(1)<p_swap){
            count[rep,j,Chain] <- 1+count[rep,j,Chain]
            Z[j,] <- Z_swap              
          } 
        }
        
        
        
         if(sum(Z[j,]!=0)){
          u6 <- array(NA,dim=c(sum(Z[j,]),sum(Z[j,]),M))
          u7 <- array(NA,dim=c(sum(Z[j,]),1,M))
          for(m in 1:M)   
          {
            u6[,,m] <- (nu[m]^(-1))*diag(1,sum(Z[j,]))-(nu[m]^(-1))*solve(((t(eta[group==m,as.logical(Z[j,])])%*%eta[group==m,as.logical(Z[j,])])/sigma[m,j])+(diag(1,sum(Z[j,]))/nu[m]))*(nu[m]^(-1))
            u7[,,m] <- ((nu[m]^(-1))*solve(((t(eta[group==m,as.logical(Z[j,])])%*%eta[group==m,as.logical(Z[j,])])/sigma[m,j])+(diag(1,sum(Z[j,]))/nu[m]))%*%((t(eta[group==m,as.logical(Z[j,])])%*%Y[group==m,j])/sigma[m,j]))
          }
          
          W_var <- solve((diag(1,sum(Z[j,]))/tau)+apply(u6,c(1,2),sum))
          W_mean <- t(apply(u7,c(1,2),sum))%*%W_var 
          
          w[j,as.logical(Z[j,])] <- mvrnorm(1,W_mean,W_var)     
          
          for (k in 1:K_model){ w[j,k] <- ifelse(Z[j,k]==0,0,w[j,k])}
        }else{
          w[j,] <- rep(0,length(Z[j,]))
        }   
        
        
        for(m in 1:M) {  
          if(sum(Z[j,]!=0)){
            lambda_var <- solve(((t(eta[group==m,as.logical(Z[j,])])%*%eta[group==m,as.logical(Z[j,])])/sigma[m,j])+(diag(1,sum(Z[j,]))/nu[m]))
            lambda_mean <- (((t(Y[group==m,j])%*%eta[group==m,as.logical(Z[j,])])/sigma[m,j])+(t(w[j,as.logical(Z[j,])])/nu[m]))%*%lambda_var ## marginal mean for lambda[j,k,m]
            lambda[j,as.logical(Z[j,]),m] <- mvrnorm(1,lambda_mean, lambda_var)
            for (k in 1:K_model){ lambda[j,k,m] <- ifelse(Z[j,k]==0,0,lambda[j,k,m])}  
          }else{
            lambda[j,,m] <- rep(0,length(Z[j,]))
          }
        }    
      }      
      
     for(m in 1:M)       
      {
        for(i in 1:n[m])  ## for each sample of each group
        {
          var_eta <- solve(diag(1,K_model)+t(lambda[,,m])%*%solve(diag(sigma[m,]))%*%lambda[,,m])
          mean_eta <- var_eta%*%(t(lambda[,,m])%*%solve(diag(sigma[m,]))%*%Y[(group==m),][i,])
          eta[(group==m),][i,] <- mvrnorm(1,mean_eta,var_eta)
        }
      }
  
      for(m in 1:M){nu[m] <- rigamma(1,0.5*sum(Z)+0.5,0.5*sum((lambda[,,m]-w)^2)+(theta*a[m])^(-1))} 
      for(m in 1:M){a[m] <- rigamma(1,1,1+((theta)^(-1)/nu[m]))} 
      theta <- rigamma(1,0.5*(M+2),1+sum((a)^(-1)/nu))
      
      sig_temp <- array(NA,dim=c(M,p))
      for(j in 1:p)
      {
        for(m in 1:M) 
        {
          sig_temp[m] <- t(Y[which(group==m),j]-eta[which(group==m),]%*%lambda[j,,m])%*%(Y[which(group==m),j]-eta[which(group==m),]%*%lambda[j,,m])
          sigma[m,j] <- rigamma(1,0.5*n[m]+c,d+(0.5*(sig_temp[m])))
        }
        }
      
      tau  <- rigamma(1,0.5*sum(Z)+tau_alpha,tau_beta+0.5*sum(w^2))

      if(rep%%thin==0){
        store_count <- store_count+1
        Indicator[,,store_count ,Chain] <- Z
        W[,,store_count ,Chain] <- w
        NU[store_count ,,Chain] <- nu
        TAU[store_count ,Chain] <- tau
        error_variance[,,store_count ,Chain] <- sigma
        ETA[,,store_count ,Chain] <- eta
        mixture_prob[store_count ,Chain] <- prob
        LAMBDA[,,,store_count ,Chain] <- lambda
      } 
    }  
  }  
  
  swap_count <- apply(count,c(2,3),sum) 
  lambda_mean <- apply(LAMBDA[,,,-(1:burn),],c(1,2,3),mean)
  sigma_mean <- apply(error_variance[,,-(1:burn),],c(1,2),mean)
  nu_est <- apply(NU[-(1:burn),,],2,mean)
  tau_est <- mean(TAU[-(1:burn),])
  eta_mean <- apply(ETA[,,-(1:burn),],c(1,2),mean)
  mcmc <- dim(LAMBDA)[4]
  
  Omega_inv <- Omega_calc <- array(NA,dim=c(p,p,M,mcmc,Nchain))                  
  Omega_calc_det <- array(NA,dim=c(M,mcmc,Nchain))
  for(Chain in 1:Nchain)
  {
    for(rep in 1:mcmc)
    {
      for(m in 1:M)
      {
        Omega_calc[,,m,rep,Chain] <- LAMBDA[,,m,rep,Chain]%*%t(LAMBDA[,,m,rep,Chain])+diag(error_variance[m,,rep,Chain])
        Omega_inv[,,m,rep,Chain] <- solve( Omega_calc[,,m,rep,Chain])
        Omega_calc_det[m,rep,Chain] <- log(det(Omega_calc[,,m,rep,Chain] ))
      }
    }
  }
  
  Omega_est_inv <- array(NA,dim=c(p,p,M))
   for(m in 1:M){Omega_est_inv[,,m] <- solve(apply(Omega_inv[,,m,-(1:burn),],c(1,2),mean))   }
  
  Loss_group <- array(NA,M)
  for(m in 1:M){Loss_group[m] <- n[m]*(tr(solve(Omega_true[,,m])%*%Omega_est_inv[,,m])-log(det(solve(Omega_true[,,m])%*%Omega_est_inv[,,m]))-p)}
  Loss_all <- sum(Loss_group)/N;Loss_all
  
  logvalue <- array(NA,dim=c(mcmc,M,Nchain))
  for(Chain in 1:Nchain)
  {
    for(m in 1:M)
    {
      for(rep in 1:mcmc)
      {
        lden <- function(x){dmvnorm(x,rep(0,p),Omega_calc[,,m,rep,Chain],log=T)}
        temp <- apply(Y[group==m,],1,lden)
        logvalue[rep,m,Chain] <- sum(temp)
      }
    }
  }
  
  (Exp_dev <- -2*mean(apply(apply(logvalue[-(1:burn),,],c(1,2),mean),1,sum)))
  Omega_est <- array(NA,dim=c(p,p,M))
  for(m in 1:M){Omega_est[,,m] <- apply(Omega_calc[,,m,-(1:burn),],c(1,2),mean)}
  dev_hat <- c()
  for(m in 1:M)
  {
    dev_fun <- function(x){dmvnorm(x,rep(0,p),Omega_est_inv[,,m],log=T)}
    dev_gr <- apply(Y[group==m,],1,dev_fun)
    dev_hat <- c(dev_hat,dev_gr)    
  }
  sum_devhat <- -2*sum(dev_hat);sum_devhat
  pD <- Exp_dev-sum_devhat;pD
  DIC <- sum_devhat+2*pD;DIC
  
  loglik <- c()      
  for(m in 1:M)
  {
    loglik_gr <- array(NA,dim=c(n[m],mcmc,Nchain))
    for(Chain in 1:Nchain)
    { 
      for(rep in 1:mcmc)
      {
        loglik_fun <- function(x){dmvnorm(x,rep(0,p),Omega_calc[,,m,rep,Chain],log=T)}
        loglik_gr[,rep,Chain] <- apply(Y[group==m,],1,loglik_fun)
      }
    }
    loglik_gr_mean <- apply(loglik_gr[,-(1:burn),],c(1,2),mean)
    loglik <- rbind(loglik,loglik_gr_mean)    
  }
  sum_loglik <- apply(loglik,2,sum)
  CPO <- 1/apply(exp(-loglik),1,mean)
  LPML <- sum(log(CPO));LPML
  
  out=list(DIC=DIC,pD=pD,Dev=sum_devhat,Expdev=Exp_dev,LPML=LPML,Loss_all=Loss_all,
           Omegaest=Omega_est_inv,Sigmaest=sigma_mean,Nuest=nu_est,Tauest=tau_est,Swap=swap_count,
             Z=Indicator,Lambda=LAMBDA,sigma=error_variance,W=W,prob=mixture_prob,NU=NU,TAU=TAU,logvalue=logvalue,count=count)
 
  return(out)
}  

## Code run and saving output files

b=BaSHF(Mod,Sim,mcmc,burn)

K_choice <- 3:8            ## Range of factors
K_model <- K_choice[Sim]   ## Sim defines the Kth factor
FileName <- paste("BaSF_Factor",K_model,"_output_File",Mod, ".rds", sep = "")
save(b,file=FileName)

