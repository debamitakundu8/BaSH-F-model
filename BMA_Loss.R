rm(list=ls())

## nSim   : Total number of simulation
## nModel : Total number of Factor Model
## This function estimates BMA loss and also LPML,DIC and Oracle Loss


Loss_estimates <- function(nSim,nModel){
  
  K_choice <- 3:8              ## Factor choices
  k_true <- K_choice[3]        ## need to define the true factor model
  
DIC_allModel <- array(NA,dim=c(nSim,nModel))
LPML_allModel <- array(NA,dim=c(nSim,nModel))
DIC_modelchoice <- array(NA,nSim)
LPML_modelchoice <- array(NA,nSim)
Loss_DIC <- array(NA,nSim)
Loss_LPML <- array(NA,nSim)
Loss_allmodel <- array(NA,dim=c(nSim,nModel))


for(model in 1:nModel)
{
  K_model <- K_choice[model]
   for( s in 1:nSim)
  {
    FileName <- paste("BaSF_Factor",K_model,"_output_File",s, ".rds", sep = "")   ## Output file 
    if (file.exists(FileName))
    {
      b=readRDS(FileName)
      Loss_allmodel[s,model] <- b$Loss_all
      DIC_allModel[s,model] <- b$DIC
      LPML_allModel[s,model] <- b$LPML
    }
    else(print(s))
  }
}

Best_DIC <- apply(DIC_allModel,1,which.min)
Best_LPML <- apply(LPML_allModel,1,which.max)

DIC_table <- table(Best_DIC)
LPML_table <- table(Best_LPML)

for(s in 1:nSim){Loss_DIC[s] <- Loss_allmodel[s,][Best_DIC[s]]}
for(s in 1:nSim){Loss_LPML[s] <- Loss_allmodel[s,][Best_LPML[s]]}

DIC_loss <- mean(Loss_DIC)
LPML_loss <- mean(Loss_LPML)
Oracle_Loss <- mean(Loss_allmodel[,k_true])
####  Bayesian model averaging estimate

model_prob <- array(NA,dim=c(nSim,nModel))
for(i in 1:nSim){
  model_prob[i,]<- exp((LPML_allModel[i,]-max(LPML_allModel[i,])))/sum(exp((LPML_allModel[i,]-max(LPML_allModel[i,]))))
}
prob_mod <-  round(apply(model_prob,2,mean),4)

M <- 3; p<- 12     ## Define number of groups and covaraite
OMEGA_bayes <- array(NA,dim=c(p,p,M,nSim))
for( s in 1:nSim)
{
  Mod <- s
  Omega_bayes_temp <- array(NA,dim=c(p,p,M,nModel))
  for(model in 1:nModel)
  {
    K_model <- K_choice[model]
    FileName <- paste("BaSF_Factor",K_model,"_new_swap_runfile",s, ".rds", sep = "")  
    if (file.exists(FileName))
    {
      b=readRDS(FileName)
      for(m in 1:M) {Omega_bayes_temp[,,m,model] <- model_prob[s,model]*b$Omegaest[,,m]}
    }
    else(print(s))
  }
  OMEGA_bayes[,,,s] <- apply(Omega_bayes_temp,c(1,2,3),sum)
}

## BMA loss
Mod <- 1 ## Choice of data
dname <- paste("Y_",Mod,".RData",sep="")
load(dname)
Y=Y
Omega_true <- Omega_true
n <- n 
N <- sum(n)  
Loss_group_bayes <- array(NA,dim=c(nSim,M))
for(s in 1:nSim)
{
  for(m in 1:M)
  {
    Loss_group_bayes[s,m] <- n[m]*(tr(solve(Omega_true[,,m])%*%OMEGA_bayes[,,m,s])-log(det(solve(Omega_true[,,m])%*%OMEGA_bayes[,,m,s]))-p)    
  }
}   
Loss_Bayes <- mean(apply(Loss_group_bayes,1,sum,na.rm=T)/N)

out= list(DIC_loss=DIC_loss,LPML_loss=LPML_loss,BMA_loss=Loss_Bayes,Oracle_Loss=Oracle_Loss,
          BMA_prob=prob_mod,DIC_table=DIC_table,LPML_table=LPML_table)
return(out)
}


