install.packages("ggplot2")
library(bbmle)
library(optimx)
library(readr)
library(readxl)
library(writexl)
library(igraph)
library(ggplot2)

#####################################################
### Import dataset #####################################

LP_distancematrix <- read_csv("usedthtis-lumpayaklangeuclideandistance.csv")
LP_distancematrix<-LP_distancematrix[,-1]# remove row name
LP_distancematrix<-LP_distancematrix/1000 #change to km

datLP <- read_excel("usedthis_kernel_lumpayaklang2.xlsx")
summary(datLP)

datBP<-read_excel("FMDboploy.xlsx")

BP_distancematrix <- read_csv("Boploy distance matrix.csv")
BP_distancematrix<-BP_distancematrix[,-1]# remove row name
BP_distancematrix<-BP_distancematrix/1000 #change to km

###############################################################
## function prepare status matrix

################################################################
#function for susceptible status matrix
status_sus_func<-function(dat){
  dat<-dat[order(dat$ID),] # sort data based on farm ID
  a <-matrix(1, nrow=nrow(dat), ncol=max(dat$stopinfectiousdate, na.rm=TRUE)) # create matrix with nrow = ID, ncol=date
  
  
  # change the susceptible status after infection date to 0, became susceptible again after waning
  for (i in 1:length(dat$ID)) {
    if(!is.na(dat$infectiondate[i]) & !is.na(dat$waningdate[i])){ a[i, dat$infectiondate[i]:dat$waningdate[i]]<-0} # susceptible =0 from infected to waning date
  if(is.na(dat$infectiondate[i]) & !is.na(dat$waningdate[i])) {a[i,1:dat$waningdate[i]]<-0} # for farm that got infected before and waning during study period
    if(!is.na(dat$infectiondate[i]) & is.na(dat$waningdate[i])) {a[i,dat$infectiondate[i]:ncol(a)]<-0} # for farm that got infection during study period and not waning during study peroid
    }
  
  suffixrow <- seq(1:nrow(a))# generate name for row and column for statussus data frame; row = FarmID; column = SusDay
  suffixcol <- seq(1:ncol(a))
  myname_row<- paste("ID", suffixrow, sep="")
  myname_col<- paste("susD", suffixcol, sep="")
  row.names(a)<-myname_row
  colnames(a)<-myname_col
  a<-data.frame(a)
  return(a)
}


#########################################################################
# function for the infected status but not yet infectious
status_infected_func<-function(dat){
  
  dat<-dat[order(dat$ID),] # sort data based on farm ID
  a <-matrix(0, nrow=nrow(dat), ncol=max(dat$stopinfectiousdate, na.rm=TRUE)) # create matrix with nrow = ID, ncol=date
  
  
  # change the status to 1 on infected day
  for (i in 1:length(dat$ID)) {
    if (!is.na(dat$infectiondate[i])){a[i,dat$infectiondate[i]] <-1}  # become infected until infectious day-1
    # if the farm did not get infected, the infection date = NA. Status infected is always remain 0
  }
  
  suffixrow <- seq(1:nrow(a))# generate name for row and column for matrix; row = FarmID; column = SusDay
  suffixcol <- seq(1:ncol(a))
  myname_row<- paste("ID", suffixrow, sep="")
  myname_col<- paste("infD", suffixcol, sep="")
  row.names(a)<-myname_row
  colnames(a)<-myname_col
  a<-data.frame(a)
  return(a)
}


###########################################################
# function for the infectious status
status_infectious_func<-function(dat){
  
  dat<-dat[order(dat$ID),] # sort data based on farm ID
  a <-matrix(0, nrow=nrow(dat), ncol=max(dat$stopinfectiousdate, na.rm=TRUE)) # create matrix with nrow = ID, ncol=date
  
  
  # change the status to 1 on infectious day
  for (i in 1:length(dat$ID)) {
    if(!is.na(dat$infectiondate[i])){a[i, dat$startinfectiousdate[i]:dat$stopinfectiousdate[i]]<-1}   
  }
  
  suffixrow <- seq(1:nrow(a))# generate name for row and column for  matrix; row = FarmID; column = infday
  suffixcol <- seq(1:ncol(a))
  myname_row<- paste("ID", suffixrow, sep="")
  myname_col<- paste("infiD", suffixcol, sep="")
  row.names(a)<-myname_row
  colnames(a)<-myname_col
  a<-data.frame(a)
  return(a)
}



####try estimate without trade ###########
kernel_estimate<- function( k0,r0,alpha){
  
  
  kernelmatrix<-(k0/(1+((distancematrix/r0)^alpha))) 
  diag( kernelmatrix) <- 0
  
  
  #create blank matrix to store lambda
  lamnda_inf<-matrix(NA , ncol(status_sus), nrow(status_sus)) 
  lamnda_esc<-matrix(NA , ncol(status_sus), nrow(status_sus))
  

  for (i in 1:ncol(status_sus)){
    #lambda for escape
    lamnda_esc[i,] = rowSums(kernelmatrix*(as.matrix(status_sus[,i])%*% as.matrix(t(status_infectious[,i])))) 
    
    
    #lambda for infected
    lamnda_inf[i,] = rowSums(kernelmatrix*(as.matrix(status_infected[,i])%*% as.matrix(t(status_infectious[,i]))))
    
  }
  
  lambdaesc<-colSums(lamnda_esc)# sum force of infection from each ID
  logPesc<-sum(lambdaesc)
  lambdinf<-colSums(lamnda_inf) # if lambda =0 , not calculate prob inf
  logPinf<-sum(log(1-(exp(-1*lambdinf[lambdinf>0]))))
  return((-1)*(logPinf-logPesc))# * -1 to minimize loglikelihood
}



####################################################################################
################### Estimate Lumpayaklang kernel without trade #####################
#prepare the data first
status_sus_LP<-status_sus_func(dat=datLP)
status_infected_LP<-status_infected_func(dat=datLP)
status_infectious_LP<-status_infectious_func(dat=datLP)
distancematrix_LP<- as.matrix(LP_distancematrix)

status_sus<-status_sus_LP
status_infected<-status_infected_LP
status_infectious<-status_infectious_LP
distancematrix<-distancematrix_LP

kernel_estimate(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289)
kernel_estimate_new(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289)

start_time <- Sys.time()#check running time
fitbaseline_LP<-mle2(kernel_estimate, start = list(k0 =  0.005 ,r0 = 0.19,alpha= 1.56), skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx", lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001))
fitbaseline_LP2<-mle2(kernel_estimate_new, start = list(k0 =  0.005 ,r0 = 0.19,alpha= 1.56), skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx", lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001))
#Call:
#  mle2(minuslogl = kernel_estimate, start = list(k0 = 0.005, r0 = 0.19, alpha = 1.56), method = "nlminb", optimizer = "optimx", 
#       skip.hessian = FALSE, lower = c(k0 = 1e-04, r0 = 1e-04, alpha = 1e-04))

#Coefficients:
#  Estimate Std. Error z value     Pr(z)    
#k0    0.0053784  0.0018946  2.8387  0.004529 ** 
#  r0    0.1878576  0.0759498  2.4734  0.013382 *  
#  alpha 1.5625178  0.1482510 10.5397 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-2 log L: 3421.981 


end_time <- Sys.time()
end_time - start_time

summary(fitbaseline_LP)
####################################################################################
################### Estimate Boploy kernel without trade #####################
#prepare the data first
status_sus<-status_sus_func(dat=datBP)
status_infected<-status_infected_func(dat=datBP)
status_infectious<-status_infectious_func(dat=datBP)
distancematrix = as.matrix(BP_distancematrix)

kernel_estimate(k0 = 0.003241679, r0 = 0.1822809, alpha = 1.003901)

start_time <- Sys.time()#check running time
fitbaseline_BP<-mle2(kernel_estimate, start = list(k0 = 0.003241679, r0 = 0.1822809, alpha = 1.003901), skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx", lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001))
end_time <- Sys.time()
end_time - start_time



####################################################################################
################### Estimate Lumpayaklang kernel with trade #####################
dc_inet2 <- components(inet2)

#get component member,select farm node in same component, exclude trader and NA used ID, keep list element as numeric

LP_df_cl<-lapply(seq_along(dc_inet2$csize)[dc_inet2$csize > 1], function(x)  as.numeric(V(inet2)$Used_ID[dc_inet2$membership %in% x & V(inet2)$status != "trader"& (V(inet2)$Used_ID)!= "NA"])) 
LP_df_cl# component from Saraburi trade with each component members in the same column

# remove component with only one member
LP_df_cl<-LP_df_cl[sapply(LP_df_cl, length)>1]

#First create blank trade matrix 
LP_component_matrix<-matrix(0, nrow=nrow(LP_distancematrix), ncol=ncol(LP_distancematrix))

# the farms from same component is 1 , while farms that are not in the same component is 0
for(i in 1:length(LP_df_cl)){
  LP_component_matrix[LP_df_cl[[i]],LP_df_cl[[i]]]<-1   }

diag(LP_component_matrix)<-0 # diagonal to 0


###########################################################################################
######## function to estimate kernel with trade
kernel_estimate_trade<- function( k0,r0,alpha,delta){
  
  
  kernelmatrix<-(k0/(1+((distancematrix/r0)^alpha))) 
  diag( kernelmatrix) <- 0
  
  component_matrix_delta<-component_matrix*delta
  diag(component_matrix_delta) <- 0
  
  
  #create blank matrix to store lambda
  lamnda_inf<-matrix(NA , ncol(status_sus), nrow(status_sus)) 
  lamnda_esc<-matrix(NA , ncol(status_sus), nrow(status_sus))
  lamnda_esc_trade<-matrix(NA , ncol(status_sus), nrow(status_sus))
  lamnda_inf_trade<-matrix(NA , ncol(status_sus), nrow(status_sus))
  
   for (i in 1:ncol(status_sus)){
    #lambda for escape
    lamnda_esc[i,] = rowSums(kernelmatrix*(as.matrix(status_sus[,i])%*% as.matrix(t(status_infectious[,i])))) 
    
    
    #lambda for infected
    lamnda_inf[i,] = rowSums(kernelmatrix*(as.matrix(status_infected[,i])%*% as.matrix(t(status_infectious[,i]))))
    
    #lambda for escape trade
   
    lamnda_esc_trade[i,] = rowSums(component_matrix_delta*(as.matrix(status_sus[,i])%*% as.matrix(t(status_infectious[,i])))) 
    
    
    #lambda for infected trade
    lamnda_inf_trade[i,] = rowSums(component_matrix_delta*(as.matrix(status_infected[,i])%*% as.matrix(t(status_infectious[,i]))))
    
  }
  
  lambdaesc<-colSums(lamnda_esc)# sum force of infection from each ID
  logPesc<-sum(lambdaesc)
  
  lambdinf<-colSums(lamnda_inf) # if lambda =0 , not calculate prob inf
  logPinf<-sum(log(1-(exp(-1*lambdinf[lambdinf>0]))))
  
  lamnda_esc_trade<-colSums(lamnda_esc_trade)# sum force of infection from each ID
  logPesc_trade<-sum(lamnda_esc_trade)
  
  lamnda_inf_trade<-colSums(lamnda_inf_trade) # if lambda =0 , not calculate prob inf
  logPinf_trade<-sum(log(1-(exp(-1*lamnda_inf_trade[lamnda_inf_trade>0]))))
  
  return((-1)*(logPinf+logPinf_trade-logPesc-logPesc_trade))# * -1 to minimize loglikelihood
}
###################################################################
###estimate kernel with trade
status_sus<-status_sus_func(dat=datLP)
status_infected<-status_infected_func(dat=datLP)
status_infectious<-status_infectious_func(dat=datLP)
distancematrix = as.matrix(LP_distancematrix)
component_matrix<-LP_component_matrix


kernel_estimate_trade(k0 = 0.005375230, r0 = 0.188186489, alpha = 1.563913534, delta = 0.003212226 )

fitbaseline_LP_trade2<-mle2(kernel_estimate_trade, start = list(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289, delta = 0.002), skip.hessian = FALSE)

start_time <- Sys.time()#check running time
fitbaseline_LP_trade<-mle2(kernel_estimate_trade, start = list(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289, delta = 0.002), skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx", lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001, delta =0.0001))
end_time <- Sys.time()
end_time - start_time
pro_fitbaseline_LP_trade<-profile(fitbaseline_LP_trade)
confint(pro_fitbaseline_LP_trade)
#2.5 %      97.5 %
#  k0    0.002842059 0.011869501
#r0    0.070430290 0.369963443
#alpha 1.293751473 1.877298796
#delta 0.001471229 0.005979808

#plot kernel
a<-function(k0,r0,alpha, distance){ k0/(1+((distance/r0)^alpha))} 
distance<-seq(0,5,0.01)
upper<-a(distance=t, k0=0.002842059 , r0= 0.070430290, alpha = 1.293751473)
middle<-a(distance=t, k0=0.005375230 , r0=0.188186489 , alpha = 1.563913534)
lower<-a(distance=t, k0=0.011869501 , r0=0.369963443 , alpha =  1.877298796)
kernelCI<-data.frame(distance,upper,middle,lower)
labelx<-seq(0,5,1)

plotk<-ggplot(data = kernelCI, aes(x = distance)) + 
  geom_line(aes(y = upper), color = "blue", linetype="dotted", size=0.5) + 
  geom_line(aes(y = middle), color="blue")+
  geom_line(aes(y = lower), color = "blue", linetype="dotted", size=0.5)
plotk  
min(kernel$middle)



###########################################################################################
######## function to estimate kernel with trade *** NEW fixed pinf
kernel_estimate_trade2 <- function( k0,r0,alpha,delta){
  
  
  kernelmatrix<-(k0/(1+((distancematrix/r0)^alpha))) 
  diag( kernelmatrix) <- 0
  
  component_matrix_delta<-component_matrix*delta
  diag(component_matrix_delta) <- 0
  
  
  #create blank matrix to store lambda
  lamda_inf<-matrix(NA , ncol(status_sus), nrow(status_sus)) 
  lamda_esc<-matrix(NA , ncol(status_sus), nrow(status_sus))
  lamda_esc_trade<-matrix(NA , ncol(status_sus), nrow(status_sus))
  lamda_inf_trade<-matrix(NA , ncol(status_sus), nrow(status_sus))
  
  for (i in 1:ncol(status_sus)){
    #lambda for escape
    lamda_esc[i,] = rowSums(kernelmatrix*(as.matrix(status_sus[,i])%*% as.matrix(t(status_infectious[,i])))) 
    
    
    #lambda for infected
    lamda_inf[i,] = rowSums(kernelmatrix*(as.matrix(status_infected[,i])%*% as.matrix(t(status_infectious[,i]))))
    
    #lambda for escape trade
    
    lamda_esc_trade[i,] = rowSums(component_matrix_delta*(as.matrix(status_sus[,i])%*% as.matrix(t(status_infectious[,i])))) 
    
    
    #lambda for infected trade
    lamda_inf_trade[i,] = rowSums(component_matrix_delta*(as.matrix(status_infected[,i])%*% as.matrix(t(status_infectious[,i]))))
    
  }
  
  lambdinf_total<-lamda_inf+lamda_inf_trade # sum lambda from distance independent and distance dependent
  lambdinf_ID<-colSums(lambdinf_total) 
  logPinf<-sum(log(1-(exp(-1*lambdinf_ID[lambdinf_ID>0]))))# if lambda =0 , not calculate prob inf
  
  
  lambdaesc_ID<-colSums(lamda_esc)# sum force of infection from each ID
  logPesc<-sum(lambdaesc_ID)
  
  
  
  lamnda_esc_trad_ID<-colSums(lamda_inf_trade)# sum force of infection from each ID
  logPesc_trade<-sum(lamnda_esc_trad_ID)
  
 
 return((-1)*(logPinf-logPesc-logPesc_trade))# * -1 to minimize loglikelihood
}
###################################################################

kernel_estimate_trade2(k0 = 0.005375230, r0 = 0.188186489, alpha = 1.563913534, delta = 0.003212226 )


status_sus<-status_sus_func(dat=datLP)
status_infected<-status_infected_func(dat=datLP)
status_infectious<-status_infectious_func(dat=datLP)
distancematrix = as.matrix(LP_distancematrix)
component_matrix<-LP_component_matrix


fitbaseline_LP_trade_new<-mle2(kernel_estimate_trade2, start = list(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289, delta = 0.002), skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx", lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001, delta =0.0001))


