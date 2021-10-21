#######################################################################################
######## This is the demo code to show how to parameter estimate from outbreak data

# required library
library(bbmle)
library(optimx)
library(readr)
library(readxl)

####### download demo data from excel #########
datafordemo <- read_excel("C:/Users/chanc004/OneDrive - Wageningen University & Research/Phd thesis/R code/R-FMD simulate/demo/datafordemo.xlsx")

distancefordemo <- read_excel("C:/Users/chanc004/OneDrive - Wageningen University & Research/Phd thesis/R code/R-FMD simulate/demo/distancefordemo.xlsx")
distancefordemo <- distancefordemo[,-1]# remove row name

networkdemo <- read_excel("C:/Users/chanc004/OneDrive - Wageningen University & Research/Phd thesis/R code/R-FMD simulate/demo/networkdemo.xlsx")
networkdemo<- networkdemo[,-1]# remove row name
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
## function for the infectious status
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

###################################################################
## prepare data #############
demo_sus<- status_sus_func(datafordemo)
demo_infected <- status_infected_func (datafordemo)
demo_infectious <- status_infectious_func (datafordemo)

#export data 
#writexl::write_xlsx(demo_sus, "demo_sus.xlsx")
#writexl::write_xlsx(demo_infected, "demo_infected.xlsx")
#writexl::write_xlsx(demo_infectious, "demo_infectious.xlsx")

###########################################################################################
######## function to estimate kernel with trade 
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
status_sus<-demo_sus
status_infected<-demo_infected
status_infectious<-demo_infectious
distancematrix<-distancefordemo
component_matrix<-networkdemo

kernel_estimate_trade2(k0 = 0.005, r0 = 0.2, alpha = 1.563913534, delta = 0.3 )

fit <-mle2(kernel_estimate_trade2, start = list(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289, delta = 0.002),
           skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx",
           lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001, delta =0.0001))
  