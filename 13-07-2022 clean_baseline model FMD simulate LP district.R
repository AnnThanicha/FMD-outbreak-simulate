# FMD outbreak simulation for LP district ----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------
## install library and files  -----------------
packages <- c("dplyr", "data.table","progress","profvis", "ggplot2", "reshape2")
lapply(packages, library, character.only = TRUE)

# import required files 
# this is a file for between-farm distance matrix
Mod_distance <- readRDS("ModLP_distance.rds")
# this is a file for trade network matrix
component_matrix <- readRDS("LP_component_matrix.rds")
# this is a file that store all the farm data such as farm type, number of animals
dat_kernel <- readRDS("datkernel_LP.rds")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------

## Set the baseline parameters for simulation -----------------------

### Baseline parameters setting -------
baseline_cond <- function(){
  
  
  # set index cases
  index_farm <- c(409,422,84) # this is index cases from kernel_LP data
  
  # transmission parameters
  k0 =   0.00537620
  r0 =  0.17187575
  alpha = 1.49815422
  delta =  0.00061939
  
  # set rate and shape for infectious duration
  shape_inf =3.0262817
  rate_inf =0.1377417
  
  #latent period
  latent = 3
  
  # calculate the kernel transmission matrix
  Mod_kernel <- (k0/(1+((Mod_distance/r0)^alpha))) 
  diag(Mod_kernel) <- 0 # diagonal = 0
  
  
  # calculate the trade network transmission matrix  
  Mod_component_matrix_delta <- component_matrix * delta
  diag(Mod_component_matrix_delta) <- 0 # diagonal = 0
  
  
  # detection time
  detection_time_index = 7
  detectiontime_cattle = 5 # for dairy cattle
  detectiontime_beef = 7 # detection time after day since infection
  detectiontime_goat = 18
  detectiontime_pig = 8
  
  # set detection time in dat_kernel dataframe
  dat_kernel$detection_time <-NA
  dat_kernel$detection_time[dat_kernel$type3=="index"]<- detection_time_index
  dat_kernel$detection_time[dat_kernel$type3=="goat"]<-  detectiontime_goat
  dat_kernel$detection_time[dat_kernel$type3=="cattle"]<-  detectiontime_cattle
  dat_kernel$detection_time[dat_kernel$type3=="pig"]<-  detectiontime_pig
 
  list( index_farm = index_farm , k0= k0, r0 = r0, alpha = alpha, delta = delta, 
        shape_inf = shape_inf, rate_inf = rate_inf, Mod_kernel = Mod_kernel,
        Mod_component_matrix_delta = Mod_component_matrix_delta, latent = latent,
        detection_time_index =  detection_time_index,
        detectiontime_cattle = detectiontime_cattle,
        detectiontime_beef = detectiontime_beef,
        detectiontime_goat = detectiontime_goat,
        detectiontime_pig = detectiontime_pig,
        dat_kernel=dat_kernel
   
       )
}

# save all the parameters as a list 
baseline_cond_list <- baseline_cond()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------

## Function to summarize number of affected farms and duration of outbreak for each iteration ---------
sum_outbreak <- function (dat) {
  
  outbreak_duration <- sapply(dat , ncol)-6
  n_susceptible <- sapply(dat , function(x) sum(x[,ncol(x)] == "S" ))
  n_latent <- sapply(dat , function(x) sum(x[,ncol(x)] == "L" ))
  n_infectious <- sapply(dat , function(x) sum(x[,ncol(x)] == "I" ))
  n_recover <- sapply(dat , function(x) sum(x[,ncol(x)] == "R" ))
  n_VI <- sapply(dat , function(x) sum(x[,ncol(x)] == "VI" ))
  n_VS <- sapply(dat , function(x) sum(x[,ncol(x)] == "VS" ))
  EMvaccinated <- sapply(dat , function(x) sum(x$day_sinceEV !=0 ))
  Quarantine <- sapply(dat , function(x) sum(x$day_sinceQ >0 ))
  n_culling <- sapply(dat , function(x) sum(x[,ncol(x)] == "C" ))
  n_UI <- sapply(dat , function(x) sum(x[,ncol(x)] == "UI" ))
  
  outbreak <- data.frame (outbreak_duration, n_susceptible, n_latent, n_infectious,
                          n_recover,  n_VI, n_VS,EMvaccinated,n_culling,Quarantine,n_UI)
  return(outbreak)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------

## Function to simulate the outbreak for baseline scenario --------------
simulate_outbreak<- function (iter) {
  # add profvis profiling
  profvis({
    # create list to keep results
    Mod_outbreaksim <- list()
    set.seed(1234)
    # for loop for iterated simulation
    for (n in 1:iter) {
      
      # create dataframe to store the status of farm in each day
      # 
      # Each row is the farm 
      # column 1 is ID, 2 is type of farm, 3 is day since infection, 
      # 4 is detection time, 5 is day since emergency vaccination, 6 is day since quarantine, 
      # and 7 is status of farm on day 1
      # in each forloop new column will be added to this dataframe to update the status of farm on that day
      # the possible status of farm S = susceptible
      # L = latent, UI = undetected infectious, I = detected infectious, R = recovered, Q = quarantine, 
      # C = culled, VI = vaccinated infectious, VS = vaccinated susceptible
      Mod_status_df <- data.frame(ID = c(1:nrow(dat_kernel)), type = dat_kernel$type3, day_sinceinf = 0,detection = dat_kernel$detection_time,
                                  day_sinceEV = 0,day_sinceQ = 0, status_D1 = "S")
      
      
      # randomly generate infectious duration
      inf_duration <- round(rgamma(nrow(dat_kernel),shape= shape_inf, rate = rate_inf))
      # plus latent period
      inf_duration <- inf_duration+latent 
   
      # set index cases
      index_farm =  index_farm
      # status of index cases on day 1 is L
      Mod_status_df$status_D1 [index_farm] <- "L"
      # day since infection of index cases is equal 1 on day 1
      Mod_status_df$day_sinceinf [index_farm] <- 1
      
      # sum number of infectious farms at the beginning
       n_infectious = sum(Mod_status_df$status_D1 %in% c("L","UI", "I", "Q","VI"))
      
      # simulated until n_infectious = 0 
      while(n_infectious > 0){
        
        # Check the infectious and susceptible status of farms on the last day
        status_infectious <- as.matrix(as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("I","UI")))
        status_susceptible <- t(as.matrix(as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("S"))))
       
        # calculate the probability of infection from transmission kernel on susceptible farm
        # by multiplying farm status matrix and transmission kernel matrix
        # row is ID of infectious farms, column is ID of susceptible farms
        dist_trans <- colSums(Mod_kernel*(status_infectious%*% status_susceptible)) 
        
        # calculate the probability of infection from trade network on susceptible farm
        # by multiplying farm status matrix and trade network transmission matrix
        trade_trans <- colSums(Mod_component_matrix_delta*(status_infectious%*% status_susceptible)) 
        
      
        # calculate probability of infection
        p_inf <- 1 - exp(-1*(dist_trans+trade_trans))
        # random infection 0 and 1 with prob. of infection
        inf <- rbinom(length(p_inf), size = 1, prob=p_inf)
        
        # create new vector to store update farm status for next day
        Mod_status_df$status <- Mod_status_df[,ncol(Mod_status_df)]  
        
        # if farm got infection. The status change to latent
        Mod_status_df$status[which(inf>0)] = "L"
        
        ## Update status of farms on the current day
        
        # +1 day since infection for if the farms are not susceptible
        Mod_status_df$day_sinceinf <- if_else( Mod_status_df$status %in% c("L","UI", "I", "Q", "R","VI"), Mod_status_df$day_sinceinf+1, Mod_status_df$day_sinceinf)
        # day since EV +1 in EM vaccinated farm
        Mod_status_df$day_sinceEV <-  if_else( Mod_status_df$day_sinceEV > 0, Mod_status_df$day_sinceEV+1, Mod_status_df$day_sinceEV)
        
        # if day_sinceinf > latent, farm status changes from latent to undetected infectious
        Mod_status_df$status[Mod_status_df$day_sinceinf > latent & Mod_status_df$day_sinceinf < Mod_status_df$detection  & Mod_status_df$status =="L"] <- "UI"
        
        # if day_sinceinf == detection time, farm status changes from undetected infectious to detected infectious
        Mod_status_df$status[Mod_status_df$day_sinceinf== Mod_status_df$detection & Mod_status_df$day_sinceinf <= inf_duration & Mod_status_df$status =="UI"] <- "I"
        
        # if day_sinceinf > inf_duration, farm status changes to recovered
        # All farm can be infected even before the recovered
        Mod_status_df$status [Mod_status_df$day_sinceinf > inf_duration  & Mod_status_df$status %in% c("I", "Q", "UI", "L","VI")] <-  "R"
        
        # after update all status from the last day, we will name the new column
        names(Mod_status_df)[ncol(Mod_status_df)] <- paste("status_D",ncol(Mod_status_df)-6,sep = "")
        
        # check if there are still infectious farms. If not, this for-loop will end and start new iterations
        n_infectious = sum(Mod_status_df[, ncol(Mod_status_df)] %in% c("L","UI", "I", "Q", "VI"))
        
      }
      # summary the result from this iteration and save the summary
      Mod_outbreaksim [[n]] <- sum_outbreak(dat = list(Mod_status_df))
      
    }
  
  return ( Mod_outbreaksim)
  })
  }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------
## Run simulation -----------
# set the condition parameters before running simulation
list2env(baseline_cond_list, envir = .GlobalEnv)

# run the simulation 
start_time <- Sys.time()
simulate_outbreak(iter=10)
end_time <- Sys.time()
end_time - start_time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------
## Profiling --------------
start_time <- Sys.time()
iter = 10
list2env(baseline_cond_list, envir = .GlobalEnv)
# add profvis for profiling
profvis({
  # create list to keep results
  Mod_outbreaksim <- list()
  set.seed(1234)
  # for loop for iterated simulation
  for (n in 1:iter) {
    
    # create dataframe to store the status of farm in each day
    # 
    # Each row is the farm 
    # column 1 is ID, 2 is type of farm, 3 is day since infection, 
    # 4 is detection time, 5 is day since emergency vaccination, 6 is day since quarantine, 
    # and 7 is status of farm on day 1
    # in each forloop new column will be added to this dataframe to update the status of farm on that day
    # the possible status of farm S = susceptible
    # L = latent, UI = undetected infectious, I = detected infectious, R = recovered, Q = quarantine, 
    # C = culled, VI = vaccinated infectious, VS = vaccinated susceptible
    Mod_status_df <- data.frame(ID = c(1:nrow(dat_kernel)), type = dat_kernel$type3, day_sinceinf = 0,detection = dat_kernel$detection_time,
                                day_sinceEV = 0,day_sinceQ = 0, status_D1 = "S")
    
    
    # randomly generate infectious duration
    inf_duration <- round(rgamma(nrow(dat_kernel),shape= shape_inf, rate = rate_inf))
    # plus latent period
    inf_duration <- inf_duration+latent 
    
    # set index cases
    index_farm =  index_farm
    # status of index cases on day 1 is L
    Mod_status_df$status_D1 [index_farm] <- "L"
    # day since infection of index cases is equal 1 on day 1
    Mod_status_df$day_sinceinf [index_farm] <- 1
    
    # sum number of infectious farms at the beginning
    n_infectious = sum(Mod_status_df$status_D1 %in% c("L","UI", "I", "Q","VI"))
    
    # simulated until n_infectious = 0 
    while(n_infectious > 0){
      
      # Check the infectious and susceptible status of farms on the last day
      status_infectious <- as.matrix(as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("I","UI")))
      status_susceptible <- t(as.matrix(as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("S"))))
      
      # calculate the probability of infection from transmission kernel on susceptible farm
      # by multiplying farm status matrix and transmission kernel matrix
      # row is ID of infectious farms, column is ID of susceptible farms
      
      dist_trans <- Mod_kernel*(status_infectious%*% status_susceptible)
      
      # calculate the probability of infection from trade network on susceptible farm
      # by multiplying farm status matrix and trade network transmission matrix
      trade_trans <- Mod_component_matrix_delta*(status_infectious%*% status_susceptible) 
      
      
      # calculate probability of infection
      p_inf <- 1 - exp(-1*(colSums(dist_trans+trade_trans)))
      # random infection 0 and 1 with prob. of infection
      inf <- rbinom(length(p_inf), size = 1, prob=p_inf)
      
      # create new vector to store update farm status for next day
      Mod_status_df$status <- Mod_status_df[,ncol(Mod_status_df)]  
      
      # if farm got infection. The status change to latent
      Mod_status_df$status[which(inf>0)] = "L"
      
      ## Update status of farms on the current day
      
      # +1 day since infection for if the farms are not susceptible
      Mod_status_df$day_sinceinf <- if_else( Mod_status_df$status %in% c("L","UI", "I", "Q", "R","VI"), Mod_status_df$day_sinceinf+1, Mod_status_df$day_sinceinf)
      # day since EV +1 in EM vaccinated farm
      Mod_status_df$day_sinceEV <-  if_else( Mod_status_df$day_sinceEV > 0, Mod_status_df$day_sinceEV+1, Mod_status_df$day_sinceEV)
      
      # if day_sinceinf > latent, farm status changes from latent to undetected infectious
      Mod_status_df$status[Mod_status_df$day_sinceinf > latent & Mod_status_df$day_sinceinf < Mod_status_df$detection  & Mod_status_df$status =="L"] <- "UI"
      
      # if day_sinceinf == detection time, farm status changes from undetected infectious to detected infectious
      Mod_status_df$status[Mod_status_df$day_sinceinf== Mod_status_df$detection & Mod_status_df$day_sinceinf <= inf_duration & Mod_status_df$status =="UI"] <- "I"
      
      # if day_sinceinf > inf_duration, farm status changes to recovered
      # All farm can be infected even before the recovered
      Mod_status_df$status [Mod_status_df$day_sinceinf > inf_duration  & Mod_status_df$status %in% c("I", "Q", "UI", "L","VI")] <-  "R"
      
      # after update all status from the last day, we will name the new column
      names(Mod_status_df)[ncol(Mod_status_df)] <- paste("status_D",ncol(Mod_status_df)-6,sep = "")
      
      # check if there are still infectious farms. If not, this for-loop will end and start new iterations
      n_infectious = sum(Mod_status_df[, ncol(Mod_status_df)] %in% c("L","UI", "I", "Q", "VI"))
      
    }
    # summary the result from this iteration and save the summary
    Mod_outbreaksim [[n]] <- sum_outbreak(dat = list(Mod_status_df))
    
  }
  
})
end_time <- Sys.time()
end_time - start_time

#' From profiling. the most time-consuming is from matrix calculation and colsums
#' Do you have any suggestion to improve that code?