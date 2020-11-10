##############################################################################################
##### within herd model  #####################################
###### SEIR model without vaccination #######################################
# susceptible-beta - > exposed - latent period -> infection -infectious period-> recovered
seir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  E = state_values [2]        # exposed
  I = state_values [3]        # infectious
  R = state_values [4]        # recovered
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      dS = (-beta * S * I)/N
      dE = (beta * S * I)/N - (delta * E)
      dI = (delta * E) - (gamma * I)
      dR = (gamma * I)
      
      # combine results
      results = c (dS, dE, dI, dR)
      list (results)
    }
  )
}

#Parameters

contact_rate = 10                     # number of contacts per day
transmission_probability = 0.07       # transmission probability
infectious_period = 14                 # infectious period
latent_period = 2                     # latent period

#Compute values of beta (tranmission rate) and gamma (recovery rate).

beta_value = 1.5 #contact_rate * transmission_probability
gamma_value = 1 / infectious_period
delta_value = 1 / latent_period

#Compute Ro - Reproductive number.

Ro = beta_value / gamma_value

#Disease dynamics parameters.

parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value)

#Initial values for sub-populations.

W = 29       # susceptible hosts
X = 0           # exposed hosts
Y = 1           # infectious hosts
Z = 0           # recovered hosts

#Compute total population.

N = W + X + Y + Z

#Initial state values for the differential equations.

initial_values = c (S = W, E = X, I = Y, R = Z)

#Output timepoints.

timepoints = seq (0, 50, by=1)

#Simulate the SEIR epidemic.

output = ode(initial_values, timepoints, seir_model, parameter_list)

##############################################################################################
########## prepare individual data ##########################################
t<-20
N0<-500

farm <- vector(mode="list", N0)
for(i in seq(farm)){
  farm[[i]]$ID<-FarmID$ID[i]
  farm[[i]]$pop <- FarmID$pop[i]
  farm[[i]]$infduration <- 0
  farm[[i]]$status <- FarmID$status[i]
  farm[[i]]$pop_sus<- FarmID$pop_sus[i]
  farm[[i]]$pop_exposed<- FarmID$pop_exposed[i]
  farm[[i]]$pop_inf<- FarmID$pop_inf[i]
  farm[[i]]$pop_recover<- FarmID$pop_recover[i]
  farm[[i]]$pop_vaccine<- FarmID$pop_vaccine[i]
}

# create dataframe to store results
result <- data.frame(matrix(nrow = t, ncol = 5))
colnames(result) <- c("day", "sus", "exposed", "inf", "recover")


######################################################################
############# Link within herd and between herd model ################################
t=10
result <- data.frame(matrix(nrow = t, ncol = 5))
colnames(result) <- c("day", "sus", "exposed", "inf", "recover")
IDinfection<-list()

outbreaksim<-function (dat){
  
  for(i in seq(t)){ # loop for each time increment
    
    is.inf <- which(sapply(dat, function(x) x$pop_inf) >0) #infection farm 
    ID.inf <-sapply(dat[c(is.inf)], function(x) x$ID ) # infection ID
    is.sus <- which(sapply(dat, function(x) x$pop_sus == x$pop)) #susceptible farm 
    ID.sus <-sapply(dat[c(is.sus)], function(x) x$ID ) # susceptible ID
    
    for(j in ID.sus){ #loop for each susceptible individual
      probinf <-rbinom(n=1, size =1, prob = 1-exp(-sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(ID.sus[j]) & kernelmatrixlong$inf %in% c(ID.inf) ])))  # calculate pob of infection for each sus ceptible individual
      dat[[j]]$status <-ifelse(probinf ==1, "exposed",dat[[j]]$status)} # if infection from binomial distribution farm change status
    
    dat<-lapply(dat, transform, infduration = ifelse(status!="sus", infduration+1, infduration))# advance infectious duration for infectious farm
    dat<-lapply(dat, transform, status = ifelse(infduration >2 , "inf", status))
    dat<-lapply(dat, transform, status = ifelse(infduration >14, "recover", status))
    
    
    #Population stats
    
    result$day[i] <- i
    result$sus[i] <- length(which(sapply(dat, function(x) x$status) == "sus")) #susceptible farm
    result$exposed[i]<-length(which(sapply(dat, function(x) x$status) == "exposed")) # exposed farm
    result$inf[i]<-length(which(sapply(dat, function(x) x$status) == "inf")) #infection farm
    result$recover[i] <-length(which(sapply(dat, function(x) x$status) == "recover")) #recover farm
    
    IDinfection[i]<-list(sapply(dat[which(sapply(dat, function(x) x$status) == "inf")], function(x) x$ID ))
    
  }
  
  list(result,IDinfection) #keep results in list
  
}

outbreakresult<-outbreaksim(dat=farm)
out = replicate(n = 3, expr = outbreaksim(dat=farm))# replicate MonteCarlo


# check sum kernel
dat=farm
is.inf <- which(sapply(dat, function(x) x$pop_inf) >0) #infection farm 
ID.inf <-sapply(dat[c(is.inf)], function(x) x$ID )
is.sus <- which(sapply(dat, function(x) x$pop_sus == x$pop)) #susceptible farm 
ID.sus <-sapply(dat[c(is.sus)], function(x) x$ID )


which(ID.sus==381)
for(j in ID.sus){ #loop for each susceptible individual
  probinf <-rbinom(n=1, size =1, prob = 1-exp(-sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(ID.sus[140]) & kernelmatrixlong$inf %in% c(ID.inf)])))  # calculate pob of infection for each sus ceptible individual
  dat[[j]]$status <-ifelse(probinf ==1, "exposed",dat[[j]]$status)} # if infection from binomial distribution farm change status

1-exp(-sum(kernelmatrixlong$value[which(kernelmatrixlong$sus %in% c(ID.sus[369]) & kernelmatrixlong$inf %in% c(ID.inf))]))
