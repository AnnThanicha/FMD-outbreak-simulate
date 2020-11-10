install.packages("deSolve")
library(deSolve) # using the "ode" function


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

output = data.frame(ode(initial_values, timepoints, seir_model, parameter_list))
output[1,2]

#Plot dynamics of Susceptibles sub-population.

plot (S ~ time, data = output, type='b', col = 'blue')       


######################################################################
######### link model  ################################
linkmodel<-function (dat){
  
  for(i in seq(t)){ # loop for each time increment
    
    is.inf <- which(sapply(dat, function(x) x$pop_inf) >0) #infection farm 
    ID.inf <-sapply(dat[c(is.inf)], function(x) x$ID ) # infection ID
    is.sus <- which(sapply(dat, function(x) x$pop_sus == x$pop)) #susceptible farm 
    ID.sus <-sapply(dat[c(is.sus)], function(x) x$ID ) # susceptible ID
    
    for(k in is.inf){ #loop for each infection individual
      N = dat[[k]]$pop
      timepoints = seq (0, i, by=1)
      initial_values = c (S = dat[[k]]$pop_sus, E = dat[[k]]$pop_exposed, I = dat[[k]]$pop_inf, R = dat[[k]]$pop_recover)
      within<-data.frame (ode(initial_values, timepoints, seir_model, parameter_list))
      dat[[k]]$pop_sus = within [i,2]
      dat[[k]]$pop_exposed = within [i,3]
      dat[[k]]$pop_inf = within [i,4]
      dat[[k]]$pop_recover = within [i,5]
      } 
    
    for(j in ID.sus){ #loop for each susceptible individual
      kernel<-sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(ID.sus[j]) & kernelmatrixlong$inf %in% c(ID.inf) ])
      probinf <-rbinom(n=1, size =1, prob = 1-exp(-kernel))  # calculate pob of infection for each sus ceptible individual
      dat[[j]]$pop <-ifelse(probinf ==1, "exposed",dat[[j]]$status)} # if infection from binomial distribution farm change status
    
    
  }
  
  
  
}

