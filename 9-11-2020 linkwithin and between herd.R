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
