library (dplyr)


b <- 0.14 # probability of birth
d <- 0.08 # probability of death
K <- 100 # carrying capacity
N0 <- 500 # starting number of individuals
t <- 100 # time of simulation

#create starting individual w attributes ("alive", "age", "color")
set.seed(1)
ind <- vector(mode="list", N0)
for(i in seq(ind)){
  ind[[i]]$alive <- 1
  ind[[i]]$age <- 0
  ind[[i]]$color <- c("blue", "red")[round(runif(1)+1)]
}

#make empty vectors to record population statistics
time <- seq(t+1)

pop <- NaN * time # population size
pop[1] <- N0

frac.blue <- NaN * time # fraction of population that is blue
cols <- sapply(ind, function(x) x$color)
frac.blue[1] <- sum(cols  == "blue") / length(cols)

med.age <- NaN * time
ages <- sapply(ind, function(x) x$age)
med.age[1] <- median(ages)


#simulation
save.alive.only <- TRUE # optional cropping of "ind" to include alive individuals only 
t1 <- Sys.time()
for(i in seq(t)){ # loop for each time increment
  
  is.alive <- which(sapply(ind, function(x) x$alive) == 1)
  for(j in is.alive){ #loop for each alive individual
    birth <- runif(1) <= (b * (1 - length(is.alive)/K)) # calculate a birth probability for each individual that is alive
    if(birth){
      len.ind <- length(ind)
      ind[[len.ind+1]] <- list(alive=1, age=0, color=ind[[j]]$color) # create offspring, inherits color of parent
    }
    death <- runif(1) <= d # calculate a death probability for each individual 
    if(death){
      ind[[j]]$alive <- 0 # if death, reset alive = 0
    } else { #else, advance age + 1
      ind[[j]]$age <- ind[[j]]$age + 1 # advance age of parent
    }
  }
  
  #optional cropping of list "ind"
  if(save.alive.only){
    is.dead <- which(sapply(ind, function(x) x$alive) == 0)
    if(length(is.dead) > 0) ind <- ind[-is.dead]
  }
  
  #Population stats
  is.alive <- which(sapply(ind, function(x) x$alive) == 1)
  pop[i+1] <- length(is.alive) 
  
  cols <- sapply(ind, function(x) x$color)
  frac.blue[i+1] <- sum(cols[is.alive]  == "blue") / length(is.alive)
  
  ages <- sapply(ind, function(x) x$age)
  med.age[i+1] <- median(ages[is.alive])
  
  print(paste(i, "of", t, "finished", "[", round(1/t*100), "%]"))
}
t2 <- Sys.time()
dt <- t2-t1
dt


#plot populations
png("pops_vs_time.png", width=6, height=4, units="in", res=400)
par(mar=c(4,4,1,1))
pop.blue <- pop * frac.blue
pop.red <- pop * (1-frac.blue)
ylim=range(c(pop.blue, pop.red))
plot(time, pop.blue, t="l", lwd=2, col=4, ylim=ylim, ylab="Population size")
lines(time, pop.red, lwd=2, col=2)
legend("topleft", legend=c("blue pop.", "red pop."), lwd=2, col=c(4,2), bty="n")
dev.off()

#plot median age
png("med_age_vs_time.png", width=6, height=4, units="in", res=400)
par(mar=c(4,4,1,1))
plot(time, med.age, t="l", lwd=2, ylab="Median age")
dev.off()


#####################################################################################
############# trial code for simulate outbreak ############################################################
library(readxl)
distance <- read_excel("D:/Phd thesis/R code/R-FMD simulate/lumpayaklang euclidean distance.xlsx")# distance matrix
distance<-distance[,-1]
distance<-distance/1000
t<-100 # time = 100 days
k0 <-0.0054 
r0<-0.19
alpha<-1.56
day <- vector(mode="list", t)


kernelmatrix<- as.matrix((k0/(1+((distance/r0)^alpha)))) # get kernel matrix

library(reshape2)
kernelmatrixlong<-melt(kernelmatrix) # make long kernel matrix
colnames(kernelmatrixlong)<-c("sus","inf","value")# change column name
kernelmatrixlong$value[kernelmatrixlong$sus %in% c(201) ] # select kernelvalue from ID
sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(201) ]) #sum kernel value
sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(201) & kernelmatrixlong$inf %in% c(202) ]) #sum kernel value from var1 and var2

FarmID <- read_excel("FarmID.xlsx")#import farm data

#create list with each list as individual farms
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


################# for loop simulate ######################

for(i in seq(t)){ # loop for each time increment
  
  is.inf <- which(sapply(farm, function(x) x$status) == "inf") #infection farm 
  is.sus <- which(sapply(farm, function(x) x$status) == "sus") #susceptible farm 
  
  for(j in is.sus){ #loop for each susceptible individual
    probinf <-rbinom(n=1, size =1, prob = 1-exp(-sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(is.sus[j]) & kernelmatrixlong$inf %in% c(is.inf) ])))  # calculate pob of infection for each sus ceptible individual
    farm[[j]]$status <-ifelse(probinf ==1, "exposed",farm[[j]]$status)} # if infection from binomial distribution farm change status
  
  farm<-lapply(farm, transform, infduration = ifelse(status!="sus", infduration+1, infduration))# advance infectious duration for infectious farm
  farm<-lapply(farm, transform, status = ifelse(infduration >2 , "inf", status))
  farm<-lapply(farm, transform, status = ifelse(infduration >14, "recover", status))
  
  
  #Population stats
  
  result$day[i] <- i
  result$sus[i] <- length(which(sapply(farm, function(x) x$status) == "sus")) #susceptible farm
  result$exposed[i]<-length(which(sapply(farm, function(x) x$status) == "exposed")) # exposed farm
  result$inf[i]<-length(which(sapply(farm, function(x) x$status) == "inf")) #infection farm
  result$recover[i] <-length(which(sapply(farm, function(x) x$status) == "recover")) #recover farm
}

result


### make a function to simulate ####
t<-5
N0<-500

farm <- vector(mode="list", N0)
for(i in seq(farm)){
  farm[[i]]$ID<-FarmID$ID[i]
  farm[[i]]$pop <- FarmID$pop[i]
  farm[[i]]$infduration <- 0
  farm[[i]]$status <- FarmID$status[i]
}

# create dataframe to store results
result <- data.frame(matrix(nrow = t, ncol = 5))
colnames(result) <- c("day", "sus", "exposed", "inf", "recover")
IDinfection<-list()

outbreaksim<-function (dat){
  
  for(i in seq(t)){ # loop for each time increment
    
    is.inf <- which(sapply(dat, function(x) x$pop_inf) >0) #infection farm 
    ID.inf <-sapply(dat[c(is.inf)], function(x) x$ID )
    is.sus <- which(sapply(dat, function(x) x$pop_sus == x$pop)) #susceptible farm 
    ID.sus <-sapply(dat[c(is.sus)], function(x) x$ID )
    
    for(j in is.sus){ #loop for each susceptible individual
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
    
    IDinfection[i]<-sapply(dat[which(sapply(dat, function(x) x$status) == "inf")], function(x) x$ID )
    
    }
  
  list(result) #keep results in list
  
}

outbreaksim(dat=farm)
out = replicate(n = 3, expr = outbreaksim(dat=farm))# replicate MonteCarlo

bind_rows(out, .id = "column_label")#merge list to one dataframe for plot


########################### test ##################################
t=10
result <- data.frame(matrix(nrow = t, ncol = 5))
colnames(result) <- c("day", "sus", "exposed", "inf", "recover")
IDinfection<-list()

outbreaksim<-function (dat){
  
  for(i in seq(t)){ # loop for each time increment
    
    is.inf <- which(sapply(dat, function(x) x$pop_inf) >0) #infection farm 
    ID.inf <-sapply(dat[c(is.inf)], function(x) x$ID )
    is.sus <- which(sapply(dat, function(x) x$pop_sus == x$pop)) #susceptible farm 
    #ID.sus <-sapply(dat[c(is.sus)], function(x) x$ID )
    
    for(j in is.sus){ #loop for each susceptible individual
      probinf <-rbinom(n=1, size =1, prob = 1-exp(-sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(dat[[j]]$ID) & kernelmatrixlong$inf %in% c(ID.inf) ])))  # calculate pob of infection for each sus ceptible individual
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


which(is.sus==381)
for(j in is.sus){ #loop for each susceptible individual
  probinf <-rbinom(n=1, size =1, prob = 1-exp(-sum(kernelmatrixlong$value[kernelmatrixlong$sus %in% c(which(is.sus==j)) & kernelmatrixlong$inf %in% c(ID.inf)])))  # calculate pob of infection for each sus ceptible individual
  dat[[j]]$status <-ifelse(probinf ==1, "exposed",dat[[j]]$status)} # if infection from binomial distribution farm change status

1-exp(-sum(kernelmatrixlong$value[which(kernelmatrixlong$sus %in% c(ID.sus[359]) & kernelmatrixlong$inf %in% c(ID.inf))]))
1-exp(-sum(kernelmatrixlong$value[which(kernelmatrixlong$sus %in% c(which(is.sus==369)) & kernelmatrixlong$inf %in% c(ID.inf))]))
for(j in is.inf){print (j)}

ID.inf <-sapply(dat[c(14)], function(x) x$ID )
dat[[1]]$ID