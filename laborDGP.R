library(ggplot2)

#we normalize wages and prices equal to 1 in the Chetty framework, 
#such that we work with income z = hours worked h
eps <- 0.5 #elasticity
N <- 10000 #individuals
reps <- 100 #number of monte carlo simulations
I <- 1000 #number of bootstrap draws
t0 <- 0 #initial tax
t1 <- 0.37 #tax bracket 1
K <- 20 #after 30 hours worked, tax increases to t1
D <- 4 #delta, window width for bunching

set.seed(1)

genData <- function(N, eps, t0, t1, K){
  #genderate a dataset with observed hours worked (income) under t0 and t1
  
  df <- data.frame("n" = rnorm(N, 40, 30)) #ability level that implies heterogeneity in preferences
  df$z0 <- ifelse(df$n > 0, df$n*(1-t0)^eps, 0) #hours worked under no tax system
  
  #if there is a tax system, then there are three options for z1 depending on ability n
  df$z1 <- NA
  for(i in 1:N){
    
    #option 1: if n <0, then better not work at all
    if(df$n[i] <0 ){
      df$z1[i] <- 0
    }
   
    #option 2: n < K/(1-t-0)^eps, in which case z0 = z1
    else if(df$n[i] > 0 & df$n[i] < K/((1 - t0)^eps)){ 
      df$z1[i] <- rnorm(1, df$n[i]*(1-t0)^eps, 1)
      
    #option 3: K/(1-t0)^eps < n < K/(1-t1)^eps, in which case z1 = K + u
    } else if (df$n[i] < K/((1-t1)^eps) & df$n[i] > K/((1-t0)^eps)){
      df$z1[i] <- rnorm(1, K, 1) 
    }
    
    #option 4:  n > K/(1-t1)^eps, in which case z1 =  n(1-t1)^eps
    else {
      df$z1[i] <- rnorm(1, df$n[i]*(1-t1)^eps, 1)
    }
  }
  return(df)
}

computeE <- function(df, K, D){
  #given a dataframe with hours worked and the kink value, compute the elasticity
  #D is the window width denoted delta in Saez
  
  #calculate fractions in each window, see page 188
  df$b1 <- ifelse(K - 2*D < df$z1 & df$z1 < K - D, 1, 0)
  df$b2 <- ifelse(K - D < df$z1 & df$z1 < K + D, 1, 0)
  df$b3 <- ifelse(K + D < df$z1 & df$z1 < K + 2*D, 1, 0)
  
  H1 <- sum(df$b1)/N
  H2 <- sum(df$b2)/N
  H3 <- sum(df$b3)/N
  
  #calculate h_, h+ and h* and B
  h1 <- H1/D
  h2 <- H2/D
  h3 <- H3/D
  B <- H2 - (H1 + H3)
  
  #we use the ABC formula to solve for epsilon in eq (5) in Saez
  a <- K*h1/2
  b <- K/2*(h3 - h1) - B
  c <- - K*h3/2
  
  X <- (-b + sqrt(b^2 - 4*a*c))/(2*a) #choose positive to ensure existence
  Y <- (1 - t0)/(1 - t1)
  eHat <- log(X)/log(Y)
  return(eHat)
}

#generate a dataset with hours worked
df <- genData(N, eps, t0, t1, K)

#for ease of plots transform data to grouped
dfGrouped <- data.frame("z" = c(df$z0, df$z1))
dfGrouped$tax <- rep(c("t0", "t1"), each = N)
dfGrouped$tax <- as.factor(dfGrouped$tax )

#plot the histograms of observed hours worked
ggplot(dfGrouped, aes(x = z, fill = tax)) +
  geom_histogram(alpha = 0.7, aes(y = ..density..), position = 'identity')

#to see whether the methodology works, do a monte carlo simulation with R draws
estimates <- data.frame("e" = rep(NA, reps))

for(r in 1:reps){
  df <- genData(N, eps, t0, t1, K)
  estimates$e[r] <- computeE(df, K, D)
}

summary(estimates$e) #indeed, the process works!

#we can use the bootstrap to calculate confidence intervals: 
#draw from the original sample N observations with replacement
#calculate ehat and store the values
#repeat this process I times

bootstrapE <- data.frame("e" = rep(NA, I))
for(i in 1:I){
  simDf <- df[sample(1:N, N, replace = TRUE), ]
  bootstrapE$e[i] <- computeE(simDf, K, D)
}
quantile(bootstrapE$e, probs= c(0.025, 0.975)) 
#note that the bootstrap does not remove any sampling bias, but rather estimation errors


