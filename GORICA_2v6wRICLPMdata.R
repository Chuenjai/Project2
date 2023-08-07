# Required packages
library(MASS)   
library(tsDyn) 
library(lavaan) 
library(restriktor) 

# Install and use package from GitHub 
library(devtools)
install_github("rebeccakuiper/ICweights")
library(ICweights)

set.seed(123)
nsim <- 5  # number of simulations

# make arrays for storing the preferred values from simulations
weightsRIvarianceTest1 <- array(NA, dim=c(nsim, 3))
colnames(weightsRIvarianceTest1) <- c("H0","H1","Hu")
weightsRIvarianceTest2 <- array(NA, dim=c(nsim, 2))
colnames(weightsRIvarianceTest2) <- c("H1","complement H1")
weightsRIvarianceTest3 <- array(NA, dim=c(nsim, 3))
colnames(weightsRIvarianceTest3) <- c("H1","H2","Hu")
weightsRIvarianceTest4 <- array(NA, dim=c(nsim, 5))
colnames(weightsRIvarianceTest4) <- c("H0","H1", "H2", "H3", "H4")
weightsRIvarianceTest5 <- array(NA, dim=c(nsim, 5))
colnames(weightsRIvarianceTest5) <- c("H0","H1", "H3a", "H4a", "Hu")
weightsRIvarianceTest6 <- array(NA, dim=c(nsim, 4))
colnames(weightsRIvarianceTest6) <- c("H1", "H2", "H3b", "H4b")

PT_weights_var1 <- array(NA, dim=c(nsim, 3))
colnames(PT_weights_var1) <-  c("H0","H1", "Hu")
PT_weights_var2 <- array(NA, dim=c(nsim, 2))
colnames(PT_weights_var2) <-  c("H1","complement H1")
PT_weights_var3 <- array(NA, dim=c(nsim, 3))
colnames(PT_weights_var3) <-  c("H1","H2","Hu")
PT_weights_var4 <- array(NA, dim=c(nsim, 5))
colnames(PT_weights_var4) <-  c("H0","H1", "H2", "H3", "H4")
PT_weights_var5 <- array(NA, dim=c(nsim, 5))
colnames(PT_weights_var5) <-  c("H0","H1", "H3a", "H4a", "Hu")
PT_weights_var6 <- array(NA, dim=c(nsim, 4))
colnames(PT_weights_var6) <-  c("H1", "H2", "H3b", "H4b")

# specify value for generate dataset
p <- 500  # number of persons and thus n=umber of rows in data
q <- 2   # number of variables
M <- 6   # number of measurement occasions
n <- q*M # number of columns in data
N<-100 # Number of iterations from VAR(1) model per person - and use only last M of those N iterations
B2<-matrix(c(0.22, -0.07, -0.07, 0.28),byrow = T, 2) # alpha, beta, gamma, delta from paper

#Omega_pop <- matrix(c(0, 0, 0, 1), byrow = T, ncol = q)
#Omega_pop <- matrix(c(0, 0, 0, 2), byrow = T, ncol = q)
#Omega_pop <- matrix(c(1, 0.5, 0.5, 2), byrow = T, ncol = q)
#Omega_pop <- matrix(c(1, -0.62, -0.62, 1), byrow = T, ncol = q) ###
#Omega_pop <- matrix(c(2, 0, 0, 2), byrow = T, ncol = q)
Omega_pop <- matrix(c(2, 1, 1, 3), byrow = T, ncol = q)

varxy <- array(NA, c(p, n)) # create array with nrow=p and ncol=n in order to keep simulate data
colnames(varxy)<-c("x1", "x2", "x3","x4", "x5", "x6", "y1", "y2", "y3","y4", "y5", "y6")
data <- array(NA, c(p, n))
colnames(data)<-c("x1", "x2", "x3","x4","x5", "x6", "y1", "y2", "y3","y4","y5", "y6")

# hypotheses with boundaries (please specify the boundaries for the variances of each random intercept)
H0 <- "kappaV == 0 ; omegaV == 0"
H1 <- "kappaV > 0.02 ; omegaV > 0.02" #the specified RI-CLPM
H2 <- "kappaV < 0 ; omegaV < 0" #which included equal to 0 and thus represents the CLPM
H3 <- "kappaV > 0.02 ; omegaV < 0"
H3a <- "kappaV > 0.02 ; omegaV == 0"
H3b <- "kappaV > 0.02 ; omegaV < 0.02"
H4 <- "kappaV < 0 ; omegaV > 0.02"
H4a <- "kappaV == 0 ; omegaV > 0.02"
H4b <- "kappaV < 0.02 ; omegaV > 0.02" 


nrposdef <- 0
round <-1 # set this to count number of iteration
WarningTeller <- 0
countNaN <- 0
countNotPosDef <-0
tellerOK <- 0
simteller<- 1

while(simteller <= nsim){
#for (simteller in 1:nsim) {
  
  #ignore warning
  oldw <- getOption("warn")
  options(warn = -1)
  
  print(paste0("Iteration", simteller))

  for(j in 1:p){
    var1 <- VAR.sim(B=B2, n=N, include="none")
    varxy[j,]<-var1[(N-M+1):N,]
  }

  xy <- mvrnorm(n = p, mu = rep(0,2), Sigma = Omega_pop, empirical=FALSE)
  data[, 1:M] <- varxy[, 1:M] + xy[,1]
  data[, (M+1):(2*M)] <- varxy[, (M+1):(2*M)] + xy[,2] #RICLPMdata, we use both models on the same data set
  
  
  # Fitting CLPM
  clpmModel<- #<- ' for specify clpmModel and eding with '
    '
  kappa =~ 1*x1 + 1*x2 + 1*x3+1*x4 + 1*x5 + 1*x6
  omega =~ 1*y1 + 1*y2 + 1*y3+1*y4 + 1*y5 + 1*y6

  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  x4 ~ mu4*1
  x5 ~ mu5*1
  x6 ~ mu6*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  y4 ~ pi4*1
  y5 ~ pi5*1
  y6 ~ pi6*1
  

  # RICLPM become CLPM because of this part
  kappa ~~ 0*kappa #variance
  omega ~~ 0*omega #variance
  kappa ~~ 0*omega #covariance

  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  p4 =~ 1*x4
  p5 =~ 1*x5
  p6 =~ 1*x6
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  q4 =~ 1*y4
  q5 =~ 1*y5
  q6 =~ 1*y6

  #constrain autoregression and cross-lagged effects to be the same across both lags.
  p6 ~ alpha*p5 + beta*q5
  p5 ~ alpha*p4 + beta*q4
  p4 ~ alpha*p3 + beta*q3
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1

  q6 ~ delta*q5 + gamma*p5
  q5 ~ delta*q4 + gamma*p4
  q4 ~ delta*q3 + gamma*p3
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1


  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  p4 ~~ u*p4
  p5 ~~ u*p5
  p6 ~~ u*p6
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  q4 ~~ v*q4
  q5 ~~ v*q5
  q6 ~~ v*q6

  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance
  p4 ~~ uv*q4
  p5 ~~ uv*q5
  p6 ~~ uv*q6'
  
  clpmConstrainedsim <- lavaan(clpmModel, data = data,
                               missing = 'ML', 
                               int.ov.free = F,
                               int.lv.free = F,
                               auto.fix.first = F,
                               auto.fix.single = F,
                               auto.cov.lv.x = F,
                               auto.cov.y = F,
                               auto.var = F)
  summary(clpmConstrainedsim, standardized = T)
  
  
  # Fitting RI-CLPM
  riclpmModel<-
    '
  kappa =~ 1*x1 + 1*x2 + 1*x3+1*x4 + 1*x5 + 1*x6
  omega =~ 1*y1 + 1*y2 + 1*y3+1*y4 + 1*y5 + 1*y6

  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  x4 ~ mu4*1
  x5 ~ mu5*1
  x6 ~ mu6*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  y4 ~ pi4*1
  y5 ~ pi5*1
  y6 ~ pi6*1

  # Random intercept part
  kappa ~~ kappa #variance
  omega ~~ omega #variance
  kappa ~~ omega #covariance

  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  p4 =~ 1*x4
  p5 =~ 1*x5
  p6 =~ 1*x6
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  q4 =~ 1*y4
  q5 =~ 1*y5
  q6 =~ 1*y6

  #constrain autoregression and cross-lagged effects to be the same across both lags.
  p6 ~ alpha*p5 + beta*q5
  p5 ~ alpha*p4 + beta*q4
  p4 ~ alpha*p3 + beta*q3
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1

  q6 ~ delta*q5 + gamma*p5
  q5 ~ delta*q4 + gamma*p4
  q4 ~ delta*q3 + gamma*p3
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1


  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  p4 ~~ u*p4
  p5 ~~ u*p5
  p6 ~~ u*p6
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  q4 ~~ v*q4
  q5 ~~ v*q5
  q6 ~~ v*q6

  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance
  p4 ~~ uv*q4
  p5 ~~ uv*q5
  p6 ~~ uv*q6'

  riclpmConstrainedsim <- lavaan(riclpmModel, data = data,
                                 missing = 'ML',
                                 int.ov.free = F,
                                 int.lv.free = F,
                                 auto.fix.first = F,
                                 auto.fix.single = F,
                                 auto.cov.lv.x = F,
                                 auto.cov.y = F,
                                 auto.var = F)
  
  warning <- any(is.na(parameterEstimates(riclpmConstrainedsim))[c(13:27,40:77),])
  notposdef <- any(eigen(lavInspect(riclpmConstrainedsim, "cov.lv"))$values <= 1e-8)
  
  if(warning == TRUE) {
    countNaN <- countNaN+1
    WarningTeller <- WarningTeller + 1
  } else if(notposdef == TRUE) {
    countNotPosDef <- countNotPosDef+1
    WarningTeller <- WarningTeller + 1
  } else { # if no warning, then do rest of code
    print(paste0("Model does not contain NAs and cov mx of latents is pos def"))
    
    #GORICA
    unstdRIvar <- coef(riclpmConstrainedsim)[c(7,8)]
    names(unstdRIvar) <- c("kappaV", "omegaV")
    unstdRIVcov <- vcov(riclpmConstrainedsim)[c(7,8), c(7,8)]
    
    goricaRIvar1 <- goric(unstdRIvar, VCOV = unstdRIVcov, hypotheses = list(H0, H1), type = "gorica")
    goricaRIvar2 <- goric(unstdRIvar, VCOV = unstdRIVcov, hypotheses = list(H1), comparison = "complement", type = "gorica")
    goricaRIvar3 <- goric(unstdRIvar, VCOV = unstdRIVcov, hypotheses = list(H1, H2), type = "gorica")
    goricaRIvar4 <- goric(unstdRIvar, VCOV = unstdRIVcov, hypotheses = list(H0, H1, H2, H3, H4), comparison = "none", type = "gorica")
    goricaRIvar5 <- goric(unstdRIvar, VCOV = unstdRIVcov, hypotheses = list(H0, H1, H3a, H4a), comparison = "unconstrained", type = "gorica")
    goricaRIvar6 <- goric(unstdRIvar, VCOV = unstdRIVcov, hypotheses = list(H1, H2, H3b, H4b), comparison = "none", type = "gorica")
    
    # Store output
    
    #GORICA weights
    weightsRIvarianceTest1[simteller,] <- goricaRIvar1$result[,7]
    weightsRIvarianceTest2[simteller,] <- goricaRIvar2$result[,7]
    weightsRIvarianceTest3[simteller,] <- goricaRIvar3$result[,7]
    weightsRIvarianceTest4[simteller,] <- goricaRIvar4$result[,7]
    weightsRIvarianceTest5[simteller,] <- goricaRIvar5$result[,7]
    weightsRIvarianceTest6[simteller,] <- goricaRIvar6$result[,7]
    
    #Penalty weights
    PT_weights_var1[simteller,] <- goricaRIvar1$result[,6]
    PT_weights_var2[simteller,] <- goricaRIvar2$result[,6]
    PT_weights_var3[simteller,] <- goricaRIvar3$result[,6]
    PT_weights_var4[simteller,] <- goricaRIvar4$result[,6]
    PT_weights_var5[simteller,] <- goricaRIvar5$result[,6]
    PT_weights_var6[simteller,] <- goricaRIvar6$result[,6]
    
    simteller<- simteller+1
  }
  #ignore warnings
  on.exit(options(warn = oldw))
  
}

