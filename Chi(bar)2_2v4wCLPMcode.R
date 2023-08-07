# Required packages
library(MASS)   # to call mvrnorm
library(tsDyn) # to call VAR.sim
library(lavaan) # to call clpmModel
library(restriktor) # to call GORICA

# Install and use package from GitHub 
library(devtools)
install_github("rebeccakuiper/ChiBarSq.DiffTest")
install_github("rebeccakuiper/ICweights")
library(ChiBarSq.DiffTest)
library(ICweights)


set.seed(123)
nsim  <-  3 # number of simulations

# make arrays for storing the preferred values from simulations
# p-values from Chi-Bar-Square Difference Test
pvalueChiBar2DiffTest_clpm_riclpm <- array(NA, dim=c(nsim, 1))
pvalueChiBar2_clpm_riclpm <- array(NA, dim=c(nsim, 1)) #checking those p-values are smaller than 0.05
pvalueChiBar2DiffTest_clpm_riclpm.kappa <- array(NA, dim=c(nsim, 1))
pvalueChiBar2_clpm_riclpm.kappa <- array(NA, dim=c(nsim, 1))  #checking those p-values are smaller than 0.05
pvalueChiBar2DiffTest_clpm_riclpm.omega <- array(NA, dim=c(nsim, 1))
pvalueChiBar2_clpm_riclpm.omega <- array(NA, dim=c(nsim, 1))  #checking those p-values are smaller than 0.05
pvalueChiBar2DiffTest_riclpm_riclpm.kappa <- array(NA, dim=c(nsim, 1))
pvalueChiBar2_riclpm_riclpm.kappa <- array(NA, dim=c(nsim, 1))  #checking those p-values are smaller than 0.05
pvalueChiBar2DiffTest_riclpm_riclpm.omega <- array(NA, dim=c(nsim, 1))
pvalueChiBar2_riclpm_riclpm.omega <- array(NA, dim=c(nsim, 1))  #checking those p-values are smaller than 0.05

# p-values from Chi-Square Difference Test
pvalueChi2DiffTest_clpm_riclpm <- array(NA, dim=c(nsim, 1))
pvalueChi2_clpm_riclpm <- array(NA, dim=c(nsim, 1)) #checking those p-values are smaller than 0.05
pvalueChi2DiffTest_clpm_riclpm.kappa <- array(NA, dim=c(nsim, 1))
pvalueChi2_clpm_riclpm.kappa <- array(NA, dim=c(nsim, 1)) #checking those p-values are smaller than 0.05
pvalueChi2DiffTest_clpm_riclpm.omega <- array(NA, dim=c(nsim, 1))
pvalueChi2_clpm_riclpm.omega <- array(NA, dim=c(nsim, 1)) #checking those p-values are smaller than 0.05
pvalueChi2DiffTest_riclpm_riclpm.kappa <- array(NA, dim=c(nsim, 1))
pvalueChi2_riclpm_riclpm.kappa <- array(NA, dim=c(nsim, 1)) #checking those p-values are smaller than 0.05
pvalueChi2DiffTest_riclpm_riclpm.omega <- array(NA, dim=c(nsim, 1))
pvalueChi2_riclpm_riclpm.omega <- array(NA, dim=c(nsim, 1)) #checking those p-values are smaller than 0.05

# AIC 
AIC_clpm_riclpm <- array(NA, dim=c(nsim, 2))
colnames(AIC_clpm_riclpm)  <-  c("CLPM","RICLPM")
AICweights_clpm_riclpm  <- array(NA, dim=c(nsim, 2))
colnames(AICweights_clpm_riclpm)  <-  c("CLPM","RICLPM")

AIC_clpm_riclpm.kappa <- array(NA, dim=c(nsim, 2))
colnames(AIC_clpm_riclpm.kappa)  <-  c("CLPM","RICLPM.kappa")
AICweights_clpm_riclpm.kappa  <- array(NA, dim=c(nsim, 2))
colnames(AICweights_clpm_riclpm.kappa)  <-  c("CLPM","RICLPM.kappa")

AIC_clpm_riclpm.omega <- array(NA, dim=c(nsim, 2))
colnames(AIC_clpm_riclpm.omega)  <-  c("CLPM","RICLPM.omega")
AICweights_clpm_riclpm.omega  <- array(NA, dim=c(nsim, 2))
colnames(AICweights_clpm_riclpm.omega)  <-  c("CLPM","RICLPM.omega")

AIC_riclpm_riclpm.kappa <- array(NA, dim=c(nsim, 2))
colnames(AIC_riclpm_riclpm.kappa)  <-  c("RICLPM","RICLPM.kappa")
AICweights_riclpm_riclpm.kappa  <- array(NA, dim=c(nsim, 2))
colnames(AICweights_riclpm_riclpm.kappa)  <-  c("RICLPM","RICLPM.kappa")

AIC_riclpm_riclpm.omega <- array(NA, dim=c(nsim, 2))
colnames(AIC_riclpm_riclpm.omega)  <-  c("RICLPM","RICLPM.omega")
AICweights_riclpm_riclpm.omega  <- array(NA, dim=c(nsim, 2))
colnames(AICweights_riclpm_riclpm.omega)  <-  c("RICLPM","RICLPM.omega")

# specify value for generating data set
p <- 25  # number of persons and thus number of rows in data
q <- 2   # number of variables
M <- 4   # number of measurement occasions
n <- q*M # number of columns in data
N<-100 # Number of iterations from VAR(1) model per person - and use only last M of those N iterations
B2<-matrix(c(0.22, -0.07, -0.07, 0.28),byrow = T, 2) # alpha, beta, gamma, delta from paper
varxy <- array(NA, c(p, n)) # create array with nrow=p and ncol=n in order to keep simulate data

colnames(varxy)<-c("x1", "x2", "x3","x4", "y1", "y2", "y3","y4")
data <- array(NA, c(p, n))
colnames(data)<-c("x1", "x2", "x3","x4", "y1", "y2", "y3","y4")


nrposdef <- 0
round <-1 # set this to count number of iteration
WarningTeller <- 0
countNaN <- 0
countNotPosDef <-0
tellerOK <- 0
simteller<- 1

while(simteller<=nsim){

  #ignore warning
  oldw <- getOption("warn")
  options(warn = -1)
  
  print(paste0("Iteration", simteller))
  
  for(j in 1:p){
    var1 <- VAR.sim(B=B2, n=N, include="none")
    data[j,]<-var1[(N-M+1):N,]
  }
  
  # Fitting CLPM
  # Fitting CLPM
  clpmModel<- # <- ' specified clpmModel and ending with '
    '
  kappa =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 
  omega =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 

  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  x4 ~ mu4*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  y4 ~ pi4*1
  

  # RICLPM become CLPM because of this part
  kappa ~~ 0*kappa #variance
  omega ~~ 0*omega #variance
  kappa ~~ 0*omega #covariance

  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  p4 =~ 1*x4
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  q4 =~ 1*y4
  
  #constrain autoregression and cross-lagged effects to be the same across both lags.
  
  p4 ~ alpha*p3 + beta*q3
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1

  
  q4 ~ delta*q3 + gamma*p3
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1


  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  p4 ~~ u*p4
  
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  q4 ~~ v*q4
  

  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance
  p4 ~~ uv*q4'
  
  
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
  
  # Fitting RI-CLPM with 1 RI .kappa
  riclpmModel.kappa<- 
    '
  kappa =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 
  omega =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 

  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  x4 ~ mu4*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  y4 ~ pi4*1
  

  # Random intercept part
  kappa ~~ kappa #variance
  omega ~~ 0*omega #variance
  kappa ~~ 0*omega #covariance

  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  p4 =~ 1*x4
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  q4 =~ 1*y4
  
  #constrain autoregression and cross-lagged effects to be the same across both lags.
  
  p4 ~ alpha*p3 + beta*q3
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1

  
  q4 ~ delta*q3 + gamma*p3
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1


  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  p4 ~~ u*p4
  
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  q4 ~~ v*q4
  

  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance
  p4 ~~ uv*q4'
  
  
  riclpmModel.kappa.Constrainedsim <- lavaan(riclpmModel.kappa, data = data,
                                             missing = 'ML', 
                                             int.ov.free = F,
                                             int.lv.free = F,
                                             auto.fix.first = F,
                                             auto.fix.single = F,
                                             auto.cov.lv.x = F,
                                             auto.cov.y = F,
                                             auto.var = F)
  
  summary(riclpmModel.kappa.Constrainedsim, standardized = T)
  
  warning <- any(is.na(parameterEstimates(riclpmModel.kappa.Constrainedsim, standardized = T)[c(9:16,28:51),]))
  notposdef <- any(eigen(lavInspect(riclpmModel.kappa.Constrainedsim, "cov.lv")[c(3,5:8), c(3,5:8)])$values <= 1e-8)#*
  
  if(warning == TRUE) {
    countNaN <- countNaN+1
    WarningTeller <- WarningTeller + 1
  } else if(notposdef == TRUE) {
    countNotPosDef <- countNotPosDef+1
    WarningTeller <- WarningTeller + 1
    
  } else { # if no warning, then do rest of code
    print(paste0("Model does not contain NAs and cov mx of latents is pos def"))
    
    
    # Fitting RI-CLPM with 1 RI .omega
    riclpmModel.omega<- 
      '
  kappa =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 
  omega =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 

  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  x4 ~ mu4*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  y4 ~ pi4*1
  

  # Random intercept part
  kappa ~~ 0*kappa #variance
  omega ~~ omega #variance
  kappa ~~ 0*omega #covariance

  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  p4 =~ 1*x4
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  q4 =~ 1*y4
  
  #constrain autoregression and cross-lagged effects to be the same across both lags.
  
  p4 ~ alpha*p3 + beta*q3
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1

  
  q4 ~ delta*q3 + gamma*p3
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1


  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  p4 ~~ u*p4
  
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  q4 ~~ v*q4
  

  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance
  p4 ~~ uv*q4'
    
    riclpmModel.omega.Constrainedsim <- lavaan(riclpmModel.omega, data = data,
                                               missing = 'ML', 
                                               int.ov.free = F,
                                               int.lv.free = F,
                                               auto.fix.first = F,
                                               auto.fix.single = F,
                                               auto.cov.lv.x = F,
                                               auto.cov.y = F,
                                               auto.var = F)
    summary(riclpmModel.omega.Constrainedsim, standardized = T)
    
    warning <- any(is.na(parameterEstimates(riclpmModel.omega.Constrainedsim, standardized = T)[c(9:16,28:51),]))
    notposdef <- any(eigen(lavInspect(riclpmModel.omega.Constrainedsim, "cov.lv")[4:8, 4:8])$values <= 1e-8)#*
    
    if(warning == TRUE) {
      countNaN <- countNaN+1
      WarningTeller <- WarningTeller + 1
    } else if(notposdef == TRUE) {
      countNotPosDef <- countNotPosDef+1
      WarningTeller <- WarningTeller + 1
      
    } else { # if no warning, then do rest of code
      print(paste0("Model does not contain NAs and cov mx of latents is pos def"))
      
      
      # Fitting RI-CLPM with 2 RIs (kappa and omega)
      riclpmModel<- 
        '
   kappa =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 
   omega =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 

  x1 ~ mu1*1 #intercepts
  x2 ~ mu2*1
  x3 ~ mu3*1
  x4 ~ mu4*1
  y1 ~ pi1*1
  y2 ~ pi2*1
  y3 ~ pi3*1
  y4 ~ pi4*1
  

  # Random intercept part
  kappa ~~ kappa #variance
  omega ~~ omega #variance
  kappa ~~ omega #covariance

  #latent vars for AR and cross-lagged effects
  p1 =~ 1*x1 #each factor loading set to 1
  p2 =~ 1*x2
  p3 =~ 1*x3
  p4 =~ 1*x4
  q1 =~ 1*y1
  q2 =~ 1*y2
  q3 =~ 1*y3
  q4 =~ 1*y4
  
  #constrain autoregression and cross-lagged effects to be the same across both lags.
  
  p4 ~ alpha*p3 + beta*q3
  p3 ~ alpha*p2 + beta*q2
  p2 ~ alpha*p1 + beta*q1

  
  q4 ~ delta*q3 + gamma*p3
  q3 ~ delta*q2 + gamma*p2
  q2 ~ delta*q1 + gamma*p1


  p1 ~~ p1 #variance
  p2 ~~ u*p2
  p3 ~~ u*p3
  p4 ~~ u*p4
  
  q1 ~~ q1 #variance
  q2 ~~ v*q2
  q3 ~~ v*q3
  q4 ~~ v*q4
  

  p1 ~~ q1 #p1 and q1 covariance
  p2 ~~ uv*q2 #p2 and q2 covariance should also be constrained to be the same
  p3 ~~ uv*q3 #p2 and q2 covariance
  p4 ~~ uv*q4'
      
      riclpmConstrainedsim <- lavaan(riclpmModel, data = data,
                                     missing = 'ML', 
                                     int.ov.free = F,
                                     int.lv.free = F,
                                     auto.fix.first = F,
                                     auto.fix.single = F,
                                     auto.cov.lv.x = F,
                                     auto.cov.y = F,
                                     auto.var = F)
      
      summary(riclpmConstrainedsim, standardized = T)
      warning <- any(is.na(parameterEstimates(riclpmConstrainedsim, standardized = T)[c(9:16,28:51),]))
      notposdef <- any(eigen(lavInspect(riclpmConstrainedsim, "cov.lv")[3:8, 3:8])$values <= 1e-8)#*
      
      if(warning == TRUE) {
        countNaN <- countNaN+1
        WarningTeller <- WarningTeller + 1
      } else if(notposdef == TRUE) {
        countNotPosDef <- countNotPosDef+1
        WarningTeller <- WarningTeller + 1
        
      } else { # if no warning, then do rest of code
        print(paste0("Model does not contain NAs and cov mx of latents is pos def"))
        
        summary(riclpmConstrainedsim, standardized = T)        
        
    #Random intercept variances
    
    #AIC
    AIC_clpm <- fitMeasures(clpmConstrainedsim, "aic") 
    AIC_riclpm <- fitMeasures(riclpmConstrainedsim, "aic") 
    AIC_riclpm_kappa <- fitMeasures(riclpmModel.kappa.Constrainedsim, "aic") 
    AIC_riclpm_omega <- fitMeasures(riclpmModel.omega.Constrainedsim, "aic") 
        
    
    #ChiBar2DiffTest
    # The full covariance matrix of the random intercepts is:
    S1 <- vcov(riclpmConstrainedsim)[c(9,10), c(9,10)]
    S2 <- vcov(riclpmModel.kappa.Constrainedsim)[9,9]
    S3 <- vcov(riclpmModel.omega.Constrainedsim)[9,9]
    S4 <- vcov(riclpmConstrainedsim)[10,10]
    S5 <- vcov(riclpmConstrainedsim)[9,9]
    
    corr.icWeight <- any(eigen(S1)$values < 0)
    if(corr.icWeight == FALSE) {
    
    # Chi-square values
    Chi2_riclpm <- fitMeasures(riclpmConstrainedsim, "chisq") 
    Chi2_clpm <- fitMeasures(clpmConstrainedsim, "chisq") 
    Chi2_riclpm_kappa <- fitMeasures(riclpmModel.kappa.Constrainedsim, "chisq") 
    Chi2_riclpm_omega <- fitMeasures(riclpmModel.omega.Constrainedsim, "chisq")
      
    #Degrees of freedom (df)
    df_clpm <- fitMeasures(clpmConstrainedsim, "df")
    df_riclpm <- fitMeasures(riclpmConstrainedsim, "df")
    df_riclpm_kappa <- fitMeasures(riclpmModel.kappa.Constrainedsim, "df") 
    df_riclpm_omega <- fitMeasures(riclpmModel.omega.Constrainedsim, "df")
      
    # Run function to do Chi-bar-square test
    ChiBar2DiffTest_clpm_riclpm <- ChiBarSq.DiffTest(2, S1, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm)  # clpm and riclpm, 2= number of variances from riclpm, S = vcov(2x2) from riclpm
    ChiBar2DiffTest_clpm_riclpm.kappa <- ChiBarSq.DiffTest(1, as.matrix(S2), Chi2_clpm, Chi2_riclpm_kappa, df_clpm, df_riclpm_kappa) # clpm and riclpm.kappa, 1= number of variance from riclpm(kappa), S = variance of kappa (1 element only) from riclpm(kappa)
    ChiBar2DiffTest_clpm_riclpm.omega <- ChiBarSq.DiffTest(1, as.matrix(S3), Chi2_clpm, Chi2_riclpm_omega, df_clpm, df_riclpm_omega) # clpm and riclpm.omega, 1= number of variance from riclpm(omega), S = variance of omega (1 element only) from riclpm(omega)
    ChiBar2DiffTest_riclpm_riclpm.kappa <- ChiBarSq.DiffTest(2, as.matrix(S4), Chi2_riclpm_kappa, Chi2_riclpm, df_riclpm_kappa, df_riclpm) # riclpm and riclpm.kappa, 2= number of variances from riclpm, S = variance of omega (1x1) from riclpm
    ChiBar2DiffTest_riclpm_riclpm.omega <- ChiBarSq.DiffTest(2, as.matrix(S5), Chi2_riclpm_omega, Chi2_riclpm, df_riclpm_omega, df_riclpm) # riclpm and riclpm.omega, 2= number of variances from riclpm, S = variance of kappa (1x1) from riclpm
    
    
    #p-values
    pValue_clpm_riclpm <- ChiBar2DiffTest_clpm_riclpm$p_value
    pValue_clpm_riclpm.kappa <- ChiBar2DiffTest_clpm_riclpm.kappa$p_value
    pValue_clpm_riclpm.omega <- ChiBar2DiffTest_clpm_riclpm.omega$p_value
    pValue_riclpm_riclpm.kappa <- ChiBar2DiffTest_riclpm_riclpm.kappa$p_value
    pValue_riclpm_riclpm.omega <- ChiBar2DiffTest_riclpm_riclpm.omega$p_value
    
    
    # Store output
    
    #p-value for Chi-Bar-Square Difference Test
    pvalueChiBar2DiffTest_clpm_riclpm[simteller,] <- pValue_clpm_riclpm
    pvalueChiBar2_clpm_riclpm[simteller,] <- pValue_clpm_riclpm < 0.05
    
    pvalueChiBar2DiffTest_clpm_riclpm.kappa[simteller,] <- pValue_clpm_riclpm.kappa
    pvalueChiBar2_clpm_riclpm.kappa[simteller,] <- pValue_clpm_riclpm.kappa < 0.05
    
    pvalueChiBar2DiffTest_clpm_riclpm.omega[simteller,] <- pValue_clpm_riclpm.omega
    pvalueChiBar2_clpm_riclpm.omega[simteller,] <- pValue_clpm_riclpm.omega < 0.05
    
    pvalueChiBar2DiffTest_riclpm_riclpm.kappa[simteller,] <- pValue_riclpm_riclpm.kappa
    pvalueChiBar2_riclpm_riclpm.kappa[simteller,] <- pValue_riclpm_riclpm.kappa < 0.05
    
    pvalueChiBar2DiffTest_riclpm_riclpm.omega[simteller,] <- pValue_riclpm_riclpm.omega
    pvalueChiBar2_riclpm_riclpm.omega[simteller,] <- pValue_riclpm_riclpm.omega < 0.05
    
    
    
    #p-value for Chi-Square Difference Test
    pvalueChi2DiffTest_clpm_riclpm[simteller,] <- 1-pchisq((Chi2_clpm-Chi2_riclpm),(df_clpm-df_riclpm))
    pvalueChi2_clpm_riclpm[simteller,] <- (1-pchisq((Chi2_clpm-Chi2_riclpm),(df_clpm-df_riclpm))) < 0.05
    
    pvalueChi2DiffTest_clpm_riclpm.kappa[simteller,] <- 1-pchisq((Chi2_clpm-Chi2_riclpm_kappa),(df_clpm-df_riclpm_kappa))
    pvalueChi2_clpm_riclpm.kappa[simteller,] <- (1-pchisq((Chi2_clpm-Chi2_riclpm_kappa),(df_clpm-df_riclpm_kappa))) < 0.05 
    
    pvalueChi2DiffTest_clpm_riclpm.omega[simteller,] <- 1-pchisq((Chi2_clpm-Chi2_riclpm_omega),(df_clpm-df_riclpm_omega))
    pvalueChi2_clpm_riclpm.omega[simteller,] <- (1-pchisq((Chi2_clpm-Chi2_riclpm_omega),(df_clpm-df_riclpm_omega))) < 0.05
    
    pvalueChi2DiffTest_riclpm_riclpm.kappa[simteller,] <- 1-pchisq((Chi2_riclpm_kappa-Chi2_riclpm),(df_riclpm_kappa-df_riclpm))
    pvalueChi2_riclpm_riclpm.kappa[simteller,] <- (1-pchisq((Chi2_riclpm_kappa-Chi2_riclpm),(df_riclpm_kappa-df_riclpm))) < 0.05 
    
    pvalueChi2DiffTest_riclpm_riclpm.omega[simteller,] <- 1-pchisq((Chi2_riclpm_omega-Chi2_riclpm),(df_riclpm_omega-df_riclpm))
    pvalueChi2_riclpm_riclpm.omega[simteller,] <- (1-pchisq((Chi2_riclpm_omega-Chi2_riclpm),(df_riclpm_omega-df_riclpm))) < 0.05
    
    
    #AIC 
    # Akaike weights (based on AIC)
    AIC_clpm_riclpm[simteller,] <- c(AIC_clpm, AIC_riclpm)
    Name_Hypo1 <- c("CLPM","RICLPM")
    AICweightsCAl_clpm_riclpm  <-  IC.weights(AIC_clpm_riclpm[simteller,], Name_Hypo1)$IC.weights
    AICweights_clpm_riclpm[simteller,] <- AICweightsCAl_clpm_riclpm
    
    AIC_clpm_riclpm.kappa[simteller,] <- c(AIC_clpm, AIC_riclpm_kappa)
    Name_Hypo2 <- c("CLPM","RICLPM.kappa")
    AICweightsCAl_clpm_riclpm.kappa  <-  IC.weights(AIC_clpm_riclpm.kappa[simteller,], Name_Hypo2)$IC.weights
    AICweights_clpm_riclpm.kappa[simteller,] <- AICweightsCAl_clpm_riclpm.kappa
    
    AIC_clpm_riclpm.omega[simteller,] <- c(AIC_clpm, AIC_riclpm_omega)
    Name_Hypo3 <- c("CLPM","RICLPM.omega")
    AICweightsCAl_clpm_riclpm.omega  <-  IC.weights(AIC_clpm_riclpm.omega[simteller,], Name_Hypo3)$IC.weights
    AICweights_clpm_riclpm.omega[simteller,] <- AICweightsCAl_clpm_riclpm.omega
    
    AIC_riclpm_riclpm.kappa[simteller,] <- c(AIC_riclpm, AIC_riclpm_kappa)
    Name_Hypo3 <- c("RICLPM","RICLPM.kappa")
    AICweightsCAl_riclpm_riclpm.kappa  <-  IC.weights(AIC_riclpm_riclpm.kappa[simteller,], Name_Hypo3)$IC.weights
    AICweights_riclpm_riclpm.kappa[simteller,] <- AICweightsCAl_riclpm_riclpm.kappa
    
    AIC_riclpm_riclpm.omega[simteller,] <- c(AIC_riclpm, AIC_riclpm_omega)
    Name_Hypo4 <- c("RICLPM","RICLPM.omega")
    AICweightsCAl_riclpm_riclpm.omega  <-  IC.weights(AIC_riclpm_riclpm.omega[simteller,], Name_Hypo4)$IC.weights
    AICweights_riclpm_riclpm.omega[simteller,] <- AICweightsCAl_riclpm_riclpm.omega
    
    simteller<- simteller+1
  }
  
  #ignore warnings
  on.exit(options(warn = oldw))
  }
  }
  }
}
  round <- round + 1

