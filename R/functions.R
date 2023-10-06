library("dplyr");library("data.table");library("stats");library("base")
library("graphics");library("grDevices");library("copula");library("VineCopula")

# Indicator function
ind <- function(Yi, y)
  ifelse(Yi >= y, 1, 0)

# Normalization function
Fn <- function(X, x) {
  sum <- 0
  for (i in c(1:length(X))) {
    sum <- sum + ind(X[i], x)
  }
  return(sum / (length(X) + 1))
}

normalize <-function(variable){
  variable <- lapply(variable, function(x) {Fn(variable, x)}) %>% unlist
}

# Power link function
g<-function(x,c) x^(c)

# Quantile function
yn<-function(Y,alpha){
  return(unname(quantile(Y, alpha)))
}

# Simulation of Pareto distribution
pareto.simu<-function(pareto.param,a,n){
  u<-runif(n,0,1)
  y<- a*u^(-1/pareto.param)
  return(y)}

# Kendall's tau from copula parameter (archimedean)
tau_arch_copula <- function(eta, copula){

  if (tolower(copula) == "copula2"){
    alpha <- eta[1]
    kappa <- eta[2]
    output <- 1 - 2 * alpha * kappa/(1 + 2 * kappa) }

  if (tolower(copula) != "copula2") {
    output <- tau(archmCopula(tolower(copula), param = eta, dim = 2)) }
  return(output)}

# Kendall's tau from copula parameter (Normal)
tau_normal_copula<-function(r){return(2/pi*(asin(r)))}

#__________________________ 3.1 Simulation of X and Y function ________________________________________________________#

# Description: this function returns a list containing the vector of simulated Y and simulated X according to model M_1

# Inputs:
#n: integer, sample size
#p: integer, dimension
#beta: numeric, vector of dimension direction
#power: numeric, vector of link function power
#sigma: numeric, matrix of standard deviaton of noise
#r: numeric, signal noise
#dist: charachter, name of Y distribution  (pareto or student)
#dist.param: numeric, value or vector of Y distribution parameters
#copula.fam: integer, defining the bivariate copula family:
#    1 = Gaussian copula
#    5 = Frank copula
#    3 = Clayton copula
#copula.param: numeric, copula parameter
# signe_tau:
#    If == 1  : positive association
#    If == -1 : negative association

# Output: the first index of the list returns the simulated Y, then the other indices return the simulated X according to the link function. For example if the size of the power vector is 4,
#that is to say we will have 4 matrices X. See example below

#X.Y <- X_Y_simu(n,p,beta, power...) where power=c(3/2,1,1/2,1/4)
#X.Y[[1]] return the Y sampling according the chosen distribution Pareto or Student
#X.Y[[2]] return the matrix X1 corresponding to link function g:x->x^{3/2}
#X.Y[[3]] return the matrix X2 corresponding to link function g:x->x^1
#X.Y[[3]] return the matrix X3 corresponding to link function g:x->sqrt(x)
#X.Y[[4]] return the matrix X4 corresponding to link function g:x->x^{1/4}


# Simulation of X and Y function
X_Y_simu <- function(n,p,beta,power,sigma,r,
                     dist,dist.param,copula.fam,
                     copula.param,signe_tau){
  # Simulation of U0
  U0 <- runif(n)
  # Simulation of Y
  Y <- dist.param[2]*(1-U0)^(-1/dist.param[1])
  # Simulation on the conditional quantile
  obj <- BiCop(family = copula.fam, par = copula.param)
  V1<-matrix(0,n,p)
  if(signe_tau < 0)
  {
    U0 <- 1 - U0
    # r <- 5*r ## Make a SNR more important
  }
  for(j in 1:p)
  {
    V1[,j] <- BiCopCondSim(n, cond.val = U0, cond.var = 1, obj)
  }
  List_X<-list()
  List_X[[1]]<-Y
  if(dist=="pareto"){
    # sd<-vector()
    for(l in 1:length(power)){
      sd_l<- g(dist.param[2]*n^{1/dist.param[1]},power[l])/(2*r)
      c_s_l <- c_s[l]
      list_x_p <- matrix(NA,n,p)
      for(j in 1:p){
        # noise <- qnorm(V1[,j],mean=0,sd=sd_l)
        noise <- abs(qnorm(V1[,j],mean=0,sd=sd_l))
        list_x_p[,j] <- g(Y,c_s_l)*beta[j] + noise
      }
      # plot(Y,list_x_p%*%beta) ; abline(0,1,col="red")
      # browser()
      List_X[[l+1]] <- list_x_p
      # sd[l] <- sd_l
      # List_X[[l+1]]<-g(Y,c_s[l])%*%t(beta) +
      #   qnorm(V1,mean=0,sd=sd[l]) #scale(qnorm(V1,mean=0,sd=sd[l]),scale = FALSE)
    }
  }
  return(List_X)
}

#__________________________ 3.2 Simulation process function ________________________________________________________#

# Description: this function returns a table containing the mean proximity criterion for each number of exceedances and link function power

# Inputs:
#n: integer, sample size
#p: integer, dimension
#beta: numeric, vector of dimension direction
#power: numeric, vector of link function power
#sigma: numeric, matrix of standard deviaton of noise
#r: numeric, signal noise
#dist: charachter, name of Y distribution  (pareto or student)
#dist.param: numeric, value or vector of Y distribution parameters
#copula.fam: integer, defining the bivariate copula family:
#    1 = Gaussian copula
#    5 = Frank copula
#    3 = Clayton copula
#copula.param: numeric, copula parameter
#N: integer, number of replications
#k.threshold<- integer, sample size coresponding to the minimum threshold


# Output: Dataframe containing the mean proximity criterion for each number of exceedances and link function power

# Simulation process function
simu_process <-function(mu0=NULL,kappa0=NULL,n,p,beta,c_s,sigma,r,
                        dist,pareto.params,copula.fam,copula.param,N,
                        k.threshold,type="vMF",signe_tau=1){
  #Cos: results storage table containing: simu = number of replications,
  #k=number of exceedances, cos^2 for each link function power
  n_powers <- length(c_s)
  n_kappas <- length(kappa0)
  n_exceedances <- k.threshold#length(k.threshold:length(Y))
  OUT <- matrix(NA,N*n_powers*n_kappas*n_exceedances,6)
  iter <- 1

  ##
  values_replications <- 1:N
  values_exceedances <- k.threshold:n
  values_powers <- 1:n_powers
  values_kappas <- 1:n_kappas
  parameters <- as.matrix(expand.grid(values_exceedances,values_powers,values_kappas))
  N_parameters <- nrow(parameters)
  all_values_output <- matrix(0,N_parameters,N)
  #Begining of the algorithm
  for(j in values_replications)
  {
    # cat(paste("Replication",j,"\n"))
    #Simulation of X and Y
    X.Y <- X_Y_simu(n,p,beta,c_s,sigma,r,dist,pareto.params,
                    copula.fam,copula.param,signe_tau=signe_tau)
    #Estimation of beta
    Y<-sort(X.Y[[1]]) #Simulated Y
    for(ii in 1:N_parameters)
    {
      parameter_ii <- parameters[ii,]
      i <- parameter_ii[1] ; l <- parameter_ii[2] ; i_kappa <- parameter_ii[3]
      y<-Y[i]
      k<-(1-(i/n))*n
      if(type=="vMF")
      {
        w <- SEPaLS(X=X.Y[[l+1]],Y=X.Y[[1]],yn=y,type = type,
                    mu0 = mu0,kappa0 = kappa0[i_kappa])
      }
      else if (type=="Laplace")
      {
        w <- SEPaLS(X=X.Y[[l+1]],Y=X.Y[[1]],yn=y,type = type,
                    lambda = kappa0[i_kappa])
      }
      all_values_output[ii,j] <- (t(w)%*%beta)^2
    }
  }

  res <- t(apply(all_values_output,1,function(cc){
    quantile(cc,probs=c(0.5,0.05,0.95) )
  }))
  out <- cbind(parameters,res)

  if(type=="vMF")
  {
    colnames(out) <- c("nb_exceed","id_power","kappa0",
                       "median","quant5","quant95")
  }
  else if (type=="Laplace")
  {
    colnames(out) <- c("nb_exceed","id_power","lambda",
                       "median","quant5","quant95")
  }
  return(out)
}
