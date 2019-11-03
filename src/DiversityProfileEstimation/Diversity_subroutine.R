CV.Ind=function(x)
{
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1); f2 <- sum(x == 2)
  if(f2 > 0){
    A <- 2*f2/((n-1)*f1+2*f2)
  }else if(f1 > 1){
    A <- 2/((n-1)*(f1-1)+2)
  }else{
    A <- 0
  }
  C.hat <- 1 - f1/n*(1-A)
  Sobs=sum(x>0)
  S0=Sobs/C.hat
  r.square=max(S0*sum(x*(x-1))/n/(n-1)-1, 0)
  r.square^0.5
}

##################################2015.09.28-(S.W.Wei)
CV.Sam <- function(y)
{
  t <- y[1]
  y <- y[-c(1)]
  y <- y[y > 0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  if(Q2 > 0){
    A <- 2*Q2/((t-1)*Q1+2*Q2)
  }else if(Q1 > 1){
    A <- 2/((t-1)*(Q1-1)+2)
  }else{
    A <- 0
  }
  C.hat <- 1 - Q1/U*(1-A)
  Sobs <- sum(y > 0)
  S0 <- Sobs/C.hat
  r.square <- max(S0*t/(t-1)*sum(y*(y-1))/U^2 - 1, 0)
  return(sqrt(r.square))
}
##################################2015.09.28

Chao1=function(x,conf=0.95)
{
  z <--qnorm((1 - conf)/2)
  x=x[x>0]
  D=sum(x>0)
  f1=sum(x==1)
  f2=sum(x==2)
  n=sum(x)
  if (f1 > 0 & f2 > 0)
  {
    S_Chao1 <- D + (n - 1)/n*f1^2/(2*f2)
    var_Chao1 <- f2*((n - 1)/n*(f1/f2)^2/2 + 
                       ((n - 1)/n)^2*(f1/f2)^3 + ((n - 1 )/n)^2*(f1/f2)^4/4)
    
    t <- S_Chao1 - D
    K <- exp(z*sqrt(log(1 + var_Chao1/t^2)))
    CI_Chao1 <- c(D + t/K, D + t*K)
  } 
  else if (f1 > 1 & f2 == 0)
  {
    S_Chao1 <- D + (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
    var_Chao1 <- (n - 1)/n*f1*(f1 - 1)/2 + 
      ((n - 1)/n)^2*f1*(2*f1 - 1)^2/4 - ((n - 1)/n)^2*f1^4/4/S_Chao1
    
    t <- S_Chao1 - D
    K <- exp(z*sqrt(log(1 + var_Chao1/t^2)))
    CI_Chao1 <- c(D + t/K, D + t*K)
  } 
  else 
  {
    S_Chao1 <- D
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i) sum(x==i)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*sum(x==i))))^2/n
    var_Chao1 <- var_obs
    P <- sum(sapply(i, function(i) sum(x==i)*exp(-i)/D))
    CI_Chao1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  return( c( round(c(S_Chao1,var_Chao1^0.5,CI_Chao1[1],CI_Chao1[2]),1),conf)    )
}

Chao1_bc=function(x,conf=0.95)
{
  z <- -qnorm((1 - conf)/2)
  x=x[x>0]
  D=sum(x>0)
  f1=sum(x==1)
  f2=sum(x==2)
  n=sum(x)
  
  S_Chao1_bc <- D + (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
  var_Chao1_bc <- (n - 1)/n*f1*(f1 - 1)/2/(f2 + 1) + 
    ((n - 1)/n)^2*f1*(2*f1 - 1)^2/4/(f2 + 1)^2 + ((n - 1)/n)^2*f1^2*f2*(f1 - 1)^2/4/(f2 + 1)^4
  
  t <- round(S_Chao1_bc - D, 5)
  if (t != 0)
  {
    K <- exp(z*sqrt(log(1 + var_Chao1_bc/t^2)))
    CI_Chao1_bc <- c(D + t/K, D + t*K)
  } 
  if(t == 0)
  {
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)sum(x==i)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*sum(x==i))))^2/n
    var_Chao1_bc <- var_obs
    P <- sum(sapply(i, function(i)sum(x==i)*exp(-i)/D))
    CI_Chao1_bc <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  round(c(S_Chao1_bc,var_Chao1_bc^0.5,CI_Chao1_bc[1],CI_Chao1_bc[2]) ,1)  
}

SpecAbunAce<- function(data, k=10, conf=0.95)
{
  data <- as.numeric(data)  
  f <- function(i, data){length(data[which(data == i)])}
  basicAbun <- function(data, k)
  {
    x <- data[which(data != 0)]
    n <- sum(x)
    D <- length(x)
    n_rare <- sum(x[which(x <= k)])
    D_rare <- length(x[which(x <= k)])
    if (n_rare != 0)
    {
      C_rare <- 1 - f(1, x)/n_rare
    } 
    else 
    {
      C_rare = 1
    } 
    n_abun <- n - n_rare
    D_abun <- length(x[which(x > k)])
    
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
    a2 <- sum(sapply(j, function(j)j*f(j, x)))
    if (C_rare != 0 & a2 >1)
    {
      gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
      gamma_rare_1_square <- max(gamma_rare_hat_square*(1 + (1 - C_rare)/C_rare*a1/(a2 - 1)), 0)
    }
    else
    {
      gamma_rare_hat_square <- 0
      gamma_rare_1_square <- 0
    }
    CV_rare <- sqrt(gamma_rare_hat_square)
    CV1_rare <- sqrt(gamma_rare_1_square)
    
    return(c( n, D, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun))
  }
  
  z <- -qnorm((1 - conf)/2)
  n <- basicAbun(data, k)[1]
  D <- basicAbun(data, k)[2]
  n_rare <- basicAbun(data, k)[3]
  D_rare <- basicAbun(data, k)[4]
  C_rare <- basicAbun(data, k)[5] 
  CV_rare <- basicAbun(data, k)[6]
  CV1_rare <- basicAbun(data, k)[7]
  n_abun <- basicAbun(data, k)[8]
  D_abun <- basicAbun(data, k)[9]
  x <- data[which(data != 0)]
  #############################
  S_ACE <- function(x, k)
  {
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
    a2 <- sum(sapply(j, function(j)j*f(j, x)))
    if (C_rare != 0 & a2 >1)
    {
      gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
    }
    else
    {
      gamma_rare_hat_square <- 0
    }
    S_ace <- D_abun + D_rare/C_rare + f(1, x)/C_rare*gamma_rare_hat_square
    #return(list(S_ace, gamma_rare_hat_square))
    return(c(S_ace, gamma_rare_hat_square))
  }
  s_ace <- S_ACE(x, k)[1]
  gamma_rare_hat_square <- S_ACE(x, k)[2]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    if (gamma_rare_hat_square != 0){
      si <- sum(sapply(u, function(u)u*(u - 1)*f(u, x)))
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(D_rare*si + f(1, x)*si) - 
             f(1, x)*D_rare*si*(-2*(1 - f(1, x)/n_rare)*(n_rare - f(1, x))/n_rare^2*n_rare*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*(2*n_rare - 1))
          )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2 - #g2
          (1 - f(1, x)/n_rare + f(1, x)*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g3
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)*(si + D_rare*q*(q - 1)) - 
             f(1, x)*D_rare*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                                  (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
          )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2 + #g2
          (q*(f(1, x))^2/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g3
      }
      return(d)
    } else {
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1 
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1
      }
      return(d)  
    }
  }
  COV.f <- function(i,j){
    if (i == j){
      cov.f <- f(i, x)*(1 - f(i, x)/s_ace)
    } else {
      cov.f <- -f(i, x)*f(j, x)/s_ace
    }     
    return(cov.f)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ace <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.f(i, j), i, j))
  if (var_ace > 0){
    var_ace <- var_ace
  } else {
    var_ace <- NA
  }
  ######################
  t <- round(s_ace - D, 5)
  if (is.nan(t) == F){
    if (t != 0){
      C <- exp(z*sqrt(log(1 + var_ace/(s_ace - D)^2)))
      CI_ACE <- c(D + (s_ace - D)/C, D + (s_ace - D)*C)
    } else {
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
      var_ace <- var_obs
      P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
      CI_ACE <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
    }
  }else{
    CI_ACE <- c(NaN, NaN)
  }
  table <- c(s_ace, sqrt(var_ace), CI_ACE)
  return(table)
}

SpecAbunAce1 <-function(data ,k=10, conf=0.95)
{
  data <- as.numeric(data)
  f <- function(i, data){length(data[which(data == i)])}
  basicAbun <- function(data, k){
    x <- data[which(data != 0)]
    n <- sum(x)
    D <- length(x)
    n_rare <- sum(x[which(x <= k)])
    D_rare <- length(x[which(x <= k)])
    if (n_rare != 0){
      C_rare <- 1 - f(1, x)/n_rare
    } else {
      C_rare = 1
    } 
    n_abun <- n - n_rare
    D_abun <- length(x[which(x > k)])
    
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
    a2 <- sum(sapply(j, function(j)j*f(j, x)))
    if (C_rare != 0 & a2 >1){
      gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
      gamma_rare_1_square <- max(gamma_rare_hat_square*(1 + (1 - C_rare)/C_rare*a1/(a2 - 1)), 0)
    }else{
      gamma_rare_hat_square <- 0
      gamma_rare_1_square <- 0
    }
    CV_rare <- sqrt(gamma_rare_hat_square)
    CV1_rare <- sqrt(gamma_rare_1_square)
    return(c( n, D, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun))
  }
  
  z <- -qnorm((1 - conf)/2)
  
  n <- basicAbun(data, k)[1]
  D <- basicAbun(data, k)[2]
  n_rare <- basicAbun(data, k)[3]
  D_rare <- basicAbun(data, k)[4]
  C_rare <- basicAbun(data, k)[5] 
  CV_rare <- basicAbun(data, k)[6]
  CV1_rare <- basicAbun(data, k)[7]
  n_abun <- basicAbun(data, k)[8]
  D_abun <- basicAbun(data, k)[9]
  x <- data[which(data != 0)]
  #############################
  S_ACE1 <- function(x, k){
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
    a2 <- sum(sapply(j, function(j)j*f(j, x)))
    if (C_rare != 0 & a2 >1){
      gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
      gamma_rare_1_square <- max(gamma_rare_hat_square*(1 + (1 - C_rare)/C_rare*a1/(a2 - 1)), 0)
    }else{
      gamma_rare_hat_square <- 0
      gamma_rare_1_square <- 0
    }
    s_ace1 <- D_abun + D_rare/C_rare + f(1, x)/C_rare*gamma_rare_1_square
    
    return(c(s_ace1, gamma_rare_1_square))
  }
  s_ace1 <- S_ACE1(x, k)[1]
  gamma_rare_1_square <- S_ACE1(x, k)[2]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    if (gamma_rare_1_square != 0){
      u <- c(1:k)
      si <- sum(sapply(u, function(u)u*(u-1)*f(u, x)))
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(D_rare*si + f(1, x)*si) - 
             f(1, x)*D_rare*si*(-2*(1 - f(1, x)/n_rare)*(n_rare - f(1, x))/n_rare^2*n_rare*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*(2*n_rare - 1))
          )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2 - #g2
          (1 - f(1, x)/n_rare + f(1, x)*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g3
          ((1 - f(1, x)/n_rare)^3*(n_rare*(n_rare - 1))^2*(2*f(1, x)*D_rare*si^2 + f(1, x)^2*si^2) - #g4
             f(1, x)^2*D_rare*si^2*(3*(1 - f(1, x)/n_rare)^2*(f(1, x) - n_rare)/(n_rare)^2*(n_rare*(n_rare - 1))^2 + 
                                      (1 - f(1, x)/n_rare)^3*2*n_rare*(n_rare - 1)^2 + (1 - f(1, x)/n_rare)^3*n_rare^2*2*(n_rare - 1)) 
          )/(1 - f(1, x)/n_rare)^6/n_rare^4/(n_rare - 1)^4 - 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(2*f(1, x)*si) - #g5
             f(1, x)^2*si*(2*(1 - f(1, x)/n_rare)*(f(1, x) - n_rare)/n_rare^2*n_rare*(n_rare - 1) + 
                             (1 - f(1, x)/n_rare)^2*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare) 
          )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)*(si + D_rare*q*(q - 1)) - 
             f(1, x)*D_rare*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                                  (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
          )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2 + #g2
          (q*(f(1, x))^2/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g3
          ((1 - f(1, x)/n_rare)^3*n_rare^2*(n_rare - 1)^2*f(1, x)^2*(si^2 + 2*D_rare*si*q*(q - 1)) - #g4
             f(1, x)^2*D_rare*si^2*(3*(1 - f(1, x)/n_rare)^2*(f(1, x)*q/n_rare^2)*(n_rare*(n_rare - 1))^2 + 
                                      2*(1 - f(1, x)/n_rare)^3*n_rare*q*(n_rare - 1)^2 + 2*(1 - f(1, x)/n_rare)^3*n_rare^2*(n_rare - 1)*q)   
          )/(1 - f(1, x)/n_rare)^6/(n_rare)^4/(n_rare - 1)^4 - 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)^2*q*(q - 1) - #g5
             f(1, x)^2*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                             (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
          )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2
      }
      return(d)
    } else {
      u <- c(1:k)
      si <- sum(sapply(u, function(u)u*(u-1)*f(u, x)))
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1 
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1
      }
      return(d)
    }
  }
  
  COV.f <- function(i,j){
    if (i == j){
      cov.f <- f(i, x)*(1 - f(i, x)/s_ace1)
    } else {
      cov.f <- -f(i, x)*f(j, x)/s_ace1
    }     
    return(cov.f)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ace1 <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.f(i, j), i, j))
  if (var_ace1 > 0){
    var_ace1 <- var_ace1
  } else {
    var_ace1 <- NA
  }
  ######################
  t <- round(s_ace1 - D, 5)
  if (is.nan(t) == F){
    if (t != 0){
      C <- exp(z*sqrt(log(1 + var_ace1/(s_ace1 - D)^2)))
      CI_ACE1 <- c(D + (s_ace1 - D)/C, D + (s_ace1 - D)*C)
    } else {
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
      var_ace1 <- var_obs
      P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
      CI_ACE1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
    }
  }else{
    CI_ACE1 <- c(NaN, NaN)
  }
  table <-c(s_ace1, sqrt(var_ace1), CI_ACE1)
  return(table)
}


EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0) 	#observed species
  n <- sum(Spec)		  	#sample size
  f1 <- sum(Spec == 1) 	#singleton 
  f2 <- sum(Spec == 2) 	#doubleton
  a <- ifelse(f1 == 0, 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n)
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  w <- a / b  			#adjusted factor for rare species in the sample
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))	#estimation of unseen species via Chao1
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(2 * f2/((n - 1) * f1 + 2 * f2), f0.hat)		#estimation of relative abundance of unseen species in the sample
  return(c(Prob.hat, Prob.hat.Unse))				#Output: a vector of estimated relative abundance
}


entropy_MEE_equ=function(X)
{
  x=X
  x=x[x>0]
  n=sum(x)
  UE <- sum(x/n*(digamma(n)-digamma(x)))
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if(f1>0)
  {
    A <-1-ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
    B=sum(x==1)/n*(1-A)^(-n+1)*(-log(A)-sum(sapply(1:(n-1),function(k){1/k*(1-A)^k})))
  }
  if(f1==0){B=0}
  UE+B
}
entropy_HT_equ<-function(X)
{
  x=X
  x=x[x>0]
  n=sum(x)
  f1=sum(x==1)
  C_head=1-f1/n
  a=-sum(C_head*(x/n)*log(C_head*(x/n))/(1-(1-C_head*(x/n))^n))
  a
}
entropy_J1_equ=function(X)
{
  X=X[X>0]
  Y=X[X>1]
  n=sum(X)
  -n*sum(X/n*log(X/n))-(n-1)/n*sum( (n-X)*(-X/(n-1)*log(X/(n-1))) )-(n-1)/n*sum(-Y*(Y-1)/(n-1)*log((Y-1)/(n-1)))   
}
entropy_MLE_equ=function(X)
{
  X=X[X>0]
  n=sum(X)
  -sum(X/n*log(X/n))
}
# entropy_MLE_bc_equ=function(X)
# {
#   entropy_MLE_equ(X)+(SpecAbunAce(X)[1]-1)/2/sum(X)
# }
entropy_MLE_bc_equ=function(X)
{
  entropy_MLE_equ(X)+(SpecAbunChao1(X,k=10,conf=0.95)[1]-1)/2/sum(X)
}
Shannon_index=function(x,boot=50)
{
  x=x[x>0]
  n=sum(x)
  MLE=entropy_MLE_equ(x) 
  MLE_bc=entropy_MLE_bc_equ(x)
  J1=entropy_J1_equ(x)
  HT=entropy_HT_equ(x)
  MEE=entropy_MEE_equ(x)
  p_hat=EstiBootComm.Ind(x)
  if(boot<=1){
    MLE_sd=NA
    MLE_bc_sd=NA
    J1_sd=NA
    HT_sd=NA
    MEE_sd=NA
    MLE_exp_sd=NA
    MLE_bc_exp_sd=NA
    J1_exp_sd=NA
    HT_exp_sd=NA
    MEE_exp_sd=NA
  }else{
    Boot.X=rmultinom(boot,n,p_hat)
    temp1=apply(Boot.X,2,entropy_MLE_equ)
    temp2=apply(Boot.X,2,entropy_MLE_bc_equ)
    temp3=apply(Boot.X,2,entropy_J1_equ)
    temp4=apply(Boot.X,2,entropy_HT_equ)
    temp5=apply(Boot.X,2,entropy_MEE_equ)
    MLE_sd=sd(temp1)
    MLE_bc_sd=sd(temp2)
    J1_sd=sd(temp3)
    HT_sd=sd(temp4)
    MEE_sd=sd(temp5)
    
    MLE_exp_sd=sd(exp(temp1))
    MLE_bc_exp_sd=sd(exp(temp2))
    J1_exp_sd=sd(exp(temp3))
    HT_exp_sd=sd(exp(temp4))  
    MEE_exp_sd=sd(exp(temp5))
  }
  
  a=matrix(0,10,4)
  a[1,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  a[2,]=c(MLE_bc,MLE_bc_sd,MLE_bc-1.96*MLE_bc_sd,MLE_bc+1.96*MLE_bc_sd)
  a[3,]=c(J1,J1_sd,J1-1.96*J1_sd,J1+1.96*J1_sd)
  a[4,]=c(HT,HT_sd,HT-1.96*HT_sd,HT+1.96*HT_sd)
  a[5,]=c(MEE,MEE_sd,MEE-1.96*MEE_sd,MEE+1.96*MEE_sd)
  a[6,]=c(exp(MLE),MLE_exp_sd,exp(MLE)-1.96*MLE_exp_sd,exp(MLE)+1.96*MLE_exp_sd)
  a[7,]=c(exp(MLE_bc),MLE_bc_exp_sd,exp(MLE_bc)-1.96*MLE_bc_exp_sd,exp(MLE_bc)+1.96*MLE_bc_exp_sd)
  a[8,]=c(exp(J1),J1_exp_sd,exp(J1)-1.96*J1_exp_sd,exp(J1)+1.96*J1_exp_sd)
  a[9,]=c(exp(HT),HT_exp_sd,exp(HT)-1.96*HT_exp_sd,exp(HT)+1.96*HT_exp_sd)
  a[10,]=c(exp(MEE),MEE_exp_sd,exp(MEE)-1.96*MEE_exp_sd,exp(MEE)+1.96*MEE_exp_sd)
  return(a)
}
simpson_MLE_equ=function(X)
{
  X=X[X>0]
  n=sum(X)
  a=sum((X/n)^2)
  a 
}
simpson_MVUE_equ=function(X)
{
  X=X[X>0]
  n=sum(X)
  a=sum(X*(X-1))/n/(n-1)
  a 
}
Simpson_index=function(x,boot=200)
{
  x=x[x>0]
  n=sum(x)
  MVUE=simpson_MVUE_equ(x)
  MLE=simpson_MLE_equ(x)
  
  ###############################################2015.10.04-(S.W.Wei)
  #ACE=SpecAbunAce(x)[1]
  #AA=sum(  ( x*(x-1)/n/(n-1)-x*(2*n-1)/n/(n-1)*MVUE  )^2  )
  #BB=sum( x*(x-1)/n/(n-1)-x*(2*n-1)/n/(n-1)*MVUE   )
  #MVUE_sd=(AA-BB^2/ACE)^0.5 
  #AA=sum(  ( (x/n)^2-2*x/n*MLE  )^2  )
  #BB=sum( (x/n)^2-2*x/n*MLE   )
  #MLE_sd=(AA-BB^2/ACE)^0.5 
  #MVUE_recip_sd=MVUE_sd/MVUE
  #MLE_recip_sd=MLE_sd/MLE
  
  p_hat=EstiBootComm.Ind(x)
  if(boot<=1){
    MVUE_sd = NA
    MLE_sd = NA
    
    MVUE_recip_sd = NA
    MLE_recip_sd = NA
  } else {
    Boot.X=rmultinom(boot,n,p_hat)
    temp1 = apply(Boot.X, 2, simpson_MVUE_equ)
    temp2 = apply(Boot.X, 2, simpson_MLE_equ)
    
    MVUE_sd = sd(temp1)
    MLE_sd = sd(temp2)
    
    MVUE_recip_sd = sd(temp1^(-1))
    MLE_recip_sd = sd(temp2^(-1))
  }
  ###############################################2015.10.04
  
  a=matrix(0,4,4)
  a[1,]=c(MVUE,MVUE_sd,MVUE-1.96*MVUE_sd,MVUE+1.96*MVUE_sd)
  a[2,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  a[3,]=c(1/MVUE,MVUE_recip_sd,1/MVUE-1.96*MVUE_recip_sd,1/MVUE+1.96*MVUE_recip_sd)
  a[4,]=c(1/MLE,MLE_recip_sd,1/MLE-1.96*MLE_recip_sd,1/MLE+1.96*MLE_recip_sd)
  return(a)
}

###############################################2015.09.15(H.W.Hsu, Y.R.Chen, S.W.Wei)
SpecInci <- function(data, k=10, conf=0.95)
{
  Chao2   <- SpecInciChao2(data, k = k, conf = conf)
  Chao2bc <- SpecInciChao2bc(data, k = k, conf = conf)
  iChao2  <- SpecInciiChao2(data, k = k, conf = conf)
  Modelh  <- SpecInciModelh(data, k = k, conf = conf)[-c(5)]
  Modelh1 <- SpecInciModelh1(data, k = k, conf = conf)[-c(5)]
  table   <- rbind(Chao2, Chao2bc, iChao2, Modelh, Modelh1)
  table   <- round(table,1)
  colnames(table) <- c("Estimate", "   s.e.", "95%Lower", "95%Upper")
  rownames(table) <- c("Chao2 (Chao, 1987)", "Chao2-bc", "iChao2", "ICE (Lee & Chao, 1994)", "ICE-1 (Lee & Chao, 1994)")
  return(table)
}

SpecInci_raw <- function(data, k=10, conf=0.95)
{
  Chao2   <- SpecInciChao2(data[-1], k = k, conf = conf)
  Chao2bc <- SpecInciChao2bc(data[-1], k = k, conf = conf)
  iChao2  <- SpecInciiChao2(data[-1], k = k, conf = conf)
  Modelh  <- SpecInciModelh(data[-2], k = k, conf = conf)[-c(5)]
  Modelh1 <- SpecInciModelh1(data[-2], k = k, conf = conf)[-c(5)]
  table   <- rbind(Chao2, Chao2bc, iChao2, Modelh, Modelh1)
  table   <- round(table,1)
  colnames(table) <- c("Estimate", "   s.e.", "95%Lower", "95%Upper")
  rownames(table) <- c("Chao2 (Chao, 1987)", "Chao2-bc", "iChao2", "ICE (Lee & Chao, 1994)", "ICE-1 (Lee & Chao, 1994)")
  return(table)
}

entropy_MLE_Inci_equ <- function(X)
{
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  H_MLE <- -sum(X/U*log(X/U))
  return(H_MLE)
}
entropy_MLE_bc_Inci_equ <- function(X)
{
  t <- X[1]
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  X_freq <- X[X > 10]
  X_infreq <- X[X <= 10]
  D_freq <- length(X_freq)
  D_infreq <- length(X_infreq)
  Q1 <- sum(X == 1)
  Q2 <- sum(X == 2)
  if(Q1 > 0 & Q2 > 0)
  {
    A <- 2*Q2/((t-1)*Q1 + 2*Q2)
  } 
  else if (Q1 > 0 & Q2 == 0)
  {
    A <- 2/((t-1)*(Q1 - 1) + 2)
  } 
  else 
  {
    A <- 1
  }
  C_infreq <- 1 - Q1/sum(X_infreq)*(1-A)
  
  j <- c(1:10)
  b1 <- sum(sapply(j, function(j){j*(j-1)*sum(X == j)}))
  b2 <- sum(sapply(j, function(j){j*sum(X == j)}))
  gamma_infreq_square <- max(D_infreq/C_infreq*t/(t-1)*b1/b2/(b2-1) - 1, 0)
  
  ICE <- D_freq + D_infreq/C_infreq + Q1/C_infreq*gamma_infreq_square
  
  H_MLE <- -sum(X/U*log(X/U))
  H_MLE_bc <- H_MLE + (ICE/U + 1/t)/1
  
  return(H_MLE_bc)
}
entropy_HT_Inci_equ <- function(X)
{
  t <- X[1]
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  Q1 <- sum(X == 1)
  Q2 <- sum(X == 2)
  if(Q1 > 0 & Q2 > 0){
    A <- 2*Q2/((t-1)*Q1 + 2*Q2)
  } else if (Q1 > 0 & Q2 == 0){
    A <- 2/((t-1)*(Q1 - 1) + 2)
  } else {
    A <- 1
  }
  C <- 1 - Q1/U*(1-A)
  H_HT <- t/U*(-sum(C*X/t*log(C*X/t)/(1-(1-C*X/t)^t))) + log(U/t)
  return(H_HT)
}
entropy_MEE_Inci_equ <- function(X)
{
  t <- X[1]
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  Q1 <- sum(X == 1)
  Q2 <- sum(X == 2)
  if(Q1 > 0 & Q2 > 0){
    A <- 2*Q2/((t-1)*Q1 + 2*Q2)
  } else if (Q1 > 0 & Q2 == 0){
    A <- 2/((t-1)*(Q1 - 1) + 2)
  } else {
    A <- 1
  }
  
  UE <- sum(X/t*(digamma(t)-digamma(X)))
  if(Q1 > 0 & A!=1){
    B <- Q1/t*(1-A)^(-t+1)*(-log(A)-sum(sapply(1:(t-1), function(k){1/k*(1-A)^k})))
    H_MEE <- t/U*(UE + B) + log(U/t)
  }else{
    H_MEE <- t/U*UE + log(U/t)
  }
  return(H_MEE)
}
Shannon_Inci_index=function(x,boot=50)
{
  x = unlist(x)
  t = x[1]
  MLE=entropy_MLE_Inci_equ(x) 
  MLE_bc=entropy_MLE_bc_Inci_equ(x)
  HT=entropy_HT_Inci_equ(x)
  MEE=entropy_MEE_Inci_equ(x)
  p_hat=EstiBootComm.Sam(x)
  if(boot>1){
    Boot.X = sapply(1:length(p_hat), function(i){
      rbinom(boot,t,p_hat[i])})
    Boot.X = cbind(rep(t,boot), Boot.X)
    temp1=apply(Boot.X,1,entropy_MLE_Inci_equ)
    temp2=apply(Boot.X,1,entropy_MLE_bc_Inci_equ)
    temp4=apply(Boot.X,1,entropy_HT_Inci_equ)
    temp5=apply(Boot.X,1,entropy_MEE_Inci_equ)
    MLE_sd=sd(temp1)
    MLE_bc_sd=sd(temp2)
    HT_sd=sd(temp4)
    MEE_sd=sd(temp5)
    
    MLE_exp_sd=sd(exp(temp1))
    MLE_bc_exp_sd=sd(exp(temp2))
    HT_exp_sd=sd(exp(temp4))  
    MEE_exp_sd=sd(exp(temp5))
  } else {
    MLE_sd=MLE_bc_sd=HT_sd=MEE_sd=MLE_exp_sd=MLE_bc_exp_sd=HT_exp_sd=MEE_exp_sd=NA
  }
  a=matrix(0,8,4)
  a[1,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  a[2,]=c(MLE_bc,MLE_bc_sd,MLE_bc-1.96*MLE_bc_sd,MLE_bc+1.96*MLE_bc_sd)
  a[3,]=c(HT,HT_sd,HT-1.96*HT_sd,HT+1.96*HT_sd)
  a[4,]=c(MEE,MEE_sd,MEE-1.96*MEE_sd,MEE+1.96*MEE_sd)
  a[5,]=c(exp(MLE),MLE_exp_sd,exp(MLE)-1.96*MLE_exp_sd,exp(MLE)+1.96*MLE_exp_sd)
  a[6,]=c(exp(MLE_bc),MLE_bc_exp_sd,exp(MLE_bc)-1.96*MLE_bc_exp_sd,exp(MLE_bc)+1.96*MLE_bc_exp_sd)
  a[7,]=c(exp(HT),HT_exp_sd,exp(HT)-1.96*HT_exp_sd,exp(HT)+1.96*HT_exp_sd)
  a[8,]=c(exp(MEE),MEE_exp_sd,exp(MEE)-1.96*MEE_exp_sd,exp(MEE)+1.96*MEE_exp_sd)
  return(a)
}

simpson_Inci_MVUE_equ=function(Y)
{ 
  t=Y[1] 
  Y=Y[-1] 
  Y=Y[Y>0]
  U=sum(Y)
  a=(sum(Y*(Y-1))/U^2/(1-1/t)) 
}
simpson_Inci_MLE_equ=function(Y)
{
  t=Y[1] 
  Y=Y[-1] 
  Y=Y[Y>0]
  a=(sum(Y^2)/sum(Y)^2)
}
Simpson_Inci_index=function(x,boot=200)
{
  x=x[x>0]
  t = x[1]
  MVUE=simpson_Inci_MVUE_equ(x)
  MLE=simpson_Inci_MLE_equ(x)
  
  p_hat=EstiBootComm.Sam(x)
  if(boot>1){
    #set.seed(1)
    Boot.X = sapply(1:length(p_hat), function(i){
      rbinom(boot,t,p_hat[i])})
    Boot.X = cbind(rep(t,boot), Boot.X)
    temp1=apply(Boot.X,1,simpson_Inci_MVUE_equ)
    temp2=apply(Boot.X,1,simpson_Inci_MLE_equ)
    
    MVUE_sd=sd(temp1)
    MLE_sd=sd(temp2)
    
    #MVUE_recip_sd=MVUE_sd/MVUE
    #MLE_recip_sd=MLE_sd/MLE
    
    MVUE_recip_sd = sd(temp1^(-1))
    MLE_recip_sd = sd(temp2^(-1))
  } else {
    MVUE_sd=NA
    MLE_sd=NA
    
    #MVUE_recip_sd=MVUE_sd/MVUE
    #MLE_recip_sd=MLE_sd/MLE
    
    MVUE_recip_sd = NA
    MLE_recip_sd = NA
  }

  
  a=matrix(0,4,4)
  a[1,]=c(MVUE,MVUE_sd,MVUE-1.96*MVUE_sd,MVUE+1.96*MVUE_sd)
  a[2,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  a[3,]=c(1/MVUE,MVUE_recip_sd,1/MVUE-1.96*MVUE_recip_sd,1/MVUE+1.96*MVUE_recip_sd)
  a[4,]=c(1/MLE,MLE_recip_sd,1/MLE-1.96*MLE_recip_sd,1/MLE+1.96*MLE_recip_sd)
  return(a)
}

###############################################2015.09.15
##from vegan(fisherfit)
alpha=function(x)
{
  Dev.logseries <- function(n.r, p, N) 
  {
    r <- as.numeric(names(n.r))
    x <- N/(N + p)
    logmu <- log(p) + log(x) * r - log(r)
    lhood <- -sum(n.r * (logmu - log(n.r)) + 1) - p * log(1 - x)
    lhood 
  }
  tmp <- as.rad(x)
  N <- sum(x)
  tmp <- tmp/N
  p <- 1/sum(tmp^2)
  n.r <- as.fisher(x)
  a <- nlm(Dev.logseries, n.r = n.r, p = p, N = N, hessian = TRUE)
  mean=a$estimate
  se=sqrt(diag(solve(a$hessian)))
  b=matrix(0,1,4)
  b[1,]=c(mean,se,mean-1.96*se,mean+1.96*se)
  b 
}

Diversity=function(X, datatype=c("abundance","incidence","incidence_raw"), q=NULL, B)
{
  if (B <= 1)
    cat("Warning: When \"nboot\" <2, the bootstrap s.e. and confidence interval can't be calculated.", 
        "\n\n") 
  X=X[,1]
  if(datatype=="abundance"){
    type="abundance"
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    
    BASIC.DATA <- matrix(round(c(sum(X), sum(X>0), 1-sum(X==1)/sum(X), CV.Ind(X)),3), ncol = 1)
    nickname <- matrix(c("n", "D", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("     Sample size", "     Number of observed species",
                              "     Estimated sample coverage",
                              "     Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    
    table0 <- matrix(0,5,4)
    table0[1,]=c(Chao1(X)[-5])
    table0[2,]=c(Chao1_bc(X))
    table0[3,]=round(SpecAbuniChao1(X, k=10, conf=0.95)[1,],1)
    table0[4,]=round(c(SpecAbunAce(X)),1)
    table0[5,]=round(c(SpecAbunAce1(X)),1)
    colnames(table0) <- c("Estimate", "s.e.", paste(Chao1(X)[5]*100,"%Lower", sep=""), paste(Chao1(X)[5]*100,"%Upper", sep=""))
    rownames(table0) <- c("     Chao1 (Chao, 1984)","     Chao1-bc ", "     iChao1","     ACE (Chao & Lee, 1992)",
                          "     ACE-1 (Chao & Lee, 1992)")
    
    SHANNON=Shannon_index(X, B)
    table1=round(SHANNON[c(1:5),],3)
    table1=table1[-2,]              ##2016.05.09
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Jackknife",
    #                      " Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("      MLE","      Jackknife",
                          "      Chao & Shen","      Chao et al. (2013)")
    
    table1_exp=round(SHANNON[c(6:10),],3)
    table1_exp=table1_exp[-2,]      ##2016.05.09
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Jackknife",
    #                         " Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("      MLE","      Jackknife",
                              "      Chao & Shen","      Chao et al. (2013)")
    
    table2=round(Simpson_index(X, B)[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("      MVUE","      MLE")
    
    table2_recip=round(Simpson_index(X)[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("      MVUE","      MLE")
    
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=NULL, from=0, to=3, interval=0.25, B, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=q, from=0, to=3, interval=0.25, B, conf=0.95))}
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2
    
    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- length(q)
    rownames(Hill) <- paste("     ",1:q_hill)
    
    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0, 
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }else if(datatype=="incidence"){
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    type="incidence"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c( T, D, U, C, CV),3), ncol = 1)
    nickname <- matrix(c( "T", "D", "U","C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c( "     Number of sampling units",
                               "     Number of observed species",
                               "     Total number of incidences",
                               "     Estimated sample coverage",
                               "     Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    #BASIC.DATA <- basicInci(X, k=10)[[1]]
    ############################################################
    table0=SpecInci(X, k=10, conf=0.95)
    rownames(table0) <- c("     Chao2 (Chao, 1987)","     Chao2-bc ", "     iChao2","     ICE (Lee & Chao, 1994)",
                          "     ICE-1 (Lee & Chao, 1994)")
    SHANNON=Shannon_Inci_index(X, B)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("      MLE","      Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("      MLE","      Chao et al. (2013)")
    
    SIMPSON=Simpson_Inci_index(X, B)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("      MVUE","      MLE")
    
    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("      MVUE","      MLE")
    
    
    ############################################################
    #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=NULL, from=0, to=3, interval=0.25, B, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=q, from=0, to=3, interval=0.25, B, conf=0.95))}
    
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2
    
    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- length(q)
    rownames(Hill) <- paste("     ",1:q_hill)
    #z <- list("BASIC.DATA"=BASIC.DATA,"HILL.NUMBERS"= Hill)
    
    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0,
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }else if(datatype=="incidence_raw"){
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    Y <- X
    X <- X[-1]
    type="incidence"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c( T, D, U, C, CV),3), ncol = 1)
    nickname <- matrix(c( "T", "D", "U","C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c( "     Number of sampling units",
                               "     Number of observed species",
                               "     Total number of incidences",
                               "     Estimated sample coverage",
                               "     Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    #BASIC.DATA <- basicInci(X, k=10)[[1]]
    ############################################################
    table0=SpecInci_raw(Y, k=10, conf=0.95)
    rownames(table0) <- c("     Chao2 (Chao, 1987)","     Chao2-bc ", "     iChao2","     ICE (Lee & Chao, 1994)",
                          "     ICE-1 (Lee & Chao, 1994)")
    SHANNON=Shannon_Inci_index(X, B)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("      MLE","      Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("      MLE","      Chao et al. (2013)")
    
    SIMPSON=Simpson_Inci_index(X,B)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("      MVUE","      MLE")
    
    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("      MVUE","      MLE")
    
    
    ############################################################
    #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=NULL, from=0, to=3, interval=0.25, B, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=q, from=0, to=3, interval=0.25, B, conf=0.95))}
    
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2
    
    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- length(q)
    rownames(Hill) <- paste("     ",1:q_hill)
    #z <- list("BASIC.DATA"=BASIC.DATA,"HILL.NUMBERS"= Hill)
    
    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0,
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }
  class(z) <- c("spadeDiv")
  return(z) 
}

#############################old##########################################
# Diversity=function(X, datatype=c("Abundance","Frequencies_of_Frequencies"), q=NULL)
# {
#   X=X[,1]
#   if(datatype=="Frequencies_of_Frequencies")
#   {
#     length(X)
#     seed=(1:length(X) %% 2 )
#     i=which(seed==1)
#     fi=which(seed==0)
#     i=X[i]
#     fi=X[fi]
#     X=rep(i,fi)
#   }
#   
#   #BASIC.DATA <- matrix(paste(c("n", "D", "C","CV"),
#   #                              c(sum(X),sum(X>0),round(c(1-sum(X==1)/sum(X),CV.Ind(X)),3)),
#   #                              sep = "="), ncol = 1)
#   #colnames(BASIC.DATA) <- c("Value")
#   #rownames(BASIC.DATA) <- c(" (Number of observed individuals)", " (Number of observed species)",
#   #                          " (Estimated sample coverage)"," (Estimated CV)")
#   
#   #############################################################2015.09.29-(S.W.Wei)
# #   BASIC.DATA <- matrix(round(c(sum(X), sum(X>0), 1-sum(X==1)/sum(X), CV.Ind(X)),3), ncol = 1)
# #   nickname <- matrix(c("n", "D", "C", "CV"), ncol = 1)
# #   BASIC.DATA <- cbind(nickname, BASIC.DATA)
# #   
# #   colnames(BASIC.DATA) <- c("Variable", "Value")
# #   rownames(BASIC.DATA) <- c("Number of observed individuals", "Number of observed species",
# #                             "Estimated sample coverage",
# #                             "Estimated CV")
# #   BASIC.DATA <- data.frame(BASIC.DATA)
#   BASIC.DATA <- basicAbun(X, k=10)[[1]][1:4,]
#   rownames(BASIC.DATA) <- c("Number of observed individuals", "Number of observed species",
#                             "Estimated sample coverage", "Estimated CV")
# 
#   #############################################################2015.09.29
#   
#   table0 <- matrix(0,5,4)
#   table0[1,]=c(Chao1(X)[-5])
#   table0[2,]=c(Chao1_bc(X))
#   table0[3,]=round(c(SpecAbuniChao1(X, k = 10, conf = 0.95)), 1)
#   table0[4,]=round(c(SpecAbunAce(X)),1)
#   table0[5,]=round(c(SpecAbunAce1(X)),1)
#   colnames(table0) <- c("Estimate", "   s.e.", paste(Chao1(X)[5]*100,"%Lower", sep=""), paste(Chao1(X)[5]*100,"%Upper", sep=""))
#   rownames(table0) <- c(" Chao1 (Chao, 1984)"," Chao1-bc "," iChao1"," ACE (Chao & Lee, 1992)",
#                         " ACE-1 (Chao & Lee, 1992)")
#   
#   SHANNON=Shannon_index(X)
#   table1=round(SHANNON[c(1:5),],3)
#   colnames(table1) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
#   rownames(table1) <- c(" MLE"," MLE_bc"," Jackknife",
#                         " Chao & Shen"," Chao (2013)")
#   
#   table1_exp=round(SHANNON[c(6:10),],3)
#   colnames(table1_exp) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
#   rownames(table1_exp) <- c(" MLE"," MLE_bc"," Jackknife",
#                             " Chao & Shen"," Chao (2013)")
#   
#   table2=round(Simpson_index(X)[c(1:2),],5)
#   colnames(table2) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
#   rownames(table2) <- c(" MVUE"," MLE")
#   
#   table2_recip=round(Simpson_index(X)[c(3:4),],5)
#   colnames(table2_recip) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
#   rownames(table2_recip) <- c(" MVUE"," MLE")
#   
#   ###############################################2015.09.22(H.W.Hsu, Y.R.Chen, S.W.Wei)
#   if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
#   if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}
#   #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
#   #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
#   #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
#   #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
#   #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
#   #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
#   #Hill<-round(Hill,3)
#   #Hill <- data.frame(Hill)
#   q_length<-length(Hill[,1])/2
#   
#   Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
#   Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
#   Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
#   Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
#   #Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Hill[1:q_length,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
#   Hill <- cbind(Hill[1:q_length,1], Hill[(q_length+1):(2*q_length),3], Chao.LCL, Chao.UCL, Hill[1:q_length,3], Emperical.LCL, Emperical.UCL)
#   Hill <- round(Hill,3)
#   Hill <- data.frame(Hill)
#   #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
#   #colnames(Hill)<-c("q","Chao","Empirical","Chao(95% Lower)","Chao(95% Upper)","Empirical(95% Lower)","Empirical(95% Upper)")
#   colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
#   ###############################################2015.09.22
#   
#   z <- list("BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0, 
#             "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
#             "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
#             "HILL.NUMBERS"= Hill)
#   class(z) <- c("spadeDiv")
#   return(z)
#   
#   #cat("\n")
#   #cat("(5)  FISHER ALPHA INDEX:\n\n")
#   #table_alpha=round(alpha(X),3)
#   #colnames(table_alpha)<-c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
#   #rownames(table_alpha)<-c(" alpha")
#   #print( table_alpha)
#   #cat("\n")
#   #cat(" See Eq. (2.9) of Magurran (1988) for a definition of Fisher's alpha index.\n")
# }
#############################old##########################################
#X=read.table("Data4a.txt")
#Y=read.table("Data4b1_t.txt")
#Diversity(datatype="Abundance",X)
#Diversity(datatype="Frequencies_of_Frequencies",Y)

#############################old##########################################
# print.spadeDiv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
#   
#   cat("\n(1)  BASIC DATA INFORMATION:\n")
#   print(x$BASIC.DATA)
#   cat("\n(2)  ESTIMATION OF SPECIES RICHNESS (DIVERSITY OF ORDER 0):\n\n")
#   print(x$SPECIES.RICHNESS)
#   cat("
#         Desriptions (See Species Part)
#         
#         Chao1 (Chao, 1984): This approach uses the numbers of singletons and doubletons to
#         estimate the number of undetected species because undetected species information is
#         mostly concentrated on those low frequency counts; see Chao (1984), and Chao and Chiu (2012).
#         
#         Chao1-bc: A bias-corrected form for the Chao1 estimator; see Chao (2005).
#       
#         iChao1: An improved Chao1 estimator; See Chiu et al. (2014).
#         
#         ACE (Abundance-based Coverage Estimator): A non-parametric estimator proposed by Chao and Lee (1992)
#         and Chao, Ma and Yang (1993). The observed species are separated as rare and abundant groups;
#         only the rare group is used to estimate the number of undetected species.
#         The estimated CV is used to characterize the degree of heterogeneity among species
#         discovery probabilities. See Eq. (2.14) in Chao and Lee (1992) or Eq. (2.2) of Chao et al. (2000).
#         
#         ACE-1: A modified ACE for highly heterogeneous communities. See Eq. (2.15) of Chao and Lee (1992).
#         
#         95% Confidence interval: A log-transformation is used for all estimators so that the lower bound 
#         of the resulting interval is at least the number of observed species. See Chao (1987).
#         ")
#   cat("\n(3a)  SHANNON INDEX:\n\n")
#   print(x$SHANNON.INDEX)
#   cat("\n")
#   cat(" For a review of the four estimators, see Chao and Shen (2003).\n")
#   cat("
#       MLE: empirical or observed entropy.
#       MLE_bc: bias-corrected empirical estimator.
#       Jackknife: see Zahl (1977).
#       Chao & Shen: based on Horvitz-Thompson estimator and sample coverage method;
#       see Chao and Shen (2003). 
#       Chao et al. (2013): A low-bias estimator of entropy; see Chao et al. (2013).
#       Estimated standard error is based on a bootstrap method.
#       \n")
#   
#   cat("(3b)  EXPONENTIAL OF SHANNON INDEX (DIVERSITY OF ORDER 1):\n\n")
#   print(x$EXPONENTIAL.OF.SHANNON.INDEX)
#   
#   cat("\n(4a)  SIMPSON INDEX:\n\n")
#   print(x$SIMPSON.INDEX)
#   
#   cat("
#       MVUE: minimum variance unbiased estimator; see Eq. (2.27) of Magurran (1988).
#       MLE: maximum likelihood estimator or empirical index; see Eq. (2.26) of Magurran (1988).
#       ")
#   
#   cat("\n(4b)  INVERSE OF SIMPSON INDEX (DIVERSITY OF ORDER 2):\n\n")
#   print(x$INVERSE.OF.SIMPSON.INDEX)
#   
#   cat("\n(5)  Chao and Jost (2015) estimates of Hill numbers of order q from 0 to 3\n\n")
#   print(x$HILL.NUMBERS)
#   
#   cat("
#       ChaoJost: see Chao and Jost (2015).
#       Empirical: maximum likelihood estimator (observed index).
#       ")
#   
# }
#############################old##########################################

print.spadeDiv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  

  
  if(x$datatype=="abundance"){
    
    cat("\n(1)  BASIC DATA INFORMATION:\n")
    print(x$BASIC.DATA)
    cat("\n(2)  ESTIMATION OF SPECIES RICHNESS (DIVERSITY OF ORDER 0):\n\n")
    print(x$SPECIES.RICHNESS)
    ################################################################2016.07.06-(P.L.Lin)
    cat("
        Descriptions of richness estimators (See Species Part)
        ")
    ################################################################
    cat("\n(3a)  SHANNON ENTROPY:\n\n")
    print(x$SHANNON.INDEX)
    #cat("\n")
    #cat(" For a review of the four estimators, see Chao and Shen (2003).\n")
    #MLE_bc: bias-corrected empirical estimator.
    cat("
        MLE: empirical or observed entropy.
        Jackknife: see Zahl (1977).
        Chao & Shen: based on the Horvitz-Thompson estimator and sample coverage method; see Chao and Shen (2003). 
        Chao et al. (2013): A nearly optimal estimator of Shannon entropy; see Chao et al. (2013).
        Estimated standard error is computed based on a bootstrap method.
        \n")
    
    cat("(3b)  SHANNON DIVERSITY (EXPONENTIAL OF SHANNON ENTROPY):\n\n")
    print(x$EXPONENTIAL.OF.SHANNON.INDEX)
    
    cat("\n(4a)  SIMPSON CONCENTRATION INDEX:\n\n")
    print(x$SIMPSON.INDEX)
    
    cat("
        MVUE: minimum variance unbiased estimator; see Eq. (2.27) of Magurran (1988).
        MLE: maximum likelihood estimator or empirical index; see Eq. (2.26) of Magurran (1988).
        ")
    
    cat("\n(4b)  SIMPSON DIVERSITY (INVERSE OF SIMPSON CONCENTRATION):\n\n")
    print(x$INVERSE.OF.SIMPSON.INDEX)
    
    cat("\n(5)  CHAO AND JOST (2015) ESTIMATES OF HILL NUMBERS\n\n")
    print(x$HILL.NUMBERS)
    
    cat("
        ChaoJost: diversity profile estimator derived by Chao and Jost (2015).
        Empirical: maximum likelihood estimator (observed index).
        ")
  }else{
    cat("\n(1)  BASIC DATA INFORMATION:\n")
    print(x$BASIC.DATA)
    cat("\n(2)  ESTIMATION OF SPECIES RICHNESS (DIVERSITY OF ORDER 0):\n\n")
    print(x$SPECIES.RICHNESS)
    cat("
        Descriptions (see Species Part)
        ")
    cat("\n(3a)  SHANNON ENTROPY:\n\n")
    print(x$SHANNON.INDEX)
    
    cat("\n(3b)  SHANNON DIVERSITY (EXPONENTIAL OF SHANNON ENTROPY):\n\n")
    print(x$EXPONENTIAL.OF.SHANNON.INDEX)
    cat("\n(4a)  SIMPSON  CONCENTRATION INDEX:\n\n")
    print(x$SIMPSON.INDEX)
    cat("\n(4b)  SIMPSON DIVERSITY (INVERSE OF SIMPSON CONCENTRATION):\n\n")
    print(x$INVERSE.OF.SIMPSON.INDEX)
    cat("\n(5)  CHAO AND JOST (2015) ESTIMATES OF HILL NUMBERS\n\n")
    print(x$HILL.NUMBERS)
    
    cat("
        ChaoJost: diversity profile estimator derived by Chao and Jost (2015).
        Empirical: maximum likelihood estimator (observed index).
        ")
    
  }
  Lower=min(x$HILL.NUMBERS[,3],x$HILL.NUMBERS[,6])
  Upper=max(x$HILL.NUMBERS[,4],x$HILL.NUMBERS[,7])
  plot(0,type="n",xlim=c(min(x$HILL.NUMBERS[,1]),max(x$HILL.NUMBERS[,1])),xlab="Order  q",ylab="Hill  numbers")
  conf.reg(x$HILL.NUMBERS[,1],x$HILL.NUMBERS[,3],x$HILL.NUMBERS[,4], col=adjustcolor(2, 0.2), border=NA)
  conf.reg(x$HILL.NUMBERS[,1],x$HILL.NUMBERS[,6],x$HILL.NUMBERS[,7], col=adjustcolor(4, 0.2), border=NA)
  lines(x$HILL.NUMBERS[,1],x$HILL.NUMBERS[,2],col=2,lwd=3)
  lines(x$HILL.NUMBERS[,1],x$HILL.NUMBERS[,5],col=4,lty=3,lwd=3)
  legend("topright", c("ChaoJost","Empirical"),col=c(2,4),lwd=c(3,3),lty=c(1,3),bty="n",cex=0.8) 
}


########################################
Diversity_Inc=function(X, q=NULL)
{ X=X[,1]
  ###############################################2015.09.15(H.W.Hsu, Y.R.Chen and S.W.Wei)
  BASIC.DATA <- basicInci(X, k=10)[[1]][c(1,2,3,4,5),]
  rownames(BASIC.DATA) <- c("Number of observed species", "Number of sample/quadrats","Total number of incidences",
                            "Estimated sample coverage", "Estimated CV")
  table0=SpecInci(X, k=10, conf=0.95)
  SHANNON=Shannon_Inci_index(X,B)
  table1=round(SHANNON[c(1,4),],3)
  colnames(table1) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
  rownames(table1) <- c(" MLE"," Chao (2013)")
  
  table1_exp=round(SHANNON[c(5,8),],3)
  colnames(table1_exp) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
  rownames(table1_exp) <- c(" MLE"," Chao (2013)")
  
  SIMPSON=Simpson_Inci_index(X,B)
  table2=round(SIMPSON[c(1:2),],5)
  colnames(table2) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
  rownames(table2) <- c(" MVUE"," MLE")
  
  table2_recip=round(SIMPSON[c(3:4),],5)
  colnames(table2_recip) <- c("Estimate", "   s.e.", paste("95%Lower"), paste("95%Upper"))
  rownames(table2_recip) <- c(" MVUE"," MLE")
  ###############################################2015.09.15
  
  ###############################################2015.09.22(H.W.Hsu, Y.R.Chen, S.W.Wei)
  #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
  if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=NULL, from=0, to=3, interval=0.25, B, conf=0.95))}
  if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=q, from=0, to=3, interval=0.25, B, conf=0.95))}
  
  #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
  #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
  #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
  #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
  #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
  #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
  #Hill<-round(Hill,3)
  #Hill <- data.frame(Hill)
  q_length<-length(Hill[,1])/2
  
  Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
  Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
  Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
  Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
  #Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Hill[1:q_length,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
  Hill <- cbind(Hill[1:q_length,1], Hill[(q_length+1):(2*q_length),3], Chao.LCL, Chao.UCL, Hill[1:q_length,3], Emperical.LCL, Emperical.UCL)
  Hill<-round(Hill,3)
  Hill <- data.frame(Hill)
  #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
  #colnames(Hill)<-c("q","Chao","Empirical","Chao(95% Lower)","Chao(95% Upper)","Empirical(95% Lower)","Empirical(95% Upper)")
  colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
  ###############################################2015.09.22
  
  ###############################################2015.09.15(H.W.Hsu, Y.R.Chen and S.W.Wei)
  z <- list("BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0,
            "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
            "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
            "HILL.NUMBERS"= Hill)
  
  
  ###############################################2015.09.15
  class(z) <- c("spadeDiv_Inc")
  return(z)
  
  #cat("\n")
  #cat("(5)  FISHER ALPHA INDEX:\n\n")
  #table_alpha=round(alpha(X),3)
  #colnames(table_alpha)<-c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
  #rownames(table_alpha)<-c(" alpha")
  #print( table_alpha)
  #cat("\n")
  #cat(" See Eq. (2.9) of Magurran (1988) for a definition of Fisher's alpha index.\n")
}
#X=read.table("Data4a.txt")
#Y=read.table("Data4b1_t.txt")
#Diversity(datatype="Abundance",X)
#Diversity(datatype="Frequencies_of_Frequencies",Y)

print.spadeDiv_Inc <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  ###############################################2015.09.15(H.W.Hsu, Y.R.Chen and S.W.Wei)
  cat("\n(1)  BASIC DATA INFORMATION:\n")
  print(x$BASIC.DATA)
  cat("\n(2)  ESTIMATION OF SPECIES RICHNESS (DIVERSITY OF ORDER 0):\n\n")
  print(x$SPECIES.RICHNESS)
  cat("
      Descriptions (see Species Part)
      ")
  cat("\n(3a)  SHANNON INDEX:\n\n")
  print(x$SHANNON.INDEX)
  cat("\n(3b)  SHANNON DIVERSITY (EXPONENTIAL OF SHANNON INDEX):\n\n")
  print(x$EXPONENTIAL.OF.SHANNON.INDEX)
  cat("\n(4a)  SIMPSON INDEX:\n\n")
  print(x$SIMPSON.INDEX)
  cat("\n(4b)  SIMPSON DIVERSITY (INVERSE OF SIMPSON):\n\n")
  print(x$INVERSE.OF.SIMPSON.INDEX)
  ###############################################2015.09.15
  cat("\n(5)  Chao and Jost (2015) estimates of Hill numbers of order q from 0 to 3\n\n")
  print(x$HILL.NUMBERS)
  
  cat("
      ChaoJost: diversity profile estimator derived by Chao and Jost (2015).
      Empirical: maximum likelihood estimator (observed index).
      ")
  
}

