SpecInciModelth1 <-
function(data, k, conf){
  data <- as.numeric(data)
  z <- -qnorm((1 - conf)/2)
  t <- length(data[1, ])
  dat <- apply(data, 1, sum) 
  x <- dat[which(dat != 0)]
  Q <- function(i, data){length(data[which(data == i)])}
  
  basicInci <- function(data, k){
    data <- as.numeric(data)
    t <- data[1]
    dat <- data[-1]
    x <- dat[which(dat != 0)]
    Q <- function(i, data){length(data[which(data == i)])}
    
    D <- length(x)
    D_infreq <- length(x[which(x <= k)])
    
    if (Q(1, x) > 0 & Q(2, x) > 0){
      A <- 2*Q(2, x)/((t-1)*Q(1, x) + 2*Q(2, x))
    } else if (Q(1, x) > 0 & Q(2, x) == 0){
      A <- 2/((t-1)*(Q(1, x) - 1) + 2)
    } else {
      A <- 1
    }
    C_infreq <- 1 - Q(1, x)/sum(x[which(x <= k)])*(1-A)
    
    j <- c(1:k)
    b1 <- sum(sapply(j, function(j)j*(j-1)*Q(j, x)))
    b2 <- sum(sapply(j, function(j)j*Q(j, x)))
    gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*b1/b2/(b2) - 1, 0)
    CV_infreq <- sqrt(gamma_infreq_square)
    D_freq <- length(x[which(x > k)])
    
    BASIC.DATA <- matrix(paste(c("D", "t", "k", "D_infreq", "C_infreq", "CV_infreq", "D_freq"),
                               c(D,t,k,D_infreq,C_infreq,CV_infreq,D_freq),
                               sep = "="), ncol=1)
    colnames(BASIC.DATA)=c("Value")
    rownames(BASIC.DATA)=c("Number of observed species","Number of sample/quadrats","Cut-off point",
                           "Number of observed species for infrequent species","Estimated sample coverage for infrequent species",
                           "Estimated CV for infrequent species",
                           "Number of observed species for frequent species")
    return(list(BASIC.DATA, D, t, D_infreq, C_infreq, CV_infreq, D_freq))
  }
  D <- basicInci(data, k)[[2]]
  D_infreq <- basicInci(data, k)[[4]]
  C_infreq <- basicInci(data, k)[[5]]
  CV_infreq <- basicInci(data, k)[[6]]
  D_freq <- basicInci(data, k)[[7]]
  
  u <- c(1:k)   
  n_infreq <- sum(x[which(x <= k)])
  si <- sum(sapply(u, function(u)u*(u-1)*Q(u, x)))
  data.f <- data[which(x <= k), ]
  nd <- as.numeric(apply(data.f, 2, sum))
  nd_c <- nd[which(nd != 0)]
  e <- c(1:length(nd_c))
  o <- rep(e ,each = length(e))
  p <- rep(e ,length(e))    
  s <- sum(mapply(function(o, p)nd_c[o]*nd_c[p], o, p)) - sum(sapply(e, function(e)nd_c[e]*nd_c[e]))
  S_Model_TH1 <- function(x, k){
    si <- sum(sapply(u, function(u)u*(u - 1)*Q(u, x)))
    gamma_infreq_square_th <- max(D_infreq/C_infreq*si/s - 1, 0)
    gamma_infreq_square_th1 <- max(gamma_infreq_square_th*(1 + Q(1, x)/C_infreq*si/s), 0)
    s_Model_th1 <- D_freq + D_infreq/C_infreq + Q(1, x)/C_infreq*gamma_infreq_square_th1
    CV_infreq_th1 <- sqrt(gamma_infreq_square_th1)
    return(unlist(list(s_Model_th1, CV_infreq_th1)))
  }   
  s_Model_th1 <- S_Model_TH1(x, k)[1]
  CV_infreq_th1 <- S_Model_TH1(x, k)[2]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    Q1 = ifelse( (Q(1, x)==0) & (Q(2, x)==0), 1 , Q(1, x))
    Q2 = ifelse( (Q(1, x)==0) & (Q(2, x)==0), 1 , Q(2, x))
    if (CV_infreq_th1 != 0){
      if ( q == 1) {
        dc_infreq <-  - (n_infreq*((t - 1)*Q1 + 2*Q2)*2*Q1*(t - 1) - 
                           (t - 1)*Q1^2*((t - 1)*(Q1 + n_infreq) + 2*Q2))/(n_infreq*((t - 1)*Q1 + 2*Q2))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          (C_infreq^2*(D_infreq*si + Q1*si) - #g3
             Q1*D_infreq*si*(2*C_infreq*dc_infreq)             
          )/C_infreq^4/s - 
          (C_infreq - Q1*dc_infreq)/C_infreq^2 + #g4
          (C_infreq^3*(Q1^2*si^2 + D_infreq*2*Q1*si^2) - D_infreq*Q1^2*si^2*(3*C_infreq^2*dc_infreq) #g5
          )/C_infreq^6/s^2 - 
          (C_infreq^2*2*Q1*si - Q1^2*si*2*C_infreq*dc_infreq #g6
          )/C_infreq^4/s
      } else if (q == 2) {
        dc_infreq <-  - ( - (t - 1)*Q1^2*(2*(t - 1)*Q1 + 2*(n_infreq + 2*Q2)))/(n_infreq*((t - 1)*Q1 + 2*Q2))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          (C_infreq^2*Q1*(si + 2*D_infreq) - Q1*D_infreq*si*( #g3
            2*C_infreq*dc_infreq)
          )/C_infreq^4/s - 
          ( - Q1*dc_infreq)/C_infreq^2 + #g4
          (C_infreq^3*Q1^2*(si^2 + D_infreq*2*si*2) - D_infreq*Q1^2*si^2*3*C_infreq^2*dc_infreq #g5
          )/C_infreq^6/s^2 - 
          (C_infreq^2*Q1^2*2 - Q1^2*si*2*C_infreq*dc_infreq #g6
          )/C_infreq^4/s
      } else if(q > k) {
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Q1^2*((t - 1)*Q1*q + 2*Q2*q))/(n_infreq*((t - 1)*Q1 + 2*Q2))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          (C_infreq^2*Q1*(si + q*(q - 1)*D_infreq) - Q1*D_infreq*si*( #g3
            2*C_infreq*dc_infreq)
          )/C_infreq^4/s - 
          ( - Q1*dc_infreq)/C_infreq^2 + #g4
          (C_infreq^3*Q1^2*(si^2 + D_infreq*2*si*q*(q - 1)) - D_infreq*Q1^2*si^2*3*C_infreq^2*dc_infreq #g5
          )/C_infreq^6/s^2 - 
          (C_infreq^2*Q1^2*q*(q - 1) - Q1^2*si*2*C_infreq*dc_infreq
          )/C_infreq^4/s
      }
      return(d)
    } else {
      if ( q == 1) {
        dc_infreq <-  - (n_infreq*((t - 1)*Q1 + 2*Q2)*2*Q1*(t - 1) - 
                           (t - 1)*Q1^2*((t - 1)*(Q1 + n_infreq) + 2*Q2))/(n_infreq*((t - 1)*Q1 + 2*Q2))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      } else if (q == 2) {
        dc_infreq <-  - ( - (t - 1)*Q1^2*(2*(t - 1)*Q1 + 2*(n_infreq + 2*Q2)))/(n_infreq*((t - 1)*Q1 + 2*Q2))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      } else if(q > k) {
        d <- 1
      } else {
        w <- c(q:k)
        ss <- sum(sapply(w, function(w)nd[w]))
        dc_infreq <-  - ( - (t - 1)*Q1^2*((t - 1)*Q1*q + 2*Q2*q))/(n_infreq*((t - 1)*Q1 + 2*Q2))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }
      return(d)
    }
  }
  
  COV.q <- function(i,j){
    if (i == j){
      cov.q <- Q(i, x)*(1 - Q(i, x)/s_Model_th1)
    } else {
      cov.q <- -Q(i, x)*Q(j, x)/s_Model_th1
    }     
    return(cov.q)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_th1 <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.q(i, j), i, j))
  ######################
  if (round(s_Model_th1 - D, 5) != 0){
    C <- exp(z*sqrt(log(1 + var_th1/(s_Model_th1-D)^2)))
    CI_Model_th1 <- c(D + (s_Model_th1 - D)/C, D + (s_Model_th1 - D)*C)
  }else{
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
    P <- sum(sapply(i, function(i)Q(i, x)*exp(-i)/D))
    CI_Model_th1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(s_Model_th1, sqrt(var_th1), CI_Model_th1), ncol = 4)
  colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  rownames(table) <- "Model(th)-1 (Lee & Chao, 1994)"
  return(list(table, CV_infreq_th1))
}
