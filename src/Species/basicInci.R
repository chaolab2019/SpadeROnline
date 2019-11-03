basicInci <-
  function(data, k){
    data <- as.numeric(data)
    t <- data[1]
    dat <- data[-1]
    x <- dat[which(dat != 0)]
    Q <- function(i, data){length(data[which(data == i)])}
    
    D <- length(x)
    D_infreq <- length(x[which(x <= k)])
    U <- sum(x)
    U_infreq <- sum(x[which(x <= k)])
    
    if (Q(1, x) > 0 & Q(2, x) > 0){
      A <- 2*Q(2, x)/((t-1)*Q(1, x) + 2*Q(2, x))
    } else if (Q(1, x) > 1 & Q(2, x) == 0){
      A <- 2/((t-1)*(Q(1, x) - 1) + 2)
    } else {
      A <- 0
    }
    C <- 1 - Q(1, x)/U*(1-A)
    C_infreq <- 1 - Q(1, x)/sum(x[which(x <= k)])*(1-A)
    CV <- CV.Sam(data)
    
    j <- c(1:k)
    b1 <- sum(sapply(j, function(j)j*(j-1)*Q(j, x)))
    b2 <- sum(sapply(j, function(j)j*Q(j, x)))
    gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*b1/b2/(b2 -1) - 1, 0)
    gamma_infreq_square_1 <- max(gamma_infreq_square*(1 + Q(1, x)/C_infreq*t/(t - 1)*b1/b2/(b2 - 1)), 0)
    CV_infreq <- sqrt(gamma_infreq_square)
    CV1_infreq <- sqrt(gamma_infreq_square_1)
    D_freq <- length(x[which(x > k)])
    
    
    #   BASIC.DATA <- matrix(paste(c("D", "t", "k", "D_infreq", "C_infreq", "CV_infreq", "D_freq"),
    #                              round(c(D,t,k,D_infreq,C_infreq,CV_infreq,D_freq), 3),
    #                              sep = "="), ncol = 1)
    BASIC.DATA <- matrix(round(c(t,D,U,C,CV,k,U_infreq,D_infreq,C_infreq,CV_infreq,CV1_infreq,D_freq), 3), ncol = 1)
    nickname <- c("T", "D",  "U", "C", "CV", "k", "U_infreq", "D_infreq", "C_infreq", "CV_infreq","CV1_infreq","D_freq")
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    colnames(BASIC.DATA)=c("Variable", "Value")
    rownames(BASIC.DATA)=c("    Number of sampling units","    Number of observed species",
                           "    Total number of incidences","    Coverage estimate for entire dataset",
                           "    CV for entire dataset","    Cut-off point","    Total number for incidences in infrequent group",
                           "    Number of observed species for infrequent group","    Estimated sample coverage for infrequent group",
                           "    Estimated CV for infrequent group in ICE",
                           "    Estimated CV1 for infrequent group in ICE-1",
                           "    Number of observed species for frequent group")
    BASIC.DATA <- data.frame(BASIC.DATA)
    return(list(BASIC.DATA, D, t, D_infreq, C_infreq, CV_infreq, CV1_infreq ,D_freq))
  }
