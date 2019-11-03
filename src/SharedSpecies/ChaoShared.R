ChaoShared <-
function(data, datatype = c("abundance", "incidence"), 
                       se = TRUE, nboot = 200, conf = 0.95) {
  
  method <- "all"
  if (se == TRUE) {
    if (nboot < 1)
      nboot <- 1
    if (nboot == 1)
      cat("Warning: When \"nboot\" <2, the bootstrap s.e. and confidence interval can't be calculated.", 
          "\n\n")  
  }
  se <- ifelse(nboot < 1, F, T)
  
  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
    cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
        "\n")
    cat("          We use \"conf\" = 0.95 to calculate!", 
        "\n\n")
    conf <- 0.95
  }
  
  datatype <- match.arg(datatype)
  if (datatype == "abundance") {
    x1 <- data[, 1]
    x2 <- data[, 2]
    Basic <- BasicFun(x1, x2, nboot)
    #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
    output <- ChaoShared.Ind(x1, x2, method, nboot, conf, se)
	colnames(output) <- c("Estimate", "   s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
  }
  if (datatype == "incidence") {
    y1 <- data[, 1]
    y2 <- data[, 2]
    Basic <- BasicFun.Sam(y1, y2, B=nboot)
    #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
    output <- ChaoShared.Sam(y1, y2, method, conf, se)
	colnames(output) <- c("Estimate", "   s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
  }
  out <- list(BASIC_DATA_INFORMATION=Basic, 
              ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES=output)
  class(out) <- c("sharedspecies")
  return(out)
}

print.sharedspecies <- function(x, ...){
	cat('\n(1) BASIC DATA INFORMATION:\n\n')
	print(x$BASIC_DATA_INFORMATION)
	cat('\n')
	cat('\n(2) ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES:\n\n')
	print(round(x$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES,3))
	cat('\n')
	if(nrow(x[[2]])==4){
	cat('
(3) DESCRIPTION OF MODELS FOR ESTIMATING SHARED SPECIES RICHNESS:
Homogeneous: This model assumes that the shared species in each community have the same discovery probabilities; see Eq. (3.11a) of Chao et al. (2000).

Heterogeneous (ACE-shared): This model allows for heterogeneous discovery probabilities among shared species; see Eq. (3.11b) of Chao et al. (2000). An extension of the ACE estimator to two communities; it is replaced by Chao1-shared when the estimated sample coverage for the rare shared species group (C12_rare in the output) is zero.

Chao1-shared: An extension of the Chao1 estimator to estimate shared species richness between two communities. It provides a lower bound of shared species richness. See Eq. (3.6) of Pan et al. (2009). It is replaced by Chao1-shared-bc for the case f[2+]=0 or f[+2]=0.
   
Chao1-shared-bc: A bias-corrected form of Chao1-shared estimator; see Pan et al. (2009).
	')
	}else{
	cat('
(3) DESCRIPTION OF MODELS FOR ESTIMATING SHARED SPECIES RICHNESS:

Chao2-shared: An extension of the Chao2 estimator to estimate shared species richness between two communities. It provides a lower bound of shared species richness. See Eq. (3.6) of Pan et al. (2009). It is replaced by Chao2-shared-bc for the case Q[2+]=0 or Q[+2]=0.

Chao2-shared-bc: A bias-corrected form of the Chao2-shared. See Pan et al. (2009).
	')}
}

