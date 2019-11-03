print.BasicFun <-
function(x) {
  #   cat("(1)  BASIC DATA INFORMATION:", "\n\n")
  cat("    Sample size  in Community 1                     n1  = ", x$n1, "\n")
  cat("    Sample size  in Community 2                     n2  = ", x$n2, "\n")
  cat("    Number of observed species in Community 1       D1  = ", x$D1, "\n")
  cat("    Number of observed species in Community 2       D2  = ", x$D2, "\n")
  cat("    Number of observed shared species               D12 = ", x$D12, "\n")
  cat("    Bootstrap replications for s.e. estimate              ", x$B,   "\n\n")
  cat("     \"Entire\" Shared Species Group:", "\n")
  cat("         Some Statistics:", "\n")
  cat("         --------------------------------------------------------------------------", "\n")
  cat("         f[11] =", x$f11, ";", "f[1+] =", x$f1.plus, ";", "f[+1] =", x$fplus.1, "; ", "f[2+] =", x$f2.plus, "; ", "f[+2] =", x$fplus.2, "; ", "f[22] =", x$f22,"\n")
  cat("         --------------------------------------------------------------------------", "\n\n")
  cat("     \"Rare\" Shared Species Group: (Both frequencies can only up to 10)", "\n")
  #####################################################################2016.07.06  J.H.Lin##########
  cat("         Some Statistics:", "\n")
  cat("         -------------------------------------------------------------------", "\n")
  cat("         f[1+]_rare =", x$f1.plus_rare, ";", "f[+1]_rare =", x$fplus.1_rare, "; ", "f[2+]_rare =", x$f2.plus_rare, "; ", "f[+2]_rare =", x$fplus.2_rare,"\n")
  cat("         -------------------------------------------------------------------", "\n")
  cat("    Number of observed individuals in Community 1     n1_rare   = ", x$n1_rare, "\n")
  cat("    Number of observed individuals in Community 2     n2_rare   = ", x$n2_rare, "\n")
  cat("    Number of observed shared species                 D12_rare  = ", x$D12_rare, "\n")
  cat("    Estimated sample coverage                         C12_rare  = ", x$C12_rare, "\n")
  cat("    Estimated CCVs                                    CCV_1     = ", x$CCV_1, "\n")
  cat("                                                      CCV_2     = ", x$CCV_2, "\n")
  cat("                                                      CCV_12    = ", x$CCV_12, "\n")
}

