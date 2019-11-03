print.BasicFun.Sam <-
function(x) {
  #cat("(1)  BASIC DATA INFORMATION:", "\n\n")
  cat("    Number of sampling units of Community 1               T1  = ", x$t1, "\n")
  cat("    Number of sampling units of Community 2               T2  = ", x$t2, "\n")
  cat("    Number of total incidences in Community 1             U1  = ", x$U1, "\n")
  cat("    Number of total incidences in Community 2             U2  = ", x$U2, "\n")
  cat("    Number of observed species in Community 1             D1  = ", x$D1, "\n")
  cat("    Number of observed species in Community 2             D2  = ", x$D2, "\n")
  cat("    Number of observed shared species in two communities  D12 = ", x$D12, "\n")
  cat("    Bootstrap replications for s.e. estimate                    ", x$B,   "\n\n")
  cat("     Some statistics:", "\n")
  cat("         --------------------------------------------------------------------------", "\n")
  cat("         Q[11] =", x$Q11, "; ", "Q[1+] =", x$Q1.plus, ";", "Q[+1] =", x$Qplus.1, "; ", "Q[2+] =", x$Q2.plus, "; ", "Q[+2] =", x$Qplus.2, "; ", "Q[22] =", x$Q22, "\n")
  cat("         --------------------------------------------------------------------------", "\n")
}
