BasicFun.Sam <-
function(y1, y2, B) {
  t1 <- y1[1]
  t2 <- y2[1]
  x1 <- y1[-1]
  x2 <- y2[-1]
  D1 <- sum(x1 > 0)
  D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)  
  U1 <- sum(x1)
  U2 <- sum(x2)
  Q11 <- sum(x1 == 1 & x2 == 1)
  Q1.plus <- sum(x1 == 1 & x2 >= 1)
  Qplus.1 <- sum(x2 == 1 & x1 >= 1)
  Q22 <- sum(x1 == 2 & x2 == 2)
  Q2.plus <- sum(x1 == 2 & x2 >= 1)
  Qplus.2 <- sum(x2 == 2 & x1 >= 1)
  
  out <- list(t1=t1, t2=t2, D1=D1, D2=D2, D12=D12, B=B, Q11=Q11, Q1.plus=Q1.plus, U1=U1, U2=U2,
              Q1.plus=Q1.plus, Qplus.1=Qplus.1, Q2.plus=Q2.plus, Qplus.2=Qplus.2, Q22=Q22)
  class(out) <- "BasicFun.Sam"
  return(out)
}
