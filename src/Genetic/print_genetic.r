GST_equ=function(X,method=c("est","mle"))
{ 
  if(method=="est"){
    X=as.matrix(X)
    no.community=length(X[1,])
    n=colSums(X)
    temp1=X%*%matrix(c(1/n),no.community,1)
    temp2=sum(sapply(1:no.community,function(k) sum((X[,k]/n[k])^2) ))
    
    Hs_hat=1-1/no.community*sum(sapply(1:no.community,function(k)  sum(X[,k]*(X[,k]-1)/n[k]/(n[k]-1))))
    Ht_hat=1-1/no.community^2*sum(sapply(1:no.community,function(k)  sum(X[,k]*(X[,k]-1)/n[k]/(n[k]-1))))-1/no.community^2*(sum(temp1^2)-temp2) 
    GST=1-Hs_hat/Ht_hat
  }else{
    X = as.matrix(X)
    N <- length(X[1,])
    n <- colSums(X)
    ps <- sapply(1:N, FUN = function(i){
      (X[,i] /n[i] )^2
    })
    Hs <- 1 - sum(ps)/N
    Ht <- 1 - sum((apply(sapply(1:N, FUN = function(i){
      (X[,i] /n[i] )
    }),MARGIN = 1, FUN = sum)/N)^2)
    GST=1-Hs/Ht
  }
return(GST)
} 
GST_se_equ <- function(X, nboot=50){
  plus_CI <-function(x){
    if(x[1] >= 1) x[1] <- 1
    c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
  }
  boot.gst.mle=rep(0,nboot)
  boot.gst.est=rep(0,nboot)
  gst.mle=GST_equ(X, method="mle")
  gst.est=GST_equ(X, method="est")
  for(i in 1:nboot)
  {
    p <- Boots.pop(X)
    boot.X=sapply(1:dim(X)[2],function(j)   rmultinom(1,sum(X[,j]),p[,j] ) )
    boot.gst.mle[i]=GST_equ(boot.X, method="mle")
    boot.gst.est[i]=GST_equ(boot.X, method="est")
  }
  se_mle=sd(boot.gst.mle)
  se_est=sd(boot.gst.est)
  out1= plus_CI(c(gst.mle, se_mle))
  out2= plus_CI(c(gst.est, se_est))
  out <- rbind(out1, out2)
  return(out)
}

print.spadeGenetic <- function(x, ...)
{
  cat('\n(1) BASIC DATA INFORMATION:\n\n')
  cat('    The loaded set includes abundance (or frequency) data from',x$info[1],'subpopulations\n')
  cat('    and a total of',x$info[2],'distinct alleles are found.\n\n')
  cat('    Sample size in each subpopulation                           n1   =', x$info[3],'\n')
  N <- x$info[1]
  q <- x$q
  for(j in 2:N){
    cat('                                                               ','n')
    cat(j,'  =',x$info[2+j],'\n')   
  }
  cat('\n')
  cat('    Number of observed alleles in one subpopulation             D1   =', x$info[N+3],'\n')
  
  for(j in 2:N){
    cat('                                                               ','D')
    cat(j,'  =',x$info[N+2+j],'\n')
  }
  cat('\n')
  cat('    Number of observed shared alleles in two subpopulations     D12  =', x$info[3+2*N], '\n')
  
  if(N>2){
    k <- 1
    for(i in 1:(N-1)){     
      for(j in (i+1):N){
        if(i==1 & j==2) next
        cat('                                                               ','D')
        cat(i,j,'  = ', x$info[3+2*N+k], '\n', sep="")
        k <- k + 1
      }
    }
  }
  cat('\n')
  if(N==3)
  {
    cat('    Number of shared alleles in three subpopulations            D123 =',rev(x$info)[2],'\n\n') 
  }
  cat('    Number of bootstrap replications for s.e. estimate                ',rev(x$info)[1],'\n\n')
  cat('(2) EMPIRICAL DIS-SIMILARITY INDICES: \n\n')
  cat('                                         Estimate       s.e.       95%Lower     95%Upper\n')
  cat('    (a) Classic richness-based dis-similarity\n\n')
  temp <- apply(as.matrix(x$Empirical_richness), 2, as.numeric)
  cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')
  cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')
  cat('    (b) Measures for comparing alleles relative abundances\n\n')
  temp <- apply(as.matrix(x$Empirical_relative), 2, as.numeric)
  cat('        1-C1');cat(N);cat('=1-U1');cat(N);cat(' (q=1, Horn)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('        1-C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
  cat('        1-U2');cat(N);cat(' (q=2, Regional diff.)      ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
  cat('        Gst                              ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  cat('    (c) Measures for comparing size-weighted alleles relative abundances\n\n')
  temp <- x$Empirical_WtRelative
  cat('        Horn size-weighted(q=1)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('(3) ESTIMATED DIS-SIMILARITY INDICES: \n\n')
  cat('                                         Estimate       s.e.       95%Lower     95%Upper\n')
  cat('    (a) Classic richness-based dis-similarity\n\n')
  temp <- apply(as.matrix(x$estimated_richness), 2, as.numeric)
  if(temp[1,1]>1) {cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",1),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
  if(temp[1,1]<=1){cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
  if(temp[2,1]>1) {cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",1) ,'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
  if(temp[2,1]<=1){cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
  cat('    (b) Measures for comparing alleles relative abundances\n\n')
  temp <- apply(as.matrix(x$estimated_relative), 2, as.numeric)
  cat('        1-C1');cat(N);cat('=1-U1');cat(N);cat(' (q=1, Horn)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('        1-C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
  cat('        1-U2');cat(N);cat(' (q=2, Regional diff.)      ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
  cat('        Gst                              ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  cat('    (c) Measures for comparing size-weighted alleles relative abundances\n\n')
  temp <- x$estimated_WtRelative
  cat('        Horn size-weighted (q=1)         ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  #if( N != 2){
    cat('(4) ESTIMATED PAIRWISE DIS-SIMILARITY:\n\n')
    if(q == 0){
      cat('    -----------------------Measure 1-C02------------------------\n\n')
      cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
      ###################################################CqN_ Equal weight
      Cqn_PC <- x$pairwise$C02
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          cat('    1-C02(')
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','   ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','   ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      cat('\n')
      cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise dis-similarity matrix: \n\n')
      C_SM=x$dissimilarity_matrix$C02
      cat('    1-C02(i,j)')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i==j){cat(round(0,0),'  ')}
          if(i<j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
      cat('    -----------------------Measure 1-U02------------------------\n\n')
      cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
      Cqn_PC <- x$pairwise$U02
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          cat('    1-U02(')
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','   ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','   ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      cat('\n')
      cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise dis-similarity matrix: \n\n')
      C_SM <- x$dissimilarity_matrix$U02
      
      cat('    1-U02(i,j)')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i==j){cat(round(0,0),'  ')}
          if(i<j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
      cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
      cat('\n')
    }
    if( q == 1 ){
      cat('    --------------------Measure 1-C12 (=1-U12)----------------------\n\n')
      cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
      ###################################################CqN_ Equal weight
      Cqn_PC <- x$pairwise$C12
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          cat('    1-C12(')
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','   ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','   ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      
      cat('\n')
      cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise dis-similarity matrix: \n\n')
      C_SM=x$dissimilarity_matrix$C12
      cat('    1-C12(i,j)')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i==j){cat(round(0,0),'  ')}
          if(i<j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
      ###################################################UqN_ Equal weight
      cat('    ------------------Measure Horn size-weighted--------------------\n\n')
      cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
      Cqn_PC <- x$pairwise$Horn
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          if(q==0){cat('    1-U02(')}
          if(q==1){cat('    Horn(')}
          #if(q==1 & x$method=="relative"){cat('    C12(')}
          #f(q==1 & x$method=="absolute"){cat('    C12*(')}
          if(q==2){cat('    1-U22(')}
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','    ',sprintf("%.3f",1)      ,'#       ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','    ',sprintf("%.3f",temp[1]),'        ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      
      cat('\n')
      cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise dis-similarity matrix: \n\n')
      C_SM <- x$dissimilarity_matrix$Horn
      
      cat('    Horn(i,j) ')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i==j){cat(round(0,0),'  ')}
          if(i<j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
      cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
      cat('\n')
    }
    if(q == 2){
      cat('    -----------------------Measure 1-C22------------------------\n\n')
      cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
      ###################################################CqN_ Equal weight
      Cqn_PC <- x$pairwise$C22
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          cat('    1-C22(')
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','   ',sprintf("%.3f",1)      ,'#        ',sprintf("%.3f",temp[2]),'         (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','   ',sprintf("%.3f",temp[1]),'         ',sprintf("%.3f",temp[2]),'         (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      cat('\n')
      cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise dis-similarity matrix: \n\n')
      C_SM=x$dissimilarity_matrix$C22
      cat('    1-C22(i,j)')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i==j){cat(round(0,0),'  ')}
          if(i<j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
      cat('    -----------------------Measure 1-U22------------------------\n\n')
      cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
      Cqn_PC <- x$pairwise$U22
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          cat('    1-U22(')
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','   ',sprintf("%.3f",1)      ,'#        ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','   ',sprintf("%.3f",temp[1]),'         ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      cat('\n')
      cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise dis-similarity matrix: \n\n')
      C_SM <- x$dissimilarity_matrix$U22
      
      cat('    1-U22(i,j)')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i==j){cat(round(0,0),'  ')}
          if(i<j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
      cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
      cat('\n')
    } 
  #}
  
      

  #cat('(2) NEARLY UNBIASED ESTIMATION OF ALLELIC DIFFERENTIATION OR MORISITA DISSIMILARITY IN ',N,'SUBPOPULATIONS:\n\n')
  #cat('    Estimator','     Estimate','      s.e.','      95% Confidence Interval\n\n')
  #if ( q==0 ) temp0n=x$overlap[1,]
  #if ( q==1 ) temp0n=x$overlap[3,]
  #if ( q==2 ) temp0n=x$overlap[4,]
  #if ( temp0n[1]<=1 ){
  #  cat('    1-C')
  #  cat(q)
  #  cat(N,sprintf("         %.3f",1-temp0n[1]),'       ',sprintf("%.3f",temp0n[2]),'        (',
  #      sprintf("%.3f",1-temp0n[4]),',',sprintf("%.3f",1-temp0n[3]),')\n')
  #}
  #if ( temp0n[1]>1 ){
  #  cat('    1-C')
  #  cat(q)
  #  cat(N,sprintf("         %.3f#",0),'      ',sprintf("%.3f",temp0n[2]),'        (',
  #      sprintf("%.3f",1-temp0n[4]),',',sprintf("%.3f",1-temp0n[3]),')\n')
  #}
  # 
  #cat('\n')
  #cat('    1-C')
  #cat(q)
  #cat(N,':This is the genetic diversity measure defined in Jost (2008) for comparing subpopulations based 
  #    on allele shared information between any two subpopulations.\n')
  #cat('    
  #    Confidence Interval: Based on an improved bootstrap percentile method. (recommend for use in the case when 
  #    similarity is close to 0 or 1 ) \n\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('    Pairwise Comparison:\n\n')
  #cat('    Estimator','     Estimate','          s.e.','      95% Confidence Interval\n\n')
  #Cqn_PC <- x$pairwise  
  #no.temp=1
  #share_index = sapply(1:choose(N,2), function(i){x$info[2+2*N+i]})
  #no.temp2=1
  #for(i in 1:(N-1))
  #{
  #  for(j in (i+1):N)
  #  {
  #    temp=Cqn_PC[no.temp,]
  #    cat('    1-C')
  #    cat(q)
  #    cat('2(')
  #    cat(i)
  #    cat(',')
  #    cat(j)
  #    if ( share_index[no.temp2]!=0 )
  #    {
  #      if ( temp[1]>1 ){
  #        cat(')','        ',sprintf("%.3f #",0),'     ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",1-temp[4]),',',sprintf("%.3f",1-temp[3]),')\n')
  #      }
  #      if ( temp[1]<=1 ){
  #        cat(')','        ',sprintf("%.3f",1-temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",1-temp[4]),',',sprintf("%.3f",1-temp[3]),')\n')
  #      }
  #    }
  #    else
  #    {
  #      cat(')','        ',sprintf("%.3f",1),'       ',sprintf("%.3f",0),'        (',
  #          sprintf("%.3f",1),',',sprintf("%.3f",1),')\n')
  #      #         cat(')','        ',sprintf("%.3f##",1),'     ',sprintf("%.3f##",0),'      (',
  #      #             sprintf("%.3f",1),',',sprintf("%.3f",1),')##\n')
  #    }
  #    no.temp=no.temp+1
  #    no.temp2=no.temp2+1
  #  }
  #}
  #cat('\n')
  #Cqn_PC <- x$pairwise
  #cat('    Average Pairwise =',sprintf("%.3f",1-mean(Cqn_PC[,1])),'\n')
  #cat('    ## There are no shared species, thus estimated similarity is zero and should be used for caution.\n\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('    1-C')
  #cat(q)
  #cat('2: This is the genetic diversity measure D defined in Jost (2008) for comparing 2 subpopulations.')
  #cat('\n\n')
  #cat('    Dissimilarity Matrix: \n\n')
  #C_SM=x$similarity.matrix
  #cat('    1-C')
  #cat(q)
  #cat('2(i,j)\t')
  #for(i in 1:N)
  #{
  #  cat(i,'\t')
  #}
  #cat('\n')
  #for(i in 1:N)
  #{
  #  cat('       ',i,'\t')
  #  for(j in 1:N)
  #  {
  #    if(i>j){cat('\t')}
  #    if(i<=j){
  #      if (1-C_SM[i,j]>=0) cat(sprintf("%.3f",abs(1-C_SM[i,j])),'\t')
  #      else cat(sprintf("%.3f#",0),'\t')
  #    }
  #  }
  #  cat('\n')
  #}
  #cat('\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #  
  #cat('(3)  NEARLY UNBIASED ESTIMATION OF MORISITA SIMILARITY IN ',N,'SUBPOPULATIONS:\n\n')
  #cat('    Estimator','     Estimate','      s.e.','      95% Confidence Interval\n\n')
  #  
  #if ( q==0 ) temp2n=x$overlap[1,]
  #if ( q==1 ) temp2n=x$overlap[3,]
  #if ( q==2 ) temp2n=x$overlap[4,]
  #if ( temp2n[1]>=1 ){
  #  cat('    C')
  #  cat(q)
  #  cat(N,sprintf("           %.3f#",1),'      ',sprintf("%.3f",temp2n[2]),'       (',
  #      sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
  #}
  #if ( temp2n[1]<1 ){
  #  cat('    C')
  #  cat(q)
  #  cat(N,sprintf("           %.3f",temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'       (',
  #      sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
  #}
  #cat('\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('\n')
  #cat('    C')
  #cat(q)
  #cat(N,': A similarity measure of comparing 3 subpopulations based on allele shared information between any two subpopulations.\n')
  #cat('\n')
  #cat('    Pairwise Comparison:\n\n')
  #cat('    Estimator','     Estimate','      s.e.','      95% Confidence Interval\n\n')
  #Cqn_PC <- x$pairwise
  #no.temp=1
  #share_index = sapply(1:choose(N,2), function(i){x$info[2+2*N+i]})
  #no.temp2=1
  #for(i in 1:(N-1))
  #{
  #  for(j in (i+1):N)
  #  {
  #    temp=Cqn_PC[no.temp,]
  #    if(q==0){cat('    C02(')}
  #    if(q==1){cat('    C12(')}
  #    if(q==2){cat('    C22(')}
  #    cat(i)
  #    cat(',')
  #    cat(j)
  #    if ( share_index[no.temp2]!=0 ){
  #      if ( temp[1]<=1 ){
  #        cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
  #      }
  #      else {
  #        cat(')','     ',sprintf("%.3f #", 1),'     ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
  #      }
  #    }
  #    else
  #    {
  #      cat(')','     ',sprintf("%.3f",0),'       ',sprintf("%.3f",0),'        (',
  #          sprintf("%.3f",0),',',sprintf("%.3f",0),')\n')
  #    }
  #    # Cqn_PC[i,1]=0
  #    no.temp=no.temp+1
  #    no.temp2=no.temp2+1
  #  }
  #}
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('\n')
  #Cqn_PC <- x$pairwise
  #cat('    Average Pairwise =',sprintf("%.3f",mean(Cqn_PC[1:choose(N,2),1])),'\n')
  #cat('\n')
  #cat('    ## There are no shared species, thus estimated similarity is zero and should be used for caution.\n\n')
  #
  #cat('    Similarity Matrix: \n\n')
  #C_SM=x$similarity.matrix
  #if(q==0){cat('    C02(i,j)   \t')}
  #if(q==1){cat('    C12(i,j)   \t')}
  #if(q==2){cat('    C22(i,j)   \t')}
  #for(i in 1:N)
  #{
  #  cat(i,'\t')
  #}
  #cat('\n')
  #for(i in 1:N)
  #{
  #  cat('       ',i,'\t')
  #  for(j in 1:N)
  #  {
  #    if(i>j){cat('\t')}
  #    if(i<=j){
  #      if (C_SM[i,j]<=1) cat(sprintf("%.3f",abs(C_SM[i,j])),'\t')
  #      else cat(sprintf("%.3f#",1),'\t')
  #    }
  #  }
  #  cat('\n')
  #}
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('\n')
  #cat('    References:\n')
  #cat('    Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-
  #    stage probabilistic approach to multiple-community similarity indices. 
  #    Biometrics, 64, 1178-1186.\n')
  #cat('    Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular
  #    Ecology, 17, 4015-4026.')

}
Genetic=function(X,q=2,nboot=200)
{ 
  if (nboot <= 1)
    cat("Warning: When \"nboot\" <2, the bootstrap s.e. and confidence interval can't be calculated.", 
        "\n\n") 
  nboot <- ifelse(nboot<=1, 1, nboot)
  type <- "abundance"
  N <- no.community <- ncol(X)
  temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
  n <- apply(X,2,sum)
  D <- apply(X,2,function(x)sum(x>0))
  
  if(N >= 2){
    temp1 <- temp2 <- rep(0, N*(N-1)/2)
    k <- 1
    for(i in 1:(N-1)){     
      for(j in (i+1):N){
        temp1[k] <- paste('D',i,j,sep="")
        temp2[k] <- sum(X[,i]>0 & X[,j]>0)
        k <- k + 1
      }
    }
  }
  names(temp2) <- temp1
  names(n) <- paste('n',1:N, sep="")
  names(D) <- paste('D',1:N, sep="")
  info <- c(temp, n, D, temp2)
  if(N == 3) info <- c(temp, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
  info <- c(info, nboot=nboot)
  ################################################################2016.07.11-(P.L.Lin)
  temp <- list()
  n <- apply(X = X, MARGIN = 2, FUN = sum)
  weight <- n/sum(n)
  weight <- - sum(weight*log(weight)) / log(N)
  plus_CI <-function(x){
    if(x[1] >= 1) x[1] <- 1
    if(x[1] <= 0) x[1] <- 0
    c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
  }
  mat2 <- GST_se_equ(X,nboot)
  MLE_ew_Gst <- mat2[1, ]
  Est_ew_Gst <- mat2[2, ]
  mat <- SimilarityMul(X,0,nboot,method="unequal weight")
  MLE_Jaccard <- plus_CI(c(1-mat$UqN[1, 1],mat$UqN[1, 2]))
  Est_Jaccard <- plus_CI(c(1-mat$UqN[2, 1],mat$UqN[2, 2]))
  MLE_Sorensen <- plus_CI(c(1-mat$CqN[1, 1],mat$CqN[1, 2]))
  Est_Sorensen <- plus_CI(c(1-mat$CqN[2, 1],mat$CqN[2, 2]))
  mat3 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("unequal"))
  MLE_Ee_Horn <- mat3$mle
  MLE_Ee_Horn <- plus_CI(c(1-MLE_Ee_Horn[1],MLE_Ee_Horn[2])) 
  Est_Ee_Horn <- mat3$est
  Est_Ee_Horn <- plus_CI(c(1-Est_Ee_Horn[1],Est_Ee_Horn[2])) 
  mat4 <- SimilarityMul(X,2,nboot,method="equal weight")
  mat5 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("equal"))
  MLE_ew_Horn <- mat5$mle
  Est_ew_Horn <- mat5$est
  MLE_ew_Horn <- plus_CI(c(1-MLE_ew_Horn[1],MLE_ew_Horn[2]))
  Est_ew_Horn <- plus_CI(c(1-Est_ew_Horn[1],Est_ew_Horn[2]))
  MLE_ew_C22 <- plus_CI(c(1-mat4$CqN[1, 1],mat4$CqN[1, 2]))
  Est_ew_C22 <- plus_CI(c(1-mat4$CqN[2, 1],mat4$CqN[2, 2]))
  MLE_ew_U22 <- plus_CI(c(1-mat4$UqN[1, 1],mat4$UqN[1, 2]))
  Est_ew_U22 <- plus_CI(c(1-mat4$UqN[2, 1],mat4$UqN[2, 2]))
  temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
  rownames(temp[[1]]) <- c("1-C0N(q=0,Sorensen)","1-U0N(q=0,Jaccard)") 
  temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22,MLE_ew_Gst)
  rownames(temp[[2]]) <- c("1-C1N=1-U1N(q=1,Horn)","1-C2N(q=2,Morisita)","1-U2N(q=2,Regional overlap)","Gst")  
  temp[[3]] <- t(as.matrix(MLE_Ee_Horn))
  rownames(temp[[3]]) <- c("Horn size weighted(q=1)")  
  temp[[4]] <- rbind(Est_Sorensen, Est_Jaccard)
  rownames(temp[[4]]) <- c("1-C0N(q=0,Sorensen)","1-U0N(q=0,Jaccard)") 
  temp[[5]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22, Est_ew_Gst)
  rownames(temp[[5]]) <- c("1-C1N=1-U1N(q=1,Horn)","1-C2N(q=2,Morisita)","1-U2N(q=2,Regional overlap)","Gst")  
  temp[[6]] <- t(as.matrix(Est_Ee_Horn))
  rownames(temp[[6]]) <- c("Horn size weighted(q=1)")  
  temp <- lapply(temp, FUN = function(x){
    colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") 
    return(x)
  })
  ################################################################
  
  #Cqn=rbind(Cqn_se_equ(X,q=0,nboot),
  #          #C1n_equ(method="relative",X,nboot),
  #          NA,
  #          C1n_equ(method="absolute",X,nboot), 
  #          Cqn_se_equ(X,q=2,nboot)[1:4])
  #if(N==3){Cqn <- rbind(Cqn, C33_se_equ(X,nboot)[1:4])}
  #colnames(Cqn) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
  #rownames(Cqn) <- c(paste("C0",N," (Sorensen)",sep=""),paste("C1",N,"(Horn)",sep=""),paste("C1",N,"*","(Horn)",sep=""),paste("C2",N," (Morisita)",sep=""),if(N==3) "C33")
  #Cqn <- Cqn[-2,]
  
    if(q == 0){
      temp_PC <- rep(0, N*(N-1)/2)
      C02=matrix(0,choose(no.community,2),4)
      U02=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(1,N,N)
      C_SM_2=matrix(1,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          if(sum( X[,i]>0 & X[,j]>0)==0){
            mat <- rbind(c(0, 0), c(0 ,0))
          }else{
            mat <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal effort')
          }
          C02[k,] <- plus_CI(c(1-mat[1, 1],mat[1, 2]))
          U02[k,] <- plus_CI(c(1-mat[2, 1],mat[2, 2]))
          temp_PC[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C02[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- U02[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C02"=C02, "U02"=U02)
      C_SM <- list("C02"=C_SM_1, "U02"=C_SM_2)
    }
    if(q == 1){
      temp_PC <- rep(0, N*(N-1)/2)
      C12=matrix(0,choose(no.community,2),4)
      Horn=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(0,N,N)
      C_SM_2=matrix(0,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          mat <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal weight')
          mat2 <-  Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal effort')
          C12[k,] <- plus_CI(c(1-mat[1, 1],mat[1, 2]))
          Horn[k,] <- plus_CI(c(1-mat2[2, 1],mat2[2, 2]))
          temp_PC[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C12[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- Horn[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C12"=C12, "Horn"=Horn)
      C_SM <- list("C12"=C_SM_1, "Horn"=C_SM_2)
    }
    if(q == 2){
      temp_PC <- rep(0, N*(N-1)/2)
      C22=matrix(0,choose(no.community,2),4)
      U22=matrix(0,choose(no.community,2),4)
      C_SM_1=matrix(0,N,N)
      C_SM_2=matrix(0,N,N)
      k=1
      for(i in 1:(N-1)){  
        for(j in (i+1):N){
          mat <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal weight')
          C22[k,] <- plus_CI(c(1-mat[1, 1],mat[1, 2]))
          U22[k,] <- plus_CI(c(1-mat[2, 1],mat[2, 2]))
          temp_PC[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
          C_SM_1[i,j] <- C_SM_1[j,i] <- C22[k,1]
          C_SM_2[i,j] <- C_SM_2[j,i] <- U22[k,1]
          k <- k+1
        }
      }
      Cqn_PC <- list("C22"=C22, "U22"=U22)
      C_SM <- list("C22"=C_SM_1,"U22"=C_SM_2)
    }
    Cqn_PC <- lapply(Cqn_PC, function(x){
      colnames(x) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL") ; rownames(x) <- temp_PC
      return(x)
    })
    
    #z <- list("datatype"=type, "info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
    
    z <- list("info"=info, "Empirical_richness"=temp[[1]], "Empirical_relative"=temp[[2]], "Empirical_WtRelative"=temp[[3]],
              "estimated_richness"=temp[[4]], "estimated_relative"=temp[[5]], "estimated_WtRelative"=temp[[6]], "pairwise"=Cqn_PC, "dissimilarity_matrix"=C_SM, "q"=q)
  # }else{
  #  z <- list("info"=info, "Empirical_richness"=temp[[1]], "Empirical_relative"=temp[[2]], "Empirical_WtRelative"=temp[[3]],
  #             "estimated_richness"=temp[[4]], "estimated_relative"=temp[[5]], "estimated_WtRelative"=temp[[6]], "q"=q)
  #}
 
  #Cqn=rbind(Cqn_se_equ(X,q=0,nboot),
  #          #C1n_equ(method="relative",X,boot),
  #          NA,
  #          C1n_equ(method="absolute",X,nboot), 
  #          Cqn_se_equ(X,q=2,nboot)[1:4])
  ##if(N==3){Cqn <- rbind(Cqn, C33_se_equ(X,boot))}
  #colnames(Cqn) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
  #rownames(Cqn) <- c(paste("C0",N,sep=""),paste("C1",N,sep=""),paste("C1",N,"*",sep=""),paste("C2",N,sep=""))#,if(N==3) "C33")
  #  
  #
  #if(q==0 || q==1){Cqn_PC=matrix(0,choose(no.community,2),4)}
  #if(q==2)        {Cqn_PC=matrix(0,choose(no.community,2),6)}
  #k=1
  #temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
  #for(i in 1:(N-1)){  
  #  for(j in (i+1):N){
  #    Cqn_PC[k,] <- Cqn_se_equ(X[,c(i,j)],q,nboot,method="absolute")
  #    temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
  #    temp_PD[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
  #    k <- k+1
  #  }
  #}
  #if(q==0 || q==1){
  #  colnames(Cqn_PC) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
  #  rownames(Cqn_PC) <- temp_PC
  #}
  #if(q==2){
  #  colnames(Cqn_PC) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL", "D.95%.LCL", "D.95%.UCL")
  #  rownames(Cqn_PC) <- temp_PC
  #}
  #  
  #C_SM=matrix(1,N,N)
  #k <- 1
  #for(i in 1:(N-1)){
  #  for(j in (i+1):N){
  #    C_SM[i,j] <- C_SM[j,i] <- Cqn_PC[k,1]
  #    k <- k+1
  #  }
  #}
  #z <- list("datatype"=type, "info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  #z <- list("info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "q"=q)
  class(z) <- c("spadeGenetic")
  z
}

