sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("src/Species//")
sourceDir("src/SharedSpecies/")
sourceDir("src/InterpolationAndExtrapolation/")
sourceDir("src/DiversityProfileEstimation/")
sourceDir("src/TwoCommunitySimilarity/")
sourceDir("src/MultipleCommunityMeasure/")
sourceDir("src/Genetic/")
source("src/sub.R")
library(ggplot2)
library(markdown)
library(igraph)
library(gridExtra)
library(shiny)

shinyServer(function(input, output) {
  
  mydata <- reactive({
    if(input$goButton==0) return(NULL)
    isolate({
      if (input$source == "demo") {
        # User has not uploaded a file yet
        if(input$type == "abun"){
          #if(input$index == "Species") dat <- data.1a
          #if(input$index == "Diversity") dat <- data.4a
          if(input$index == "IE") dat <- data.sp
          if(input$index == "Shared Species") dat <- data.2a
          if(input$index == "Two-Community Measures") dat <- data.5a
          if(input$index == "Multiple-Community Measures") dat <- data.6a
          #if(input$index == "Genetics Measures") dat <- data.7a
        }
        if(input$type == "inci"){
          #if(input$index == "Species") dat <- data.1e
          #if(input$index == "Diversity") dat <- data.4e
          if(input$index == "IE") dat <- data.an
          if(input$index == "Shared Species") dat <- data.2c
          if(input$index == "Two-Community Measures") dat <- data.5b
          if(input$index == "Multiple-Community Measures") dat <- data.6b
          #if(input$index == "Genetics Measures") dat <- NULL
        }
        if(input$type == "inci_raw"){
          if(input$index == "Shared Species") dat <- data.2c
          #if(input$index == "Diversity") dat <- data.frame("a" = data.5c[,1])
          if(input$index == "Two-Community Measures") dat <- data.5c
          if(input$index == "Multiple-Community Measures") dat <- data.6c
          if(input$index == "Genetics Measures") dat <- NULL
        }
        if(input$typeG == "abun1") { if(input$index == "Genetics Measures") dat <- data.7a }
        if(input$typeG == "abun2") { if(input$index == "Genetics Measures") dat <- data.7b }
        if(input$type2 == "abun_infr"){ 
          if(input$index == "Species") dat <- data.1b
          if(input$index == "Diversity") dat <- data.4b
          }
        if(input$type2 == "inci_count"){ 
          if(input$index == "Species") dat <- data.1d 
          if(input$index == "Diversity") dat <- data.1e 
          }
        if(input$type2 == "abun_speci"){ 
          if(input$index == "Species") dat <- data.1a 
          if(input$index == "Diversity") dat <- data.4a
          }
        if(input$type2 == "inci_speci"){ 
          if(input$index == "Species") dat <- data.1e 
          if(input$index == "Diversity") dat <- data.4e
          }
        if(input$type2 == "inci_raw"){ 
          if(input$index == "Species") {
            data.1c <- read.table("Data/Data1c(II).txt")
            k <- input$cut
            t <- ncol(data.1c)
            data <- rowSums(data.1c)
            data <- as.integer(data)
            t_infreq <- sum(colSums(data.1c[which(data<k),])>=1)
            data.1c <- data
            data.1c <- data.frame("data"=c(t_infreq, t , data.1c))
            rownames(data.1c) <- c("T_infreq" , "T" , paste("y",1:(nrow(data.1c)-2),sep=""))
            dat <- data.1c
          }
          if(input$index == "Diversity") {
            data.1c <- read.table("Data/Data1c.txt")
            k <- 10
            t <- ncol(data.1c)
            data <- rowSums(data.1c)
            data <- as.integer(data)
            t_infreq <- sum(colSums(data.1c[which(data<k),])>=1)
            data.1c <- data
            data.1c <- data.frame("data"=c(t_infreq, t , data.1c))
            rownames(data.1c) <- c("T_infreq" , "T" , paste("y",1:(nrow(data.1c)-2),sep=""))
            dat <- data.1c
          }
        }
        
      } else {
        ##dat <- read.table(input$files$datapath, header=TRUE)
         dat <- read.table(input$files$datapath) ## by KH MA 2015.01.26
         if(input$index == "Species" | input$index == "Diversity"){
          if(nrow(dat) <= 2){
          if(input$type2 == "abun_infr"){
            dat <- as.integer(dat)
            lengthdat <- length(dat)
            dat <- data.frame(rep(dat[seq(1,lengthdat,2)],dat[seq(2,lengthdat,2)]))
            rownames(dat) <- paste("x",1:nrow(dat),sep="")
            colnames(dat) <- c("data")
          }
          if(input$type2 == "inci_count"){
            t <- dat[,1]
            dat <- dat[,-c(1)]
            dat <- as.integer(dat)
            lengthdat <- length(dat)
            dat <- data.frame(rep(dat[seq(1,lengthdat,2)],dat[seq(2,lengthdat,2)]))
            dat <- rbind(t,dat)
            rownames(dat) <- c("T", paste("y",1:(nrow(dat)-1),sep=""))
            colnames(dat) <- c("data")
          }
        }
        if(input$type2 == "inci_raw"){
          
             k <- input$cut
             t <- ncol(dat)
             data <- rowSums(dat)
             data <- as.integer(data)
             t_infreq <- sum(colSums(dat[which(data<k),])>=1)
             dat <- data
             dat <- data.frame("data"=c(t_infreq, t , dat))
             rownames(dat) <- c("T_infreq" , "T" , paste("y",1:(nrow(dat)-2),sep=""))
          
        }
        }else{
          if(input$type == "inci_raw"){
            t=as.numeric(sapply(readLines(textConnection(input$keyin_t)), function(x) scan(text = x, what = 'char')))
            if(ncol(dat) != sum(t)) stop("Number of columns does not euqal to the sum of key in sampling units")
            data <- matrix(0, ncol = length(t), nrow = nrow(dat))
            n <- 0
            for(i in 1:length(t)){
              data[, i] <- as.integer(rowSums(dat[,(n+1):(n+t[i])] ) )
              n <- n+t[i]
            }
            t <- as.integer(t)
            data <- apply(data, MARGIN = 2, as.integer)
            dat <- data.frame(rbind(t, data),row.names = NULL)
            rownames(dat) <- c("T" , paste("y",1:(nrow(dat)-1),sep = ""))
            colnames(dat) <- 1:length(t)
          }else{
            dat
          }
        }
      }
      dat
    })
  })
  
  
  rawdat <- reactive({
    if(input$goButton==0) return(NULL)
    isolate({
    if (input$source == "demo"){
        if(input$type == "inci_raw" | input$type2 == "inci_raw"){
           if(input$index == "Species") dat <- read.table("Data/Data1c(II).txt")
           if(input$index == "Diversity"){
           dat <- read.table("Data/Data1c.txt")
           dat <- sapply(1:ncol(dat) , function(i) as.integer(dat[,i]))
           dat <- data.frame(dat, row.names = NULL)
           }
           if(input$index == "Shared Species") {
           dat <-  read.csv("Data/Data2c_raw.csv")[,-1]
           #dat1 <- get(dat)[[1]] ; dat2 <- get(dat)[[2]]
           #dat1 <- data.frame(dat1, row.names = NULL) ; dat2 <- data.frame(dat2, row.names = NULL)
           #dat1 <- sapply(1:ncol(dat1) , function(i) as.integer(dat1[,i]))
           #dat2 <- sapply(1:ncol(dat2) , function(i) as.integer(dat2[,i]))
           #dat <- cbind(dat1,dat2)
           dat <- data.frame(dat)
           }
           if(input$index == "Two-Community Measures"){
           dat <- load("Data/ciliates.rda")
           dat1 <- get(dat)[[1]] ; dat2 <- get(dat)[[2]]
           dat1 <- data.frame(dat1, row.names = NULL) ; dat2 <- data.frame(dat2, row.names = NULL)
           dat1 <- sapply(1:ncol(dat1) , function(i) as.integer(dat1[,i]))
           dat2 <- sapply(1:ncol(dat2) , function(i) as.integer(dat2[,i]))
           dat <- cbind(dat1,dat2)
           }
           if(input$index == "Multiple-Community Measures"){
           dat <- load("Data/ciliates.rda")
           dat1 <- get(dat)[[1]] ; dat2 <- get(dat)[[2]] ; dat3 <- get(dat)[[3]]
           dat1 <- data.frame(dat1, row.names = NULL) 
           dat2 <- data.frame(dat2, row.names = NULL)
           dat3 <- data.frame(dat3, row.names = NULL)
           dat1 <- sapply(1:ncol(dat1) , function(i) as.integer(dat1[,i]))
           dat2 <- sapply(1:ncol(dat2) , function(i) as.integer(dat2[,i]))
           dat3 <- sapply(1:ncol(dat3) , function(i) as.integer(dat3[,i]))
           dat <- cbind(dat1,dat2,dat3)
           }   
        }
    }else{
    if(input$type == "inci_raw" | input$type2 == "inci_raw") dat <- read.table(input$files$datapath)
    }  
    return(dat)
    })
  })
  
  output$rawdata <- renderUI({
    if(input$goButton==0) return(NULL)
    isolate({
      type1 <- c("abun","inci")
      type2 <- c("abun_speci","abun_infr","inci_speci","inci_count")
      section1 <- c("Species","Diversity")
      section2 <- c("Shared Species","Two-Community Measures","Multiple-Community Measures")
      if( (input$index %in% section1) & !(input$index %in% section2)){
        if( (input$type2 == "inci_raw" & !(input$type2 %in% type2) ) ){
          downloadLink("dlrawtab", "Download demo/uploaded raw data")
        }
      }else if( !(input$index %in% section1) & (input$index %in% section2)){
        if( (input$type == "inci_raw" & !(input$type %in% type1) ) ){
          downloadLink("dlrawtab", "Download demo/uploaded raw data")
        }
      }
        })
    })
  

  #output$rawdata <- renderUI({
  #  if(input$goButton==0) return(NULL)
  #  isolate({
  #    if(input$type == "inci_raw" | input$type2 == "inci_raw"){
  #    if(input$index == "Diversity" | input$index == "Species"){
  #      tagList(
  #        h4("Raw Data (Demo for part of original raw data)"),
  #        tableOutput("table_raw"),
  #        downloadLink("dlrawtab", "Download raw dat as txt file")
  #      )
  #    }else if(input$index == "Two-Community Measures" | input$index == "Shared Species"){
  #        tagList(
  #          h4("Raw Data (Demo for part of original raw data)"),
  #          tableOutput("table_raw2"),
  #          tableOutput("table_raw3"),
  #          downloadLink("dlrawtab", "Download raw data as txt file")
  #        )
  #    }else {
  #        return(NULL)
  #      }
  #    }else {
  #      return(NULL)
  #    }
  #  })
  #})
  #  
  #output$table_raw <- renderTable({
  #  if(input$goButton==0) return(NULL)
  #  isolate({
  #    if(input$type == "inci_raw" | input$type2 == "inci_raw"){
  #      if(input$index == "Diversity" | input$index == "Species"){
  #        rawdat()[1:5,1:5]
  #      }else{
  #        return(NULL)
  #      }
  #    }else{
  #      return(NULL)
  #    }
  #  })
  #})
  #
  #output$table_raw2 <- renderTable({
  #  if(input$goButton==0) return(NULL)
  #  isolate({
  #    if(input$type == "inci_raw" | input$type2 == "inci_raw"){
  #      if(input$index == "Two-Community Measures" | input$index == "Shared Species"){
  #        rawdat()[[1]][1:5,1:5]
  #      }else{
  #        return(NULL)
  #      }
  #    }else{
  #      return(NULL)
  #    }
  #  })
  #})
  #
  #output$table_raw3 <- renderTable({
  #  if(input$goButton==0) return(NULL)
  #  isolate({
  #    if(input$type == "inci_raw" | input$type2 == "inci_raw"){
  #      if(input$index == "Two-Community Measures" | input$index == "Shared Species"){
  #    rawdat()[[2]][1:5,1:5]
  #      }else{
  #        return(NULL)
  #      }
  #    }else{
  #      return(NULL)
  #    }
  #  })
  #})
  
  output$table <- renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      if(input$index == "Species"| input$index == "Diversity"){
        if(input$type2 == "inci_count" | input$type2 == "inci_speci"){
        dat <- as.data.frame(mydata()[1:11,])
        rownames(dat) <- c("T", paste("y",1:10,sep=""))
        colnames(dat) <-  colnames(mydata())
        dat
        }else if(input$type2 == "inci_raw"){
        
         dat <- as.data.frame(mydata()[1:12,])
         rownames(dat) <- rownames(mydata())[1:12]
         colnames(dat) <-  colnames(mydata())
         dat
       }else{
         dat <- as.matrix(mydata()[1:10,])
         rownames(dat) <- c(paste("x",1:10,sep=""))
         colnames(dat) <-  colnames(mydata())
         dat
       }
      }else if(input$index== "Genetics Measures"){
          dat <- as.matrix(mydata()[1:10,])
          rownames(dat) <- c(paste("x",1:10,sep=""))
          colnames(dat) <-  colnames(mydata())
          dat
      }else{
        if(input$type == "inci" | input$type == "inci_raw"){
          dat <- as.data.frame(mydata()[1:11,])
          rownames(dat) <- c("T", paste("y",1:10,sep=""))
          colnames(dat) <-  colnames(mydata())
          dat
        }else{
          dat <- as.matrix(mydata()[1:10,])
          rownames(dat) <- c(paste("x",1:10,sep=""))
          colnames(dat) <-  colnames(mydata())
          dat
        }
      }
    })
  })
  
  
  computation <- reactive({
    if(input$goButton==0) return(cat(""))
    isolate({
      if(input$type2=="abun_speci" | input$type2=="abun_infr" ){
        if(input$index=="Species"){
          out <- ChaoSpecies_all(data=unlist(mydata()), datatype="abundance", k=input$cut)
        }  
        if(input$index=="Diversity"){
          out <- Diversity(X=mydata(),B = input$nboot, datatype="abundance",q=as.numeric(sapply(readLines(textConnection(input$orderq)), function(x) scan(text = x, what = 'char'))))
        }
      }
      if(input$typeG=="abun1" | input$typeG=="abun2"){
        if(input$index=="Genetics Measures"){
          out <- Genetic(mydata(), q=as.numeric(input$order), nboot=input$nboot)
        }  
      }
      if(input$type=="abun"){        
        if(input$index=="IE"){
          x <- as.numeric(c(mydata())[[1]])
          endpoint <- ifelse(is.na(input$endpt), 2*sum(x), input$endpt)
          out <- iNEXT(x, q=as.numeric(input$order), datatype="abundance", endpoint=endpoint,
                       knot=input$knot, se=(input$nboot>0), nboot=input$nboot)
        }
        if(input$index=="Shared Species"){
          out <- ChaoShared(data=mydata(), datatype="abundance", nboot=input$nboot, conf=0.95)
        }        
        if(input$index=="Two-Community Measures"){
          out <- Two_Community_Similarity(X=mydata(), datatype="abundance", boot=input$nboot)
        }
        if(input$index=="Multiple-Community Measures"){
          if(as.numeric(input$order) == 0 )  method <- "absolute"
          if(as.numeric(input$order) != 0 )  method <- input$method
          out <- Multiple_Community_Measure(X=mydata(), datatype="abundance", q=as.numeric(input$order), nboot=input$nboot, method)
        }
      } 
      if(input$type2=="inci_speci"){
        if(input$index=="Species"){
          out <- ChaoSpecies_all(data=mydata(), datatype="incidence", k=input$cut)
        } 
        if(input$index=="Diversity"){
          out <- Diversity(X=mydata(),B = input$nboot, datatype="incidence",q=as.numeric(sapply(readLines(textConnection(input$orderq)), function(x) scan(text = x, what = 'char'))))
        }
      }
      if(input$type2=="inci_raw"){
        if(input$index=="Species"){
          out <- ChaoSpecies_all(data=mydata(), datatype="incidence_raw", k=input$cut)
        } 
        if(input$index=="Diversity"){
          out <- Diversity(X=mydata(),B = input$nboot, datatype="incidence_raw",q=as.numeric(sapply(readLines(textConnection(input$orderq)), function(x) scan(text = x, what = 'char'))))
        }
      }
      if(input$type2=="inci_count"){
        if(input$index=="Species"){
          out <- ChaoSpecies_all(data=mydata(), datatype="incidence", k=input$cut)
        } 
        if(input$index=="Diversity"){
          out <- Diversity(X=mydata(),B = input$nboot, datatype="incidence",q=as.numeric(sapply(readLines(textConnection(input$orderq)), function(x) scan(text = x, what = 'char'))))
        }
      }
      if(input$type=="inci" | input$type=="inci_raw"){
        if(input$index=="IE"){
          x <- unlist(mydata())
          endpoint <- ifelse(is.na(input$endpt), 2*max(x), input$endpt)
          out <- iNEXT(x, q=as.numeric(input$order), datatype="incidence", 
                       knot=input$knot, endpoint=endpoint, se=(input$nboot>0), nboot=input$nboot)
        }
        if(input$index=="Shared Species"){
          out <- ChaoShared(data=mydata(), datatype="incidence", nboot=input$nboot, conf=0.95)
        }        
        if(input$index=="Two-Community Measures"){
          out <- Two_Community_Similarity(X=mydata(), datatype="incidence", boot=input$nboot)
        }
        if(input$index=="Multiple-Community Measures"){
          if(as.numeric(input$order) == 0 )  method <- "absolute"
          if(as.numeric(input$order) != 0 )  method <- input$method
          out <- Multiple_Community_Measure(X=mydata(), datatype="incidence", q=as.numeric(input$order), nboot=input$nboot, method)
        }
      }
      out
    })
  })
  
  output$estimation <- renderPrint({
    computation()
  })
  
  myplot <- reactive({
    if(input$goButton==0) return(NULL)
    isolate({
      if(input$index=="Species"){
        if(input$type2=="abun_speci"){
          name.ag <- c("Homo","MLE","Chao1","Chao1-bc","iChao1","ACE","ACE-1","Jack1","Jack2")
        }
        if(input$type2=="abun_infr"){
          name.ag <- c("Homo","MLE","Chao1","Chao1-bc","iChao1","ACE","ACE-1","Jack1","Jack2")
        } 
        if(input$type2=="inci_speci"){
          name.ag <- c("Homo","Chao2","Chao2-bc","iChao2","ICE","ICE-1","Jack1","Jack2")
        }
        if(input$type2=="inci_count"){
          name.ag <- c("Homo","Chao2","Chao2-bc","iChao2","ICE","ICE-1","Jack1","Jack2")
        }
        if(input$type2=="inci_raw"){
          name.ag <- c("Homo","Chao2","Chao2-bc","iChao2","ICE","ICE-1","Jack1","Jack2")
        }
        tab <- computation()$Species.Table        
        #bp <- barplot(tab[,1], beside=T, names.arg=name.ag,
        #              ylim=c(0,max(tab[,4])),
        #              main="Comparison of estimators",ylab="Species richness", las=2)
        plot(tab[,1], ylim=c(min(tab[,3]),max(tab[,4])), pch=2, cex=1.5, xlab="",
             main="Comparison of estimators",ylab="Species richness", las=2, xaxt="n")
        axis(1, at=1:nrow(tab), labels=name.ag)
        arrows(1:nrow(tab), tab[,4], 1:nrow(tab), tab[,3], angle=90, code=3, length=.1)
        
        
      }
      if(input$index=="IE"){
        g1 <- ggiNEXT(computation(), type = 1)
        g2 <- ggiNEXT(computation(), type = 2)
        g3 <- ggiNEXT(computation(), type = 3)
        grid.arrange(g1, g2, g3, ncol=3)
        #par(mfrow=c(1,3), cex=1.5)
        #plot(computation(),style="N2D")
        #plot(computation(),style="N2SC")
        #plot(computation(),style="SC2D")
      }
      
      if(input$index=="Diversity"){
        out <- computation()[["HILL.NUMBERS"]]
        out$ChaoL <- out[,3]
        out$ChaoU <- out[,4]
        out$MleL <- out[,6]
        out$MleU <- out[,7]
        
        ################################################################2015.09.15-(S.W.Wei)
        ymin = min(c(out$MleL, out$ChaoL))
        ymax = max(c(out$MleU, out$ChaoU))
        if(!is.na(ymin)) plot(out$q, out$ChaoJost, type = "l", xlab = "Order q", ylab = "Hill numbers", col = 2, ylim = c(ymin, ymax))
        else plot(out$q, out$ChaoJost, type = "l", xlab = "Order q", ylab = "Hill numbers", col = 2)
        points(out$q, out$Empirical, type = "l", col = 4, lty = 2)
        conf.reg(out$q, out$ChaoL, out$ChaoU, border = NA, col = adjustcolor(2, 0.2))
        conf.reg(out$q, out$MleL, out$MleU, border = NA, col = adjustcolor(4,0.2))
        legend("topright", c("ChaoJost", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
        #plot(out$q,out$Chao,type="l",xlab="Order q",ylab = "Hill numbers",lty=2,col=4,ylim = c(ymin,ymax))
        #points(out$q,out$Empirical,type="l",col=2)
        #conf.reg(out$q,out$ChaoL,out$ChaoU,border=NA,col=adjustcolor(4,0.2))
        #conf.reg(out$q,out$MleL,out$MleU,border=NA,col=adjustcolor(2,0.2))
        #legend("topright",c("Chao","Empirical"),col=c(2,4),lty=c(1,2),bty="n")
        ################################################################2015.09.15
        
      }
      
      
      
      if( input$index=="Multiple-Community Measures" | input$index=="Two-Community Measures" ){
        out <- computation()
        N <- out$info[1]
        if(input$index=="Two-Community Measures") N <- 2
        if(input$method2 == "relative" ){
          mle <- rbind(out$Empirical_richness[1,],out$Empirical_relative[1,], out$Empirical_relative[2,] )
          est <- rbind(out$estimated_richness[1,],out$estimated_relative[1,], out$estimated_relative[2,] )
          
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          if(is.na(ymin)){
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(C[q][N]~(relative~abundance), list(N = N)), col = 2)
          } else {
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(C[q][N]~(relative~abundance), list(N = N)), col = 2, ylim = c(ymin, ymax))
            
          }
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          
        }else{
          mle <- rbind(out$Empirical_richness[1,],out$Empirical_absolute[1,], out$Empirical_absolute[2,] )
          est <- rbind(out$estimated_richness[1,],out$estimated_absolute[1,], out$estimated_absolute[2,] )
          
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          if(is.na(ymin)){
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(C[q][N]~(relative~abundance), list(N = N)), col = 2)
          } else {
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(C[q][N]~(relative~abundance), list(N = N)), col = 2, ylim = c(ymin, ymax))
            
          }
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
        }
      }
      
      if(input$index=="Shared Species"){
        if(input$type=="abun"){
          name.ag <- c("Homo","ACE","Chao1-shared","Chao1-shared-bc")
        } else{
          name.ag <- c("Chao2-shared","Chao2-shared-bc")
        }
        tab <- computation()$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES
        
        #bp <- barplot(tab[,1], beside=T, names.arg=name.ag,
        #              ylim=c(0,max(tab[,4])),
        #              main="Comparison of estimators", ylab="Shared species", las=2)
        #arrows(bp, tab[,4], bp, tab[,3], angle=90, code=3, length=.1)
        ymin = min(tab[,3])
        ymax = max(tab[,4])
        # plot(tab[,1], ylim=c(ymin, ymax),
        #      pch=2, cex=1.5, xlab="",main="Comparison of estimators",ylab="Shared species", las=2, xaxt="n")
        
        if(!is.na(ymin)) plot(tab[,1], ylim=c(min(tab[,3]),max(tab[,4])),
                              pch=2, cex=1.5, xlab="",main="Comparison of estimators",ylab="Shared species", las=2, xaxt="n")
        else plot(tab[,1],ylim = c(min(tab[ ,1], na.rm = T),max(tab[,4], na.rm = T)), pch=2, cex=1.5, xlab="",main="Comparison of estimators",ylab="Shared species", las=2, xaxt="n")
        axis(1, at=1:nrow(tab), labels=name.ag)
        arrows(1:nrow(tab), tab[,4], 1:nrow(tab), tab[,3], angle=90, code=3, length=.1)
        
      }
      #       if(input$index=="Multiple-Community Measures"){
      #         if(input$type=="abun"){
      #           termMatrix <- computation()$similarity.matrix
      #           rownames(termMatrix) <- colnames(termMatrix) <- names(mydata())
      #           g <- graph.adjacency(termMatrix, weighted=T, mode = "undirected")
      #           g <- simplify(g)
      #           set.seed(1016)
      #           layout1 <- layout.fruchterman.reingold(g)
      #           w <- 10^(E(g)$weight)
      #           E(g)$width <- 10*((w-min(w))/(max(w)-min(w))) + 0.1
      #           E(g)$color <-  sample(colors()[30:300], length(E(g)$weight), replace=TRUE)
      #           V(g)$label <- V(g)$name
      #           plot.igraph(g, layout=layout1, vertex.label.family="STHeiti", vertex.label.cex=1.4)
      #         }else {
      #           cat("Currently we only provide Abundance data for Multiple-Community Diversity Measure.")
      #         }
      #       }
      
    })
    
  })
  
  myplot2 <- reactive({
    if(input$goButton==0) return(NULL)
    isolate({
      if( input$index=="Multiple-Community Measures" | input$index=="Two-Community Measures" ){
        out <- computation()
        N <- out$info[1]
        if(input$index=="Two-Community Measures") N <- 2
        if(input$method2 == "relative" ){
          mle <- rbind(out$Empirical_richness[2,],out$Empirical_relative[1,], out$Empirical_relative[3,] )
          est <- rbind(out$estimated_richness[2,],out$estimated_relative[1,], out$estimated_relative[3,] )
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          if(is.na(ymin)){
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(U[q][N]~(relative~abundance), list(N = N)), col = 2, ylim = c(min(c(mle[,1], est[,1]), na.rm = T), max(c(mle[,1], est[,1]), na.rm = T)))
          } else {
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(U[q][N]~(relative~abundance), list(N = N)), col = 2, ylim = c(ymin, ymax))
          }
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          
        }else{
          mle <- rbind(out$Empirical_richness[2,],out$Empirical_absolute[1,], out$Empirical_absolute[3,] )
          est <- rbind(out$estimated_richness[2,],out$estimated_absolute[1,], out$estimated_absolute[3,] )
          
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          if(is.na(ymin)){
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(U[q][N]~(relative~abundance), list(N = N)), col = 2,
                 ylim = c(min(c(mle[,1], est[,1]), na.rm = T), max(c(mle[,1], est[,1]), na.rm = T)))
          } else {
            plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(U[q][N]~(relative~abundance), list(N = N)), col = 2, ylim = c(ymin, ymax))
          }
          
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
        }
        
      }
    })
  })
  
  output$visualization <- renderPlot({
    myplot()
  })
  
  output$visualization2 <- renderPlot({
    myplot2()
  })
  
  output$dlrawtab <- downloadHandler(
    filename = function() { paste('data_', Sys.Date(), '_[spadeR].csv', sep='') },
    content = function(file) { 
      out <- rawdat()
      write.csv(out, file=file, row.names=FALSE)
    }
  )
  
  output$dltab <- downloadHandler(
    filename = function() { paste('data_', Sys.Date(), '_[spadeR].txt', sep='') },
    content = function(file) { 
      out <- mydata()
      write.table(out, file=file, sep="\t", dec=".", 
                  quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
    }
  )
  
  output$dlest <- downloadHandler(
    filename = function() { paste('output_', Sys.Date(), '_[spadeR].txt', sep='') },
    content = function(file) { 
      out <- capture.output(computation(), file=file, append=TRUE)
      write.table(out, file=file, sep="\t", dec=".", 
                  quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
    }
  )
  
  ########################################################2015.09.28-(S.W.Wei)
  output$dlvis <- downloadHandler(
    filename = function() { paste('figure_', Sys.Date(), '_[spadeR].png', sep='') },
    content <- function(file){
      if(input$index=="Species"){
        png(file)
        if(input$type2=="abun_speci"){
          name.ag <- c("Homo","MLE","Chao1","Chao1-bc","iChao1","ACE","ACE-1","Jack1","Jack2")
        }
        if(input$type2=="abun_infr"){
          name.ag <- c("Homo","MLE","Chao1","Chao1-bc","iChao1","ACE","ACE-1","Jack1","Jack2")
        }
        if(input$type2=="inci_speci"){
          name.ag <- c("Homo","Chao2","Chao2bc","iChao2","ICE","ICE-1","Jack1","Jack2")
        }
        if(input$type2=="inci_count"){
          name.ag <- c("Homo","Chao2","Chao2bc","iChao2","ICE","ICE-1","Jack1","Jack2")
        }
        if(input$type2=="inci_raw"){
          name.ag <- c("Homo","Chao2","Chao2bc","iChao2","ICE","ICE-1","Jack1","Jack2")
        }
        tab <- computation()$Species.Table
        #bp <- barplot(tab[,1], beside=T, names.arg=name.ag,
        #              ylim=c(0,max(tab[,4])),
        #              main="Comparison of estimators",ylab="Species richness", las=2)
        plot(tab[,1], ylim=c(min(tab[,3]),max(tab[,4])), pch=2, cex=1.5, xlab="",
             main="Comparison of estimators",ylab="Species richness", las=2, xaxt="n")
        axis(1, at=1:length(name.ag), labels=name.ag)
        arrows(1:nrow(tab), tab[,4], 1:nrow(tab), tab[,3], angle=90, code=3, length=.1)
        dev.off()
        
      }
      #       if(input$index=="IE"){
      #         pdf(file, width=7, height=21)
      #         g1 <- ggiNEXT(computation(), type = 1)
      #         g2 <- ggiNEXT(computation(), type = 2)
      #         g3 <- ggiNEXT(computation(), type = 3)
      #         grid.arrange(g1, g2, g3, ncol=1)
      #         dev.off()
      #       }
      if(input$index=="Shared Species"){
        png(file)
        if(input$type=="abun"){
          name.ag <- c("Homo","ACE","Chao1-shared","Chao1-shared-bc")
        } else{
          name.ag <- c("Chao2-shared","Chao2-shared-bc")
        }
        tab <- computation()$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES
        
        #bp <- barplot(tab[,1], beside=T, names.arg=name.ag,
        #              ylim=c(0,max(tab[,4])),
        #              main="Comparison of estimators", ylab="Shared species", las=2)
        #arrows(bp, tab[,4], bp, tab[,3], angle=90, code=3, length=.1)
        plot(tab[,1], ylim=c(min(tab[,3]),max(tab[,4])), pch=2, cex=1.5, xlab="",
             main="Comparison of estimators",ylab="Shared species", las=2, xaxt="n")
        axis(1, at=1:nrow(tab), labels=name.ag)
        arrows(1:nrow(tab), tab[,4], 1:nrow(tab), tab[,3], angle=90, code=3, length=.1)
        dev.off()
        
      }
      if(input$index == "Diversity"){
        png(file)
        out <- computation()[["HILL.NUMBERS"]]
        out$ChaoL <- out[,3]
        out$ChaoU <- out[,4]
        out$MleL <- out[,6]
        out$MleU <- out[,7]
        
        ymin = min(c(out$MleL, out$ChaoL))
        ymax = max(c(out$MleU, out$ChaoU))
        
        plot(out$q, out$ChaoJost, type = "l", xlab = "Order q", ylab = "Hill numbers", col = 2, ylim = c(ymin, ymax))
        points(out$q, out$Empirical, type = "l", col = 4, lty = 2)
        conf.reg(out$q, out$ChaoL, out$ChaoU, border = NA, col = adjustcolor(2, 0.2))
        conf.reg(out$q, out$MleL, out$MleU, border = NA, col = adjustcolor(4,0.2))
        legend("topright", c("ChaoJost", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
        dev.off()
      }
      if( input$index=="Multiple-Community Measures" | input$index=="Two-Community Measures" ){
        png(file)
        out <- computation()
        N <- out$info[1]
        if(input$index=="Two-Community Measures") N <- 2
        if(input$method2 == "relative" ){
          mle <- rbind(out$Empirical_richness[1,],out$Empirical_relative[1,], out$Empirical_relative[2,] )
          est <- rbind(out$estimated_richness[1,],out$estimated_relative[1,], out$estimated_relative[2,] )
          
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          
          plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(C[q][N], list(N = N)), col = 2, ylim = c(ymin, ymax))
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          
        }else{
          mle <- rbind(out$Empirical_richness[1,],out$Empirical_absolute[1,], out$Empirical_absolute[2,] )
          est <- rbind(out$estimated_richness[1,],out$estimated_absolute[1,], out$estimated_absolute[2,] )
          
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          
          plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(C[q][N], list(N = N)), col = 2, ylim = c(ymin, ymax))
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
        }
        dev.off()
      }
    }
  )
  output$dlvis2 <- downloadHandler(
    filename = function() { paste('figure_', Sys.Date(), '_[spadeR].png', sep='') },
    content <- function(file){
      if(input$index=="Multiple-Community Measures" | input$index=="Two-Community Measures" ){
        png(file)
        out <- computation()
        N <- out$info[1]
        if(input$index=="Two-Community Measures") N <- 2
        if(input$method2 == "relative" ){
          mle <- rbind(out$Empirical_richness[2,],out$Empirical_relative[1,], out$Empirical_relative[3,] )
          est <- rbind(out$estimated_richness[2,],out$estimated_relative[1,], out$estimated_relative[3,] )
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          
          plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(U[q][N], list(N = N)), col = 2, ylim = c(ymin, ymax))
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          
        }else{
          mle <- rbind(out$Empirical_richness[2,],out$Empirical_absolute[1,], out$Empirical_absolute[3,] )
          est <- rbind(out$estimated_richness[2,],out$estimated_absolute[1,], out$estimated_absolute[3,] )
          
          ymin = min(c(mle[,3], est[,3]))
          ymax = max(c(mle[,4], est[,4]))
          
          plot(c(0,1,2), est[,1], type = "l", xlab = "Order q", ylab = substitute(U[q][N], list(N = N)), col = 2, ylim = c(ymin, ymax))
          points(c(0,1,2), mle[,1], type = "l", col = 4, lty = 2)
          conf.reg(c(0,1,2), est[,3], est[,4], border = NA, col = adjustcolor(2, 0.2))
          conf.reg(c(0,1,2), mle[,3], mle[,4], border = NA, col = adjustcolor(4,0.2))
          if(mle[1,1]>mle[3,1]) legend("topright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
          if(mle[1,1]<mle[3,1]) legend("bottomright", c("Estimated", "Empirical"), col = c(2, 4), lty = c(1, 2), cex = 1.5, bty = "n")
        }
        dev.off()
      }
      
    }
  )
  ########################################################2015.09.28
  #   output$dlvis <- downloadHandler(
  #     filename = function() { paste('figure_', Sys.Date(), '_[spadeR].pdf', sep='') },
  #     content = function(file) {
  #       if(input$index=="Species"){
  #         pdf(file, width=9, height=6)
  #         if(input$type=="abun"){
  #           name.ag <- c("Homo","MLE","Chao1","Chao1-bc","iChao1","ACE","ACE-1","Jack1","Jack2")
  #         } else{
  #           name.ag <- c("Homo","Chao2","Chao2bc","iChao2","ICE","ICE-1","Jack1","Jack2")
  #         }
  #         tab <- computation()$Species.Table
  #         #bp <- barplot(tab[,1], beside=T, names.arg=name.ag,
  #         #              ylim=c(0,max(tab[,4])),
  #         #              main="Comparison of estimators",ylab="Species richness", las=2)
  #         plot(tab[,1], ylim=c(min(tab[,3]),max(tab[,4])), pch=2, cex=1.5, xlab="",
  #              main="Comparison of estimators",ylab="Species richness", las=2, xaxt="n")
  #         axis(1, at=1:length(name.ag), labels=name.ag)
  #         arrows(1:nrow(tab), tab[,4], 1:nrow(tab), tab[,3], angle=90, code=3, length=.1)
  #         dev.off()
  #         
  #       }
  #       #       if(input$index=="IE"){
  #       #         pdf(file, width=7, height=21)
  #       #         g1 <- ggiNEXT(computation(), type = 1)
  #       #         g2 <- ggiNEXT(computation(), type = 2)
  #       #         g3 <- ggiNEXT(computation(), type = 3)
  #       #         grid.arrange(g1, g2, g3, ncol=1)
  #       #         dev.off()
  #       #       }
  #       if(input$index=="Shared Species"){
  #         pdf(file)
  #         if(input$type=="abun"){
  #           name.ag <- c("Homo","ACE","Chao1-shared","Chao1-shared-bc")
  #         } else{
  #           name.ag <- c("Chao2-shared","Chao2-shared-bc")
  #         }
  #         tab <- computation()$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES
  #         
  #         #bp <- barplot(tab[,1], beside=T, names.arg=name.ag,
  #         #              ylim=c(0,max(tab[,4])),
  #         #              main="Comparison of estimators", ylab="Shared species", las=2)
  #         #arrows(bp, tab[,4], bp, tab[,3], angle=90, code=3, length=.1)
  #         plot(tab[,1], ylim=c(min(tab[,3]),max(tab[,4])), pch=2, cex=1.5, xlab="",
  #              main="Comparison of estimators",ylab="Shared species", las=2, xaxt="n")
  #         axis(1, at=1:nrow(tab), labels=name.ag)
  #         arrows(1:nrow(tab), tab[,4], 1:nrow(tab), tab[,3], angle=90, code=3, length=.1)
  #         dev.off()
  #         
  #       }
  #       if(input$index == "Diversity"){
  #         pdf(file)
  #         out <- computation()[["HILL.NUMBERS"]]
  #         out$ChaoL <- out[,4]
  #         out$ChaoU <- out[,5]
  #         out$MleL <- out[,6]
  #         out$MleU <- out[,7]
  #         
  #         ymin = min(c(out$MleL, out$ChaoL))
  #         ymax = max(c(out$MleU, out$ChaoU))
  #         
  #         plot(out$q, out$Chao, type = "l", xlab = "Order q", ylab = "Hill numbers", col = 2, ylim = c(ymin, ymax))
  #         points(out$q, out$Empirical, type = "l", col = 4, lty = 2)
  #         conf.reg(out$q, out$ChaoL, out$ChaoU, border = NA, col = adjustcolor(2, 0.2))
  #         conf.reg(out$q, out$MleL, out$MleU, border = NA, col = adjustcolor(4,0.2))
  #         legend("topright", c("Chao", "Empirical"), col = c(2, 4), lty = c(1, 2), bty = "n")
  #         dev.off()
  #       }
  #       #       if(input$index=="Multiple-Community Measures"){
  #       #         pdf(file)
  #       #         if(input$type=="abun"){
  #       #           termMatrix <- computation()$similarity.matrix
  #       #           rownames(termMatrix) <- colnames(termMatrix) <- names(mydata())
  #       #           g <- graph.adjacency(termMatrix, weighted=T, mode = "undirected")
  #       #           g <- simplify(g)
  #       #           set.seed(1016)
  #       #           layout1 <- layout.fruchterman.reingold(g)
  #       #           w <- 10^(E(g)$weight)
  #       #           E(g)$width <- 10*((w-min(w))/(max(w)-min(w))) + 0.1
  #       #           E(g)$color <-  sample(colors()[30:300], length(E(g)$weight), replace=TRUE)
  #       #           V(g)$label <- V(g)$name
  #       #           plot.igraph(g, layout=layout1, vertex.label.cex=1.4)
  #       #         }else {
  #       #           cat("Currently we only provide Abundance data for Multiple-Community Diversity Measure.")
  #       #         }
  #       #         dev.off()  
  #       #       }
  #       
  #     }
  #   )
  
})