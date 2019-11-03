MynumericInput <- function (inputId, label, value, min = NA, max = NA, step = NA, placeholder = "") 
{
  inputTag <- tags$input(id = inputId, type = "number", value = value, placeholder = placeholder)
  if (!is.na(min)) 
    inputTag$attribs$min = min
  if (!is.na(max)) 
    inputTag$attribs$max = max
  if (!is.na(step)) 
    inputTag$attribs$step = step
  tagList(tags$label(label, `for` = inputId), inputTag)
}


saveList2txt <- function(out, file) {
  for (i in seq_along(out)){
    write(names(out)[i], file=file, append=TRUE)  #writes the name of the list elements
    write(out[[i]], file=file, append=TRUE)  #writes the data.frames
  }
}

# Example data one community (abundance data)
# 
data.1a <- read.table("Data/Data1a.txt")
rownames(data.1a) <- paste("x",1:nrow(data.1a),sep="")

data.1b <- as.integer(read.table("Data/Data1b.txt"))
length1b <- length(data.1b)
data.1b <- data.frame(rep(data.1b[seq(1,length1b,2)],data.1b[seq(2,length1b,2)]))
rownames(data.1b) <- paste("x",1:nrow(data.1b),sep="")
colnames(data.1b) <- c("data")

data.1c <- read.table("Data/Data1c.txt")

data.1d <- as.integer(read.table("Data/Data1d.txt")[-c(1)])
t <- read.table("Data/Data1d.txt")[,1]
length1d <- length(data.1d)
data.1d <- data.frame(rep(data.1d[seq(1,length1d,2)],data.1d[seq(2,length1d,2)]))
data.1d <- rbind(t,data.1d)
rownames(data.1d) <- c("T", paste("y",1:(nrow(data.1d)-1),sep=""))
colnames(data.1d) <- c("data")



data.4a <- read.table("Data/Data4a.txt")
rownames(data.4a) <- paste("x",1:nrow(data.4a),sep="")

data.4b <- as.integer(read.table("Data/data4b.txt"))
length4b <- length(data.4b)
data.4b <- data.frame(rep(data.4b[seq(1,length4b,2)],data.4b[seq(2,length4b,2)]))
rownames(data.4b) <- paste("x",1:nrow(data.4b),sep="")
colnames(data.4b) <- c("data")


# Example data one community (incidence data)
#data.o2 <- read.table("Data/Data1c.txt")
#data.o2 <- data.frame(V1=c(t=ncol(data.o2), y=rowSums(data.o2)))


data.1e <- read.table("Data/Data1e.txt")
rownames(data.1e) <- c("T", paste("y",1:(nrow(data.1e)-1),sep=""))

data.4e <- read.table("Data/Data5b.txt")
data.4e <- data.frame("soil"=data.4e[,1])
rownames(data.4e) <- c("T", paste("y",1:(nrow(data.4e)-1),sep=""))

# Example data for two communities (abundance data)
#data.t1 <- read.table("Data/Data5a.txt")
#data.t1 <- data.frame("seedlings"=data.t1[,1], "trees"=data.t1[,2]+data.t1[,3])
#rownames(data.t1) <- paste("x",1:nrow(data.t1),sep="")
data.2a <- read.table("Data/Data2a.txt")
rownames(data.2a) <- paste("x",1:nrow(data.2a),sep="")

data.5a <- read.table("Data/Data5a.txt")[,c(1,3)]
rownames(data.5a) <- paste("x",1:nrow(data.5a),sep="")

data.6a <- read.table("Data/Data6a.txt")
rownames(data.6a) <- paste("x",1:nrow(data.6a),sep="")

# Example data for two communities (incidence data)
# Hong Kong Big Bird Race (BBR) Data
data.2b <- read.table("Data/Data2b.txt")
rownames(data.2b) <- c("T", paste("y",1:(nrow(data.2b)-1),sep=""))

# Example data for multiple community (abundance data)
data.sp <- read.table("Data/spider.txt", header=TRUE)
rownames(data.sp) <- paste("x",1:nrow(data.sp),sep="")
data.sp <- data.sp["Girdled"]

data.6a <- read.table("Data/Data6a.txt", header=FALSE)
rownames(data.6a) <- paste("x",1:nrow(data.6a),sep="")

#Example data for two communities (incidence raw data)
ciliates <- get(load("Data/ciliates.rda"))
a <- ciliates[[1]] ; b <- ciliates[[2]] ; c <- ciliates[[3]]
t1 <- ncol(a) ; t2 <- ncol(b) ; t3 <- ncol(c)
data1 <- as.integer(rowSums(a)) ; data2 <- as.integer(rowSums(b)) ; data3 <- as.integer(rowSums(c))
data.5c <- data.frame("a" = c(t1, data1), "b" = c(t2, data2))
rownames(data.5c) <- c("T" , paste("y",1:(nrow(data.5c)-1),sep = ""))

data.6c <- data.frame("a" = c(t1, data1), "b" = c(t2, data2), "c" = c(t3, data3))
rownames(data.6c) <- c("T" , paste("y",1:(nrow(data.6c)-1),sep = ""))

data.2c <- load("Data/Data2c.rda")
data.2c <- get(data.2c)

data.2c <- data.frame(data.2c, row.names = NULL)
for(i in 1:ncol(data.2c)){
  data.2c[,i] <- as.integer(data.2c[,i])
}



# Example data for multiple community (incidence data)
# tropical rainforest ants 
data.5b <- read.table("Data/Data5b.txt")
colnames(data.5b) <- c("soil", "fogging", "trap")
data.6b <- data.5b
rownames(data.6b) <- c("T", paste("y",1:(nrow(data.6b)-1),sep=""))
data.5b <- data.5b[,c(1,3)]
rownames(data.5b) <- c("T", paste("y",1:(nrow(data.5b)-1),sep=""))



#data.an <- as.data.frame(y500)
#rownames(data.an) <- c("T", paste("y",1:(nrow(data.an)-1),sep=""))


# Example data for genetic (abundance data)
data.7a <- read.table("Data/Data7a.txt", header=FALSE)
rownames(data.7a) <- paste("x",1:nrow(data.7a),sep="")

data.7b <- read.table("Data/Data7b.txt", header=FALSE)
rownames(data.7b) <- paste("x",1:nrow(data.7b),sep="")

