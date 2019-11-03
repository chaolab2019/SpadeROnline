set.seed(123)
library(foreach)
library(msm)
R <- 100
n <- round(seq(100, 5000, length.out = R))
df <- foreach(i=1:R, .combine = rbind, .errorhandling='pass')%do%{
  Fun <- function(j){
  rate <- rtnorm(1, mean=0.5, sd=0.2, lower=0, upper=1) # runif(1, 0, 1) # Shared rate
  total_Sp <- sample(100:1000, 1)
  S12 <- round(rate*total_Sp) + 1
  V1_Sp <- sample((total_Sp-S12), 1)
  V2_Sp <- total_Sp - S12 - V1_Sp

  mat <- cbind(V1=c(rep(1,S12), rep(1,V1_Sp), rep(0, V2_Sp)),
               V2=c(rep(1,S12), rep(0,V1_Sp), rep(1, V2_Sp)))
 
  x <- mat*runif(total_Sp*2)
  B <- 1000
  V1 <- rmultinom(B, n[i], x[,1])
  V2 <- rmultinom(B, n[i], x[,2])
  arr <- array(0, dim = c(total_Sp,2,B))
  arr[,1,] <- V1
  arr[,2,] <- V2
  out <- t(apply(arr, 3, function(x){
    Sg <- sum(rowSums(x)>0)
    Sa <- mean(colSums(x>0))
    Sb <- Sg/Sa
    c(Sg, Sa, Sb, i)
  }))
  out
  }
  #do.call(rbind, lapply(1:10, Fun))
  Fun(1)
}
plot(df[,c(1,3)])
df <- data.frame(df)
colnames(df) <- c("gamma", "alpha", "beta", "replication")

library(ggplot2)

g <- ggplot(df) + aes(x=gamma, y=beta) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#FEFEFF", blues9))(256), 
                       name="density^0.25") +
  ylim(c(1,2)) + theme(legend.position="bottom") +
  ggtitle("Model 1")
g + theme(text=element_text(size=18), legend.text=element_text(size=8)) + ggtitle("Model 1")

ggsave("DivInd-1.png", width = 7, height = 7)



head(df)
ggplot(df) + aes(x=alpha, y=beta) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#FEFEFF", blues9))(256), 
                       name="density^0.25") +
  theme(legend.position="bottom") +
  ggtitle("Model 2") + theme(text=element_text(size=18), legend.text=element_text(size=8))


ggplot(df) + aes(x=gamma, y=beta) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#FEFEFF", blues9))(256), 
                       name="density^0.25") +
  theme(legend.position="bottom") +
  ggtitle("Model 2") + theme(text=element_text(size=18), legend.text=element_text(size=8))

library(dplyr)
df.sub <- filter(df, replication%in%sample(R,3))
ggplot(df) + aes(x=gamma, y=beta) + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) + 
#  scale_fill_gradientn(colours = colorRampPalette(c("#FEFEFF", blues9))(256), 
#                       name="density^0.25") +
  geom_point(aes(colour=factor(replication)), data=df.sub, alpha=0.5)+
  theme(legend.position="bottom") + 
  ggtitle("Model 2") + theme(text=element_text(size=18), legend.text=element_text(size=8))
