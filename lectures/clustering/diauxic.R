library(GEOquery)

gse <- "GSE28"
eset <- getGEO(gse)[[1]]
save(eset,file="eset.rda")

# get the times
time <- as.character(pData(eset)$title)
m1 <- regexpr(":", time)+2
m2 <- regexpr("hr",time)-2
time <- mapply(function(x,i,j) as.numeric(substr(x,i,j)), time,m1,m2)

o <- order(time)
x <- exprs(eset)[,o]

# remove genes with any NAs
hasNas <- rowSums(is.na(x)) > 0
x <- x[!hasNas,]

geneMeans <- rowMeans(x)
geneSds <- matrixStats::rowSds(x)

sx <- x
for (i in seq(len=nrow(x))) {
    sx[i,] <- (x[i,] - geneMeans[i]) / geneSds[i]
}

matplot(sx[1,], type="b", lwd=1.4,pch=19,ylab="expression")
abline(v=2,lty=2)

matplot(t(sx[1:10,]),type="b", lwd=1.4,pch=19,ylab="expression")
abline(v=2,lty=2)

matplot(t(sx[1:100,]),type="b", lwd=1.4,pch=19,ylab="expression")
abline(v=2,lty=2)

plotgene=function(x, ...) {
  matplot(t(x),type="b", lwd=1.4,pch=19,ylab="expression", ...)
  abline(v=2,lty=2)  
}

center1 <- sx[1,]
plotgene(rbind(center1), main="center 1")

get_distance <- function(x,y) {
  sqrt(sum((x-y)^2))
}

center1Distance <- apply(sx, 1, function(x) get_distance(x, center1))
center2Index <- which.max(center1Distance)
center2 <- sx[center2Index,]
plotgene(rbind(center2), main="center 2")

center2Distance <- apply(sx, 1, function(x) get_distance(x, center2))
minDistance <- apply(cbind(center1Distance, center2Distance),1,min)
center3Index <- which.max(minDistance)
center3 <- sx[center3Index,]
plotgene(rbind(center3), main="center 3")

center3Distance <- apply(sx, 1, function(x) get_distance(x, center3))
distances <- cbind(center1Distance, center2Distance, center3Distance)
clusterIndex <- apply(distances, 1, which.min)
table(clusterIndex)

for (k in seq(len=3)) {
  indx <- which(clusterIndex == k)
  indx <- indx[order(minDistance[indx])][1:50]
  plotgene(sx[indx,], main=paste("cluster", k))
}

res=kmeans(sx, centers=6)
for (k in seq(len=6)) {
  plotgene(rbind(res$centers[k,]), main=paste("center", k))
}

for (k in seq(len=6)) {
  ind <- which(res$cluster == k)
  plotgene(sx[ind,], main=paste("center", k))
}


#####
mu1 <- 70
mu2 <- 65
sd <- 3

n <- 100
heights <- rnorm(n, mean=c(mu1,mu2), sd=sd)
hist(heights, nc=20)
rug(heights)

my_guess <- mean(heights)
abline(v=my_guess)
mean(sqrt((heights-my_guess)^2))

guesses <- c(60, 73)
abline(v=guesses, lty=2)
distances <- sapply(1:2, function(i) abs(heights-guesses[i]))

z <- exp(-distances) / rowSums(exp(-distances))
weighted_guesses <- sapply(1:2, function(i) sum(z[,i] * heights) / sum(z[,i]))
abline(v=weighted_guesses, lty=2, col="blue", lwd=1.5)

mean(sqrt((heights-rep(rev(weighted_guesses),len=n))^2))
