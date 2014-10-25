library(GEOquery)

# dowload and save dataset from GEO
gse <- "GSE28"
eset <- getGEO(gse)[[1]]
save(eset,file="eset.rda")

# get the sample times times
time <- as.character(pData(eset)$title)
m1 <- regexpr(":", time)+2
m2 <- regexpr("hr",time)-2
time <- mapply(function(x,i,j) as.numeric(substr(x,i,j)), time,m1,m2)

# order samples by time
o <- order(time)
x <- exprs(eset)[,o]

# remove genes with any NAs
hasNas <- rowSums(is.na(x)) > 0
x <- x[!hasNas,]

# center and scale gene data
geneMeans <- rowMeans(x)
geneSds <- matrixStats::rowSds(x)

sx <- x
for (i in seq(len=nrow(x))) {
    sx[i,] <- (x[i,] - geneMeans[i]) / geneSds[i]
}

# plot the first gene
matplot(sx[1,], type="b", lwd=1.4,pch=19,ylab="expression")
abline(v=2,lty=2)

# plot the first ten genes
matplot(t(sx[1:10,]),type="b", lwd=1.4,pch=19,ylab="expression")
abline(v=2,lty=2)

# plot the first 100 genes
matplot(t(sx[1:100,]),type="b", lwd=1.4,pch=19,ylab="expression")
abline(v=2,lty=2)

# make a function to plot genes
plotgene=function(x, ...) {
  matplot(t(x),type="b", lwd=1.4,pch=19,ylab="expression", ...)
  abline(v=2,lty=2)  
}

# let's illustrate furtherst first traversal

# use the first gene as first center
center1 <- sx[1,]
plotgene(rbind(center1), main="center 1")

# this function computes distance between two points
get_distance <- function(x,y) {
  sqrt(sum((x-y)^2))
}

# compute the distance between each point and first center
center1Distance <- apply(sx, 1, function(x) get_distance(x, center1))

# choose the gene that's furthest away as second center
center2Index <- which.max(center1Distance)
center2 <- sx[center2Index,]
plotgene(rbind(center2), main="center 2")

# compute distance between each gene and second center
center2Distance <- apply(sx, 1, function(x) get_distance(x, center2))

# is center1 or center2 closest center
minDistance <- apply(cbind(center1Distance, center2Distance),1,min)

# choose the gene that's furthest away from it's closest center as third center
center3Index <- which.max(minDistance)
center3 <- sx[center3Index,]
plotgene(rbind(center3), main="center 3")

# compute distance to third center
center3Distance <- apply(sx, 1, function(x) get_distance(x, center3))
distances <- cbind(center1Distance, center2Distance, center3Distance)

# assign each gene to cluster based on closest center
clusterIndex <- apply(distances, 1, which.min)
table(clusterIndex)

# plot gene clusters
for (k in seq(len=3)) {
  indx <- which(clusterIndex == k)
  indx <- indx[order(minDistance[indx])][1:50]
  plotgene(sx[indx,], main=paste("cluster", k))
}

# run k-means to find six centers
res=kmeans(sx, centers=6)

# plot the 6 "center" genes
for (k in seq(len=6)) {
  plotgene(rbind(res$centers[k,]), main=paste("center", k))
}

# plot the 6 gene clusters
for (k in seq(len=6)) {
  ind <- which(res$cluster == k)
  plotgene(sx[ind,], main=paste("center", k))
}


#####
# a different example of how EM and fuzzy kmeans works

# we're using fake human height data
mu1 <- 70
mu2 <- 65
sd <- 3

# generate fake human height data
n <- 100
heights <- rnorm(n, mean=c(mu1,mu2), sd=sd)


# plot a histogram of height data
hist(heights, nc=20)
rug(heights)

# mean height is not a good estimate 
my_guess <- mean(heights)
abline(v=my_guess,lwd=2)
mean(sqrt((heights-my_guess)^2))

# can we learn average height for male and female
# without knowing the sex of the people we measured?
# let's start by guessing average female height is smallest height we observed
# and male height is largest height we observed
guesses <- range(heights)
abline(v=guesses, lty=2,col=c("blue","red"),lwd=2)

# compute distance to our two guesses
distances <- sapply(1:2, function(i) abs(heights-guesses[i]))

# plot the distances
plot(heights, distances[,1],pch=19,col="blue")
points(heights, distances[,2], pch=19, col="red")

# this behaves better...
plot(heights, exp(-distances[,1]),pch=19,col="blue")
points(heights, exp(-distances[,2]), pch=19, col="red")

# let's turn that into a probability!
# this is called the E-step in the EM algorithm
z <- exp(-distances) / rowSums(exp(-distances))
plot(heights, z[,1],pch=19,col="blue")
points(heights, z[,2], pch=19, col="red")

# let's guess again but weigh by our probability that each subject is male or
# female (based on our current guess)
# this is called the M-step in the EM algorithm
weighted_guesses <- sapply(1:2, function(i) sum(z[,i] * heights) / sum(z[,i]))

# let's plot again
hist(heights, nc=20)
rug(heights)

my_guess <- mean(heights)
abline(v=my_guess,lwd=1)
mean(sqrt((heights-my_guess)^2))

guesses <- range(heights)
abline(v=guesses, lty=2,col=c("blue","red"),lwd=1)

abline(v=weighted_guesses, col=c("blue","red"), lwd=2, lty=3)

truth=rep(2:1, len=100)
# now we have much better estimates
mean(sqrt((heights-weighted_guesses[truth])^2))
mean(sqrt((heights-c(mu2,mu1)[truth])^2))

