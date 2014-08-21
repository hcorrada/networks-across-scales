tab=read.csv("grades_20131218.csv",stringsAsFactors=FALSE)

tab2=tab[,c(1,2,3,6)]
names(tab2)=tab2[1,]

tab2=tab2[-1,]

ii=tab[1,] %in% c("mid1","mid2","p1","p2","p3","p4","hw","participation","final","ext")
tab3=tab[,ii]

names(tab3)=tab3[1,]
tab3=tab3[-1,]

theTab=cbind(tab2,tab3)

mid1=as.integer(theTab$mid1)
smid1=((mid1-median(mid1))/mad(mid1))
gmid1=ifelse(smid1< -1, "D",
  ifelse(smid1< -0.5, "C",
         ifelse(smid1 < 0, "B-",
                ifelse(smid1< 0.5,"B",
                       ifelse(smid1<0.75,"B+",
                              ifelse(smid1<1,"A-",
                                     ifelse(smid1 < 1.2,"A", "A+")))))))

cutoffs=c(-1,-0.5,0,0.5,0.75,1,1.2)
cutoffs*mad(mid1) + median(mid1)

                             
theTab$gmid1=factor(gmid1,levels=c("D","C","B-","B","B+","A-","A","A+"))

mid2=as.integer(theTab$mid2)
smid2=((mid2-median(mid2))/mad(mid2))

gmid2=ifelse(smid2< -3, "F",
  ifelse(smid2 < -2, "D",
  ifelse(smid2< -1, "C",
         ifelse(smid2 < -0.5, "B-",
                ifelse(smid2< 0,"B",
                       ifelse(smid2<0.5,"B+",
                              ifelse(smid2<1,"A-",
                                     ifelse(smid2 < 1.5,"A", "A+"))))))))

theTab$gmid2=factor(gmid2,levels=c("F","D","C","B-","B","B+","A-","A","A+"))

cutoffs=c(-3,-2,-1,-0.5,0,0.5,1,1.5)
cutoffs*mad(mid2) + median(mid2)

scoreMap <- c(50,65,75,80,85,88,90,95,98)
names(scoreMap) <- levels(theTab$gmid2)

theTab$smid1 <- scoreMap[as.character(theTab$gmid1)]
theTab$smid2 <- scoreMap[as.character(theTab$gmid2)]

theTab$part <- as.integer(theTab$part) * 20
theTab$p3 <- as.integer(theTab$p3) / 40 * 100

final=as.integer(theTab$final)
sfinal <- ((final-median(final))/mad(final))

gfinal=ifelse(sfinal< -3, "F",
  ifelse(sfinal < -1.5, "D",
  ifelse(sfinal< -.5, "C",
         ifelse(sfinal < -0.25, "B-",
                ifelse(sfinal< 0.25,"B",
                       ifelse(sfinal<0.5,"B+",
                              ifelse(sfinal<.75,"A-",
                                     ifelse(sfinal < 1.25,"A", "A+"))))))))
theTab$gfinal=factor(gfinal,levels=c("F","D","C","B-","B","B+","A-","A","A+"))
theTab$sfinal=scoreMap[as.character(theTab$gfinal)]

# fix Sal's by hand
theTab$p4[1] <- 73

# fix P. Liu's by hand
ii <- which(theTab$Dir == "pliu416")
theTab$p4[ii] <- 58

# fix J. Wang's
ii <- which(theTab$Dir == "jxwang")
theTab$p3[ii] <- 100

ii <- c("smid1","smid2","sfinal","p1","p2","p3","p4","hw","part","ext")

theScores <- as.matrix(theTab[,ii])
theScores <- matrix(as.integer(theScores),nr=nrow(theScores))


w1 <- c(.15,.15,.2,.1,.1,.05,.05,.15,.05)
theTab$finalGrade1 <- rowSums(sweep(theScores[,-10],2,w1,"*"))/sum(w1) + theScores[,10]

w2 <- c(.125,.125,.15,.1,.1,.1,.1,.15,.05)
theTab$finalGrade2 <- rowSums(sweep(theScores[,-10],2,w2,"*"))/sum(w2) + theScores[,10]

theTab$finalGrade <- pmax(theTab$finalGrade1,theTab$finalGrade2)

cutoffs=c(100,91.9,88,85,80,75.5,70,67,60,45,0)
grades=cut(theTab$finalGrade,breaks=rev(cutoffs))
levels(grades) <- c("F","D","C","C+","B-","B","B+","A-","A","A+")

pdf("grades.pdf")
o <- order(theTab$finalGrade)
with(theTab,{
  dotchart(finalGrade[o])
  text(finalGrade[o],seq(along=finalGrade),paste(theTab[o,1],grades[o]),cex=.5,pos=2)
})
abline(v=cutoffs,lty=2)
dev.off()

finalTab=data.frame(theTab[,1:4],theTab$finalGrade,grades)
write.csv(finalTab,file="CMSC423_201301_finalGrades.csv")

bigTab=data.frame(theTab[,1:4],theTab$finalGrade,grades,theTab[,c("gmid1","gmid2","gfinal","p1","p2","p3","p4","hw","part","ext")])
write.csv(bigTab,file="CMSC423_201301_allGrades.csv")




####
hwTab=read.csv("hw_grades.csv",stringsAsFactors=FALSE)
hwCols=6:15
hwScores=as.matrix(hwTab[,hwCols])

# count the number not submitted
nMissed <- rowSums(is.na(hwScores))

for (i in 2:nrow(hwScores)) {
  isMissed <- is.na(hwScores[i,])

  if (nMissed[i] < 2) {
    # find the worst scores
    scores <- hwScores[i,!isMissed] / hwScores[1,!isMissed]
    nDrop <- 2-nMissed[i]
    idxDrop <- which(!isMissed)[order(scores)[1:nDrop]]
    hwScores[i,idxDrop] <- NA
    next
  }

  iMissed <- which(isMissed)
  # set all but the two highest weight missing to 0
  iSkip <- order(hwScores[1,iMissed],decreasing=FALSE)[1:2]
  hwScores[i,iMissed[-iSkip]] <- 0
}

hwPct <- sapply(2:nrow(hwScores), function(i) {
  sum(hwScores[i, !is.na(hwScores[i,])]) / sum(hwScores[1,!is.na(hwScores[i,])])
})

hwGrade=cbind(hwTab[-1,c(1,3,4)], hwPct=round(100*(hwPct)))

upGrade=data.frame(id=hwGrade[,3],ass="hw",hwGrade[,4],com="")
write.table(upGrade,file="hw_upload.csv",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)
