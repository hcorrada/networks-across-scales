tab=read.csv("grades_20141218.csv",stringsAsFactors=FALSE)

tab2=tab[,c(1,2,3,6)]
names(tab2)=tab2[1,]

tab2=tab2[-1,]

ii=tab[1,] %in% c("mid1","mid2","hw1","hw2","hw3","hw4","hw5","hw5_bonus", "r1","r2","r3","r4", "quiz", "part","final")
tab3=tab[,ii]

names(tab3)=tab3[1,]
tab3=tab3[-1,]

theTab=cbind(tab2,tab3)

mid1=as.integer(theTab$mid1)
gmid1=ifelse(mid1 < 27, "F",
             ifelse(mid1 < 35, "D",
             ifelse(mid1 < 37, "D+",
                    ifelse(mid1 < 41, "C-",
                    ifelse(mid1 < 44, "C",
                    ifelse(mid1 < 47, "C+",
                           ifelse(mid1 < 50, "B-",
                           ifelse(mid1 < 56, "B",
                           ifelse(mid1 < 57, "B+", 
                                  ifelse(mid1 < 60, "A-",
                                  ifelse(mid1 < 69, "A", "A+")))))))))))
                                  
theTab$gmid1=factor(gmid1,levels=c("F","D-","D","D+", 
                                   "C-", "C", "C+",
                                   "B-","B","B+",
                                   "A-","A","A+"))

mid2=as.integer(theTab$mid2)
gmid2=ifelse(mid2 < 46, "F",
             ifelse(mid2 < 51, "D",
                    ifelse(mid2 < 56, "C-",
                    ifelse(mid2 < 62, "C",
                    ifelse(mid2 < 66, "C+",
                           ifelse(mid2 < 70, "B-",
                           ifelse(mid2 < 76, "B",        
                           ifelse(mid2 < 81, "B+", 
                                  ifelse(mid2 < 89, "A-",
                                  ifelse(mid2 < 97, "A", "A+"))))))))))                                  
theTab$gmid2=factor(gmid2,levels=c("F","D-","D","D+",
                                   "C-","C","C+",
                                   "B-","B","B+",
                                   "A-","A","A+"))

scoreMap <- c(50,60,65,68,70,75,78,80,85,88,90,95,98)
names(scoreMap) <- levels(theTab$gmid2)

final=as.integer(theTab$final)
sfinal <- ((final-median(final))/mad(final))


theTab$smid1 <- scoreMap[as.character(theTab$gmid1)]
theTab$smid2 <- scoreMap[as.character(theTab$gmid2)]

gfinal=ifelse(final < 50, "F",
              ifelse(final < 60, "D",
                     ifelse(final < 65, "C-",
                     ifelse(final < 70, "C",
                     ifelse(final < 75, "C+",
                            ifelse(final < 80, "B-",
                            ifelse(final < 85, "B",
                            ifelse(final < 89, "B+",
                                   ifelse(final < 95, "A-",
                                   ifelse(final < 100, "A", "A+"))))))))))

theTab$gfinal=factor(gfinal,levels=c("F","D-","D","D+",
  "C-","C","C+",
  "B-","B","B+",
  "A-","A","A+"))

theTab$sfinal=scoreMap[as.character(theTab$gfinal)]


theTab$part <- as.integer(theTab$part) * 20

rosalindScores <- matrix(as.integer(as.matrix(theTab[,c("r1","r2","r3","r4")])),nc=4)
ii <- which(theTab$Dir %in% c("dmstein", "wdefinba", "eshi", "kwang07", "jgzamora","ipersons"))
rosalindScores[ii,1] <- rosalindScores[ii,1] + 35

rosalindPoints <- c(72, 54, 42, 6)
rosalindGrade <- round(100*rowSums(rosalindScores) / sum(rosalindPoints))
theTab$rosalind <- scoreMap[ifelse(rosalindGrade < 50, "C",
                          ifelse(rosalindGrade < 80, "B", "A"))]

hwScores <- matrix(as.integer(as.matrix(theTab[,c("hw1","hw2","hw3","hw4","hw5")])), nc=5)
hwPoints <- c(36, 26, 18, 12, 6)
hwGrade <- round(100*rowSums(hwScores) / sum(hwPoints))
theTab$hw <- hwGrade + as.integer(theTab$hw5_bonus)

theTab$squiz <- round(as.integer(theTab$quiz) / 35 * 100)

ii <- c("smid1","smid2","sfinal","rosalind", "hw","part","squiz")

theScores <- as.matrix(theTab[,ii])
theScores <- matrix(as.integer(theScores),nr=nrow(theScores))

w1 <- c(.10,.10,.15,.15,.25,.1,.15)
theTab$finalGrade <- ceiling(rowSums(sweep(theScores,2,w1,"*")))

cutoffs=c(100,91.9,88,85,80,75.5,70,67,60,45,0)
grades=cut(theTab$finalGrade,breaks=rev(cutoffs))
levels(grades) <- c("F","D","C","C+","B-","B","B+","A-","A","A+")

cutoffs2=c(100,94,87,82,79,75,70,67,60,45,0)
grades2=cut(theTab$finalGrade,breaks=rev(cutoffs2))
levels(grades2) <- c("F","D","C","C+","B-","B","B+","A-","A","A+")


pdf("grades.pdf")
o <- order(theTab$finalGrade)
with(theTab,{
  dotchart(finalGrade[o])
  text(finalGrade[o],seq(along=finalGrade),paste(theTab[o,1],grades2[o]),cex=.5,pos=2)
})
abline(v=cutoffs2,lty=2)
dev.off()

finalTab=data.frame(theTab[,1:4],theTab$finalGrade,grades)
write.csv(finalTab,file="CMSC423_201401_finalGrades.csv")




