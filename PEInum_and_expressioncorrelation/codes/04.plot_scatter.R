## library(DescTools)

setwd('/Users/hsl/Desktop/forJinplot/addHiC/finalNEEcodes/filtered_resultsmoreE/')
data1<-read.table(file='housekeeping_forplot11.SAT.all.txt',sep='\t',header=F)
data2<-read.table(file='housekeeping_forplot12.SAT.all.txt',sep='\t',header=F)

#===================== new method 1 start ====================================


pdf("/Users/hsl/Desktop/forJinplot/addHiC/finalNEEcodes/moreEnumandControl_cor.correct.filtered.pdf",w=6,h=6)
result1<-lm(data1[,2]~log(as.numeric(data1[,1])-54))

b1<-as.numeric(result1[[1]][1])
a1<-as.numeric(result1[[1]][2])

myx1<-1:(94-54)
myy1<-log(myx1)*a1+b1


myx1_new<-(1:(94-54)+54)

plot(data1[,1:2],pch=20,col="orange",ylab="Gene expression correlation",xlab="Divergence time",ylim=c(0.65,0.8))
points(data2[,1:2],pch=20,col="gray")

lines(myx1_new,myy1,col="orange",lwd=3)



result2<-lm(data2[,2]~log(as.numeric(data2[,1])-54))

b2<-as.numeric(result2[[1]][1])
a2<-as.numeric(result2[[1]][2])

myx2<-1:(94-54)
myy2<-log(myx1)*a2+b2


myx2_new<-(1:(94-54)+54)

lines(myx2_new,myy2,col="gray",lwd=3)

legend("top",legend=c("With more promoters","Control"),col=c("Orange","Gray"),lty=1,pch=-1,lwd=3,bty="n")

dev.off()


#===================== new method 1 end ====================================

#===================== new method 2 start =================================

result1<-lm(data1[,2]~log(as.numeric(data1[,1])))

b1<-as.numeric(result1[[1]][1])
a1<-as.numeric(result1[[1]][2])

myx1<-1:(94)
myy1<-log(myx1)*a1+b1


myx1_new<-(1:(94))

plot(data1[,1:2],pch=20,col="orange",ylab="Gene expression correlation",xlab="Divergence time",ylim=c(0.63,0.83))
points(data2[,1:2],pch=20,col="gray")

lines(myx1_new,myy1,col="orange",lwd=3)



result2<-lm(data2[,2]~log(as.numeric(data2[,1])-54))

b2<-as.numeric(result2[[1]][1])
a2<-as.numeric(result2[[1]][2])

myx2<-1:(94)
myy2<-log(myx1)*a2+b2


myx2_new<-(1:(94))

lines(myx2_new,myy2,col="gray",lwd=3)

#===================== new method 2 end =====================================




#============ old codes start =================================================

##data3<-read.table(file='housekeeping_forplot21.SAT.all.txt',sep='\t',header=F)
## data4<-read.table(file='housekeeping_forplot22.SAT.all.txt',sep='\t',header=F)

plot(data1[,1:2],pch=20,col="orange",ylab="Gene expression correlation",xlab="Divergence time",ylim=c(0.6,0.85))
points(data2[,1:2],pch=20,col="gray")
my1<-lowess(data1[,1],data1[,2])
lines(my1,col="orange")

## lines(data3,col="orange")

my2<-lowess(data2)
lines(my2,col="gray")
## lines(data4,col="gray")

#============ old codes end =================================================
