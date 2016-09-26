library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
res<-data.frame()
std<-data.frame()
tab<-read.table(args[1],header=TRUE)
outTab<-tab[,c(1:7,9)]
L<-as.numeric(as.character(args[2]))
Lfragm<-as.numeric(as.character(args[3]))
GS<-as.numeric(as.character(args[4]))

pdf("binding_sites.pdf",width=10,height=3.5)
op<-par(ps=10)

makeTransparent<-function(someColor, alpha=40){
newColor<-col2rgb(someColor)
apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
#col<-sample(colors(),1)
colfunc<-colorRampPalette(c("lightblue", "red"))
col<-colfunc(101)[(GS*100)+1]

h<-0.4

layout(t(1:2), widths=c(6,1))

plot(c(1,L),c(0,0),lwd=15,type="l",ylim=c(-1,1),col="grey",yaxt='n', xaxt="n", ann=FALSE,bty="n")
axis(1, at = round(seq(1, L, length.out= 10)), las=2) # rounded position given number of intervals
p <- par('usr')
text(mean(p[1:2]), p[3]+0.5, labels = 'RNA sequence',xpd=TRUE,font=2)

bindingsites<-data.frame()
bindstsTop3<-data.frame()

res<-tab[which(tab$RBDs==1 & tab$percentileRBDs>=0.979),c(1:4,9)]
if(dim(res)[1]==0){
res<-tail(tab[order(tab$normScore),c(1:4,9)],1)
}
res$V1<-rep("prot",dim(res)[1])
res$V2<-rep("rna",dim(res)[1])
tab<-res[,c(6,7,1:5)]
colnames(tab)<-c("V1","V2","V3","V4","V5","V6","V7")
coord<-unique(tab[,5:6])
coord$V7<-1:nrow(coord)
#x<-unique(sort(setDT(coord)[,V5:V6,by=V7]$V1))
df<-data.frame()
for(i in 1:dim(coord)[1]){
    s<-seq(coord[i,1],coord[i,2])
    df<-rbind(df,data.frame(s,rep(coord[i,3],length(s))))
}
colnames(df)<-c("V1","V2")
x<-unique(sort(df$V1))
coord<-data.frame(coord)
df<-data.frame("r"=coord[1,1],"V2"=coord[1,2]);r<-x[1];for(i in 1:(length(x)-1)){if(x[i+1]!=x[i]+1){df<-rbind(df,cbind(r,x[i])); r<-x[i+1]} }
coord<-df
if(nrow(coord)>1){coord<-coord[-1,]}

# keep the top 3 scoring fragment in contiguous regions

newcoord<-data.frame()
for(i in 1:nrow(coord)){
c<-tab[which(tab$V5>=coord[i,1] & tab$V6<=coord[i,2]),]
newcoord<-rbind(newcoord,c[which(c$V7%in%tail(sort(c$V7),3)),5:6])
}

# drop fragments that are more distant from the mass center (~5 times RNA fragment length)

mc<-mean(c(newcoord[1,1],newcoord[nrow(newcoord),2]))
dist<-abs(apply(newcoord,1,mean)-mc)
l<-unique(apply(newcoord,1,diff)) # RNA fragms length
n<-floor((Lfragm*5)/100)*100
std<-rbind(std,data.frame("V1"="rna","V2"=sd(dist)))

if(max(dist)>n){
newcoord<-newcoord[-which(dist>n),]
}

var<-sd(dist)*1.8
if(var>n | is.na(var)){
for(i in 1:dim(coord)[1]){
x1<-rep(coord[i,1])
x2<-rep(coord[i,2])
lines(c(x1,x2),rep(h,2),lwd=7,col=makeTransparent(col))
bindingsites<-rbind(bindingsites,data.frame("V1"="rna","V2"=x1,"V3"=x2))
}
}

coord<-unique(newcoord)

for(i in 1:dim(coord)[1]){
x1<-rep(coord[i,1])
x2<-rep(coord[i,2])
lines(c(x1,x2),rep(h,2),lwd=8,col=col)
bindingsites<-rbind(bindingsites,data.frame("V1"="rna","V2"=x1,"V3"=x2))
}

image(1, 0:1001/1000, t(seq_along(0:1001)), main=substitute(italic('Global Score')), col=colfunc(1001), axes=FALSE, ylim=c(0,1), ylab="", xlab="")
axis(2,col="grey",las=1)

dev.off()

bindingsites<-unique(bindingsites)
write.table(bindingsites,"bindingsites.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

newres<-res[res$startRNA%in%bindingsites$V2 & res$endRNA%in%bindingsites$V3,]
newres<-newres[order(newres$normScore,decreasing=TRUE),-c(6,7)]
write.table(newres,"smallTable.txt",quote=FALSE,col.names=TRUE,row.names=FALSE) # RBDs and select RNA fragments

#res<-res[order(res$normScore,decreasing=TRUE),-c(6,7)]
#write.table(res,"bigTable.txt",quote=FALSE,col.names=FALSE,row.names=FALSE) # RBDs and all RNA fragments
write.table(outTab,"bigTable.txt",quote=FALSE,col.names=TRUE,row.names=FALSE) # all the fragments with no selections at all

#bindstsTop3<-data.frame()
#mydf<-data.frame()
#	for(n in 1:dim(bindsts)[1]){
#	a<-tab[which(tab$V5>=bindsts$V2[n] & tab$V6<=bindsts$V3[n]),]
#	mydf<-rbind(mydf,tail(a[order(a$V7),],1))
#}
#bindstsTop3<-rbind(bindstsTop3,head(unique(mydf[order(mydf$V7,decreasing=TRUE),c(1,5,6)]),3)) # top3 scoring RNA regions


