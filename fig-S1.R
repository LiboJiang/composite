

load("final-par.RData")
LR1<-final_par[,1]


load("allpar_log.RData")
LR2<-allpar_log[,1]

marker_table_P <- read.table("poplar-marker.txt")[-1,1:2]

snpinfo <- as.matrix(marker_table_P)
scaffoldsum <- length(table(as.character(snpinfo[,1])))
scaffoldc <- sort(as.numeric(names(table(snpinfo[,1]))))

for(i in 1:scaffoldsum){
  snpinfo[which(snpinfo[,1]==scaffoldc[i]),1] <- i
}

snpinfo[which( as.numeric(snpinfo[,1])>19)]<-20


info_2 <- cbind(as.numeric(snpinfo[,1]),as.numeric(snpinfo[,2])/100000,as.numeric(LR1),as.numeric(LR2))


info_1<-rbind(info_2[which(info_2[,1]==1),],info_2[which(info_2[,1]==3),],info_2[which(info_2[,1]==5),],
              info_2[which(info_2[,1]==8),],info_2[which(info_2[,1]==10),],info_2[which(info_2[,1]==11),],
              info_2[which(info_2[,1]==14),],info_2[which(info_2[,1]==16),],info_2[which(info_2[,1]==18),],
              info_2[which(info_2[,1]==19),], info_2[which(info_2[,1]==20),])


AA<-names(table(info_1[,1]))


for (i in 1:length(AA)){ 
  ndx_1 <- which(info_1[, 1]==as.numeric(AA[i]))
  lstMrk_1 <- max(info_1[ndx_1, 2])
  if (i < scaffoldsum_1) ndx2_1 <-which(info_1[, 1]==as.numeric(AA[i+1]))
  if (i < scaffoldsum_1) info_1[ndx2_1, 2] <- info_1[ndx2_1, 2] + lstMrk_1
}


bpMidVec <- vector(length=length(AA))

for (i in 1:length(AA)){
  ndx_1 <- which(info_1[, 1]==as.numeric(AA[i]))
  posSub_1 <- info_1[ndx_1, 2]
  bpMidVec[i] <- ((max(posSub_1) - min(posSub_1))/2) + min(posSub_1)
}

sigsnp1 <- info_1
colnames(sigsnp1) <- c("lg","cm","lr1","lr2")
rownames(sigsnp1) <- NULL
sigsnp1 <- as.data.frame(sigsnp1)

sigsnp1






library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

pdf("fig2.pdf",width=11.2,height =8)


height<-3000
length<-85*2+10
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(-120,height+100)); 



sub_rc<-c(0,3000,85,3000-max(sigsnp1$cm)-60)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="purple4",lwd=1,col="lightyellow1")

nn<-function(LR_len,a,b){
  n1<-which(LR_len[geno_type_st$marker2]> a )
  
  n2<-which(LR_len[geno_type_st$marker3]> b)
  
  n11<-c()
  for(i in 1:length(n1)){
    n11<-c(n11,which(LR_len==LR_len[geno_type_st$marker2][n1[i]]))
  }
  
  n22<-c()
  for(i in 1:length(n2)){
    n22<-c(n22,which(LR_len==LR_len[geno_type_st$marker3][n2[i]]))
  }
  
  n12<-unique(sort(c(n11,n22)))
  n12
}


gene_mat1<-function(LR,AB){
  num<-c()
  for(i in 1:length(AB)){
    num<-c(num,which(sigsnp1$lr1==LR[AB[i]]))
  }
  unique(num)
}

gene_mat2<-function(LR,AB){
  num<-c()
  for(i in 1:length(AB)){
    num<-c(num,which(sigsnp1$lr2==LR[AB[i]]))
  }
  unique(num)
}



for(i in 1:c(length(AA)-1)){
  n1<-max(sigsnp1$cm[which(sigsnp1$lg==as.numeric(AA[i]))])
  segments(sub_rc[1],sub_rc[2]-n1,sub_rc[3],sub_rc[2]-n1,lty=2,col="grey70")
  
}

k=80/320
segments(sub_rc[3]-168*k,sub_rc[2],sub_rc[3]-168*k,sub_rc[4],lwd=2,col="green3",lty=2)
segments(sub_rc[3]-216*k,sub_rc[2],sub_rc[3]-216*k,sub_rc[4],lwd=2,col="green3")




AB1<-unique(c(nn(LR1,168,400),nn(LR2,97.51880,400)))
AB2<-unique(c(nn(LR1,400,216),nn(LR2,400,120.2317)))

plot_num<-gene_mat1(LR1,AB1)

for(i in 1:length(plot_num)){
  points(sub_rc[3]-sigsnp1$lr1[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=1,col="green3")
}

plot_n<-plot_num[3]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,13,adj=0.5)

plot_n<-plot_num[6]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,14,adj=0.5)

plot_n<-plot_num[17]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,19,adj=0.5)

plot_n<-plot_num[19]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,25,adj=0.5)

plot_n<-plot_num[21]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,27,adj=0.5)

plot_num<-gene_mat1(LR1,AB2)

for(i in 1:length(plot_num)){
  points(sub_rc[3]-sigsnp1$lr1[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=3,col="green3")
}

plot_n<-plot_num[1]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,1,adj=0.5)


plot_n<-plot_num[2]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,2,adj=0.5)


plot_n<-plot_num[5]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+1,sub_rc[2]-sigsnp1$cm[plot_n]+40,3,adj=0.5)


plot_n<-plot_num[7]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+1,sub_rc[2]-sigsnp1$cm[plot_n]-40,4,adj=0.5)

plot_n<-plot_num[8]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]+40,5,adj=0.5)

plot_n<-plot_num[11]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]-40,6,adj=0.5)

plot_n<-plot_num[12]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n],7,adj=0.5)


plot_n<-plot_num[14]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,8,adj=0.5)

plot_n<-plot_num[23]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+1,sub_rc[2]-sigsnp1$cm[plot_n]+40,9,adj=0.5)

plot_n<-plot_num[26]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,10,adj=0.5)

plot_n<-plot_num[34]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,11,adj=0.5)

plot_n<-plot_num[35]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,12,adj=0.5)

plot_n<-plot_num[36]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,15,adj=0.5)

plot_n<-plot_num[38]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,16,adj=0.5)

plot_n<-plot_num[40]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]+40,17,adj=0.5)

plot_n<-plot_num[41]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,18,adj=0.5)

plot_n<-plot_num[51]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,20,adj=0.5)

plot_n<-plot_num[52]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,21,adj=0.5)

plot_n<-plot_num[53]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,22,adj=0.5)

plot_n<-plot_num[74]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,23,adj=0.5)

plot_n<-plot_num[83]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,24,adj=0.5)

plot_n<-plot_num[97]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+2,sub_rc[2]-sigsnp1$cm[plot_n]+40,26,adj=0.5)

plot_n<-plot_num[101]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,28,adj=0.5)



for(i in 0:4){
  segments(sub_rc[3]-20*i,sub_rc[4],sub_rc[3]-20*i,sub_rc[4]+20,lwd=1)
  text(sub_rc[3]-20*i,sub_rc[4]-60,80*i,adj=0.5,family="Times New Roman",cex=1.6)
}

text(85/2,sub_rc[4]-180,adj=0.5,expression("LR"),family="Times New Roman",cex=1.6)


text(85/2,sub_rc[2]+100,adj=0.5,expression("Dissection Model"),family="Times New Roman",cex=1.6)




######################################################################################
sub_rc<-c(85,3000,95,3000-max(sigsnp1$cm)-60)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="purple4",lwd=1,col="mediumslateblue")



for(i in 1:c(length(AA)-1)){
  text(sub_rc[1]+5,sub_rc[2]-bpMidVec[i],AA[i], adj = 0.5,cex=1.6)
}


for(i in 11){
  text(sub_rc[1]+5,sub_rc[2]-bpMidVec[i]-30,"scaffold", adj = 0.5,cex=0.9)
}

for(i in 1:c(length(AA)-1)){
  n1<-max(sigsnp1$cm[which(sigsnp1$lg==as.numeric(AA[i]))])
  segments(sub_rc[1],sub_rc[2]-n1,sub_rc[3],sub_rc[2]-n1,col="grey70")
  
}

text(sub_rc[1]+5,sub_rc[2]+100,adj=0.5,expression("Chromosome"),family="Times New Roman",cex=1.6)

######################################################################

sub_rc<-c(95,3000,95+85,3000-max(sigsnp1$cm)-60)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="purple4",lwd=1,col="lightyellow1")


for(i in 1:c(length(AA)-1)){
  n1<-max(sigsnp1$cm[which(sigsnp1$lg==as.numeric(AA[i]))])
  segments(sub_rc[1],sub_rc[2]-n1,sub_rc[3],sub_rc[2]-n1,lty=2,col="grey70")
  
}

k=80/160
segments(sub_rc[1]+97.51880*k,sub_rc[2],sub_rc[1]+97.51880*k,sub_rc[4],lwd=2,col="red",lty=2)
segments(sub_rc[1]+120.2317*k,sub_rc[2],sub_rc[1]+120.2317*k,sub_rc[4],lwd=2,col="red")



AB1<-unique(c(nn(LR1,168,400),nn(LR2,97.51880,400)))
AB2<-unique(c(nn(LR1,400,216),nn(LR2,400,120.2317)))

plot_num<-gene_mat2(LR2,AB1)

for(i in 1:length(plot_num)){
  points(sub_rc[1]+sigsnp1$lr2[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=1,col="red")
}

plot_n<-plot_num[3]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,13,adj=0.5)

plot_n<-plot_num[6]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,14,adj=0.5)

plot_n<-plot_num[17]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,19,adj=0.5)

plot_n<-plot_num[19]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,25,adj=0.5)

plot_n<-plot_num[21]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,27,adj=0.5)

plot_num<-gene_mat2(LR2,AB2)

for(i in 1:length(plot_num)){
  points(sub_rc[1]+sigsnp1$lr2[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=3,col="red")
}

plot_n<-plot_num[1]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,1,adj=0.5)


plot_n<-plot_num[2]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,2,adj=0.5)


plot_n<-plot_num[5]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,3,adj=0.5)


plot_n<-plot_num[7]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,4,adj=0.5)

plot_n<-plot_num[8]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,5,adj=0.5)

plot_n<-plot_num[11]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,6,adj=0.5)

plot_n<-plot_num[12]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,7,adj=0.5)


plot_n<-plot_num[14]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,8,adj=0.5)

plot_n<-plot_num[23]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,9,adj=0.5)

plot_n<-plot_num[26]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,10,adj=0.5)

plot_n<-plot_num[34]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,11,adj=0.5)

plot_n<-plot_num[35]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,12,adj=0.5)

plot_n<-plot_num[36]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,15,adj=0.5)

plot_n<-plot_num[38]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,16,adj=0.5)

plot_n<-plot_num[40]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]-40,17,adj=0.5)

plot_n<-plot_num[41]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,18,adj=0.5)

plot_n<-plot_num[51]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,20,adj=0.5)

plot_n<-plot_num[52]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,21,adj=0.5)

plot_n<-plot_num[53]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,22,adj=0.5)

plot_n<-plot_num[74]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,23,adj=0.5)

plot_n<-plot_num[83]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,24,adj=0.5)

plot_n<-plot_num[97]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,26,adj=0.5)

plot_n<-plot_num[101]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,28,adj=0.5)







for(i in 0:4){
  segments(sub_rc[1]+20*i,sub_rc[4],sub_rc[1]+20*i,sub_rc[4]+20,lwd=1)
  text(sub_rc[1]+20*i,sub_rc[4]-60,40*i,adj=0.5,family="Times New Roman",cex=1.6)
}
text(sub_rc[1]+85/2,sub_rc[4]-180,adj=0.5,expression("LR"),family="Times New Roman",cex=1.6)
text(sub_rc[1]+85/2,sub_rc[2]+100,adj=0.5,expression("Composite Model"),family="Times New Roman",cex=1.6)



gene1<-read.table("gene.txt")
gene1<-as.character(gene1[,3])
AAA<-c(1:28)
for(i in 0:3){
  for(j in 1:7){
    sub_rc<-c(0+45*i,500-100*(j-1),8+45*i,500-100*j)
    rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],col="thistle1",border = "thistle1")
    text(4+45*i,500-100*j+50,AAA[i*7+j],cex=1.2,adj=0.5,family="Times New Roman")
    
  }
}


AAA<-c(1:28)
for(i in 0:3){
  for(j in 1:7){
    sub_rc<-c(8+45*i,500-100*(j-1),45+45*i,500-100*j)
    rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],col="lightcyan",border = "lightcyan")
    text(25+45*i,500-100*j+50,gene1[i*7+j],cex=1.2,adj=0.5,family="Times New Roman")
    
  }
}

dev.off()





library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

jpeg("fig2.jpeg",width=1120,height =720, units = "px", pointsize = 12,quality =150)


height<-3000
length<-85*2+10
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(-120,height+100)); 



sub_rc<-c(0,3000,85,3000-max(sigsnp1$cm)-60)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="purple4",lwd=1,col="lightyellow1")

nn<-function(LR_len,a,b){
  n1<-which(LR_len[geno_type_st$marker2]> a )
  
  n2<-which(LR_len[geno_type_st$marker3]> b)
  
  n11<-c()
  for(i in 1:length(n1)){
    n11<-c(n11,which(LR_len==LR_len[geno_type_st$marker2][n1[i]]))
  }
  
  n22<-c()
  for(i in 1:length(n2)){
    n22<-c(n22,which(LR_len==LR_len[geno_type_st$marker3][n2[i]]))
  }
  
  n12<-unique(sort(c(n11,n22)))
  n12
}


gene_mat1<-function(LR,AB){
  num<-c()
  for(i in 1:length(AB)){
    num<-c(num,which(sigsnp1$lr1==LR[AB[i]]))
  }
  unique(num)
}

gene_mat2<-function(LR,AB){
  num<-c()
  for(i in 1:length(AB)){
    num<-c(num,which(sigsnp1$lr2==LR[AB[i]]))
  }
  unique(num)
}



for(i in 1:c(length(AA)-1)){
  n1<-max(sigsnp1$cm[which(sigsnp1$lg==as.numeric(AA[i]))])
  segments(sub_rc[1],sub_rc[2]-n1,sub_rc[3],sub_rc[2]-n1,lty=2,col="grey70")
  
}

k=80/320
segments(sub_rc[3]-168*k,sub_rc[2],sub_rc[3]-168*k,sub_rc[4],lwd=2,col="green3",lty=2)
segments(sub_rc[3]-216*k,sub_rc[2],sub_rc[3]-216*k,sub_rc[4],lwd=2,col="green3")




AB1<-unique(c(nn(LR1,168,400),nn(LR2,97.51880,400)))
AB2<-unique(c(nn(LR1,400,216),nn(LR2,400,120.2317)))

plot_num<-gene_mat1(LR1,AB1)

for(i in 1:length(plot_num)){
  points(sub_rc[3]-sigsnp1$lr1[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=1,col="green3")
}

plot_n<-plot_num[3]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,13,adj=0.5)

plot_n<-plot_num[6]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,14,adj=0.5)

plot_n<-plot_num[17]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,19,adj=0.5)

plot_n<-plot_num[19]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,25,adj=0.5)

plot_n<-plot_num[21]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,27,adj=0.5)

plot_num<-gene_mat1(LR1,AB2)

for(i in 1:length(plot_num)){
  points(sub_rc[3]-sigsnp1$lr1[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=3,col="green3")
}

plot_n<-plot_num[1]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,1,adj=0.5)


plot_n<-plot_num[2]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,2,adj=0.5)


plot_n<-plot_num[5]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+1,sub_rc[2]-sigsnp1$cm[plot_n]+40,3,adj=0.5)


plot_n<-plot_num[7]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+1,sub_rc[2]-sigsnp1$cm[plot_n]-40,4,adj=0.5)

plot_n<-plot_num[8]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]+40,5,adj=0.5)

plot_n<-plot_num[11]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]-40,6,adj=0.5)

plot_n<-plot_num[12]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n],7,adj=0.5)


plot_n<-plot_num[14]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,8,adj=0.5)

plot_n<-plot_num[23]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+1,sub_rc[2]-sigsnp1$cm[plot_n]+40,9,adj=0.5)

plot_n<-plot_num[26]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,10,adj=0.5)

plot_n<-plot_num[34]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,11,adj=0.5)

plot_n<-plot_num[35]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,12,adj=0.5)

plot_n<-plot_num[36]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,15,adj=0.5)

plot_n<-plot_num[38]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,16,adj=0.5)

plot_n<-plot_num[40]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]+40,17,adj=0.5)

plot_n<-plot_num[41]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,18,adj=0.5)

plot_n<-plot_num[51]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,20,adj=0.5)

plot_n<-plot_num[52]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,21,adj=0.5)

plot_n<-plot_num[53]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,22,adj=0.5)

plot_n<-plot_num[74]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,23,adj=0.5)

plot_n<-plot_num[83]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,24,adj=0.5)

plot_n<-plot_num[97]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k+2,sub_rc[2]-sigsnp1$cm[plot_n]+40,26,adj=0.5)

plot_n<-plot_num[101]
text(sub_rc[3]-sigsnp1$lr1[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,28,adj=0.5)



for(i in 0:4){
  segments(sub_rc[3]-20*i,sub_rc[4],sub_rc[3]-20*i,sub_rc[4]+20,lwd=1)
  text(sub_rc[3]-20*i,sub_rc[4]-60,80*i,adj=0.5,family="Times New Roman",cex=1.6)
}

text(85/2,sub_rc[4]-180,adj=0.5,expression("LR"),family="Times New Roman",cex=1.6)


text(85/2,sub_rc[2]+100,adj=0.5,expression("Dissection Model"),family="Times New Roman",cex=1.6)




######################################################################################
sub_rc<-c(85,3000,95,3000-max(sigsnp1$cm)-60)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="purple4",lwd=1,col="mediumslateblue")



for(i in 1:c(length(AA)-1)){
  text(sub_rc[1]+5,sub_rc[2]-bpMidVec[i],AA[i], adj = 0.5,cex=1.6)
}


for(i in 11){
  text(sub_rc[1]+5,sub_rc[2]-bpMidVec[i]-30,"scaffold", adj = 0.5,cex=0.9)
}

for(i in 1:c(length(AA)-1)){
  n1<-max(sigsnp1$cm[which(sigsnp1$lg==as.numeric(AA[i]))])
  segments(sub_rc[1],sub_rc[2]-n1,sub_rc[3],sub_rc[2]-n1,col="grey70")
  
}

text(sub_rc[1]+5,sub_rc[2]+100,adj=0.5,expression("Chromosome"),family="Times New Roman",cex=1.6)

######################################################################

sub_rc<-c(95,3000,95+85,3000-max(sigsnp1$cm)-60)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="purple4",lwd=1,col="lightyellow1")


for(i in 1:c(length(AA)-1)){
  n1<-max(sigsnp1$cm[which(sigsnp1$lg==as.numeric(AA[i]))])
  segments(sub_rc[1],sub_rc[2]-n1,sub_rc[3],sub_rc[2]-n1,lty=2,col="grey70")
  
}

k=80/160
segments(sub_rc[1]+97.51880*k,sub_rc[2],sub_rc[1]+97.51880*k,sub_rc[4],lwd=2,col="red",lty=2)
segments(sub_rc[1]+120.2317*k,sub_rc[2],sub_rc[1]+120.2317*k,sub_rc[4],lwd=2,col="red")



AB1<-unique(c(nn(LR1,168,400),nn(LR2,97.51880,400)))
AB2<-unique(c(nn(LR1,400,216),nn(LR2,400,120.2317)))

plot_num<-gene_mat2(LR2,AB1)

for(i in 1:length(plot_num)){
  points(sub_rc[1]+sigsnp1$lr2[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=1,col="red")
}

plot_n<-plot_num[3]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,13,adj=0.5)

plot_n<-plot_num[6]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,14,adj=0.5)

plot_n<-plot_num[17]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,19,adj=0.5)

plot_n<-plot_num[19]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,25,adj=0.5)

plot_n<-plot_num[21]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,27,adj=0.5)

plot_num<-gene_mat2(LR2,AB2)

for(i in 1:length(plot_num)){
  points(sub_rc[1]+sigsnp1$lr2[plot_num[i]]*k,sub_rc[2]-sigsnp1$cm[plot_num[i]],pch=3,col="red")
}

plot_n<-plot_num[1]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,1,adj=0.5)


plot_n<-plot_num[2]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,2,adj=0.5)


plot_n<-plot_num[5]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,3,adj=0.5)


plot_n<-plot_num[7]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,4,adj=0.5)

plot_n<-plot_num[8]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,5,adj=0.5)

plot_n<-plot_num[11]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,6,adj=0.5)

plot_n<-plot_num[12]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,7,adj=0.5)


plot_n<-plot_num[14]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,8,adj=0.5)

plot_n<-plot_num[23]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,9,adj=0.5)

plot_n<-plot_num[26]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,10,adj=0.5)

plot_n<-plot_num[34]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,11,adj=0.5)

plot_n<-plot_num[35]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,12,adj=0.5)

plot_n<-plot_num[36]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,15,adj=0.5)

plot_n<-plot_num[38]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,16,adj=0.5)

plot_n<-plot_num[40]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k-1,sub_rc[2]-sigsnp1$cm[plot_n]-40,17,adj=0.5)

plot_n<-plot_num[41]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,18,adj=0.5)

plot_n<-plot_num[51]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,20,adj=0.5)

plot_n<-plot_num[52]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,21,adj=0.5)

plot_n<-plot_num[53]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,22,adj=0.5)

plot_n<-plot_num[74]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,23,adj=0.5)

plot_n<-plot_num[83]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]-40,24,adj=0.5)

plot_n<-plot_num[97]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,26,adj=0.5)

plot_n<-plot_num[101]
text(sub_rc[1]+sigsnp1$lr2[plot_n]*k,sub_rc[2]-sigsnp1$cm[plot_n]+40,28,adj=0.5)







for(i in 0:4){
  segments(sub_rc[1]+20*i,sub_rc[4],sub_rc[1]+20*i,sub_rc[4]+20,lwd=1)
  text(sub_rc[1]+20*i,sub_rc[4]-60,40*i,adj=0.5,family="Times New Roman",cex=1.6)
}
text(sub_rc[1]+85/2,sub_rc[4]-180,adj=0.5,expression("LR"),family="Times New Roman",cex=1.6)
text(sub_rc[1]+85/2,sub_rc[2]+100,adj=0.5,expression("Composite Model"),family="Times New Roman",cex=1.6)



gene1<-read.table("gene.txt")
gene1<-as.character(gene1[,3])
AAA<-c(1:28)
for(i in 0:3){
  for(j in 1:7){
    sub_rc<-c(0+45*i,500-100*(j-1),8+45*i,500-100*j)
    rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],col="thistle1",border = "thistle1")
    text(4+45*i,500-100*j+50,AAA[i*7+j],cex=1.2,adj=0.5,family="Times New Roman")
    
  }
}


AAA<-c(1:28)
for(i in 0:3){
  for(j in 1:7){
    sub_rc<-c(8+45*i,500-100*(j-1),45+45*i,500-100*j)
    rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],col="lightcyan",border = "lightcyan")
    text(25+45*i,500-100*j+50,gene1[i*7+j],cex=1.2,adj=0.5,family="Times New Roman")
    
  }
}

dev.off()
























