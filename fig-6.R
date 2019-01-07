



library(showtext)
showtext_auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")



pdf("3.pdf",width=10,height =7.2)




height<-250
length<-122
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(-2,length),ylim=c(0,height)); 

a<-40
b<-90


sub_rc<-c(6,70+b/2,34.5,10+b/2)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6,col="#FEE8C8")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}
#max(a_F,a_F_Y,a_F_A)
#min(a_F,a_F_Y,a_F_A)

k=0.01
seq<- 58/k/14
n<-6
for(i in -3:4){
  segments(sub_rc[1],sub_rc[4]+2+ 2*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,2*k*i, cex=1,adj=1,font=1)
}
segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)


lines(seq(1,24,0.01)+sub_rc[1]+1,a_F*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,a_F_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,a_F_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)

text(sub_rc[1]+25+0.4,a_F[2301]*seq+k*seq*n+sub_rc[4]+2,expression(a["F"]),cex=1,adj=0,family="Times New Roman")


text(-4,130,expression("Genetic Effect"~(cm^3) ),cex=1.2,srt=90,family="Times New Roman")

text(123,130,expression("Genetic Effect"~(cm^3)),cex=1.2,srt=-90,family="Times New Roman")

text(sub_rc[1]+2,sub_rc[2]-5,expression(B),cex=1,adj=0,family="Times New Roman")

######################################################
sub_rc<-c(6+a,70+b/2*4,34.5+a,10+b/2*4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6,col="#FEE8C8")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}

#max(a_H,a_H_Y,a_H_A)
#min(a_H,a_H_Y,a_H_A)

k=0.1
seq<- 58/k/8
n<-6
for(i in -3:1){
  segments(sub_rc[1],sub_rc[4]+2+ 2*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,2*k*i, cex=1,adj=1,font=1)
}
segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)

lines(seq(1,24,0.01)+sub_rc[1]+1,a_H*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,a_H_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,a_H_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)

text(sub_rc[1]+25+0.4,a_H[2301]*seq+k*seq*n+sub_rc[4]+2,expression(a["H"]),cex=1,adj=0,family="Times New Roman")

text(sub_rc[1]+2,sub_rc[2]-5,expression(A),cex=1,adj=0,family="Times New Roman")
###############################################################################
sub_rc<-c(6+2*a,70+b/2,34.5+2*a,10+b/2)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6,col="#FEE8C8")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}

#max(a_D,a_D_Y,a_D_A,d_D,d_D_Y,d_D_A)
#min(a_D,a_D_Y,a_D_A,d_D,d_D_Y,d_D_A)

k=0.1
seq<- 58/k/8
n<-1.5
for(i in 0:3){
  segments(sub_rc[3],sub_rc[4]+2+ 2*k*seq*i+k*seq*n,sub_rc[3]-0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[3]+0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,2*k*i, cex=1,adj=0,font=1)
}

segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)

lines(seq(1,24,0.01)+sub_rc[1]+1,a_D*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,a_D_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,a_D_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)

lines(seq(1,24,0.01)+sub_rc[1]+1,d_D*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2,lty=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,d_D_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2,lty=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,d_D_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2,lty=2)


text(sub_rc[1]+25+0.4,a_D[2301]*seq+k*seq*n+sub_rc[4]+2,expression(a["D"]),cex=1,adj=0,family="Times New Roman")
text(sub_rc[1]+25+0.4,d_D[2301]*seq+k*seq*n+sub_rc[4]+2,expression(d["D"]),cex=1,adj=0,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-5,expression(C),cex=1,adj=0,family="Times New Roman")
############################################################################################

sub_rc<-c(6,70+b/2*3,34.5,10+b/2*3)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}

#max(i_FH,i_FH_Y,i_FH_A)
#min(i_FH,i_FH_Y,i_FH_A)


k=0.001
seq<- 58/k/5
n<-2
for(i in -2:3){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,1*k*i, cex=1,adj=1,font=1)
}

segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)

lines(seq(1,24,0.01)+sub_rc[1]+1,i_FH*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,i_FH_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,i_FH_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)


text(sub_rc[1]+25+0.4,i_FH[2301]*seq+k*seq*n+sub_rc[4]+2,expression(i["FH"]),cex=1,adj=0,family="Times New Roman")


text(sub_rc[1]+2,sub_rc[2]-5,expression(AB),cex=1,adj=0,family="Times New Roman")

########################################################################
sub_rc<-c(6+a,70,34.5+a,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}

#max(i_FD,i_FD_Y,i_FD_A,j_FD,j_FD_Y,j_FD_A)
#min(i_FD,i_FD_Y,i_FD_A,j_FD,j_FD_Y,j_FD_A)


k=0.005
seq<- 58/k/5
n<-2.4
for(i in -1:1){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,1*k*i, cex=1,adj=1,font=1)
}

for(i in -2){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,"-0.010", cex=1,adj=1,font=1)
}
for(i in 2){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,"0.010", cex=1,adj=1,font=1)
}

segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)

lines(seq(1,24,0.01)+sub_rc[1]+1,i_FD*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,i_FD_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,i_FD_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)

lines(seq(1,24,0.01)+sub_rc[1]+1,j_FD*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2,lty=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,j_FD_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2,lty=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,j_FD_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2,lty=2)

text(sub_rc[1]+25+0.4,i_FD[2301]*seq+k*seq*n+sub_rc[4]+2,expression(i["FD"]),cex=1,adj=0,family="Times New Roman")
text(sub_rc[1]+25+0.4,j_FD[2301]*seq+k*seq*n+sub_rc[4]+2,expression(j["FD"]),cex=1,adj=0,family="Times New Roman")


text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-14,"Age (year)",cex=1.2,family="Times New Roman")


text(sub_rc[1]+2,sub_rc[2]-5,expression(BC),cex=1,adj=0,family="Times New Roman")


####################################################################3
sub_rc<-c(6+2*a,70+b/2*3,34.5+2*a,10+b/2*3)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}

#max(i_HD,i_HD_Y,i_HD_A,j_HD,j_HD_Y,j_HD_A)
#min(i_HD,i_HD_Y,i_HD_A,j_HD,j_HD_Y,j_HD_A)

k=0.01
seq<- 58/k/14
n<-6
for(i in -3:4){
  segments(sub_rc[3],sub_rc[4]+2+ 2*k*seq*i+k*seq*n,sub_rc[3]-0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[3]+0.6,sub_rc[4]+2+ 2*k*seq*i+k*seq*n,2*k*i, cex=1,adj=0,font=1)
}

segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)

lines(seq(1,24,0.01)+sub_rc[1]+1,i_HD*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,i_HD_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,i_HD_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)

lines(seq(1,24,0.01)+sub_rc[1]+1,j_HD*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2,lty=2)
lines(seq(1,t1_max,0.01)+sub_rc[1]+1,j_HD_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2,lty=2)
lines(seq(1,24,0.01)+sub_rc[1]+1,j_HD_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2,lty=2)

text(sub_rc[1]+25+0.4,i_HD[2301]*seq+k*seq*n+sub_rc[4]+2+4,expression(i["HD"]),cex=1,adj=0,family="Times New Roman")
text(sub_rc[1]+25+0.4,j_HD[2301]*seq+k*seq*n+sub_rc[4]+2-2,expression(j["HD"]),cex=1,adj=0,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-5,expression(AC),cex=1,adj=0,family="Times New Roman")
######################################################################################

sub_rc<-c(34.5+7.2+3,10+b/2*4-15,6+2*a-7.2+3,70+15)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=0.6)


for(i in 0:5){
  segments(sub_rc[1]+1+5*1.25*i,sub_rc[4],sub_rc[1]+1+5*1.25*i,sub_rc[4]+1.2,font=2,lwd=0.6)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1,font=1)
}

#max(j_FHD,j_FHD_Y,j_FHD_A)
#min(j_FHD,j_FHD_Y,j_FHD_A)

k=0.005
seq<- 88/k/5
n<-2.4
for(i in -1:1){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,1*k*i, cex=1,adj=1,font=1)
}
for(i in -2){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,"-0.010", cex=1,adj=1,font=1)
}
for(i in 2){
  segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,"0.010", cex=1,adj=1,font=1)
}
segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)


lines(seq(1,24,0.01)*1.25+sub_rc[1]+1,j_FHD*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2,lty=2)
lines(seq(1,t1_max,0.01)*1.25+sub_rc[1]+1,j_FHD_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2,lty=2)
lines(seq(1,24,0.01)*1.25+sub_rc[1]+1,j_FHD_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2,lty=2)


text(sub_rc[1]+25*1.25+0.3,j_FHD[2301]*seq+k*seq*n+sub_rc[4]+2+3,expression(j["FHD"]),cex=1,adj=0,family="Times New Roman")


#max(i_FHD,i_FHD_Y,i_FHD_A)
#min(i_FHD,i_FHD_Y,i_FHD_A)


#k=0.001
#seq<- 88/k/5.2
#n<-3
#for(i in -3:2){
#segments(sub_rc[1],sub_rc[4]+2+ 1*k*seq*i+k*seq*n,sub_rc[1]+0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,font=2,lwd=0.6)
#text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*k*seq*i+k*seq*n,1*k*i, cex=1,adj=1,font=1)
#}

#segments(sub_rc[1],k*seq*n+sub_rc[4]+2,sub_rc[3],k*seq*n+sub_rc[4]+2,lty=3,col="gray60",lwd=1.2)

lines(seq(1,24,0.01)*1.25+sub_rc[1]+1,i_FHD*seq+k*seq*n+sub_rc[4]+2,col="green3",lwd=2)
lines(seq(1,t1_max,0.01)*1.25+sub_rc[1]+1,i_FHD_Y*seq+k*seq*n+sub_rc[4]+2,col="red",lwd=2)
lines(seq(1,24,0.01)*1.25+sub_rc[1]+1,i_FHD_A*seq+k*seq*n+sub_rc[4]+2,col="blue",lwd=2)

text(sub_rc[1]+25*1.25+0.3,i_FHD[2301]*seq+k*seq*n+sub_rc[4]+2-3,expression(i["FHD"]),cex=1,adj=0,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-5,expression(ABC),cex=1,adj=0,family="Times New Roman")

ang<-  atan(6/15)
arrows(35,55,41-1.8*tan(ang),40+1.8,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)

for(i in 1:30){
  arrows(41-1.8*tan(ang),40+1.8,41,40,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}

arrows(20,116.5,20,141.5,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)

for(i in 1:30){
  arrows(20,141.5,20,143.5,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}


ang<-  atan(9.2/15)
arrows(35,70,44.2-1.2*tan(ang),85-1.2,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)


for(i in 1:30){
  arrows(44.2-1.2*tan(ang),85-1.2,44.2,85,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}


ang<-  atan(5/15)
arrows(41,220,36-1.8*tan(ang),205-1.8,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)

for(i in 1:30){
  arrows(36-1.8*tan(ang),205-1.8,36,205,lwd=2,col="grey60",length=0.1,angle=i,code=1,lty=1)
}


ang<-  atan(10.5/15)
arrows(75,220,85.5-1.8*tan(ang),205+1.8,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)
for(i in 1:30){
  arrows(85.5-1.8*tan(ang),205+1.8,85.5,205,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}



arrows(60.25,188.5,60.25,178.5,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)
for(i in 1:30){
  arrows(60.25,178.5,60.25,176.5,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}


arrows(100,116.5,100,141.5,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)
for(i in 1:30){
  arrows(100,141.5,100,143.5,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}



ang<-  atan(10.5/15)
arrows(85.5,55,75+1.2*tan(ang),40+1.2,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)
for(i in 1:30){
  arrows(75+1.2*tan(ang),40+1.2,75,40,lwd=2,col="grey60",length=0.1,angle=i,code=2,lty=1)
}


ang<-  atan(3.2/12)
arrows(85.5,70,82.3-1.2*tan(ang),82+1.2,lwd=4,col="grey60",length=0.3,angle=15,code=0,lty=1)
for(i in 1:30){
  arrows(82.3-1.2*tan(ang),82+1.2,82.3,82,lwd=2,col="grey60",length=0.1,angle=i,code=1,lty=1)
}


dev.off()