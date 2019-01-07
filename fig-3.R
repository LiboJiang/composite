

par<-final_par[145255,-1]




library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")


pdf("fig-3.pdf",width=12,height =10)
height<-192
length<-74
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(-1,length),ylim=c(0,height)); 



#############################################################################################

sub_rc<-c(5,70,32.5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}
seq<-58/5/0.02

para_plot1<- par[1:8]
para_plot2<- par[21:28]
para_plot3<- par[41:48]

F1<-Legendre.model(seq(1,24,0.01),para_plot1)
F2<-Legendre.model(seq(1,24,0.01),para_plot2)
F3<-Legendre.model(seq(1,24,0.01),para_plot3)

a<-1/2*(F1-F3)
m<-1/2*(F1+F3)
d<-F2-m

lines(seq(1,24,0.01)+sub_rc[1]+1,a*seq+0.02*seq*2+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(1,24,0.01)+sub_rc[1]+1,d*seq+0.02*seq*2+sub_rc[4]+2,col="green3",lwd=3,lty=2)

for(i in -2:3){
  segments(sub_rc[1],sub_rc[4]+2+ 0.02*seq*i+0.02*seq*2,sub_rc[1]+0.3,sub_rc[4]+2+ 0.02*seq*i+0.02*seq*2,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.02*seq*i+0.02*seq*2,0.02*i, cex=1.6,adj=1,font=1)
}


segments(sub_rc[1],0.02*seq*2+sub_rc[4]+2,sub_rc[3],0.02*seq*2+sub_rc[4]+2,lty=3,col="gray60",lwd=2)

text(sub_rc[1]-6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-2,expression("Genetic Effect" ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-5,sub_rc[2],expression(C),cex=2,family="Times New Roman")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##########################################################################################


sub_rc<-c(5,70+63,32.5,10+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}



para_plot1<- par[15:20]
para_plot2<- par[35:40]
para_plot3<- par[55:60]

F1<-fn_P(seq(0,24,0.01),para_plot1)
F2<-fn_P(seq(0,24,0.01),para_plot2)
F3<-fn_P(seq(0,24,0.01),para_plot3)

a_F<-1/2*(F1-F3)
m_F<-1/2*(F1+F3)
d_F<-F2-m_F

t1_max_1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max_2<- -(log(1/0.96-1)-log(para_plot2[2]))/para_plot2[3]
t1_max_3<- -(log(1/0.96-1)-log(para_plot3[2]))/para_plot3[3]
t1_max<-max(t1_max_1,t1_max_2,t1_max_3)


Y1<-fn_P_1(seq(0,t1_max,0.01),para_plot1[1:3])
Y2<-fn_P_1(seq(0,t1_max,0.01),para_plot2[1:3])
Y3<-fn_P_1(seq(0,t1_max,0.01),para_plot3[1:3])

a_Y<-1/2*(Y1-Y3)
m_Y<-1/2*(Y1+Y3)
d_Y<-Y2-m_Y


A1<-fn_P_1(seq(0,24,0.01),para_plot1[4:6])
A2<-fn_P_1(seq(0,24,0.01),para_plot2[4:6])
A3<-fn_P_1(seq(0,24,0.01),para_plot3[4:6])

a_A<-1/2*(A1-A3)
m_A<-1/2*(A1+A3)
d_A<-A2-m_A


#max(a_F,d_F,a_Y,d_Y,a_A,d_A)

seq<-23.2

for(i in c(-3,-1,0,1)){
  segments(sub_rc[1],sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,sub_rc[1]+0.3,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,0.5*i, cex=1.6,adj=1,font=1)
}
for(i in -2){
  segments(sub_rc[1],sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,sub_rc[1]+0.3,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,"-1.0", cex=1.6,adj=1,font=1)
}
for(i in 2){
  segments(sub_rc[1],sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,sub_rc[1]+0.3,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,"1.0", cex=1.6,adj=1,font=1)
}




lines(seq(0,24,0.01)+sub_rc[1]+1,a_F*seq+0.5*seq*3+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_F*seq+0.5*seq*3+sub_rc[4]+2,col="green3",lwd=3,lty=2)

lines(seq(0,t1_max,0.01)+sub_rc[1]+1,a_Y*seq+0.5*seq*3+sub_rc[4]+2,col="red",lwd=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,d_Y*seq+0.5*seq*3+sub_rc[4]+2,col="red",lwd=3,lty=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,a_A*seq+0.5*seq*3+sub_rc[4]+2,col="blue",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_A*seq+0.5*seq*3+sub_rc[4]+2,col="blue",lwd=3,lty=2)


segments(sub_rc[1],0.5*seq*3+sub_rc[4]+2,sub_rc[3],0.5*seq*3+sub_rc[4]+2,lty=3,col="gray60",lwd=2)


text(sub_rc[1]-6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Genetic Effect (cm)"),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Age (year)",cex=2,font=1)
text(sub_rc[1]-5,sub_rc[2],expression(B),cex=2,family="Times New Roman")

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

###############################################################################

#par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 


sub_rc<-c(5,70+63+63,32.5,10+63+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


para_plot1<- H0_P_V1$par[1:6]
para_plot2<- H0_P_V2$par[1:6]
para_plot3<- H0_P_V3$par[1:6]

F1<-fn_P(seq(0,24,0.01),para_plot1)
F2<-fn_P(seq(0,24,0.01),para_plot2)
F3<-fn_P(seq(0,24,0.01),para_plot3)

a_F<-1/2*(F1-F3)
m_F<-1/2*(F1+F3)
d_F<-F2-m_F

t1_max_1<- -(log(1/0.9-1)-log(para_plot1[5]))/para_plot1[6]
t1_max_2<- -(log(1/0.9-1)-log(para_plot2[5]))/para_plot2[6]
t1_max_3<- -(log(1/0.9-1)-log(para_plot3[5]))/para_plot3[6]
t1_max<-max(t1_max_1,t1_max_2,t1_max_3)


Y1<-fn_P_1(seq(0,t1_max,0.01),para_plot1[1:3])
Y2<-fn_P_1(seq(0,t1_max,0.01),para_plot2[1:3])
Y3<-fn_P_1(seq(0,t1_max,0.01),para_plot3[1:3])

a_Y<-1/2*(Y1-Y3)
m_Y<-1/2*(Y1+Y3)
d_Y<-Y2-m_Y


A1<-fn_P_1(seq(0,24,0.01),para_plot1[4:6])
A2<-fn_P_1(seq(0,24,0.01),para_plot2[4:6])
A3<-fn_P_1(seq(0,24,0.01),para_plot3[4:6])

a_A<-1/2*(A1-A3)
m_A<-1/2*(A1+A3)
d_A<-A2-m_A


#max(a_F,d_F,a_Y,d_Y,a_A,d_A)

seq<-7.8


for(i in -3:4){
  segments(sub_rc[1],sub_rc[4]+2+ 1*seq*i+1*seq*3.2,sub_rc[1]+0.3,sub_rc[4]+2+ 1*seq*i+1*seq*3.2,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*seq*i+1*seq*3.2,1*i, cex=1.6,adj=1,font=1)
}

lines(seq(0,24,0.01)+sub_rc[1]+1,a_F*seq+1*seq*3.2+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_F*seq+1*seq*3.2+sub_rc[4]+2,col="green3",lwd=3,lty=2)

lines(seq(0,t1_max,0.01)+sub_rc[1]+1,a_Y*seq+1*seq*3.2+sub_rc[4]+2,col="red",lwd=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,d_Y*seq+1*seq*3.2+sub_rc[4]+2,col="red",lwd=3,lty=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,a_A*seq+1*seq*3.2+sub_rc[4]+2,col="blue",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_A*seq+1*seq*3.2+sub_rc[4]+2,col="blue",lwd=3,lty=2)


segments(sub_rc[1],1*seq*3.2+sub_rc[4]+2,sub_rc[3],1*seq*3.2+sub_rc[4]+2,lty=3,col="gray60",lwd=2)

text(sub_rc[1]-6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Genetic Effect (m)"),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]-5.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height"),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Time (year)",cex=2,font=1)
text(sub_rc[1]-5,sub_rc[2],expression(A),cex=2,family="Times New Roman")

##################################################################

sub_rc<-c(40.5,height/2+45,68,height/2-45)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}



para_plot1<- g_V1$par
para_plot2<- g_V2$par
para_plot3<- g_V3$par

F1<-fn_P(seq(0,24,0.01),para_plot1)/1000000
F2<-fn_P(seq(0,24,0.01),para_plot2)/1000000
F3<-fn_P(seq(0,24,0.01),para_plot3)/1000000

a_F<-1/2*(F1-F3)
m_F<-1/2*(F1+F3)
d_F<-F2-m_F

t1_max_1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t1_max_2<- -(log(1/0.96-1)-log(para_plot2[5]))/para_plot2[6]
t1_max_3<- -(log(1/0.96-1)-log(para_plot3[5]))/para_plot3[6]
t1_max<-max(t1_max_1,t1_max_2,t1_max_3)


Y1<-fn_P_1(seq(0,t1_max,0.01),para_plot1[1:3])/1000000
Y2<-fn_P_1(seq(0,t1_max,0.01),para_plot2[1:3])/1000000
Y3<-fn_P_1(seq(0,t1_max,0.01),para_plot3[1:3])/1000000

a_Y<-1/2*(Y1-Y3)
m_Y<-1/2*(Y1+Y3)
d_Y<-Y2-m_Y


A1<-fn_P_1(seq(0,24,0.01),para_plot1[4:6])/1000000
A2<-fn_P_1(seq(0,24,0.01),para_plot2[4:6])/1000000
A3<-fn_P_1(seq(0,24,0.01),para_plot3[4:6])/1000000

a_A<-1/2*(A1-A3)
m_A<-1/2*(A1+A3)
d_A<-A2-m_A


#max(a_F,d_F,a_Y,d_Y,a_A,d_A)
#min(a_F,d_F,a_Y,d_Y,a_A,d_A)
k<-0.01
seq<-88/8/k

for(i in -4:4){
  segments(sub_rc[3],sub_rc[4]+2+ 1*k*seq*i+1*k*seq*4,sub_rc[3]-0.3,sub_rc[4]+2+ 1*k*seq*i+1*k*seq*4,font=2,lwd=1)
  text(sub_rc[3]+0.6,sub_rc[4]+2+ 1*k*seq*i+1*k*seq*4,1*0.01*i, cex=1.6,adj=0,font=1)
}


lines(seq(0,24,0.01)+sub_rc[1]+1,a_F*seq+1*k*seq*4+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_F*seq+1*k*seq*4+sub_rc[4]+2,col="green3",lwd=3,lty=2)

lines(seq(0,t1_max,0.01)+sub_rc[1]+1,a_Y*seq+1*k*seq*4+sub_rc[4]+2,col="red",lwd=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,d_Y*seq+1*k*seq*4+sub_rc[4]+2,col="red",lwd=3,lty=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,a_A*seq+1*k*seq*4+sub_rc[4]+2,col="blue",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_A*seq+1*k*seq*4+sub_rc[4]+2,col="blue",lwd=3,lty=2)


segments(sub_rc[1],1*k*seq*4+sub_rc[4]+2,sub_rc[3],1*k*seq*4+sub_rc[4]+2,lty=3,col="gray60",lwd=2)



text(sub_rc[3]+6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Genetic Effect"~(m^3) ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-1.5,sub_rc[2],expression(D),cex=2,family="Times New Roman")

arrows(32.7,96,40,96,lwd=1,col="black",length=0.3,angle=15,code=0,lty=2)

arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)

arrows(32.7,36,40,95,lwd=1.5,col="black",length=0.3,angle=15,code=0,lty=2)
for(i in 1:15){
  arrows(40- 7.3/59* 4.3,95- 4.3,40,95,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}



arrows(32.7,162,40,97,lwd=1.5,col="black",length=0.3,angle=15,code=0,lty=2)

for(i in 1:15){
  arrows(40- 7.3/65* 4.3,97+ 4.3,40,97,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}

dev.off()





jpeg("fig-3.jpeg",width=1200,height =1000, units = "px", pointsize = 12,quality =150)
height<-192
length<-74
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(-1,length),ylim=c(0,height)); 



#############################################################################################

sub_rc<-c(5,70,32.5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}
seq<-58/5/0.02

para_plot1<- par[1:8]
para_plot2<- par[21:28]
para_plot3<- par[41:48]

F1<-Legendre.model(seq(1,24,0.01),para_plot1)
F2<-Legendre.model(seq(1,24,0.01),para_plot2)
F3<-Legendre.model(seq(1,24,0.01),para_plot3)

a<-1/2*(F1-F3)
m<-1/2*(F1+F3)
d<-F2-m

lines(seq(1,24,0.01)+sub_rc[1]+1,a*seq+0.02*seq*2+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(1,24,0.01)+sub_rc[1]+1,d*seq+0.02*seq*2+sub_rc[4]+2,col="green3",lwd=3,lty=2)

for(i in -2:3){
  segments(sub_rc[1],sub_rc[4]+2+ 0.02*seq*i+0.02*seq*2,sub_rc[1]+0.3,sub_rc[4]+2+ 0.02*seq*i+0.02*seq*2,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.02*seq*i+0.02*seq*2,0.02*i, cex=1.6,adj=1,font=1)
}


segments(sub_rc[1],0.02*seq*2+sub_rc[4]+2,sub_rc[3],0.02*seq*2+sub_rc[4]+2,lty=3,col="gray60",lwd=2)

text(sub_rc[1]-6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-2,expression("Genetic Effect" ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-5,sub_rc[2],expression(C),cex=2,family="Times New Roman")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##########################################################################################


sub_rc<-c(5,70+63,32.5,10+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}



para_plot1<- par[15:20]
para_plot2<- par[35:40]
para_plot3<- par[55:60]

F1<-fn_P(seq(0,24,0.01),para_plot1)
F2<-fn_P(seq(0,24,0.01),para_plot2)
F3<-fn_P(seq(0,24,0.01),para_plot3)

a_F<-1/2*(F1-F3)
m_F<-1/2*(F1+F3)
d_F<-F2-m_F

t1_max_1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max_2<- -(log(1/0.96-1)-log(para_plot2[2]))/para_plot2[3]
t1_max_3<- -(log(1/0.96-1)-log(para_plot3[2]))/para_plot3[3]
t1_max<-max(t1_max_1,t1_max_2,t1_max_3)


Y1<-fn_P_1(seq(0,t1_max,0.01),para_plot1[1:3])
Y2<-fn_P_1(seq(0,t1_max,0.01),para_plot2[1:3])
Y3<-fn_P_1(seq(0,t1_max,0.01),para_plot3[1:3])

a_Y<-1/2*(Y1-Y3)
m_Y<-1/2*(Y1+Y3)
d_Y<-Y2-m_Y


A1<-fn_P_1(seq(0,24,0.01),para_plot1[4:6])
A2<-fn_P_1(seq(0,24,0.01),para_plot2[4:6])
A3<-fn_P_1(seq(0,24,0.01),para_plot3[4:6])

a_A<-1/2*(A1-A3)
m_A<-1/2*(A1+A3)
d_A<-A2-m_A


#max(a_F,d_F,a_Y,d_Y,a_A,d_A)

seq<-23.2

for(i in c(-3,-1,0,1)){
  segments(sub_rc[1],sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,sub_rc[1]+0.3,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,0.5*i, cex=1.6,adj=1,font=1)
}
for(i in -2){
  segments(sub_rc[1],sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,sub_rc[1]+0.3,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,"-1.0", cex=1.6,adj=1,font=1)
}
for(i in 2){
  segments(sub_rc[1],sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,sub_rc[1]+0.3,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 0.5*seq*i+0.5*seq*3,"1.0", cex=1.6,adj=1,font=1)
}




lines(seq(0,24,0.01)+sub_rc[1]+1,a_F*seq+0.5*seq*3+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_F*seq+0.5*seq*3+sub_rc[4]+2,col="green3",lwd=3,lty=2)

lines(seq(0,t1_max,0.01)+sub_rc[1]+1,a_Y*seq+0.5*seq*3+sub_rc[4]+2,col="red",lwd=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,d_Y*seq+0.5*seq*3+sub_rc[4]+2,col="red",lwd=3,lty=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,a_A*seq+0.5*seq*3+sub_rc[4]+2,col="blue",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_A*seq+0.5*seq*3+sub_rc[4]+2,col="blue",lwd=3,lty=2)


segments(sub_rc[1],0.5*seq*3+sub_rc[4]+2,sub_rc[3],0.5*seq*3+sub_rc[4]+2,lty=3,col="gray60",lwd=2)


text(sub_rc[1]-6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Genetic Effect (cm)"),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Age (year)",cex=2,font=1)
text(sub_rc[1]-5,sub_rc[2],expression(B),cex=2,family="Times New Roman")

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

###############################################################################

#par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 


sub_rc<-c(5,70+63+63,32.5,10+63+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


para_plot1<- H0_P_V1$par[1:6]
para_plot2<- H0_P_V2$par[1:6]
para_plot3<- H0_P_V3$par[1:6]

F1<-fn_P(seq(0,24,0.01),para_plot1)
F2<-fn_P(seq(0,24,0.01),para_plot2)
F3<-fn_P(seq(0,24,0.01),para_plot3)

a_F<-1/2*(F1-F3)
m_F<-1/2*(F1+F3)
d_F<-F2-m_F

t1_max_1<- -(log(1/0.9-1)-log(para_plot1[5]))/para_plot1[6]
t1_max_2<- -(log(1/0.9-1)-log(para_plot2[5]))/para_plot2[6]
t1_max_3<- -(log(1/0.9-1)-log(para_plot3[5]))/para_plot3[6]
t1_max<-max(t1_max_1,t1_max_2,t1_max_3)


Y1<-fn_P_1(seq(0,t1_max,0.01),para_plot1[1:3])
Y2<-fn_P_1(seq(0,t1_max,0.01),para_plot2[1:3])
Y3<-fn_P_1(seq(0,t1_max,0.01),para_plot3[1:3])

a_Y<-1/2*(Y1-Y3)
m_Y<-1/2*(Y1+Y3)
d_Y<-Y2-m_Y


A1<-fn_P_1(seq(0,24,0.01),para_plot1[4:6])
A2<-fn_P_1(seq(0,24,0.01),para_plot2[4:6])
A3<-fn_P_1(seq(0,24,0.01),para_plot3[4:6])

a_A<-1/2*(A1-A3)
m_A<-1/2*(A1+A3)
d_A<-A2-m_A


#max(a_F,d_F,a_Y,d_Y,a_A,d_A)

seq<-7.8


for(i in -3:4){
  segments(sub_rc[1],sub_rc[4]+2+ 1*seq*i+1*seq*3.2,sub_rc[1]+0.3,sub_rc[4]+2+ 1*seq*i+1*seq*3.2,font=2,lwd=1)
  text(sub_rc[1]-0.6,sub_rc[4]+2+ 1*seq*i+1*seq*3.2,1*i, cex=1.6,adj=1,font=1)
}

lines(seq(0,24,0.01)+sub_rc[1]+1,a_F*seq+1*seq*3.2+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_F*seq+1*seq*3.2+sub_rc[4]+2,col="green3",lwd=3,lty=2)

lines(seq(0,t1_max,0.01)+sub_rc[1]+1,a_Y*seq+1*seq*3.2+sub_rc[4]+2,col="red",lwd=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,d_Y*seq+1*seq*3.2+sub_rc[4]+2,col="red",lwd=3,lty=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,a_A*seq+1*seq*3.2+sub_rc[4]+2,col="blue",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_A*seq+1*seq*3.2+sub_rc[4]+2,col="blue",lwd=3,lty=2)


segments(sub_rc[1],1*seq*3.2+sub_rc[4]+2,sub_rc[3],1*seq*3.2+sub_rc[4]+2,lty=3,col="gray60",lwd=2)

text(sub_rc[1]-6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Genetic Effect (m)"),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]-5.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height"),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Time (year)",cex=2,font=1)
text(sub_rc[1]-5,sub_rc[2],expression(A),cex=2,family="Times New Roman")

##################################################################

sub_rc<-c(40.5,height/2+45,68,height/2-45)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}



para_plot1<- g_V1$par
para_plot2<- g_V2$par
para_plot3<- g_V3$par

F1<-fn_P(seq(0,24,0.01),para_plot1)/1000000
F2<-fn_P(seq(0,24,0.01),para_plot2)/1000000
F3<-fn_P(seq(0,24,0.01),para_plot3)/1000000

a_F<-1/2*(F1-F3)
m_F<-1/2*(F1+F3)
d_F<-F2-m_F

t1_max_1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t1_max_2<- -(log(1/0.96-1)-log(para_plot2[5]))/para_plot2[6]
t1_max_3<- -(log(1/0.96-1)-log(para_plot3[5]))/para_plot3[6]
t1_max<-max(t1_max_1,t1_max_2,t1_max_3)


Y1<-fn_P_1(seq(0,t1_max,0.01),para_plot1[1:3])/1000000
Y2<-fn_P_1(seq(0,t1_max,0.01),para_plot2[1:3])/1000000
Y3<-fn_P_1(seq(0,t1_max,0.01),para_plot3[1:3])/1000000

a_Y<-1/2*(Y1-Y3)
m_Y<-1/2*(Y1+Y3)
d_Y<-Y2-m_Y


A1<-fn_P_1(seq(0,24,0.01),para_plot1[4:6])/1000000
A2<-fn_P_1(seq(0,24,0.01),para_plot2[4:6])/1000000
A3<-fn_P_1(seq(0,24,0.01),para_plot3[4:6])/1000000

a_A<-1/2*(A1-A3)
m_A<-1/2*(A1+A3)
d_A<-A2-m_A


#max(a_F,d_F,a_Y,d_Y,a_A,d_A)
#min(a_F,d_F,a_Y,d_Y,a_A,d_A)
k<-0.01
seq<-88/8/k

for(i in -4:4){
  segments(sub_rc[3],sub_rc[4]+2+ 1*k*seq*i+1*k*seq*4,sub_rc[3]-0.3,sub_rc[4]+2+ 1*k*seq*i+1*k*seq*4,font=2,lwd=1)
  text(sub_rc[3]+0.6,sub_rc[4]+2+ 1*k*seq*i+1*k*seq*4,1*0.01*i, cex=1.6,adj=0,font=1)
}


lines(seq(0,24,0.01)+sub_rc[1]+1,a_F*seq+1*k*seq*4+sub_rc[4]+2,col="green3",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_F*seq+1*k*seq*4+sub_rc[4]+2,col="green3",lwd=3,lty=2)

lines(seq(0,t1_max,0.01)+sub_rc[1]+1,a_Y*seq+1*k*seq*4+sub_rc[4]+2,col="red",lwd=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,d_Y*seq+1*k*seq*4+sub_rc[4]+2,col="red",lwd=3,lty=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,a_A*seq+1*k*seq*4+sub_rc[4]+2,col="blue",lwd=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,d_A*seq+1*k*seq*4+sub_rc[4]+2,col="blue",lwd=3,lty=2)


segments(sub_rc[1],1*k*seq*4+sub_rc[4]+2,sub_rc[3],1*k*seq*4+sub_rc[4]+2,lty=3,col="gray60",lwd=2)



text(sub_rc[3]+6.6,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Genetic Effect"~(m^3) ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-1.5,sub_rc[2],expression(D),cex=2,family="Times New Roman")

arrows(32.7,96,40,96,lwd=1,col="black",length=0.3,angle=15,code=0,lty=2)

arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(39.5,96,40,96,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)

arrows(32.7,36,40,95,lwd=1.5,col="black",length=0.3,angle=15,code=0,lty=2)
for(i in 1:15){
  arrows(40- 7.3/59* 4.3,95- 4.3,40,95,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}



arrows(32.7,162,40,97,lwd=1.5,col="black",length=0.3,angle=15,code=0,lty=2)

for(i in 1:15){
  arrows(40- 7.3/65* 4.3,97+ 4.3,40,97,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}



dev.off()


