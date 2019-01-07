setwd("D:/Rdocument/p-pleiotropic model")


fn_P_1<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))
}

library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

pdf("fig1.pdf",width=12,height =10)

showtext_begin()

height<-192
length<-71
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 

sub_rc<-c(5,height/2+45,32.5,height/2-45)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
para_plot<- H0_P_V$par


ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


#lines(t+sub_rc[1]+1-0.5,fn_P(t,para_plot)/12000+sub_rc[4]+5,col="red",lwd=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+5,lwd=2.4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+5,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+5,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+5,lwd=2.4,lty=2,col="green3")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+5,lwd=2,lty=2,col="blue")

#segments(ti_2+sub_rc[1]+1-0.5,gi_2/12000+sub_rc[4]+5,ti_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2/12000+sub_rc[4]+5,pch=17)
#segments(ti_1+sub_rc[1]+1-0.5,gi_1/12000+sub_rc[4]+5,ti_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1/12000+sub_rc[4]+5,pch=17)
text(ti_1+sub_rc[1]+1.5,gi_1/12000+sub_rc[4]+2,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")
text(ti_2+sub_rc[1]+2,gi_2/12000+sub_rc[4]+2,expression(italic(T["I"])),cex=1.4,family="Times New Roman1")

#segments(ta_2+sub_rc[1]+1-0.5,ga_2/12000+sub_rc[4]+5,ta_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2/12000+sub_rc[4]+5,pch=17)
#segments(ta_1+sub_rc[1]+1-0.5,ga_1/12000+sub_rc[4]+5,ta_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1/12000+sub_rc[4]+5,pch=17)
text(ta_1+sub_rc[1]+0.3,ga_1/12000+sub_rc[4]+8,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2/12000+sub_rc[4]+2,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

#segments(td_2+sub_rc[1]+1-0.5,gd_2/12000+sub_rc[4]+5,td_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2/12000+sub_rc[4]+5,pch=17)
#segments(td_1+sub_rc[1]+1-0.5,gd_1/12000+sub_rc[4]+5,td_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1/12000+sub_rc[4]+5,pch=17)
text(td_1+sub_rc[1]+1.5,gd_1/12000+sub_rc[4]+2,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1],gd_2/12000+sub_rc[4]+8,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+5
#region_y<-rep(sub_rc[2]-sub_rc[4],length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)
#polygon(region_x+sub_rc[1]+1-0.5,region_y+sub_rc[4],density = 40)

#segments(t2_min+sub_rc[1]+1-0.5,fn_P(t2_min,para_plot[1:6])/12000+sub_rc[4]+5,t2_min+sub_rc[1]+1-0.5,sub_rc[4],lty=1,lwd=2,col="black")
#segments(t1_max+sub_rc[1]+1-0.5,fn_P(t1_max,para_plot[1:6])/12000+sub_rc[4]+5,t1_max+sub_rc[1]+1-0.5,sub_rc[4],lty=1,lwd=2,col="black")

arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+0.5,sub_rc[4]+5+40+2.5,pch=2,cex=1)
text(ta_1+sub_rc[1]+1+1,sub_rc[4]+5+80-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+5+80-3.5,pch=2,cex=1)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+5+80-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+5+60+2.5,pch=2,cex=1)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+5+70+3,expression(italic(L)),family="Times New Roman1",cex=1.4)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+5+200000/12000*i,sub_rc[1]+0.30,sub_rc[4]+5+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+5+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+5+200000/12000*i,sub_rc[1]+0.30,sub_rc[4]+5+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+5+200000/12000*i,"1.0", cex=1.6,font=1)
}

text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stemwood Volume " (m^3)),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-4.5,sub_rc[2],expression(A),cex=2,family="Times New Roman")

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
#showtext_end()
#dev.off()


#############################################################################################

sub_rc<-c(45.5,70,73,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot<- H0_P_psi$par[1:8]
#Legendre.model(0:24,para_plot)
#Legendre.model(0:26,para_plot)
#Legendre.model(24:26,para_plot)

for(i in 1:dim(dat_P_st$psi)[1]){
  lines(c(1:24)+sub_rc[1]+1,(dat_P_st$psi[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

#lines(seq(0,26,0.01)+sub_rc[1]+1-0.5,Legendre.model(seq(0,26,0.01),para_plot)*45*1.2+sub_rc[4]+2,col="black",lwd=2,lty=3)
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

para_plot<- H0_P_DIA$par
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)

para_plot<- H0_P_HT$par
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)


text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-2,expression("Stem Form Factor" ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-4.5,sub_rc[2],expression(D),cex=2,family="Times New Roman")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##########################################################################################
#par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 

sub_rc<-c(45.5,70+63,73,10+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot<- H0_P_DIA$par

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.3,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm)),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Time (year)",cex=2,font=1)
text(sub_rc[1]-4.5,sub_rc[2],expression(C),cex=2,family="Times New Roman")

ti_2<-log(para_plot[5])/para_plot[6]
gi_2<-para_plot[4]/2*1.4

ta_2<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_2<-para_plot[4]*(3-sqrt(3))/6*1.4

td_2<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_2<-para_plot[4]*(3+sqrt(3))/6*1.4


ti_1<-log(para_plot[2])/para_plot[3]
gi_1<-para_plot[1]/2*1.4

ta_1<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_1<-para_plot[1]*(3-sqrt(3))/6*1.4

td_1<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_1<-para_plot[1]*(3+sqrt(3))/6*1.4

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,lty=2,col="green3")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,lty=2,col="blue")


#segments(ti_2+sub_rc[1]+1-0.5,gi_2/12000+sub_rc[4]+5,ti_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
#segments(ti_1+sub_rc[1]+1-0.5,gi_1/12000+sub_rc[4]+5,ti_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+1.7,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

#segments(ta_2+sub_rc[1]+1-0.5,ga_2/12000+sub_rc[4]+5,ta_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
#segments(ta_1+sub_rc[1]+1-0.5,ga_1/12000+sub_rc[4]+5,ta_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+5,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+0.1,gd_1+sub_rc[4]+4,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2,gd_2+sub_rc[4]-2,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)



region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])*1.4+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,angle = 135,col="gray60",lty=1)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+4,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43+3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

###############################################################################

#par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 


sub_rc<-c(45.5,70+63+63,73,10+63+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot<- H0_P_HT$par

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}


text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height " (m)),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Time (year)",cex=2,font=1)
text(sub_rc[1]-4.5,sub_rc[2],expression(B),cex=2,family="Times New Roman")

ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2*2.2

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6*2.2

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6*2.2


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2*2.2

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6*2.2

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6*2.2

t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,lty=2,col="green3")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,lty=2,col="blue")


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+1.7,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2.2,gd_1+sub_rc[4],expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.2,gd_2+sub_rc[4]-1.2,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)



region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])*2.2+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.2,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+2.5,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+35+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+1,sub_rc[4]+2+43+3,expression(italic(L)),family="Times New Roman1",cex=1.4)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
#arrows(32.7,96,38.5,96,lwd=2,col="black",length=0.2,angle=20,code=2,lty=2)
arrows(32.7,96,38.5,96,lwd=1,col="black",length=0.3,angle=15,code=0,lty=2)
#arrows(32.7,96,38.5,96,lwd=2,col="black",length=0.3,angle=10,code=2,lty=2)
#arrows(32.7,96,38.5,96,lwd=2,col="black",length=0.3,angle=5,code=2,lty=2)

#arrows(38.3,96,38.5,96,lwd=2,col="black",length=0.2,angle=20,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)
#segments(34,36,34,162,lwd=2,col="black")

#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=20,code=2)
#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=15,code=2)
#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=10,code=2)
#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=5,code=2)

#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=20,code=2)
#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=15,code=2)
#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=10,code=2)
#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=5,code=2)


#arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=20,code=2,lty=2)
arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=15,code=0,lty=2)
#arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=10,code=2,lty=2)
#arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=5,code=2,lty=2)

arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)
#arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=20,code=2,lty=2)
arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=15,code=0,lty=2)
#arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=10,code=2,lty=2)
#arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=5,code=2,lty=2)

arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)

showtext_end()
dev.off()







library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

jpeg("fig-1.jpeg",width=1200,height =1000, units = "px", pointsize = 12,quality =150)

showtext_begin()

height<-192
length<-71
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 

sub_rc<-c(5,height/2+45,32.5,height/2-45)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
para_plot<- H0_P_V$par


ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


#lines(t+sub_rc[1]+1-0.5,fn_P(t,para_plot)/12000+sub_rc[4]+5,col="red",lwd=2)

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+5,lwd=2.4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+5,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+5,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+5,lwd=2.4,lty=2,col="green3")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+5,lwd=2,lty=2,col="blue")

#segments(ti_2+sub_rc[1]+1-0.5,gi_2/12000+sub_rc[4]+5,ti_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2/12000+sub_rc[4]+5,pch=17)
#segments(ti_1+sub_rc[1]+1-0.5,gi_1/12000+sub_rc[4]+5,ti_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1/12000+sub_rc[4]+5,pch=17)
text(ti_1+sub_rc[1]+1.5,gi_1/12000+sub_rc[4]+2,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")
text(ti_2+sub_rc[1]+2,gi_2/12000+sub_rc[4]+2,expression(italic(T["I"])),cex=1.4,family="Times New Roman1")

#segments(ta_2+sub_rc[1]+1-0.5,ga_2/12000+sub_rc[4]+5,ta_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2/12000+sub_rc[4]+5,pch=17)
#segments(ta_1+sub_rc[1]+1-0.5,ga_1/12000+sub_rc[4]+5,ta_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1/12000+sub_rc[4]+5,pch=17)
text(ta_1+sub_rc[1]+0.3,ga_1/12000+sub_rc[4]+8,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2/12000+sub_rc[4]+2,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

#segments(td_2+sub_rc[1]+1-0.5,gd_2/12000+sub_rc[4]+5,td_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2/12000+sub_rc[4]+5,pch=17)
#segments(td_1+sub_rc[1]+1-0.5,gd_1/12000+sub_rc[4]+5,td_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1/12000+sub_rc[4]+5,pch=17)
text(td_1+sub_rc[1]+1.5,gd_1/12000+sub_rc[4]+2,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1],gd_2/12000+sub_rc[4]+8,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+5
#region_y<-rep(sub_rc[2]-sub_rc[4],length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)
#polygon(region_x+sub_rc[1]+1-0.5,region_y+sub_rc[4],density = 40)

#segments(t2_min+sub_rc[1]+1-0.5,fn_P(t2_min,para_plot[1:6])/12000+sub_rc[4]+5,t2_min+sub_rc[1]+1-0.5,sub_rc[4],lty=1,lwd=2,col="black")
#segments(t1_max+sub_rc[1]+1-0.5,fn_P(t1_max,para_plot[1:6])/12000+sub_rc[4]+5,t1_max+sub_rc[1]+1-0.5,sub_rc[4],lty=1,lwd=2,col="black")

arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+5+80,td_1+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+0.5,sub_rc[4]+5+40+2.5,pch=2,cex=1)
text(ta_1+sub_rc[1]+1+1,sub_rc[4]+5+80-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+5+80,td_2+sub_rc[1]+1,sub_rc[4]+5+80,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+5+80-3.5,pch=2,cex=1)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+5+80-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+5+70,ti_2+sub_rc[1]+1,sub_rc[4]+5+70,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+5+60+2.5,pch=2,cex=1)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+5+70+3,expression(italic(L)),family="Times New Roman1",cex=1.4)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+5+200000/12000*i,sub_rc[1]+0.30,sub_rc[4]+5+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+5+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+5+200000/12000*i,sub_rc[1]+0.30,sub_rc[4]+5+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+5+200000/12000*i,"1.0", cex=1.6,font=1)
}

text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stemwood Volume " (m^3)),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-4.5,sub_rc[2],expression(A),cex=2,family="Times New Roman")

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
#showtext_end()
#dev.off()


#############################################################################################

sub_rc<-c(45.5,70,73,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot<- H0_P_psi$par[1:8]
#Legendre.model(0:24,para_plot)
#Legendre.model(0:26,para_plot)
#Legendre.model(24:26,para_plot)

for(i in 1:dim(dat_P_st$psi)[1]){
  lines(c(1:24)+sub_rc[1]+1,(dat_P_st$psi[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

#lines(seq(0,26,0.01)+sub_rc[1]+1-0.5,Legendre.model(seq(0,26,0.01),para_plot)*45*1.2+sub_rc[4]+2,col="black",lwd=2,lty=3)
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

para_plot<- H0_P_DIA$par
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)

para_plot<- H0_P_HT$par
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)


text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-2,expression("Stem Form Factor" ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-12,"Age (year)",cex=2,family="Times New Roman")
text(sub_rc[1]-4.5,sub_rc[2],expression(D),cex=2,family="Times New Roman")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##########################################################################################
#par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 

sub_rc<-c(45.5,70+63,73,10+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot<- H0_P_DIA$par

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.3,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm)),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Time (year)",cex=2,font=1)
text(sub_rc[1]-4.5,sub_rc[2],expression(C),cex=2,family="Times New Roman")

ti_2<-log(para_plot[5])/para_plot[6]
gi_2<-para_plot[4]/2*1.4

ta_2<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_2<-para_plot[4]*(3-sqrt(3))/6*1.4

td_2<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_2<-para_plot[4]*(3+sqrt(3))/6*1.4


ti_1<-log(para_plot[2])/para_plot[3]
gi_1<-para_plot[1]/2*1.4

ta_1<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_1<-para_plot[1]*(3-sqrt(3))/6*1.4

td_1<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_1<-para_plot[1]*(3+sqrt(3))/6*1.4

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,lty=2,col="green3")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,lty=2,col="blue")


#segments(ti_2+sub_rc[1]+1-0.5,gi_2/12000+sub_rc[4]+5,ti_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
#segments(ti_1+sub_rc[1]+1-0.5,gi_1/12000+sub_rc[4]+5,ti_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+1.7,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

#segments(ta_2+sub_rc[1]+1-0.5,ga_2/12000+sub_rc[4]+5,ta_2+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="blue")
segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
#segments(ta_1+sub_rc[1]+1-0.5,ga_1/12000+sub_rc[4]+5,ta_1+sub_rc[1]+1-0.5,sub_rc[2],lty=3,lwd=1.6,col="red")
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+5,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+0.1,gd_1+sub_rc[4]+4,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2,gd_2+sub_rc[4]-2,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)



region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])*1.4+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,angle = 135,col="gray60",lty=1)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+4,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43+3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

###############################################################################

#par(mar=c(0,0,0,0),oma=c(0,0,0,0))
#plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 


sub_rc<-c(45.5,70+63+63,73,10+63+63)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot<- H0_P_HT$par

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}


text(sub_rc[1]-5.2,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height " (m)),cex=2,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-6,"Time (year)",cex=2,font=1)
text(sub_rc[1]-4.5,sub_rc[2],expression(B),cex=2,family="Times New Roman")

ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2*2.2

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6*2.2

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6*2.2


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2*2.2

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6*2.2

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6*2.2

t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,lty=2,col="green3")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,lty=2,col="blue")


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+1.7,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2.2,gd_1+sub_rc[4],expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.2,gd_2+sub_rc[4]-1.2,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)



region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])*2.2+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.2,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+2.5,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+35+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+1,sub_rc[4]+2+43+3,expression(italic(L)),family="Times New Roman1",cex=1.4)

rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
#arrows(32.7,96,38.5,96,lwd=2,col="black",length=0.2,angle=20,code=2,lty=2)
arrows(32.7,96,38.5,96,lwd=1,col="black",length=0.3,angle=15,code=0,lty=2)
#arrows(32.7,96,38.5,96,lwd=2,col="black",length=0.3,angle=10,code=2,lty=2)
#arrows(32.7,96,38.5,96,lwd=2,col="black",length=0.3,angle=5,code=2,lty=2)

#arrows(38.3,96,38.5,96,lwd=2,col="black",length=0.2,angle=20,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(38.2,96,38.5,96,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)
#segments(34,36,34,162,lwd=2,col="black")

#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=20,code=2)
#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=15,code=2)
#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=10,code=2)
#arrows(34,36,37,36,lwd=2,col="black",length=0.08,angle=5,code=2)

#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=20,code=2)
#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=15,code=2)
#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=10,code=2)
#arrows(34,162,37,162,lwd=2,col="black",length=0.08,angle=5,code=2)


#arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=20,code=2,lty=2)
arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=15,code=0,lty=2)
#arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=10,code=2,lty=2)
#arrows(32.7,96,38.5,36,lwd=1.5,col="black",length=0.3,angle=5,code=2,lty=2)

arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(38.5-5.8/60*5.8,36+5.8,38.5,36,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)
#arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=20,code=2,lty=2)
arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=15,code=0,lty=2)
#arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=10,code=2,lty=2)
#arrows(32.7,96,38.5,162,lwd=1.5,col="black",length=0.2,angle=5,code=2,lty=2)

arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=15,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=13,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=11,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=9,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=7,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=5,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=3,code=2,lty=1)
arrows(38.5-5.8/60*5.8,162-5.8,38.5,162,lwd=2,col="black",length=0.3,angle=1,code=2,lty=1)

showtext_end()
dev.off()




