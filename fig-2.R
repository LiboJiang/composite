

par<-final_par[145255,-1]
V_1<-Legendre.model(1:24,par[1:8])*fn_P(1:24,NH1_P_V_HT$par[1:6])*(fn_P(1:24,par[15:20])/100)^2
V_2<-Legendre.model(1:24,par[21:28])*fn_P(1:24,NH1_P_V_HT$par[7:12])*(fn_P(1:24,par[35:40])/100)^2
V_3<-Legendre.model(1:24,par[41:48])*fn_P(1:24,NH1_P_V_HT$par[13:18])*(fn_P(1:24,par[55:60])/100)^2




fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
}

loss1_P<-function(par,t,y_P){
  par.covar1_P <- par[7:8]
  SAD1_P <- AR1.get_mat(par.covar1_P,t)
  fy1_P<-dmvnorm(y_P,fn_P(t,par[1:6]),SAD1_P)
  LL_P  <- -sum(log(fy1_P))
  LL_P
} 


data_pheno_P<-V_1[1:8]

plot(1:8,data_pheno_P,pch=16,xlab="Time",ylab="HT",xlim=c(1,24),ylim=c(0,1),type="l")
par0<-c(0.1,0.1,0.1)
for (k in 1:10){
  fn_P<-function(t,par){
    par[1]/(1+ par[2]*exp(-par[3]*t))
  }
  loss_P<-function(par,t,y){
    sum((y-fn_P(t,par))^2)
  }
  g_P_HT<-optim(par0,loss_P,t=1:8,y=data_pheno_P,method="SANN",control=list(maxit=10000))
  par0<-c(g_P_HT$par)
}
lines(1:8,fn_P(1:8,g_P_HT$par),col="red") 


data_pheno_P<-V_1*1000000
#data_pheno_P<-V_1*10
T_P<-c(1:24)

plot(1:24,data_pheno_P,pch=16,xlab="Time",ylab="HT",xlim=c(1,24),ylim=c(0,1000000),type="l")
par0<-c( 7.729790e+05, 4.913101e+02, 3.074683e-01, 2.436702e+05, 1.416980e+02, 7.941173e-01)
#par0<-c(8.2639660, 403.3804382,   0.2952012,   2.4330644, 113.4960326,   0.7638803)
for (k in 1:10){
  fn_P<-function(t,par){
    par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
  }
  loss_P<-function(par,t,y){
    sum((y-fn_P(t,par))^2)
  }
  g_V1<-optim(par0,loss_P,t=1:24,y=data_pheno_P,method="BFGS",control=list(maxit=10000))#,trace = TRUE
  par0<-c(g_V1$par)
}
lines(1:24,fn_P(1:24,g_V1$par),col="red") 
lines(1:24,fn_P_1(1:24,g_V1$par[1:3]),col="red",lty=2)

lines(1:24,fn_P_1(1:24,g_V1$par[4:6]),col="red",lty=3)


data_pheno_P<-V_2*1000000
#data_pheno_P<-V_2*10
T_P<-c(1:24)

plot(1:24,data_pheno_P,pch=16,xlab="Time",ylab="HT",xlim=c(1,24),ylim=c(0,1000000),type="l")
par0<-c( 7.729790e+05, 4.913101e+02, 3.074683e-01, 2.436702e+05, 1.416980e+02, 7.941173e-01)
#par0<-c(8.2639660, 403.3804382,   0.2952012,   2.4330644, 113.4960326,   0.7638803)
for (k in 1:10){
  fn_P<-function(t,par){
    par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
  }
  loss_P<-function(par,t,y){
    sum((y-fn_P(t,par))^2)
  }
  g_V2<-optim(par0,loss_P,t=1:24,y=data_pheno_P,method="BFGS",control=list(maxit=10000))#,trace = TRUE
  par0<-c(g_V2$par)
}
lines(1:24,fn_P(1:24,g_V2$par),col="red") 

lines(1:24,fn_P_1(1:24,g_V2$par[1:3]),col="red",lty=2)

lines(1:24,fn_P_1(1:24,g_V2$par[4:6]),col="red",lty=3)



data_pheno_P<-V_3*1000000
#data_pheno_P<-V_3*10
T_P<-c(1:24)

plot(1:24,data_pheno_P,pch=16,xlab="Time",ylab="HT",xlim=c(1,24),ylim=c(0,1000000),type="l")
par0<-c( 7.729790e+05, 4.913101e+02, 3.074683e-01, 2.436702e+05, 1.416980e+02, 7.941173e-01)
#par0<-c(8.2639660, 403.3804382,   0.2952012,   2.4330644, 113.4960326,   0.7638803)
for (k in 1:10){
  fn_P<-function(t,par){
    par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
  }
  loss_P<-function(par,t,y){
    sum((y-fn_P(t,par))^2)
  }
  g_V3<-optim(par0,loss_P,t=1:24,y=data_pheno_P,method="BFGS",control=list(maxit=10000))#,trace = TRUE
  par0<-c(g_V3$par)
}
lines(1:24,fn_P(1:24,g_V3$par),col="red") 


lines(1:24,fn_P_1(1:24,g_V3$par[1:3]),col="red",lty=2)

lines(1:24,fn_P_1(1:24,g_V3$par[4:6]),col="red",lty=3)


library(showtext)
showtext_auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")


pdf("fig3.pdf",width=13.6,height =15)

height<-298
length<-90

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 


sub_rc<-c(6,102,33.5,12)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

par<-final_par[145255,-1]
#V_11<-Legendre.model(seq(1,26,0.01),par[1:8])*fn_P(seq(1,26,0.01),NH1_P_V_HT$par[1:6])*(fn_P(seq(1,26,0.01),par[15:20])/100)^2
#lines(seq(1,24,0.01)+sub_rc[1]+1-0.5,V_11[1:2301]*1000000/10000+sub_rc[4]+2,lwd=2,col="green3")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.30,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,"1.0", cex=1.6,font=1)
}
#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Wood Volume " (m^3) ),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-3,expression("Stemwood Volume " (m^3) ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]-6,sub_rc[2]+3,expression("D" ),cex=2,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=2,family="Times New Roman")
#text(sub_rc[1]+2,sub_rc[2]-4,expression(D),cex=1.6,font=1)
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,305,expression("CC"),cex=2,family="Times New Roman")

para_plot<-g_V1$par
ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2/12000

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6/12000

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6/12000


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2/12000

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6/12000

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6/12000

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+0.1,gi_1+sub_rc[4]+6,expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)


segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+0.1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4]-1,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]-0.5,gd_2+sub_rc[4]+5,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+3+84-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+2+84-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+73+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

P1<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[4:6]))/12000+sub_rc[4]+2
P2<-max(fn_P_1(seq(0,24,0.01),para_plot[1:3]))/12000+sub_rc[4]+2
P3<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))/12000+sub_rc[4]+2


segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##############################################################################################

sub_rc<-c(6+29.5,102,33.5+29.5,12)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

par<-final_par[145255,-1]
#V_11<-Legendre.model(seq(1,26,0.01),par[1:8])*fn_P(seq(1,26,0.01),NH1_P_V_HT$par[1:6])*(fn_P(seq(1,26,0.01),par[15:20])/100)^2
#lines(seq(1,24,0.01)+sub_rc[1]+1-0.5,V_11[1:2301]*1000000/10000+sub_rc[4]+2,lwd=2,col="green3")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}


#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Wood Volume " (m^3) ),cex=2,srt=90,font=1)

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=2,family="Times New Roman")
#text(sub_rc[1]+2,sub_rc[2]-4,expression(D),cex=1.6,font=1)
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,305,expression("CT"),cex=2,family="Times New Roman")

para_plot<-g_V2$par

ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2/12000

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6/12000

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6/12000


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2/12000

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6/12000

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6/12000

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+0.1,gi_1+sub_rc[4]+6,expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)


segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+0.1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4]-1,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+1.5,gd_2+sub_rc[4]-3,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+3+84-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+2+84-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+73+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
#############################################################################################

sub_rc<-c(6+29.5+29.5,102,33.5+29.5+29.5,12)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}


#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Wood Volume " (m^3) ),cex=2,srt=90,font=1)

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=2,family="Times New Roman")
#text(sub_rc[1]+2,sub_rc[2]-4,expression(D),cex=1.6,font=1)
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,305,expression("TT"),cex=2,family="Times New Roman")

para_plot<-g_V3$par
ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2/12000

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6/12000

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6/12000


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2/12000

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6/12000

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6/12000

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+0.1,gi_1+sub_rc[4]+6,expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)


segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+0.1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4]-1,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+1.5,gd_2+sub_rc[4]-3,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+3+84-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+2+84-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+73+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
###################################################################################################
SNP_P <- (dat_P_st$geno_table_P)[,145255]
SNP_P<-as.character(SNP_P)
snp.type_P <- names(table(SNP_P))

miss.type_P <- grep("\\.",snp.type_P)
if(length(miss.type_P)>0){
  snp.type_P <- snp.type_P[-miss.type_P]
}else{
  snp.type_P <- snp.type_P
}

sub_rc<-c(6,102+65,33.5,12+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot1<- par[1:8]

for(i in 1:dim(dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[1]){
  lines(c(1:24)+sub_rc[1]+1,((dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Form Factor" ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]-6,sub_rc[2]+3,expression("C" ),cex=2,family="Times New Roman")

para_plot<- par[15:20]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)

para_plot<- H0_P_V1$par[1:6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)


rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

###############################################################################################################3

sub_rc<-c(6+29.5,102+65,33.5+29.5,12+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


para_plot1<- par[21:28]
for(i in 1:dim(dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[1]){
  lines(c(1:24)+sub_rc[1]+1,((dat_P_st$psi[which(SNP_P==snp.type_P[2]),])[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

para_plot<- par[35:40]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)


para_plot<- H0_P_V2$par[1:6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)


rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
########################################################################################################

sub_rc<-c(6+29.5+29.5,102+65,33.5+29.5+29.5,12+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot1<- par[41:48]
for(i in 1:dim(dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[1]){
  lines(c(1:24)+sub_rc[1]+1,((dat_P_st$psi[which(SNP_P==snp.type_P[3]),])[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}
para_plot<- par[55:60]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)

para_plot<- H0_P_V3$par[1:6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##################################################################################################################


sub_rc<-c(6,102+65+65,33.5,12+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2]+3,expression("B" ),cex=2,family="Times New Roman")
para_plot<- par[15:20]

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


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+2,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+6,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

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
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)


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
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43-3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

P1<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[1:3]))*1.4+sub_rc[4]+2
P2<-max(fn_P_1(seq(0,24,0.01),para_plot[4:6]))*1.4+sub_rc[4]+2
P3<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))*1.4+sub_rc[4]+2


segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")


##################################################################################################################


sub_rc<-c(6+29.5,102+65+65,33.5+29.5,12+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


para_plot<- par[35:40]


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


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+2,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+6,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+0.1,gd_1+sub_rc[4]+4,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)


region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])*1.4+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)



arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+2.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43-3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
########################################################################################################


sub_rc<-c(6+29.5+29.5,102+65+65,33.5+29.5+29.5,12+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}

para_plot<- par[55:60]


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


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+2,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+6,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

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
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+2.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43-3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
############################################################################################################


sub_rc<-c(6,102+65+65+65,33.5,12+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height " (m)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2]+3,expression("A" ),cex=2,family="Times New Roman")
para_plot<- H0_P_V1$par[1:6]

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
text(ti_2+sub_rc[1]+2.5,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2.5,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4],expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

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
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4,sub_rc[4]+2+43+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

P1<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[4:6]))*2.2+sub_rc[4]+2
P2<-max(fn_P_1(seq(0,24,0.01),para_plot[1:3]))*2.2+sub_rc[4]+2
P3<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))*2.2+sub_rc[4]+2


segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")

###########################################################################################################


sub_rc<-c(6+29.5,102+65+65+65,33.5+29.5,12+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

para_plot<- H0_P_V2$par[1:6]


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
text(ti_2+sub_rc[1]+2.5,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2.5,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4],expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

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
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4,sub_rc[4]+2+43+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
#############################################################################################################

sub_rc<-c(6+29.5+29.5,102+65+65+65,33.5+29.5+29.5,12+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}


para_plot<- H0_P_V3$par[1:6]


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
text(ti_2+sub_rc[1]+2.5,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+8,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2.5,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+0.2,gd_1+sub_rc[4]+6,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

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
text(ta_1+sub_rc[1]+1+0.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+3,sub_rc[4]+2+43+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
dev.off()







jpeg("fig-2.jpeg",width=1360,height =1500, units = "px", pointsize = 12,quality =150)

height<-298
length<-90

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 


sub_rc<-c(6,102,33.5,12)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

par<-final_par[145255,-1]
#V_11<-Legendre.model(seq(1,26,0.01),par[1:8])*fn_P(seq(1,26,0.01),NH1_P_V_HT$par[1:6])*(fn_P(seq(1,26,0.01),par[15:20])/100)^2
#lines(seq(1,24,0.01)+sub_rc[1]+1-0.5,V_11[1:2301]*1000000/10000+sub_rc[4]+2,lwd=2,col="green3")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.30,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,"1.0", cex=1.6,font=1)
}
#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Wood Volume " (m^3) ),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-3,expression("Stemwood Volume " (m^3) ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]-6,sub_rc[2]+3,expression("D" ),cex=2,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=2,family="Times New Roman")
#text(sub_rc[1]+2,sub_rc[2]-4,expression(D),cex=1.6,font=1)
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,305,expression("CC"),cex=2,family="Times New Roman")

para_plot<-g_V1$par
ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2/12000

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6/12000

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6/12000


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2/12000

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6/12000

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6/12000

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+0.1,gi_1+sub_rc[4]+6,expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)


segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+0.1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4]-1,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]-0.5,gd_2+sub_rc[4]+5,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+3+84-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+2+84-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+73+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

P1<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[4:6]))/12000+sub_rc[4]+2
P2<-max(fn_P_1(seq(0,24,0.01),para_plot[1:3]))/12000+sub_rc[4]+2
P3<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))/12000+sub_rc[4]+2


segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##############################################################################################

sub_rc<-c(6+29.5,102,33.5+29.5,12)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

par<-final_par[145255,-1]
#V_11<-Legendre.model(seq(1,26,0.01),par[1:8])*fn_P(seq(1,26,0.01),NH1_P_V_HT$par[1:6])*(fn_P(seq(1,26,0.01),par[15:20])/100)^2
#lines(seq(1,24,0.01)+sub_rc[1]+1-0.5,V_11[1:2301]*1000000/10000+sub_rc[4]+2,lwd=2,col="green3")

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}


#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Wood Volume " (m^3) ),cex=2,srt=90,font=1)

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=2,family="Times New Roman")
#text(sub_rc[1]+2,sub_rc[2]-4,expression(D),cex=1.6,font=1)
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,305,expression("CT"),cex=2,family="Times New Roman")

para_plot<-g_V2$par

ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2/12000

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6/12000

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6/12000


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2/12000

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6/12000

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6/12000

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+0.1,gi_1+sub_rc[4]+6,expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)


segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+0.1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4]-1,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+1.5,gd_2+sub_rc[4]-3,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+3+84-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+2+84-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+73+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
#############################################################################################

sub_rc<-c(6+29.5+29.5,102,33.5+29.5+29.5,12)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}


#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Wood Volume " (m^3) ),cex=2,srt=90,font=1)

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=2,family="Times New Roman")
#text(sub_rc[1]+2,sub_rc[2]-4,expression(D),cex=1.6,font=1)
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,305,expression("TT"),cex=2,family="Times New Roman")

para_plot<-g_V3$par
ti_1<-log(para_plot[5])/para_plot[6]
gi_1<-para_plot[4]/2/12000

ta_1<-log(para_plot[5]*(2-sqrt(3)))/para_plot[6]
ga_1<-para_plot[4]*(3-sqrt(3))/6/12000

td_1<-log(para_plot[5]*(2+sqrt(3)))/para_plot[6]
gd_1<-para_plot[4]*(3+sqrt(3))/6/12000


ti_2<-log(para_plot[2])/para_plot[3]
gi_2<-para_plot[1]/2/12000

ta_2<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga_2<-para_plot[1]*(3-sqrt(3))/6/12000

td_2<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd_2<-para_plot[1]*(3+sqrt(3))/6/12000

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P(seq(24,26,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(24,26,0.01)+sub_rc[1]+1,fn_P_1(seq(24,26,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+0.1,gi_1+sub_rc[4]+6,expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)


segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+0.1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4]-1,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+1.5,gd_2+sub_rc[4]-3,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+2+84,td_1+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+3+84-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+2+84,td_2+sub_rc[1]+1,sub_rc[4]+2+84,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+2+84-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+73,ti_2+sub_rc[1]+1,sub_rc[4]+2+73,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+73+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)

region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])/12000+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
###################################################################################################
SNP_P <- (dat_P_st$geno_table_P)[,145255]
SNP_P<-as.character(SNP_P)
snp.type_P <- names(table(SNP_P))

miss.type_P <- grep("\\.",snp.type_P)
if(length(miss.type_P)>0){
  snp.type_P <- snp.type_P[-miss.type_P]
}else{
  snp.type_P <- snp.type_P
}

sub_rc<-c(6,102+65,33.5,12+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot1<- par[1:8]

for(i in 1:dim(dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[1]){
  lines(c(1:24)+sub_rc[1]+1,((dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Form Factor" ),cex=2,srt=90,family="Times New Roman")

text(sub_rc[1]-6,sub_rc[2]+3,expression("C" ),cex=2,family="Times New Roman")

para_plot<- par[15:20]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)

para_plot<- H0_P_V1$par[1:6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)


rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

###############################################################################################################3

sub_rc<-c(6+29.5,102+65,33.5+29.5,12+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


para_plot1<- par[21:28]
for(i in 1:dim(dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[1]){
  lines(c(1:24)+sub_rc[1]+1,((dat_P_st$psi[which(SNP_P==snp.type_P[2]),])[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

para_plot<- par[35:40]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)


para_plot<- H0_P_V2$par[1:6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)


rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
########################################################################################################

sub_rc<-c(6+29.5+29.5,102+65,33.5+29.5+29.5,12+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

para_plot1<- par[41:48]
for(i in 1:dim(dat_P_st$psi[which(SNP_P==snp.type_P[1]),])[1]){
  lines(c(1:24)+sub_rc[1]+1,((dat_P_st$psi[which(SNP_P==snp.type_P[3]),])[i,]-0.1)*56+sub_rc[4]+2,col="#C7E9C0",lty=1,lwd=1)
}

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.1)*56+sub_rc[4]+2,col="green3",lwd=2)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.3,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}
para_plot<- par[55:60]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)

para_plot<- H0_P_V3$par[1:6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
region_x<-seq(t2_min,t1_max,0.01)
region_y<-rep(60,length(region_x))
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20,col="gray60",lty=1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
##################################################################################################################


sub_rc<-c(6,102+65+65,33.5,12+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2]+3,expression("B" ),cex=2,family="Times New Roman")
para_plot<- par[15:20]

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


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+2,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+6,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

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
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)


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
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43-3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

P1<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[1:3]))*1.4+sub_rc[4]+2
P2<-max(fn_P_1(seq(0,24,0.01),para_plot[4:6]))*1.4+sub_rc[4]+2
P3<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))*1.4+sub_rc[4]+2


segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")


##################################################################################################################


sub_rc<-c(6+29.5,102+65+65,33.5+29.5,12+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


para_plot<- par[35:40]


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


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+2,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+6,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+0.1,gd_1+sub_rc[4]+4,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)


region_x<-seq(t2_min,t1_max,0.01)
region_y<-fn_P(seq(t2_min,t1_max,0.01),para_plot[1:6])*1.4+2
region_x<-c(region_x[1],region_x,tail(region_x,1))
region_y<-c(0,region_y,0)
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)



arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+2.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43-3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
########################################################################################################


sub_rc<-c(6+29.5+29.5,102+65+65,33.5+29.5+29.5,12+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}

para_plot<- par[55:60]


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


segments(ti_2+sub_rc[1]+1,sub_rc[4],ti_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ti_2+sub_rc[1]+1,gi_2+sub_rc[4]+2,pch=17)
segments(ti_1+sub_rc[1]+1,sub_rc[4],ti_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ti_1+sub_rc[1]+1,gi_1+sub_rc[4]+2,pch=17)
text(ti_1+sub_rc[1]+2,gi_1+sub_rc[4],expression(italic(t["I"])),family="Times New Roman1",cex=1.4)
text(ti_2+sub_rc[1]+2,gi_2+sub_rc[4]-1,expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="royalblue2")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+4,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1],ga_2+sub_rc[4]+6,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

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
polygon(region_x+sub_rc[1]+1,region_y+sub_rc[4],density = 20, angle = 135,col="gray60",lty=1)


arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=20,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=15,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=10,code=3)
arrows(ta_1+sub_rc[1]+1,sub_rc[4]+55,td_1+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="red",length=0.08,angle=5,code=3)
#points(ta_1+sub_rc[1]+1-0.5+1,sub_rc[4]+2+25+2.5,pch=2,cex=0.8)
text(ta_1+sub_rc[1]+1+2.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4+2.5,sub_rc[4]+2+43-3,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
############################################################################################################


sub_rc<-c(6,102+65+65+65,33.5,12+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height " (m)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2]+3,expression("A" ),cex=2,family="Times New Roman")
para_plot<- H0_P_V1$par[1:6]

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
text(ti_2+sub_rc[1]+2.5,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2.5,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4],expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

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
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4,sub_rc[4]+2+43+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

P1<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[4:6]))*2.2+sub_rc[4]+2
P2<-max(fn_P_1(seq(0,24,0.01),para_plot[1:3]))*2.2+sub_rc[4]+2
P3<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))*2.2+sub_rc[4]+2


segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")

###########################################################################################################


sub_rc<-c(6+29.5,102+65+65+65,33.5+29.5,12+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

para_plot<- H0_P_V2$par[1:6]


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
text(ti_2+sub_rc[1]+2.5,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+6,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2.5,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+2,gd_1+sub_rc[4],expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

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
text(ta_1+sub_rc[1]+1+1.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+4,sub_rc[4]+2+43+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
#############################################################################################################

sub_rc<-c(6+29.5+29.5,102+65+65+65,33.5+29.5+29.5,12+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}


para_plot<- H0_P_V3$par[1:6]


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
text(ti_2+sub_rc[1]+2.5,gi_2+sub_rc[4],expression(italic(T["I"])),family="Times New Roman1",cex=1.4)

segments(ta_2+sub_rc[1]+1,sub_rc[4],ta_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(ta_2+sub_rc[1]+1,ga_2+sub_rc[4]+2,pch=17)
segments(ta_1+sub_rc[1]+1,sub_rc[4],ta_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(ta_1+sub_rc[1]+1,ga_1+sub_rc[4]+2,pch=17)
text(ta_1+sub_rc[1],ga_1+sub_rc[4]+8,expression(italic(t["a"])),family="Times New Roman1",cex=1.4)
text(ta_2+sub_rc[1]+2.5,ga_2+sub_rc[4]-1,expression(italic(T["a"])),family="Times New Roman1",cex=1.4)

segments(td_2+sub_rc[1]+1,sub_rc[4],td_2+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="blue")
points(td_2+sub_rc[1]+1,gd_2+sub_rc[4]+2,pch=17)
segments(td_1+sub_rc[1]+1,sub_rc[4],td_1+sub_rc[1]+1,sub_rc[2],lty=3,lwd=1,col="red")
points(td_1+sub_rc[1]+1,gd_1+sub_rc[4]+2,pch=17)
text(td_1+sub_rc[1]+0.2,gd_1+sub_rc[4]+6,expression(italic(t["d"])),family="Times New Roman1",cex=1.4)
text(td_2+sub_rc[1]+2.5,gd_2+sub_rc[4]-1,expression(italic(T["d"])),family="Times New Roman1",cex=1.4)

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
text(ta_1+sub_rc[1]+1+0.5,sub_rc[4]+55-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=20,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=15,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=10,code=3)
arrows(ta_2+sub_rc[1]+1,sub_rc[4]+55,td_2+sub_rc[1]+1,sub_rc[4]+55,lwd=1,col="blue",length=0.08,angle=5,code=3)
#points(ta_2+sub_rc[1]+1-0.5+2+0.5,sub_rc[4]+2+45+2.5,pch=2,cex=0.8)
text(ta_2+sub_rc[1]+1+2+1,sub_rc[4]+55-3,expression(italic(Delta*T)),family="Times New Roman1",cex=1.4)

arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=20,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=15,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=10,code=3)
arrows(ti_1+sub_rc[1]+1,sub_rc[4]+2+43,ti_2+sub_rc[1]+1,sub_rc[4]+2+43,lwd=1,col="green3",length=0.08,angle=5,code=3)
#points(ti_1+sub_rc[1]+1-0.5+4+0.5,sub_rc[4]+2+38+2.5,pch=2,cex=0.8)
text(ti_1+sub_rc[1]+1+3,sub_rc[4]+2+43+2.5,expression(italic(L)),family="Times New Roman1",cex=1.4)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

segments(sub_rc[1],P1,sub_rc[3],P1,lty=2,lwd=1,col="black")
segments(sub_rc[1],P2,sub_rc[3],P2,lty=2,lwd=1,col="black")
segments(sub_rc[1],P3,sub_rc[3],P3,lty=2,lwd=1,col="black")
dev.off()






