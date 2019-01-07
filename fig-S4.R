setwd("D:/Rdocument/p-simulation")
load("par_F_100_5.RData")
load("par_HT_100_5.RData")
load("par_DIA_100_5.RData")
load("par_V_100_5_1.RData")
load("par_F_200_10.RData")
load("par_HT_200_10.RData")
load("par_DIA_200_10.RData")
load("par_V_200_10_1.RData")



par_F_1<-colMeans(par_F_100_5)
par_HT_1<-colMeans(par_HT_100_5)
par_DIA_1<-colMeans(par_DIA_100_5)
V_1<-Legendre.model(1:24,par_F_1[1:8])*fn_P(1:24,par_HT_1[1:6])*(fn_P(1:24,par_DIA_1[1:6])/100)^2
V_2<-Legendre.model(1:24,par_F_1[9:16])*fn_P(1:24,par_HT_1[7:12])*(fn_P(1:24,par_DIA_1[7:12])/100)^2
V_3<-Legendre.model(1:24,par_F_1[17:24])*fn_P(1:24,par_HT_1[13:18])*(fn_P(1:24,par_DIA_1[13:18])/100)^2


par_F_2<-colMeans(par_F_200_10)
par_HT_2<-colMeans(par_HT_200_10)
par_DIA_2<-colMeans(par_DIA_200_10)
V1_1<-Legendre.model(1:24,par_F_2[1:8])*fn_P(1:24,par_HT_2[1:6])*(fn_P(1:24,par_DIA_2[1:6])/100)^2
V1_2<-Legendre.model(1:24,par_F_2[9:16])*fn_P(1:24,par_HT_2[7:12])*(fn_P(1:24,par_DIA_2[7:12])/100)^2
V1_3<-Legendre.model(1:24,par_F_2[17:24])*fn_P(1:24,par_HT_2[13:18])*(fn_P(1:24,par_DIA_2[13:18])/100)^2

V_mol_1<-Legendre.model(1:24,par_F1)*fn_P(1:24,par_HT1)*(fn_P(1:24,par_DIA1)/100)^2
V_mol_2<-Legendre.model(1:24,par_F2)*fn_P(1:24,par_HT2)*(fn_P(1:24,par_DIA2)/100)^2
V_mol_3<-Legendre.model(1:24,par_F3)*fn_P(1:24,par_HT3)*(fn_P(1:24,par_DIA3)/100)^2


par_V_1<-colMeans(par_V_100_5_1)
par_V_2<-colMeans(par_V_200_10_1)


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

g<-function(V_1){
  data_pheno_P<-V_1*1000000
  T_P<-c(1:24)
  
  plot(1:24,data_pheno_P,pch=16,xlab="Time",ylab="HT",xlim=c(1,24),ylim=c(0,1000000),type="l")
  par0<-c( 7.729790e+05, 4.913101e+02, 3.074683e-01, 2.436702e+05, 1.416980e+02, 7.941173e-01)
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
  g_V1$par
}
g_V1_par<-g(V_1)
g_V2_par<-g(V_2)
g_V3_par<-g(V_3)
g1_V1_par<-g(V1_1)
g1_V2_par<-g(V1_2)
g1_V3_par<-g(V1_3)
g_mol_V1_par<-g(V_mol_1)
g_mol_V2_par<-g(V_mol_2)
g_mol_V3_par<-g(V_mol_3)

fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
}

fn_P_1<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))
}

library(showtext)
showtext_auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

pdf("fig-8.pdf",width=13.6,height =18)

height<-385
length<-87

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height));

sub_rc<-c(6,100,32.5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,"1.0", cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-3,expression("Stemwood volume " (m^3) ),cex=1.8,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2],expression("E" ),cex=2,family="Times New Roman")
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=1.8,family="Times New Roman")

para_plot<-par_V_1[1:6]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*10/1.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*10/1.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*10/1.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)


para_plot<-par_V_2[1:6]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*10/1.2+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*10/1.2+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*10/1.2+sub_rc[4]+2,lwd=3,col="blue",lty=3)

para_plot1<-c(8.020638e+05, 3.427003e+02, 2.858530e-01, 2.417763e+05, 1.028924e+02, 6.619138e-01)

t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")

##########################################################################################################
sub_rc<-c(6+28.5,100,32.5+28.5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,"1.0", cex=1.6,font=1)
}

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=1.8,family="Times New Roman")

para_plot<-par_V_1[7:12]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*10/1.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*10/1.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*10/1.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)


para_plot<-par_V_2[7:12]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*10/1.2+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*10/1.2+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*10/1.2+sub_rc[4]+2,lwd=3,col="blue",lty=3)

para_plot1<-c(7.133049e+05, 3.549579e+02, 2.979209e-01, 2.029372e+05, 1.399852e+02, 8.026018e-01)

t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")

#####################################################################################################
sub_rc<-c(6+28.5+28.5,100,32.5+28.5+28.5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,"1.0", cex=1.6,font=1)
}

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-16,expression("Age (year)"),cex=1.8,family="Times New Roman")

para_plot<-par_V_1[13:18]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*10/1.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*10/1.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*10/1.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<-par_V_2[13:18]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*10/1.2+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*10/1.2+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*10/1.2+sub_rc[4]+2,lwd=3,col="blue",lty=3)


para_plot1<-c(7.273785e+05, 1.527081e+02, 2.557606e-01, 2.090475e+05, 4.766965e+02, 9.627867e-01)

t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")

################################################################################################
sub_rc<-c(6,100+95,32.5,10+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+200000/12000*i,"1.0", cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-3,expression("Stemwood volume " (m^3) ),cex=1.8,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2],expression("D" ),cex=2,family="Times New Roman")

para_plot<-g_V1_par

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<-g1_V1_par

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=3)

para_plot1<-g_mol_V1_par

t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])/12000+sub_rc[4]+2,lwd=3,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])/12000+sub_rc[4]+2,lwd=3,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])/12000+sub_rc[4]+2,lwd=3,col="blue")


###############################################################################################3
sub_rc<-c(6+28.5,100+95,32.5+28.5,10+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
}

para_plot<-g_V2_par

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<-g1_V2_par

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=3,col="blue",lty=3)


para_plot1<-g_mol_V2_par

t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")

###############################################################################
sub_rc<-c(6+28.5+28.5,100+95,32.5+28.5+28.5,10+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+200000/12000*i,sub_rc[1]+0.40,sub_rc[4]+2+200000/12000*i,font=2,lwd=1)
}

para_plot<-g_V3_par

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<-g1_V3_par

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max


lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])/12000+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])/12000+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])/12000+sub_rc[4]+2,lwd=3,col="blue",lty=3)


para_plot1<-g_mol_V3_par

t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])/12000+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])/12000+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])/12000+sub_rc[4]+2,lwd=2,col="blue")

###########################################################################################
sub_rc<-c(6,100+95+65,32.5,10+95+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+2+ 18*i,sub_rc[1]+0.4,sub_rc[4]+2+ 18*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+ 18*i,0.2*i+0.2, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem form factor" ),cex=1.8,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2]+3,expression("C" ),cex=2,family="Times New Roman")
para_plot<-par_F_1[1:8] #par_1[1:8]
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.2)*85+sub_rc[4]+2,col="green3",lwd=2,lty=2)

para_plot<-par_F_2[1:8] #par_1[1:8]
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.2)*85+sub_rc[4]+2,col="green3",lwd=3,lty=3)

para_plot1<- c(0.343763043, -0.107975321,  0.081683157, -0.003651312,  0.007474049, -0.031597981,  0.022869350, -0.004370809)
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.2)*85+sub_rc[4]+2,col="green3",lwd=2)

#################################################################################################
sub_rc<-c(6+28.5,100+95+65,32.5+28.5,10+95+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)
for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+2+ 18*i,sub_rc[1]+0.4,sub_rc[4]+2+ 18*i,font=2,lwd=1)
}

para_plot<- par_F_1[9:16]#par_1[21:28]

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.2)*85+sub_rc[4]+2,col="green3",lwd=2,lty=2)

para_plot<- par_F_2[9:16]#par_1[21:28]

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.2)*85+sub_rc[4]+2,col="green3",lwd=3,lty=3)

para_plot1<- c(0.336629708, -0.080396811,  0.093624247, -0.039648608,  0.034869980, -0.036736270,  0.016930957, -0.004751406)
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.2)*85+sub_rc[4]+2,col="green3",lwd=2)

#################################################################################################
sub_rc<-c(6+28.5+28.5,100+95+65,32.5+28.5+28.5,10+95+95)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1,font=1)
}

for(i in 0:3){
  segments(sub_rc[1],sub_rc[4]+2+ 18*i,sub_rc[1]+0.4,sub_rc[4]+2+ 18*i,font=2,lwd=1)
}

para_plot<- par_F_1[17:24]#par_1[41:48]

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.2)*85+sub_rc[4]+2,col="green3",lwd=2,lty=2)
para_plot<- par_F_2[17:24]#par_1[41:48]

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.2)*85+sub_rc[4]+2,col="green3",lwd=3,lty=3)

para_plot1<- c(0.339233974, -0.058399805,  0.061893012, -0.010961303,  0.002420895, -0.027216529,  0.046103764, -0.031198423)
lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot1)-0.2)*85+sub_rc[4]+2,col="green3",lwd=2)

################################################################################################

sub_rc<-c(6,100+95+65+65,32.5,10+95+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}


text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm)),cex=1.8,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2],expression("B" ),cex=2,family="Times New Roman")

para_plot<- par_DIA_1[1:6]#par_1[15:20]

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<- par_DIA_2[1:6]#par_1[15:20]

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=3,col="blue",lty=3)

para_plot1<- c(23.1952453,   7.6817903,   0.4415837,  12.5459971, 195.6246114,   0.2741480)

t1_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t2_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[5])/para_plot1[6]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue")


###############################################################################################
sub_rc<-c(6+28.5,100+95+65+65,32.5+28.5,10+95+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
}

para_plot<- par_DIA_1[7:12]#par_1[35:40]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<- par_DIA_2[7:12]#par_1[35:40]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=3,col="blue",lty=3)


para_plot1<- c(20.8669738,  12.4717218,   0.6150039,  12.4033256, 200.8813002,   0.2878006)

t1_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t2_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[5])/para_plot1[6]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue")

###############################################################################################

sub_rc<-c(6+28.5+28.5,100+95+65+65,32.5+28.5+28.5,10+95+95+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.4,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
}

para_plot<-par_DIA_1[13:18]#par_1[55:60]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<-par_DIA_2[13:18]#par_1[55:60]
t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=3,col="blue",lty=3)


para_plot1<- c( 21.0672222,  15.4764854 ,  0.6127538,  14.3016875, 203.5803395 ,  0.2757311)

t1_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t2_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_min1<-2*log(para_plot1[5])/para_plot1[6]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])*1.4+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[1:3])*1.4+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[4:6])*1.4+sub_rc[4]+2,lwd=2,col="blue")

############################################################################################
sub_rc<-c(6,100+95+65+65+65,32.5,10+95+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.4,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4]-3,expression("Stem height " (m)),cex=1.8,srt=90,family="Times New Roman")
text(sub_rc[1]-6,sub_rc[2],expression("A" ),cex=2,family="Times New Roman")
para_plot<- par_HT_1[1:6]#par_1[9:14]

t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<- par_HT_2[1:6]#par_1[9:14]

t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=3,col="blue",lty=3)

para_plot1<- c(15.2541915, 32.5000220,  0.2771275, 10.1702390,  4.5673914,  0.5111841)

t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue")

####################################################################################################
sub_rc<-c(6+28.5,100+95+65+65+65,32.5+28.5,10+95+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.4,sub_rc[4]+2+11*i,font=2,lwd=1)
}

para_plot<- par_HT_1[7:12]#par_1[29:34]

t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<- par_HT_2[7:12]#par_1[29:34]

t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=3,col="blue",lty=3)

para_plot1<- c(13.3489054, 33.1604386,  0.2845272, 11.4044683,  6.8261126,  0.6386955)

t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue")


########################################################################################
sub_rc<-c(6+28.5+28.5,100+95+65+65+65,32.5+28.5+28.5,10+95+95+65+65)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.4,sub_rc[4]+2+11*i,font=2,lwd=1)
}

para_plot<- par_HT_1[13:18]#par_1[49:54]

t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)

para_plot<- par_HT_2[13:18]#par_1[49:54]

t1_max<- -(log(1/0.9-1)-log(para_plot[5]))/para_plot[6]
t2_max<- -(log(1/0.9-1)-log(para_plot[2]))/para_plot[3]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=3,col="green3",lty=3)
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=3,col="red",lty=3)
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=3,col="blue",lty=3)
para_plot1<- c(18.2498348, 11.0692102,  0.3207969,  4.5062786, 20.2115823,  1.6530021)

t1_max1<- -(log(1/0.96-1)-log(para_plot1[5]))/para_plot1[6]
t2_max1<- -(log(1/0.96-1)-log(para_plot1[2]))/para_plot1[3]
t2_min1<-2*log(para_plot1[2])/para_plot1[3]-t2_max1

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot1[1:6])*2.2+sub_rc[4]+2,lwd=2,col="green3")
lines(seq(0,t1_max1,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max1,0.01),para_plot1[4:6])*2.2+sub_rc[4]+2,lwd=2,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot1[1:3])*2.2+sub_rc[4]+2,lwd=2,col="blue")

dev.off()