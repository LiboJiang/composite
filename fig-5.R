

#DIA 3
#HT 8 F 60


SNP_P_2 <- (marker_st$marker)[,60]
SNP_P_2<-as.character(SNP_P_2)
snp.type_P_2 <- names(table(SNP_P_2))
snp.type_P_2



par_F_1<-c(0.341185178,  -0.089094437,   0.092596674,  -0.036587244,   0.033783474,  -0.043541718,   0.037673553,  -0.021033835)
par_F_2<-c(0.334439352,  -0.070404165,   0.072595764,  -0.013735935,   0.009810583,  -0.025486838,   0.017950544,  -0.005197664)



par_HT_1<-c(16.216266368, 14.866226115,   0.323050281,   6.409429522,   7.689820668 ,  0.984480829)
par_HT_2<-c(17.815207433,  14.073109663,   0.223720333,   8.371496194,  8.235158926,   0.708394290)



par_DIA_1<-c(22.380732811,  13.635399119 ,  0.538552438,  12.814511498, 260.825901557,   0.274356991)
par_DIA_2<-c(20.813399258,  12.799217654,   0.643287252,  13.028334204, 342.053298356,   0.321299928)
par_DIA_3<-c(21.690352918,  10.834594390,   0.514430742,  11.266581847, 326.202369212,   0.303955134)





par_V<-function(V_1_1_1){
  data_pheno_P<-V_1_1_1*1000000
  T_P<-c(1:24)
  
  plot(1:24,data_pheno_P,pch=16,xlab="Time",ylab="HT",xlim=c(1,24),ylim=c(0,1000000),type="l")
  par0<-c( 7.729790e+05, 4.913101e+02, 3.074683e-01, 2.436702e+05, 1.416980e+02, 7.941173e-01)
  #for (k in 1:10){
  fn_P<-function(t,par){
    par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
  }
  loss_P<-function(par,t,y){
    sum((y-fn_P(t,par))^2)
  }
  g_V1<-optim(par0,loss_P,t=1:24,y=data_pheno_P,method="BFGS",control=list(maxit=10000))#,trace = TRUE
  #par0<-c(g_V1$par)
  #}
  lines(1:24,fn_P(1:24,g_V1$par),col="red") 
  lines(1:24,fn_P_1(1:24,g_V1$par[1:3]),col="red",lty=2)
  lines(1:24,fn_P_1(1:24,g_V1$par[4:6]),col="red",lty=3)
  g_V1$par
}



V_1_1_1<-Legendre.model(1:24,par_F_1)*fn_P(1:24,par_HT_1)*(fn_P(1:24,par_DIA_1)/100)^2
V_1_1_2<-Legendre.model(1:24,par_F_1)*fn_P(1:24,par_HT_1)*(fn_P(1:24,par_DIA_2)/100)^2
V_1_1_3<-Legendre.model(1:24,par_F_1)*fn_P(1:24,par_HT_1)*(fn_P(1:24,par_DIA_3)/100)^2
V_1_2_1<-Legendre.model(1:24,par_F_1)*fn_P(1:24,par_HT_2)*(fn_P(1:24,par_DIA_1)/100)^2
V_1_2_2<-Legendre.model(1:24,par_F_1)*fn_P(1:24,par_HT_2)*(fn_P(1:24,par_DIA_2)/100)^2
V_1_2_3<-Legendre.model(1:24,par_F_1)*fn_P(1:24,par_HT_2)*(fn_P(1:24,par_DIA_3)/100)^2


V_2_1_1<-Legendre.model(1:24,par_F_2)*fn_P(1:24,par_HT_1)*(fn_P(1:24,par_DIA_1)/100)^2
V_2_1_2<-Legendre.model(1:24,par_F_2)*fn_P(1:24,par_HT_1)*(fn_P(1:24,par_DIA_2)/100)^2
V_2_1_3<-Legendre.model(1:24,par_F_2)*fn_P(1:24,par_HT_1)*(fn_P(1:24,par_DIA_3)/100)^2
V_2_2_1<-Legendre.model(1:24,par_F_2)*fn_P(1:24,par_HT_2)*(fn_P(1:24,par_DIA_1)/100)^2
V_2_2_2<-Legendre.model(1:24,par_F_2)*fn_P(1:24,par_HT_2)*(fn_P(1:24,par_DIA_2)/100)^2
V_2_2_3<-Legendre.model(1:24,par_F_2)*fn_P(1:24,par_HT_2)*(fn_P(1:24,par_DIA_3)/100)^2





par_1_1_1<-par_V(V_1_1_1)
par_1_1_2<-par_V(V_1_1_2)
par_1_1_3<-par_V(V_1_1_3)
par_1_2_1<-par_V(V_1_2_1)
par_1_2_2<-par_V(V_1_2_2)
par_1_2_3<-par_V(V_1_2_3)



par_2_1_1<-par_V(V_2_1_1)
par_2_1_2<-par_V(V_2_1_2)
par_2_1_3<-par_V(V_2_1_3)
par_2_2_1<-par_V(V_2_2_1)
par_2_2_2<-par_V(V_2_2_2)
par_2_2_3<-par_V(V_2_2_3)





max_value_DIA<-function(par_DIA_1){
  para_plot<- par_DIA_1#par_1[15:20]
  
  t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
  t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
  t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max
  
  v1<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))
  v2<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[1:3]))
  v3<-max(fn_P_1(seq(0,24,0.01),para_plot[4:6]))
  
  v<-cbind(v1,v2,v3)
  v
}


max_DIA_1<-max_value_DIA(par_DIA_1)
max_DIA_2<-max_value_DIA(par_DIA_2)
max_DIA_3<-max_value_DIA(par_DIA_3)

#32.11123#21.4844#11.29788


max_value_HT<-function(par_HT_1){
  
  para_plot<- par_HT_1#par_1[15:20]
  
  t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
  t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
  t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
  
  v1<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))
  v2<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[4:6]))
  v3<-max(fn_P_1(seq(0,24,0.01),para_plot[1:3]))
  
  v<-cbind(v1,v2,v3)
  v
}

max_HT_1<-max_value_HT(par_HT_1)
max_HT_2<-max_value_HT(par_HT_2)
#25.09084 #8.036043 #16.71934

max_value_V<-function(par_1_1_1){
  
  para_plot<- par_1_1_1#par_1[15:20]
  
  t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
  t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
  t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max
  
  v1<-max(fn_P(seq(0,24,0.01),para_plot[1:6]))
  v2<-max(fn_P_1(seq(0,t1_max,0.01),para_plot[4:6]))
  v3<-max(fn_P_1(seq(0,24,0.01),para_plot[1:3]))
  
  v<-cbind(v1,v2,v3)
  v
}
max_1<-max_value_V(par_1_1_1)
max_2<-max_value_V(par_1_1_2)
max_3<-max_value_V(par_1_1_3)
max_4<-max_value_V(par_1_2_1)
max_5<-max_value_V(par_1_2_2)
max_6<-max_value_V(par_1_2_3)

max_7<-max_value_V(par_2_1_1)
max_8<-max_value_V(par_2_1_2)
max_9<-max_value_V(par_2_1_3)
max_10<-max_value_V(par_2_2_1)
max_11<-max_value_V(par_2_2_2)
max_12<-max_value_V(par_2_2_3)

max(max_1[1],max_2[1],max_3[1],max_4[1],max_5[1],max_6[1],
    max_7[1],max_8[1],max_9[1],max_10[1],max_11[1],max_12[1])

max(max_1[2],max_2[2],max_3[2],max_4[2],max_5[2],max_6[2],
    max_7[2],max_8[2],max_9[2],max_10[2],max_11[2],max_12[2])

max(max_1[3],max_2[3],max_3[3],max_4[3],max_5[3],max_6[3],
    max_7[3],max_8[3],max_9[3],max_10[3],max_11[3],max_12[3])

#828758.1#247928.9#636403.8



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


pdf("1.pdf",width=20,height =12)

height<-267
length<-218

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length),ylim=c(0,height)); 

a<-29.5
b<-65
d<-2
###################################################################
sub_rc<-c(6,70,33.5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.8,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  text(sub_rc[1]-3.2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

text(sub_rc[1]-9,72.5,expression("Stem Form Factor" ),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-14,expression("Age (year)"),cex=2,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<-par_F_2 #par_1[1:8]

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.1)*56+sub_rc[4]+2,col="green3",lwd=4)



sub_rc<-c(6,70+b,33.5,10+b)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}

for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+ 11.2*i,sub_rc[1]+0.8,sub_rc[4]+2+ 11.2*i,font=2,lwd=1)
  text(sub_rc[1]-3.2,sub_rc[4]+2+ 11.2*i,0.2*i+0.1, cex=1.6,font=1)
}

text(sub_rc[1]+2,sub_rc[2]-3,expression("CT" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<-par_F_1 #par_1[1:8]

lines(seq(1,24,0.01)+sub_rc[1]+1,(Legendre.model(seq(1,24,0.01),para_plot)-0.1)*56+sub_rc[4]+2,col="green3",lwd=4)



###################################################################
sub_rc<-c(6+d+a*5/2,70+5+b*2,33.5+d+a*5/2,10+5+b*2)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.8,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm) ),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("AA" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_DIA_1#par_1[15:20]

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],32.11123*1.4+sub_rc[4]+2,sub_rc[3],32.11123*1.4+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],21.4844*1.4+sub_rc[4]+2,sub_rc[3],21.4844*1.4+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],11.29788*1.4+sub_rc[4]+2,sub_rc[3],11.29788*1.4+sub_rc[4]+2,lwd=2,col="blue",lty=2)

sub_rc<-c(6+d+a*7/2,70+5+b*2,33.5+d+a*7/2,10+5+b*2)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.8,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  #text(sub_rc[1]-3.2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}

#text(sub_rc[1]-9,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm) ),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("AG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_DIA_2#par_1[15:20]

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],32.11123*1.4+sub_rc[4]+2,sub_rc[3],32.11123*1.4+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],21.4844*1.4+sub_rc[4]+2,sub_rc[3],21.4844*1.4+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],11.29788*1.4+sub_rc[4]+2,sub_rc[3],11.29788*1.4+sub_rc[4]+2,lwd=2,col="blue",lty=2)




sub_rc<-c(6+d+a*9/2,70+5+b*2,33.5+d+a*9/2,10+5+b*2)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1.5,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+10*1.4*i,sub_rc[1]+0.8,sub_rc[4]+2+10*1.4*i,font=2,lwd=1)
  #text(sub_rc[1]-3.2,sub_rc[4]+2+10*1.4*i,10*i, cex=1.6,font=1)
}

#text(sub_rc[1]-9,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("DBH " (cm) ),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("GG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_DIA_3#par_1[15:20]

t1_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t2_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[5])/para_plot[6]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*1.4+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[1:3])*1.4+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[4:6])*1.4+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],32.11123*1.4+sub_rc[4]+2,sub_rc[3],32.11123*1.4+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],21.4844*1.4+sub_rc[4]+2,sub_rc[3],21.484*1.4+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],11.29788*1.4+sub_rc[4]+2,sub_rc[3],11.29788*1.4+sub_rc[4]+2,lwd=2,col="blue",lty=2)


#######################################################################################


sub_rc<-c(6+d+a*3,70+5+b*3,33.5+d+a*3,10+5+b*3)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height " (m)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CC" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_HT_1#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],25.09084*2.2+sub_rc[4]+2,sub_rc[3],25.09084*2.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],8.036043*2.2+sub_rc[4]+2,sub_rc[3],8.036043*2.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],16.71934*2.2+sub_rc[4]+2,sub_rc[3],16.71934*2.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)




sub_rc<-c(6+d+a*4,70+5+b*3,33.5+d+a*4,10+5+b*3)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+5*i,sub_rc[4]-3,5*i, cex=1.6,font=1)
}


for(i in 0:5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.3,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,5*i, cex=1.6,font=1)
}

#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem Height " (m)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_HT_2#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*2.2+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*2.2+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*2.2+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],25.09084*2.2+sub_rc[4]+2,sub_rc[3],25.09084*2.2+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],8.036043*2.2+sub_rc[4]+2,sub_rc[3],8.036043*2.2+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],16.71934*2.2+sub_rc[4]+2,sub_rc[3],16.71934*2.2+sub_rc[4]+2,lwd=2,col="blue",lty=2)
####################################################################################



sub_rc<-c(6+d+a*1,70+b*1,33.5+d+a*1,10+b*1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CT/CC/AA" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_1_1_1#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)




sub_rc<-c(6+d+a*2,70+b*1,33.5+d+a*2,10+b*1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CT/CC/AG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_1_1_2#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)



sub_rc<-c(6+d+a*3,70+b*1,33.5+d+a*3,10+b*1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CT/CC/GG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_1_1_3#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)






sub_rc<-c(6+d+a*4,70+b*1,33.5+d+a*4,10+b*1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CT/CG/AA" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_1_2_1#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)



sub_rc<-c(6+d+a*5,70+b*1,33.5+d+a*5,10+b*1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CT/CG/AG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_1_2_2#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)


sub_rc<-c(6+d+a*6,70+b*1,33.5+d+a*6,10+b*1)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[3],sub_rc[4]+2+11*i,sub_rc[3]-0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[3]+2.8,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[3],sub_rc[4]+2+11*i,sub_rc[3]-0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[3]+2.8,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("CT/CG/GG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_1_2_3#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")

segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)




##########################################################################################################################


sub_rc<-c(6+d+a*1,70,33.5+d+a*1,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  # text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



text(33.5+d+a*6+8,72.5,expression("Stem wood volume " (m^3)),cex=2,srt=270,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT/CC/AA" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_2_1_1#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)


sub_rc<-c(6+d+a*2,70,33.5+d+a*2,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT/CC/AG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_2_1_2#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

sub_rc<-c(6+d+a*3,70,33.5+d+a*3,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT/CC/GG" ),cex=1.6,adj=0,family="Times New Roman")
text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-14,expression("Age (year)"),cex=2,family="Times New Roman")
para_plot<- par_2_1_3#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)







sub_rc<-c(6+d+a*4,70,33.5+d+a*4,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+1+5*i,sub_rc[4]-6,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT/CG/AA" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_2_2_1#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)


sub_rc<-c(6+d+a*5,70,33.5+d+a*5,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[1],sub_rc[4]+2+11*i,sub_rc[1]+0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  #text(sub_rc[1]-2.4,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT/CG/AG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_2_2_2#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)

sub_rc<-c(6+d+a*6,70,33.5+d+a*6,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)


for(i in 0:5){
  segments(sub_rc[1]+1+5*i,sub_rc[4],sub_rc[1]+1+5*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+1+5*i,sub_rc[4]-4,5*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[3],sub_rc[4]+2+11*i,sub_rc[3]-0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[3]+2.8,sub_rc[4]+2+11*i,0.2*i, cex=1.6,font=1)
}

for(i in 5){
  segments(sub_rc[3],sub_rc[4]+2+11*i,sub_rc[3]-0.80,sub_rc[4]+2+11*i,font=2,lwd=1)
  text(sub_rc[3]+2.8,sub_rc[4]+2+11*i,"1.0", cex=1.6,font=1)
}



#text(sub_rc[1]-7,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Stem wood volume " (m^3)),cex=2,srt=90,family="Times New Roman")
text(sub_rc[1]+2,sub_rc[2]-3,expression("TT/CG/GG" ),cex=1.6,adj=0,family="Times New Roman")

para_plot<- par_2_2_3#par_1[15:20]

t2_max<- -(log(1/0.96-1)-log(para_plot[2]))/para_plot[3]
t1_max<- -(log(1/0.96-1)-log(para_plot[5]))/para_plot[6]
t2_min<-2*log(para_plot[2])/para_plot[3]-t2_max

lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P(seq(0,24,0.01),para_plot[1:6])*11/200000+sub_rc[4]+2,lwd=4,col="green3")
lines(seq(0,t1_max,0.01)+sub_rc[1]+1,fn_P_1(seq(0,t1_max,0.01),para_plot[4:6])*11/200000+sub_rc[4]+2,lwd=4,col="red")
lines(seq(0,24,0.01)+sub_rc[1]+1,fn_P_1(seq(0,24,0.01),para_plot[1:3])*11/200000+sub_rc[4]+2,lwd=4,col="blue")
segments(sub_rc[1],828758.1*11/200000+sub_rc[4]+2,sub_rc[3],828758.1*11/200000+sub_rc[4]+2,lwd=2,col="green3",lty=2)
segments(sub_rc[1],247928.9*11/200000+sub_rc[4]+2,sub_rc[3],247928.9*11/200000+sub_rc[4]+2,lwd=2,col="red",lty=2)
segments(sub_rc[1],636403.8*11/200000+sub_rc[4]+2,sub_rc[3],636403.8*11/200000+sub_rc[4]+2,lwd=2,col="blue",lty=2)


segments(35.4,-5,35.4,height+5,lwd=2)

segments(-5,140,length+5,140,lwd=2)
dev.off()
