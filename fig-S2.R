
library(mvtnorm)

SAD1<- function(par0, times_1,traits=1, options=list()){
  par<-par0;
  if (class(par0)=="list")
    par <- unlist(par0);
  
  t_len_1 <- length( times_1 );
  SAD.1 <- array(0, dim=c((t_len_1)*traits,(t_len_1)*traits));
  for (i0 in 1:traits)
    for (i1 in 1:traits)
    {
      if (i0==i1)
        for (k0 in 1:t_len_1)
          for (k1 in 1:t_len_1)
          {
            SAD.1[(i0-1)*t_len_1+k0,(i1-1)*t_len_1+k1] <- par[i0*2-1]^abs(k0-k1)*
              (1-par[i0*2-1]^(2*min(k0,k1)))/(1-par[i0*2-1]^2)*par[i0*2]^2;
          }
    }
  return(SAD.1);
}

###################################################################################################

fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))
}


loss1<-function(par,t,y_P){
  par.covar1 <- par[4:5]
  AR1_P <- SAD1(par.covar1,t,traits = 1)
  fy1_P<-c()
  for(i in 1:dim(y_P)[1]){
    y1<-y_P[i,]
    t1<-t
    t2<-t1[-which(y1==0)]
    y2<-y1[-which(y1==0)]
    
    AR2_P<-AR1_P[-which(y1==0),-which(y1==0)]
    
    t2_1<-t1[1:length(t1)][-which(y1[1:length(t1)]==0)]
    
    y3_1<-fn_P(t2_1,par[1:3])
    
    yy<-dmvnorm(y2,y3_1,AR2_P)
    fy1_P<-c(fy1_P,yy)
  }
  
  LL_P  <- -sum(log(fy1_P))
  LL_P
} 


par0<-c(21.0012640, 380.3726350 ,  0.2250795,0.8,100)
H0_len<-optim(par0,loss1,t=dat_P_st$TIME,y_P=dat_P_st$lateral_root_len,
              method="BFGS",control=list(maxit=10000))


par0<-c(15.1088471, 20.9177517,  0.1243092,  0.6966802 , 3.9718366)
H0_ave<-optim(par0,loss1,t=dat_P_st$TIME,y_P=dat_P_st$lateral_root_ave,
              method="BFGS",control=list(maxit=10000))



par0<-c(41.55090420, 117.02996102,   0.11183293,0.93537430, 3.68245873)
H0_num<-optim(par0,loss1,t=dat_P_st$TIME,y_P=dat_P_st$lateral_root_num,
              method="BFGS",control=list(maxit=10000))
###################################################################




library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

pdf("fig1.pdf",width=11.2,height =8)


height<-140
length<-118+8*2+20+5
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length-3),ylim=c(0,height)); 


sub_rc<-c(10,100,59+8,40)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

k<-1000/55


for(i in 0:5){
  segments(sub_rc[1]+0.5+10*i,sub_rc[4],sub_rc[1]+0.5+10*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+0.5+10*i,sub_rc[4]-3,10*i, cex=1.4,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+250/k*i,sub_rc[1]+0.6,sub_rc[4]+2+250/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+250/k*i,25*i, cex=1.4,font=1)
}




pheno<-dat_P_st$lateral_root_len
#pheno[which(pheno>1200)]<-NA
for(i in 1:dim(pheno)[1]){
  lines(dat_P_st$TIME[-which(pheno[i,]==0)]+sub_rc[1]+0.5,
        pheno[i,][-which(pheno[i,]==0)]/k+sub_rc[4]+2,col="#C7E9C0")
}

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),H0_len$par[1:3])/k+sub_rc[4]+2,col="green3",lwd=2) 
lines(seq(47,47+8,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+8,0.01),H0_len$par[1:3])/k+sub_rc[4]+2,col="green3",lwd=2,lty=2) 


para_plot<- H0_len$par
ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]-2,gi/k+sub_rc[4]+6,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]-2,ga/k+sub_rc[4]+6,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)

text(td+sub_rc[1]-2,gd/k+sub_rc[4]+6,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")

for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+40,td+sub_rc[1]+0.5,sub_rc[4]+5+40,lwd=1,col="red",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+40+5,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

text(sub_rc[1]-11,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Total Lateral Root Length " (cm)),cex=1.6,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-10,"Time (day)",cex=1.6,family="Times New Roman")

text(sub_rc[1]-4.5,sub_rc[2]+3,expression(A),cex=1.6,family="Times New Roman")


#########################################################################################



sub_rc<-c(77+25,70,134+25,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

pheno<-dat_P_st$lateral_root_num

k<-60/55


for(i in 0:5){
  segments(sub_rc[1]+0.5+10*i,sub_rc[4],sub_rc[1]+0.5+10*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+0.5+10*i,sub_rc[4]-3,10*i, cex=1.4,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+15/k*i,sub_rc[1]+0.6,sub_rc[4]+2+15/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+15/k*i,15*i, cex=1.4,font=1)
}


for(i in 1:dim(pheno)[1]){
  lines(dat_P_st$TIME[-which(pheno[i,]==0)]+sub_rc[1]+0.5,
        pheno[i,][-which(pheno[i,]==0)]/k+sub_rc[4]+2,col="#C7E9C0")
}

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),H0_num$par[1:3])/k+sub_rc[4]+2,col="green3",lwd=2) 
lines(seq(47,47+8,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+8,0.01),H0_num$par[1:3])/k+sub_rc[4]+2,col="green3",lwd=2,lty=2) 



para_plot<- H0_num$par
ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]-2,gi/k+sub_rc[4]+6,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]-2,ga/k+sub_rc[4]+6,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)

text(td+sub_rc[1]-2,gd/k+sub_rc[4]+6,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")


for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+40,td+sub_rc[1]+0.5,sub_rc[4]+5+40,lwd=1,col="red",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+40+5,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

text(sub_rc[1]-10,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Lateral Root Number "),cex=1.6,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-10,"Time (day)",cex=1.6,family="Times New Roman")

text(sub_rc[1]-4.5,sub_rc[2]+3,expression(C),cex=1.6,family="Times New Roman")



##################################################################################


sub_rc<-c(77+25,70+70,134+25,10+70)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

pheno<-dat_P_st$lateral_root_ave

k<-48/55

for(i in 0:5){
  segments(sub_rc[1]+0.5+10*i,sub_rc[4],sub_rc[1]+0.5+10*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+10*i,sub_rc[2]+3,30*i, cex=1.6,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+12/k*i,sub_rc[1]+0.6,sub_rc[4]+2+12/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+12/k*i,1.2*i, cex=1.4,font=1)
}


for(i in 1:dim(pheno)[1]){
  lines(dat_P_st$TIME[-which(pheno[i,]==0)]+sub_rc[1]+0.5,
        pheno[i,][-which(pheno[i,]==0)]/k+sub_rc[4]+2,col="#C7E9C0")
}

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),H0_ave$par[1:3])/k+sub_rc[4]+2,col="green3",lwd=2) 
lines(seq(47,47+8,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+8,0.01),H0_ave$par[1:3])/k+sub_rc[4]+2,col="green3",lwd=2,lty=2) 



para_plot<- H0_ave$par
ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]-2,gi/k+sub_rc[4]+6,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]-2,ga/k+sub_rc[4]+6,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)

text(td+sub_rc[1]-2,gd/k+sub_rc[4]+6,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")


for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+40,td+sub_rc[1]+0.5,sub_rc[4]+5+40,lwd=1,col="red",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+40+5,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

text(sub_rc[1]-10,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Average Lateral Length " (cm)),cex=1.6,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-10,"Time (day)",cex=2,family="Times New Roman")

text(sub_rc[1]-4.5,sub_rc[2]+3,expression(B),cex=1.6,family="Times New Roman")

arrows(67.2,70,87.2,40,lwd=2,col="black",length=0.3,angle=15,code=2,lty=2)

for(i in 1:20){
  arrows(87.2-3*20/30,43,87.2,40,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}


arrows(67.2,70,87.2,110,lwd=2,col="black",length=0.3,angle=15,code=2,lty=2)

for(i in 1:20){
  arrows(87.2-3*20/40,107,87.2,110,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}

dev.off()
