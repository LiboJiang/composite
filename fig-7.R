nn<-5739#which.max(LR)#7490

plot_len<-final_par_len_del[nn,-1]
plot_ot<-final_par_del[nn,-1]

SNP <- (dat_P_st$geno_table)[,nn]
SNP<-as.character(SNP)
snp.type <- names(table(SNP))

miss.type <- grep("-",snp.type)
if(length(miss.type)>0){
  snp.type <- snp.type[-miss.type]
}else{
  snp.type <- snp.type
}



library(showtext)
showtext.auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")

pdf("fig3.pdf",width=11.2,height =8)


height<-140
length<-120+8*2+20+5
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(0,length-3),ylim=c(0,height)); 

####################################################################################



sub_rc<-c(10,70,69,10)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

pheno<-dat_P_st$lateral_root_num

k<-48/55


for(i in 0:5){
  segments(sub_rc[1]+0.5+10*i,sub_rc[4],sub_rc[1]+0.5+10*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+0.5+10*i,sub_rc[4]-3,10*i, cex=1.4,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+12/k*i,sub_rc[1]+0.6,sub_rc[4]+2+12/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+12/k*i,12*i, cex=1.4,font=1)
}


para_plot<- plot_ot[4:6]

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),para_plot)/k+sub_rc[4]+2,col="red",lwd=2) 
lines(seq(47,47+10,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+10,0.01),para_plot)/k+sub_rc[4]+2,col="red",lwd=2,lty=2) 


ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]-2,gi/k+sub_rc[4]+5,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]-2,ga/k+sub_rc[4]+5,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)

text(td+sub_rc[1]-2,gd/k+sub_rc[4]+5,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")

for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+50,td+sub_rc[1]+0.5,sub_rc[4]+5+50,lwd=1,col="red",length=0.08,angle=i,code=3)
}
text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+50+3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)


para_plot<- plot_ot[4:6+6]

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),para_plot)/k+sub_rc[4]+2,col="green3",lwd=2) 
lines(seq(47,47+10,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+10,0.01),para_plot)/k+sub_rc[4]+2,col="green3",lwd=2,lty=2) 


ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]+3,gi/k+sub_rc[4]-1,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]+3,ga/k+sub_rc[4]-1,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)
text(td+sub_rc[1]+3,gd/k+sub_rc[4]-1,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")



for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+48,td+sub_rc[1]+0.5,sub_rc[4]+5+48,lwd=1,col="green3",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+48-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)

text(sub_rc[1]-10,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Lateral Root Number "),cex=1.6,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-10,"Time (day)",cex=1.6,family="Times New Roman")

text(sub_rc[1]-4.5,sub_rc[2]+3,expression(B),cex=1.6,family="Times New Roman")


#########################################################################################################


sub_rc<-c(10,70+70,69,10+70)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

pheno<-dat_P_st$lateral_root_ave

k<-20/55

for(i in 0:5){
  segments(sub_rc[1]+0.5+10*i,sub_rc[4],sub_rc[1]+0.5+10*i,sub_rc[4]+1,font=2,lwd=1)
  #text(sub_rc[1]+0.5+10*i,sub_rc[2]+3,30*i, cex=1.6,font=1)
}


for(i in c(0,1,3)){
  segments(sub_rc[1],sub_rc[4]+2+5/k*i,sub_rc[1]+0.6,sub_rc[4]+2+5/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+5/k*i,0.5*i, cex=1.4,font=1)
}

for(i in 2){
  segments(sub_rc[1],sub_rc[4]+2+5/k*i,sub_rc[1]+0.6,sub_rc[4]+2+5/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+5/k*i,expression("1.0"), cex=1.4,font=1)
}

for(i in 4){
  segments(sub_rc[1],sub_rc[4]+2+5/k*i,sub_rc[1]+0.6,sub_rc[4]+2+5/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+5/k*i,expression("2.0"), cex=1.4,font=1)
}


para_plot<- plot_ot[1:3]

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),para_plot)/k+sub_rc[4]+2,col="red",lwd=2) 
lines(seq(47,47+10,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+10,0.01),para_plot)/k+sub_rc[4]+2,col="red",lwd=2,lty=2) 


ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]-2,gi/k+sub_rc[4]+5,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]-2,ga/k+sub_rc[4]+5,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)

text(td+sub_rc[1]-2,gd/k+sub_rc[4]+5,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")


for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+50,td+sub_rc[1]+0.5,sub_rc[4]+5+50,lwd=1,col="red",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+50+3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)



para_plot<- plot_ot[1:3+6]

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),para_plot)/k+sub_rc[4]+2,col="green3",lwd=2) 
lines(seq(47,47+10,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+10,0.01),para_plot)/k+sub_rc[4]+2,col="green3",lwd=2,lty=2) 


ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]+3,gi/k+sub_rc[4]-1,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]+3,ga/k+sub_rc[4]-1,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)
text(td+sub_rc[1]+3,gd/k+sub_rc[4]-1,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")


for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+48,td+sub_rc[1]+0.5,sub_rc[4]+5+48,lwd=1,col="green3",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+48-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)




text(sub_rc[1]-10,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Average Lateral Length " (cm)),cex=1.6,srt=90,family="Times New Roman")

#text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-10,"Time (day)",cex=2,family="Times New Roman")

text(sub_rc[1]-4.5,sub_rc[2]+3,expression(A),cex=1.6,family="Times New Roman")


#########################################################################################


sub_rc<-c(79+25,100,138+25,40)
rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1)

k<-600/53


for(i in 0:5){
  segments(sub_rc[1]+0.5+10*i,sub_rc[4],sub_rc[1]+0.5+10*i,sub_rc[4]+1,font=2,lwd=1)
  text(sub_rc[1]+0.5+10*i,sub_rc[4]-3,10*i, cex=1.4,font=1)
}


for(i in 0:4){
  segments(sub_rc[1],sub_rc[4]+2+150/k*i,sub_rc[1]+0.6,sub_rc[4]+2+150/k*i,font=2,lwd=1)
  text(sub_rc[1]-4,sub_rc[4]+2+150/k*i,15*i, cex=1.4,font=1)
}

para_plot<- plot_len[1:3]

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),para_plot)/k+sub_rc[4]+2,col="red",lwd=2) 
lines(seq(47,47+10,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+10,0.01),para_plot)/k+sub_rc[4]+2,col="red",lwd=2,lty=2) 



ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]-2,gi/k+sub_rc[4]+5,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]-2,ga/k+sub_rc[4]+5,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="red")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)

text(td+sub_rc[1]-2,gd/k+sub_rc[4]+5,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")

for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+50,td+sub_rc[1]+0.5,sub_rc[4]+5+50,lwd=1,col="red",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+50+3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)





para_plot<- plot_len[1:3+3]

lines(seq(1,47,0.01)+sub_rc[1]+0.5,fn_P(seq(1,47,0.01),para_plot)/k+sub_rc[4]+2,col="green3",lwd=2) 
lines(seq(47,47+10,0.01)+sub_rc[1]+0.5,fn_P(seq(47,47+10,0.01),para_plot)/k+sub_rc[4]+2,col="green3",lwd=2,lty=2) 


ti<-log(para_plot[2])/para_plot[3]
gi<-para_plot[1]/2

ta<-log(para_plot[2]*(2-sqrt(3)))/para_plot[3]
ga<-para_plot[1]*(3-sqrt(3))/6

td<-log(para_plot[2]*(2+sqrt(3)))/para_plot[3]
gd<-para_plot[1]*(3+sqrt(3))/6


segments(ti+sub_rc[1]+0.5,sub_rc[4],ti+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(ti+sub_rc[1]+0.5,gi/k+sub_rc[4]+2,pch=17)
text(ti+sub_rc[1]+3,gi/k+sub_rc[4]-1,expression(italic(t["I"])),cex=1.4,family="Times New Roman1")

segments(ta+sub_rc[1]+0.5,sub_rc[4],ta+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(ta+sub_rc[1]+0.5,ga/k+sub_rc[4]+2,pch=17)
text(ta+sub_rc[1]+3,ga/k+sub_rc[4]-1,expression(italic(t["a"])),cex=1.4,family="Times New Roman1")

segments(td+sub_rc[1]+0.5,sub_rc[4],td+sub_rc[1]+0.5,sub_rc[2],lty=3,lwd=1,col="green3")
points(td+sub_rc[1]+0.5,gd/k+sub_rc[4]+2,pch=17)
text(td+sub_rc[1]+3,gd/k+sub_rc[4]-1,expression(italic(t["d"])),cex=1.4,family="Times New Roman1")

for(i in 1:20){
  arrows(ta+sub_rc[1]+0.5,sub_rc[4]+5+48,td+sub_rc[1]+0.5,sub_rc[4]+5+48,lwd=1,col="green3",length=0.08,angle=i,code=3)
}

text(ta+sub_rc[1]+0.5+5,sub_rc[4]+5+48-3,expression(italic(Delta*t)),family="Times New Roman1",cex=1.4)







text(sub_rc[1]-10,(sub_rc[2]-sub_rc[4])/2+sub_rc[4],expression("Total Lateral Root Length " (cm)),cex=1.6,srt=90,family="Times New Roman")

text(sub_rc[1]+(sub_rc[3]-sub_rc[1])/2,sub_rc[4]-10,"Time (day)",cex=1.6,family="Times New Roman")

text(sub_rc[1]-4.5,sub_rc[2]+3,expression(C),cex=1.6,family="Times New Roman")


#########################################################################################



##################################################################################


arrows(69.2,40,91.2,70,lwd=2,col="black",length=0.3,angle=15,code=2,lty=2)

for(i in 1:20){
  arrows(91.2-3*22/30,67,91.2,70,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}


arrows(69.2,110,91.2,70,lwd=2,col="black",length=0.3,angle=15,code=2,lty=2)

for(i in 1:20){
  arrows(91.2-3*22/40,73,91.2,70,lwd=2,col="black",length=0.3,angle=i,code=2,lty=1)
}

dev.off()

