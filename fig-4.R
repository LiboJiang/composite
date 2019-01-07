
draw.point.sig_1<-function(xo,yo,xc, yc,  r2,w_1, col=col,  cex=cex,pch=pch){
  theangle2 <- w_1;
  l2        <- r2;
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  segments(xo,yo,x2,y2,lty=2,lwd=2,col="#EE9A00")
  
  points(xo,yo, col=col,  cex=cex,pch=pch);
  
  
}



draw.point.sig_2<-function(xo,yo,xc, yc,  r2,w_1,w_2, w_3,  col=col,  cex=cex,pch=pch){
  theangle2 <- w_1;
  l2        <- r2;
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  segments(xo,yo,x2,y2,lty=2,lwd=2,col="forestgreen")
  
  theangle3 <- w_2;
  l3        <- r2;
  theangle3 <- (theangle3/360)*2*pi;
  x3        <- xc+l3*cos(theangle3);
  y3        <- yc-l3*sin(theangle3);
  segments(xo,yo,x3,y3,lty=2,lwd=2,col="darkviolet")
  
  points(xo,yo, col=col,  cex=cex,pch=pch);
  
  
}

p2<-function(a,b){
  p1_1<-a*1500+1
  p1_1
}


library(showtext)
showtext_auto(enable=TRUE)
font_add("Times New Roman","times.ttf")
font_add("Times New Roman1",regular = "timesi.ttf")


pdf("fig-4.pdf", 7.2, 8)



xc=400
yc=440
r=340

#col.out=1

chr.po <- ciro_x_1$chr.po
chr.num     <- nrow(chr.po)

plot(c(1,800), c(1,840), type="n", axes=FALSE, xlab="", ylab="", main="")
#redtrans <- rgb(0, 0, 0, 4, maxColorValue=255) 
col.out=rep(c("blue"),20)#,"red","green"

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  draw.arc.s(xc, yc, r=r, w1, w2, col=col.out[chr.i], lwd=8)
  w.m <- (w1+w2)/2
  r.j <- 25
  chr.t <- chr.po[chr.i,1]
  chr.t <- gsub("LG","",chr.t)
  draw.text.rt_1(xc, yc, r+r.j, w.m, chr.t, cex=1.2,LG=chr.i)#
}


for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  if(chr.i<=19){
    draw.line(xc, yc,w2, 330,350, col="black", lwd=1.6, lend=1)
  }
}



chr.po <- ciro_x_2$chr.po
chr.num     <- nrow(chr.po)

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  draw.arc.s(xc, yc, r=r, w1, w2, col=col.out[chr.i], lwd=8)
  w.m <- (w1+w2)/2
  r.j <- 25
  chr.t <- chr.po[chr.i,1]
  chr.t <- gsub("LG","",chr.t)
  draw.text.rt_2(xc, yc, r+r.j, w.m, chr.t, cex=1.2,LG=chr.i)#
}

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  if(chr.i<=19){
    draw.line(xc, yc,w2, 330,350, col="black", lwd=1.6, lend=1)
  }
}


chr.po <- ciro_x_3$chr.po
chr.num     <- nrow(chr.po)

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  draw.arc.s(xc, yc, r=r, w1, w2, col=col.out[chr.i], lwd=8)
  w.m <- (w1+w2)/2
  r.j <- 25
  chr.t <- chr.po[chr.i,1]
  chr.t <- gsub("LG","",chr.t)
  draw.text.rt_3(xc, yc, r+r.j, w.m, chr.t, cex=1.2,LG=chr.i)#
  
}

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  if(chr.i<=19){
    draw.line(xc, yc,w2, 330,350, col="black", lwd=1.6, lend=1)
  }
}
text(340,48,expression("S " ),cex=1.6,srt=-6,adj=0,family="Times New Roman")
text(340+20,44,expression("t " ),cex=1.6,srt=-4,adj=0,family="Times New Roman")
text(340+15+15,40,expression("e " ),cex=1.6,srt=-2,adj=0,family="Times New Roman")
text(340+15+15+15,38,expression("m " ),cex=1.6,srt=0,adj=0,family="Times New Roman")


text(423,42,expression("F" ),cex=1.6,srt=4,adj=0,family="Times New Roman")
text(440,36,expression("o" ),cex=1.6,srt=6,adj=0,family="Times New Roman")
text(440+17,38,expression("r" ),cex=1.6,srt=8,adj=0,family="Times New Roman")
text(440+17+12,40,expression("m" ),cex=1.6,srt=12,adj=0,family="Times New Roman")


text(495,50,expression(" F" ),cex=1.6,srt=16,adj=0,family="Times New Roman")
text(520,52,expression("a" ),cex=1.6,srt=18,adj=0,family="Times New Roman")
text(520+15,58,expression("c" ),cex=1.6,srt=20,adj=0,family="Times New Roman")
text(520+15+15,66,expression("t" ),cex=1.6,srt=22,adj=0,family="Times New Roman")
text(520+15+15+10,70,expression("o" ),cex=1.6,srt=24,adj=0,family="Times New Roman")
text(520+15+15+10+15,78,expression("r" ),cex=1.6,srt=26,adj=0,family="Times New Roman")


text(15+4,540,expression("S"),cex=1.6,srt=74,adj=0,family="Times New Roman")
text(15+10,540+20,expression("t"),cex=1.6,srt=72,adj=0,family="Times New Roman")
text(15+10+5,540+30,expression("e"),cex=1.6,srt=70,adj=0,family="Times New Roman")
text(15+10+5+5,540+46,expression("m"),cex=1.6,srt=68,adj=0,family="Times New Roman")
#text(15,540,expression("Stem"),cex=1.6,srt=65,adj=0,family="Times New Roman")
#text(50,640,expression("Height"),cex=1.6,srt=54,adj=0,family="Times New Roman")
#text(100,740,expression("(m)"),cex=1.6,srt=40,adj=0,family="Times New Roman")
text(45,640,expression("H"),cex=1.6,srt=60,adj=0,family="Times New Roman")
text(45+15,640+20,expression("e"),cex=1.6,srt=58,adj=0,family="Times New Roman")
text(45+10+10,640+20+15,expression("i"),cex=1.6,srt=56,adj=0,family="Times New Roman")
text(45+10+10+10,640+15+15+8,expression("g "),cex=1.6,srt=54,adj=0,family="Times New Roman")
text(45+10+10+10+5,640+15+15+15+15,expression("h "),cex=1.6,srt=52,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
#text(15,540,expression("( "),cex=1.6,srt=33,adj=0,family="Times New Roman")
#text(15,540,expression("m "),cex=1.6,srt=30,adj=0,family="Times New Roman")
#text(15,540,expression(")"),cex=1.6,srt=27,adj=0,family="Times New Roman")

text(720,680,expression("DBH " (cm)),cex=1.6,srt=130+180,family="Times New Roman")



draw.point.sig_2(200,640,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(68),240), 
                 col="red",  cex=2.4,pch=16)

draw.point.sig_1(200,640,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)

draw.point.sig_2(600,640,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(37),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(16),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(39),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(44),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(45),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(48),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(52),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(64),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(67),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(68),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(76),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(90),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(94),120),
                 col="red",  cex=2.4,pch=16)


draw.point.sig_2(400,640,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(34),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(400,640,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(400,640,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)



draw.point.sig_2(300,440,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(23),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(45),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(48),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(52),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(90),120),
                 col="red",  cex=2.4,pch=16)


draw.point.sig_2(200,380,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(8),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(200,380,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(200,380,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)



dev.off()









bmp("fig-4.bmp", width=720, height =800, units = "px", pointsize = 12)



xc=400
yc=440
r=340

#col.out=1

chr.po <- ciro_x_1$chr.po
chr.num     <- nrow(chr.po)

plot(c(1,800), c(1,840), type="n", axes=FALSE, xlab="", ylab="", main="")
#redtrans <- rgb(0, 0, 0, 4, maxColorValue=255) 
col.out=rep(c("blue"),20)#,"red","green"

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  draw.arc.s(xc, yc, r=r, w1, w2, col=col.out[chr.i], lwd=8)
  w.m <- (w1+w2)/2
  r.j <- 25
  chr.t <- chr.po[chr.i,1]
  chr.t <- gsub("LG","",chr.t)
  draw.text.rt_1(xc, yc, r+r.j, w.m, chr.t, cex=1.2,LG=chr.i)#
}


for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  if(chr.i<=19){
    draw.line(xc, yc,w2, 330,350, col="black", lwd=1.6, lend=1)
  }
}



chr.po <- ciro_x_2$chr.po
chr.num     <- nrow(chr.po)

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  draw.arc.s(xc, yc, r=r, w1, w2, col=col.out[chr.i], lwd=8)
  w.m <- (w1+w2)/2
  r.j <- 25
  chr.t <- chr.po[chr.i,1]
  chr.t <- gsub("LG","",chr.t)
  draw.text.rt_2(xc, yc, r+r.j, w.m, chr.t, cex=1.2,LG=chr.i)#
}

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  if(chr.i<=19){
    draw.line(xc, yc,w2, 330,350, col="black", lwd=1.6, lend=1)
  }
}


chr.po <- ciro_x_3$chr.po
chr.num     <- nrow(chr.po)

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  draw.arc.s(xc, yc, r=r, w1, w2, col=col.out[chr.i], lwd=8)
  w.m <- (w1+w2)/2
  r.j <- 25
  chr.t <- chr.po[chr.i,1]
  chr.t <- gsub("LG","",chr.t)
  draw.text.rt_3(xc, yc, r+r.j, w.m, chr.t, cex=1.2,LG=chr.i)#
  
}

for (chr.i in 1:chr.num){
  w1 <- as.numeric(chr.po[chr.i,2])
  w2 <- as.numeric(chr.po[chr.i,3])
  if(chr.i<=19){
    draw.line(xc, yc,w2, 330,350, col="black", lwd=1.6, lend=1)
  }
}
text(340,48,expression("S " ),cex=1.6,srt=-6,adj=0,family="Times New Roman")
text(340+20,44,expression("t " ),cex=1.6,srt=-4,adj=0,family="Times New Roman")
text(340+15+15,40,expression("e " ),cex=1.6,srt=-2,adj=0,family="Times New Roman")
text(340+15+15+15,38,expression("m " ),cex=1.6,srt=0,adj=0,family="Times New Roman")


text(423,42,expression("F" ),cex=1.6,srt=4,adj=0,family="Times New Roman")
text(440,36,expression("o" ),cex=1.6,srt=6,adj=0,family="Times New Roman")
text(440+17,38,expression("r" ),cex=1.6,srt=8,adj=0,family="Times New Roman")
text(440+17+12,40,expression("m" ),cex=1.6,srt=12,adj=0,family="Times New Roman")


text(495,50,expression(" F" ),cex=1.6,srt=16,adj=0,family="Times New Roman")
text(520,52,expression("a" ),cex=1.6,srt=18,adj=0,family="Times New Roman")
text(520+15,58,expression("c" ),cex=1.6,srt=20,adj=0,family="Times New Roman")
text(520+15+15,66,expression("t" ),cex=1.6,srt=22,adj=0,family="Times New Roman")
text(520+15+15+10,70,expression("o" ),cex=1.6,srt=24,adj=0,family="Times New Roman")
text(520+15+15+10+15,78,expression("r" ),cex=1.6,srt=26,adj=0,family="Times New Roman")


text(15+4,540,expression("S"),cex=1.6,srt=74,adj=0,family="Times New Roman")
text(15+10,540+20,expression("t"),cex=1.6,srt=72,adj=0,family="Times New Roman")
text(15+10+5,540+30,expression("e"),cex=1.6,srt=70,adj=0,family="Times New Roman")
text(15+10+5+5,540+46,expression("m"),cex=1.6,srt=68,adj=0,family="Times New Roman")
#text(15,540,expression("Stem"),cex=1.6,srt=65,adj=0,family="Times New Roman")
#text(50,640,expression("Height"),cex=1.6,srt=54,adj=0,family="Times New Roman")
#text(100,740,expression("(m)"),cex=1.6,srt=40,adj=0,family="Times New Roman")
text(45,640,expression("H"),cex=1.6,srt=60,adj=0,family="Times New Roman")
text(45+15,640+20,expression("e"),cex=1.6,srt=58,adj=0,family="Times New Roman")
text(45+10+10,640+20+15,expression("i"),cex=1.6,srt=56,adj=0,family="Times New Roman")
text(45+10+10+10,640+15+15+8,expression("g "),cex=1.6,srt=54,adj=0,family="Times New Roman")
text(45+10+10+10+5,640+15+15+15+15,expression("h "),cex=1.6,srt=52,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
text(45+10+10+10+10+8,640+15+15+15+15+15,expression("t "),cex=1.6,srt=50,adj=0,family="Times New Roman")
#text(15,540,expression("( "),cex=1.6,srt=33,adj=0,family="Times New Roman")
#text(15,540,expression("m "),cex=1.6,srt=30,adj=0,family="Times New Roman")
#text(15,540,expression(")"),cex=1.6,srt=27,adj=0,family="Times New Roman")

text(720,680,expression("DBH " (cm)),cex=1.6,srt=130+180,family="Times New Roman")



draw.point.sig_2(200,640,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(68),240), 
                 col="red",  cex=2.4,pch=16)

draw.point.sig_1(200,640,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)

draw.point.sig_2(600,640,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(37),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(16),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(39),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(44),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(45),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(48),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(52),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(64),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(67),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(68),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(76),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(90),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(600,640,xc, yc,  r-5,point_angle_D(p2(94),120),
                 col="red",  cex=2.4,pch=16)


draw.point.sig_2(400,640,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(34),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(400,640,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(400,640,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)



draw.point.sig_2(300,440,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(23),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(45),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(48),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(52),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(300,440,xc, yc,  r-5,point_angle_D(p2(90),120),
                 col="red",  cex=2.4,pch=16)


draw.point.sig_2(200,380,xc, yc,  r-5,point_angle_D(p2(3),0),point_angle_D(p2(8),240), 
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(200,380,xc, yc,  r-5,point_angle_D(p2(60),120),
                 col="red",  cex=2.4,pch=16)
draw.point.sig_1(200,380,xc, yc,  r-5,point_angle_D(p2(72),120),
                 col="red",  cex=2.4,pch=16)



dev.off()





