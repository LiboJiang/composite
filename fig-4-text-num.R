
draw.text.rt_1 <- function(xc, yc, r, w, n, col="black", cex=1.2, side="out",LG){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 18;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w ;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  if(LG==1)
    text(x-15, y-7, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==5)
    text(x-13, y, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==10)
    text(x-15, y+12, adj=0, offset=1, labels=n, srt=the.w+270-3, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==15)
    text(x-10, y+12, adj=0, offset=1, labels=n, srt=the.w+270-5, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==19)
    text(x-8, y+12, adj=0, offset=1, labels=n, srt=the.w+90, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  
}



draw.text.rt_2 <- function(xc, yc, r, w, n, col="black", cex=1.2, side="out",LG){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 18;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w ;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  if(LG==1)
    text(x-8, y-7, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==5)
    text(x-13, y-4, adj=0, offset=1, labels=n, srt=the.w+270+2, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==10)
    text(x+18, y-2, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==15)
    text(x+15, y-9, adj=0, offset=1, labels=n, srt=the.w+270+2, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==19)
    text(x+12, y-10, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  
}




draw.text.rt_3 <- function(xc, yc, r, w, n, col="black", cex=1.2, side="out",LG){
  
  w     <- w%%360;
  the.o <- w;
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 18;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w ;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.w <- the.w +180;
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  if(LG==1)
    text(x+15, y-15, adj=0, offset=1, labels=n, srt=the.w+270+5, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==5)
    text(x+10, y+10, adj=0, offset=1, labels=n, srt=the.w+270+5, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==10)
    text(x+18, y+5, adj=0, offset=1, labels=n, srt=the.w+270-5, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==15)
    text(x+18, y+2, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  if(LG==19)
    text(x+18, y-5, adj=0, offset=1, labels=n, srt=the.w+270, 
         pos=the.pos, col=col, cex=cex,family="Times New Roman");
  
}



