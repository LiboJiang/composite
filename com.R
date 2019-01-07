ret.1 <- c()
for(i in 1:76){
  A <- paste("par-",i,".RData",sep="")
  ret.tmp <- try(load(A))
  if((class(ret.tmp) == "try-error") )
    ret.1 <- rbind(ret.1,matrix(NA,2060,100))
  else
    ret.1 <- rbind(ret.1,ret)
}

final_par<-ret.1
save(final_par,file="final-par.RData")