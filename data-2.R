setwd("/home/workdir/mmsang/smm/cegen")


dat_P.load <- function(geno_table="Map-Genotype.csv",
                       pheno1_table="final 10.csv",
                       pheno2_table="final 11.csv",
                       pheno3_table="final 12.csv"){
  dat_P<-list(TIME=NULL,
              geno_table=NULL,
              lateral_root_len=NULL,
              lateral_root_num=NULL,
              lateral_root_ave=NULL,
              ID=NULL,
              trait=1)
  
  A_P<-read.csv("Map-Genotype.csv",header=F)
  A_P <- A_P[,-1:-3]
  A_P <- t(A_P)
  
  B_P<-read.csv("final 10.csv")
  B_P<-as.matrix(B_P)
  
  D_P<-read.csv("final 11.csv")
  D_P<-as.matrix(D_P)
  
  E_P<-read.csv("final 12.csv")
  E_P<-as.matrix(E_P)
  
  A1_P<- A_P[,1]
  B1_P<- B_P[,1]
  D1_P<- D_P[,1]
  E1_P<- E_P[,1]
  AB_P<- table(c(A1_P,B1_P,D1_P,E1_P))
  raw.id <- as.numeric(names(AB_P))
  C_P<-which(AB_P==4)
  ID_P<- sort(raw.id[C_P])
  
  
  AA_P<- c()
  BB_P<- c()
  DD_P<- c()
  EE_P<- c()
  for (i in 1:length(ID_P)){
    m<-which(A_P[,1]==ID_P[i])
    AA_P<-rbind(AA_P,A_P[m,])
    n1<-which(B_P[,1]==ID_P[i])
    BB_P<-rbind(BB_P,B_P[n1,])
    n2<-which(D_P[,1]==ID_P[i])
    DD_P<-rbind(DD_P,D_P[n2,])
    n3<-which(E_P[,1]==ID_P[i])
    EE_P<-rbind(EE_P,E_P[n3,])
  }
  
  rownames(AA_P) <- AA_P[,1]
  A_P.1 <- AA_P[,-1]
  rownames(BB_P) <- BB_P[,1]
  B_P.1 <- BB_P[,-1:-5]
  rownames(DD_P) <- DD_P[,1]
  D_P.1 <- DD_P[,-1]
  rownames(EE_P) <- EE_P[,1]
  E_P.1 <- EE_P[,-1]
  
  
  E_P2 <- EE_P[,-1]
  id2<-c()
  for(i in 1:dim( E_P2 )[1]){
    id1<-length(which(E_P2[i,]==0))
    if(id1>9){
      id2<-c(id2,i)
    }
  }
  
  dat_P$geno_table <- A_P.1[-id2,]
  dat_P$lateral_root_len <- B_P.1[-id2,]
  dat_P$lateral_root_num <- D_P.1[-id2,]
  dat_P$lateral_root_ave <- E_P.1[-id2,]
  
  dat_P$ID <- ID_P[-id2]
  dat_P$TIME <- c(1,3,5,7,9,11,13,16,19,23,28,32,39,47)
  return(dat_P)
}

dat_P_st <- dat_P.load(geno_table="Map-Genotype.csv",
                       pheno1_table="final 10.csv",
                       pheno2_table="final 11.csv",
                       pheno3_table="final 12.csv")


###################################################################


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
  par.covar1 <- par[7:10]
  AR1_P <- SAD1(par.covar1,t,traits = 2)
  fy1_P<-c()
  for(i in 1:dim(y_P)[1]){
    y1<-y_P[i,]
    t1<-t
    t2<-t1[-which(y1==0)]
    y2<-y1[-which(y1==0)]
    
    AR2_P<-AR1_P[-which(y1==0),-which(y1==0)]
    
    t2_1<-t1[1:length(t1)][-which(y1[1:length(t1)]==0)]
    t2_2<-t1[1:length(t1)][-which(y1[1:length(t1)+length(t1)]==0)]
    
    y3_1<-fn_P(t2_1,par[1:3])
    y3_2<-fn_P(t2_2,par[4:6])
    yy<-dmvnorm(y2,c(y3_1,y3_2),AR2_P)
    fy1_P<-c(fy1_P,yy)
  }
  
  LL_P  <- -sum(log(fy1_P))
  LL_P
} 






loss2<-function(par,t,y_P,snp.type,SNP){
  snp.index <- length(snp.type)
  par.covar <- par[(snp.index*6+1):(snp.index*6+4)]
  AR1_P <- SAD1(par.covar,t,traits=2)
  
  fy_P <- 0
  for(ii in 1:snp.index){
    
    curve.mean_P_1<-c()
    yy<-y_P[which(SNP==snp.type[ii]),]
    for(j in 1:dim(yy)[1]){
      y1<-yy[j,]
      t1<-t
      #t2<-t1[-which(y1==0)]
      y2<-y1[-which(y1==0)]
      
      AR2_P<-AR1_P[-which(y1==0),-which(y1==0)]
      
      t2_1<-t1[1:length(t1)][-which(y1[1:length(t1)]==0)]
      t2_2<-t1[1:length(t1)][-which(y1[1:length(t1)+length(t1)]==0)]
      
      y3_1<-fn_P(t2_1,par[((ii-1)*6+1):((ii-1)*6+3)])
      y3_2<-fn_P(t2_2,par[((ii-1)*6+4):((ii-1)*6+6)])
      fy<-dmvnorm(y2,c(y3_1,y3_2),AR2_P)
      
      curve.mean_P_1<-c(curve.mean_P_1,fy)
    }
    fy_P <- fy_P + -sum(log(curve.mean_P_1))
    
  }
  fy_P
}


########################################################################################################

get_con_param<-function(parm.id)
{
  for (e in commandArgs())
  {
    ta = strsplit(e,"=", fixed=TRUE);
    if(! is.na( ta[[1]][2]))
    {
      temp = ta[[1]][2];
      if( ta[[1]][1] == parm.id) {
        return (as.character(temp));
      }
    }
  }
  
  return(NA);
}


test1 <- get_con_param("test1")
test2 <- get_con_param("test2")
test3 <- get_con_param("test3")

test1 <- as.numeric(test1)
test2 <- as.numeric(test2)
test3 <- as.numeric(test3)

#########################################


pheno1<-dat_P_st$lateral_root_ave
pheno2<-dat_P_st$lateral_root_num
geno<-dat_P_st$geno_table

par0<-c( 15.1088471, 20.9177517,  0.1243092,  41.55090420, 117.02996102,   0.11183293,
         0.6966802,  3.9718366,   0.93537430,
         3.68245873)
H0<-optim(par0,loss1,t=dat_P_st$TIME,y_P=cbind(pheno1,pheno2),
          method="BFGS",control=list(maxit=10000))


nm <- dim(geno)[2]
n1 <- test1
n2 <- test2 
if(n2 >nm)
  n2 <- nm
res <- matrix(NA,nrow=length(c(n1:n2)),ncol=100)


for(i in n1:n2){#
  SNP <- (geno)[,i]
  SNP<-as.character(SNP)
  snp.type <- names(table(SNP))
  
  miss.type <- grep("-",snp.type)
  if(length(miss.type)>0){
    snp.type <- snp.type[-miss.type]
  }else{
    snp.type <- snp.type
  }
  
  miss.snp <- grep("-",(geno)[,i])
  if(length(miss.snp)>0){
    ny_P_1 <- pheno1[-miss.snp,]
    ny_P_2 <- pheno2[-miss.snp,]
  }else{
    ny_P_1 <- pheno1
    ny_P_2 <- pheno2
  } 
  
  NH0 <- optim(H0$par,loss1,t=dat_P_st$TIME,y_P=cbind(ny_P_1,ny_P_2),
               method="BFGS",
               control=list(maxit=20000))#,trace = TRUE
  
  
  
  if(length(snp.type)==2)
    npar<-c(NH0$par[1:6],NH0$par[1:6],NH0$par[7:10])
  
  if(length(snp.type)==3)
    npar<-c(NH0$par[1:6],NH0$par[1:6],NH0$par[1:6],NH0$par[7:10])
  
  NH1 <- optim(npar,loss2,t=dat_P_st$TIME,
               y_P=cbind(pheno1,pheno2),
               snp.type=snp.type,
               SNP=SNP,method="BFGS",
               control=list(maxit=20000))#,trace = TRUE
  
  
  LR<- 2*(NH0$value-NH1$value)
  cat("SNP=",i,"LR=",LR,"\n")
  allpar<-c(LR,NH1$par)
  res[(i-(n1-1)),(1:length(allpar))] <- allpar
}

ret_del<-res


filename <- paste("par_del-",test3,".RData",sep="")
save(ret_del,file=filename)


