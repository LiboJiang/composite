setwd("/checkpoints/home/lbjiang/smm/p")



dat_P.load <- function(geno_table_P="poplar-marker.txt",
                       pheno_table_P_DIA="poplar-DIA.txt",pheno_table_P_HT="poplar-HT.txt",
                       pheno_table_P_V="poplar-V.txt"){
  dat_P<-list(TIME=NULL,
              geno_table_P=NULL,
              pheno_table_P_HT=NULL,
              pheno_table_P_DIA=NULL,
              pheno_table_P_V=NULL,
              psi=NULL,
              ID=NULL,
              trait=1)
  
  A_P<-read.table("poplar-marker.txt")
  A_P<-t(A_P)
  A_P<-as.matrix(A_P[-1:-4,])
  
  B_P<-read.table("poplar-DIA.txt")
  B_P<-as.matrix(B_P)
  
  D_P<-read.table("poplar-HT.txt")
  D_P<-as.matrix(D_P)
  
  E_P<-read.table("poplar-V.txt")
  E_P<-as.matrix(E_P)
  
  ID_P<-sort(A_P[,1])
  ID_P<-ID_P[-67]
  AA_P<-c()
  for (i in 1:length(ID_P)){
    m<-which(A_P[,1]==ID_P[i])
    AA_P<-rbind(AA_P,A_P[m,])
  }
  rownames(AA_P) <- AA_P[,1]
  A_P.1 <- AA_P[,-1]
  
  B_P.1 <- B_P
  D_P.1 <- D_P
  E_P.1 <- E_P
  
  dat_P$geno_table_P <- A_P.1
  dat_P$pheno_table_P_DIA <- B_P.1
  dat_P$pheno_table_P_HT <- D_P.1*100
  dat_P$pheno_table_P_V <- E_P.1*1000000
  dat_P$psi<- (E_P.1*1000000)/(D_P.1*100)/(B_P.1)^2
  dat_P$ID <- ID_P
  dat_P$TIME <- c(1:24)
  return(dat_P)
}

dat_P_st <- dat_P.load(geno_table_P="poplar-marker.txt",
                       pheno_table_P_DIA="poplar-DIA.txt",pheno_table_P_HT="poplar-HT.txt",
                       pheno_table_P_V="poplar-V.txt")

###########################################################################################################

library(mvtnorm)


AR1.get_mat <- function(par0, times, traits=1, options=list()){
  par<-par0;
  if (class(par0)=="list")
    par <- unlist(par0);
  
  t_len <- length( times );
  
  Ar.1 <- array(0, dim=c(t_len*traits,t_len*traits));
  for (i0 in 1:traits)
    for (i1 in 1:traits)
    {
      if (i0==i1)
        for (k0 in 1:t_len)
          for (k1 in 1:t_len)
          {
            Ar.1[(i0-1)*t_len+k0,(i1-1)*t_len+k1] <- par[i0*2]^2 * par[i0*2-1]^abs( k0 - k1 );
          }
    }
  return(Ar.1);
}

###################################################################################################

fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
}

Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}

###########################################################################################
loss1<-function(par,t,y_P){
  par.covar1_P <- par[21:26]
  AR1_P <- AR1.get_mat(par.covar1_P,t,3)
  y_psi<-Legendre.model(t,par[1:8])
  y_psi<-as.numeric(log(ifelse(y_psi>0,y_psi, 0.1)))
  y_HT<-fn_P(t,par[9:14])
  y_HT<-as.numeric(log( ifelse(y_HT>0,y_HT, 0.1)))
  y_DIA<-fn_P(t,par[15:20])
  y_DIA<-as.numeric(2*log(ifelse(y_DIA>0,y_DIA, 0.1)))
  fy1_P<-dmvnorm(y_P,c(y_psi,y_HT,y_DIA),AR1_P)
  LL_P  <- -sum(log(fy1_P))
  LL_P
} 

loss2<-function(par,t,y_P,snp.type_P,SNP_P){
  snp.index_P <- length(snp.type_P)
  par.covar_P <- par[(snp.index_P*20+1):(snp.index_P*20+6)]
  AR1_P <- AR1.get_mat(par.covar_P,t,traits=3)
  sigma_P<-AR1_P
  fy_P <- 0
  for(ii in 1:length(snp.type_P)){
    curve.mean_P_1 <- Legendre.model(t,par[((ii-1)*20+1):((ii-1)*20+8)])
    curve.mean_P_1 <- log(ifelse(curve.mean_P_1>0,curve.mean_P_1, 0.1))
    curve.mean_P_2 <- fn_P(t,par[((ii-1)*20+9):((ii-1)*20+14)])
    curve.mean_P_2 <- log(ifelse(curve.mean_P_2>0,curve.mean_P_2, 0.1))
    curve.mean_P_3 <- fn_P(t,par[((ii-1)*20+15):((ii-1)*20+20)])
    curve.mean_P_3 <- 2*log(ifelse(curve.mean_P_3>0,curve.mean_P_3, 0.1))
    fy_P <- fy_P + -sum(log(dmvnorm(y_P[which(SNP_P==snp.type_P[ii]),],c(curve.mean_P_1,curve.mean_P_2,curve.mean_P_3),sigma_P)))
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

#######################################################################################################################
times=c(1:24)

par0<-c(0.33892278,  -0.07913492,   0.08122050,  -0.02406407,   0.01988429,  -0.03250651,   0.02573655,  -0.01180772,
        14.51368947,  25.09435099,   0.28605094,
        9.77192923,   6.45274827,   0.70016692,  
        20.96732798,
        12.31035065,   0.58920523,  12.97249578, 193.73751525,   0.28293380,
        0.82529331,   0.18942746,
        0.92861137 ,  0.18531271,   0.97237785,   0.58144749)

H0_P_V<-optim(par0,loss1,t=times,y_P=cbind(log(dat_P_st$psi),log(dat_P_st$pheno_table_P_HT/100),2*log(dat_P_st$pheno_table_P_DIA)),
              method="BFGS",control=list(maxit=20000))#,trace=TRUE

nm <- dim(dat_P_st$geno_table_P)[2]
n1 <- test1
n2 <- test2 
if(n2 >nm)
  n2 <- nm
res <- matrix(NA,nrow=length(c(n1:n2)),ncol=100)
for(i in n1:n2){
  SNP_P <- (dat_P_st$geno_table_P)[,i]
  SNP_P<-as.character(SNP_P)
  snp.type_P <- names(table(SNP_P))
  
  miss.type_P <- grep("\\.",snp.type_P)
  if(length(miss.type_P)>0){
    snp.type_P <- snp.type_P[-miss.type_P]
  }else{
    snp.type_P <- snp.type_P
  }
  
  miss.snp_P <- grep("\\.",(dat_P_st$geno_table_P)[,i])
  if(length(miss.snp_P)>0){
    ny_P_psi <- dat_P_st$psi[-miss.snp_P,]
    ny_P_HT <- (dat_P_st$pheno_table_P_HT/100)[-miss.snp_P,]
    ny_P_DIA <- dat_P_st$pheno_table_P_DIA[-miss.snp_P,]
  }else{
    ny_P_psi <- dat_P_st$psi
    ny_P_HT <- dat_P_st$pheno_table_P_HT/100
    ny_P_DIA <- dat_P_st$pheno_table_P_DIA
  } 
  
  NH0_P_V <- optim(H0_P_V$par,loss1,t=1:24,y_P=cbind(log(ny_P_psi),log(ny_P_HT),2*log(ny_P_DIA)),
                   method="BFGS",
                   control=list(maxit=20000))#,trace = TRUE
  
  
  if(length(snp.type_P)==2)
    npar<-c(NH0_P_V$par[1:20],NH0_P_V$par[1:20],NH0_P_V$par[21:26])
  
  if(length(snp.type_P)==3)
    npar<-c(NH0_P_V$par[1:20],NH0_P_V$par[1:20],NH0_P_V$par[1:20],NH0_P_V$par[21:26])
  
  
  NH1_P_V <- optim(npar,loss2,t=1:24,
                   y_P=cbind(log(dat_P_st$psi),log(dat_P_st$pheno_table_P_HT/100),2*log(dat_P_st$pheno_table_P_DIA)),
                   snp.type_P=snp.type_P,
                   SNP_P=SNP_P,method="BFGS",
                   control=list(maxit=20000))#,trace = TRUE
  
  LR<- 2*(NH0_P_V$value-NH1_P_V$value)
  cat("SNP=",i,"LR=",LR,"\n")
  allpar<-c(LR,NH1_P_V$par)
  res[(i-(n1-1)),(1:length(allpar))] <- allpar
}


ret<-res


filename <- paste("par-",test3,".RData",sep="")
save(ret,file=filename)



