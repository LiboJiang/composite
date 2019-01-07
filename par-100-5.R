setwd("D:/Rdocument/p-simulation")

par_F1<-c(0.343763043,  -0.107975321,   0.081683157,  -0.003651312,   0.007474049,  -0.031597981,   0.022869350,  -0.004370809)
par_F2<-c(0.336629708,  -0.080396811,   0.093624247,  -0.039648608,   0.034869980,  -0.036736270,   0.016930957,  -0.004751406)
par_F3<-c(0.339233974,  -0.058399805,   0.061893012,  -0.010961303,   0.002420895,  -0.027216529,   0.046103764,  -0.031198423)


par_HT1<-c(15.2541915, 32.5000220,  0.2771275, 10.1702390,  4.5673914,  0.5111841)   
par_HT2<-c(13.3489054, 33.1604386,  0.2845272, 11.4044683,  6.8261126,  0.6386955)   
par_HT3<-c(18.2498348, 11.0692102,  0.3207969,  4.5062786, 20.2115823,  1.6530021)  


par_DIA1<-c(23.195245268,   7.681790270,  0.441583672,  12.545997121, 195.624611397,   0.274148011) 
par_DIA2<-c(20.866973757,  12.471721801,  0.615003879,  12.403325627, 200.881300198,   0.287800576)             
par_DIA3<-c(21.067222182,  15.476485424,  0.612753796,  14.301687476, 203.580339517,   0.275731109)



timepoint <- 1:24


############--------------marker data simulation---------------############
pro<-function(n,mean,sd){
  A1<-rnorm(n,mean,sd) 
  num_1<-which(A1<0.05 | A1>0.95)
  
  for(i in num_1){
    A1[i]<-runif(1,0.05,0.5)
  }
  
  A2<-rnorm(n,mean,sd)
  num_2<-which(A2<0.05 | A2>0.95)
  
  for(i in num_2){
    A2[i]<-runif(1,0.05,0.5)
  }
  
  A3<-1-A1-A2
  
  num_3<-which(A3<0.05 | A3>0.95)
  for(i in num_3){
    A1[i]<-runif(1,0.05,0.45)
    A2[i]<-runif(1,0.05,0.45)
  }
  
  A3<-1-A1-A2
  
  A11<-as.matrix(A1,ncol=1)
  A21<-as.matrix(A2,ncol=1)
  A31<-as.matrix(A3,ncol=1)
  
  pro1<-cbind(A11,A21,A31)
  pro1
}


#marker<-function(num,n,mean,sd){
#marker1<- matrix(0,num,n)
#pro_1<-pro(n,mean,sd)
#for(i in 1: dim(pro_1)[1]){
# marker1[,i]<-sample(c(0,1,2), num, replace = TRUE,prob=c(1/4,1/2,1/4))

#}
# marker1
#}

marker<-function(num){
  marker1<- matrix(0,num,1)
  for(i in 1: dim( marker1)[2]){
    marker1[,i]<-sample(c(0,1,2), num, replace = TRUE,prob=c(1/4,1/2,1/4))
    
  }
  marker1
}


#########----------------growth curve----------------#########

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


fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
}


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


##########------------phenotype simulation--------------########## 

rho1 <- 0.821383805;s1<-0.188520747*3*0.95
rho2 <- 0.930985766;s2<-0.181897417*0.95 
rho3 <- 0.969354160;s3<-0.517247460*0.95 


library(mvtnorm)
log_F1<-log(Legendre.model(timepoint,par_F1))
log_F2<-log(Legendre.model(timepoint,par_F2))
log_F3<-log(Legendre.model(timepoint,par_F3))

log_HT1<-log(fn_P(timepoint,par_HT1))
log_HT2<-log(fn_P(timepoint,par_HT2))
log_HT3<-log(fn_P(timepoint,par_HT3))

log_DIA1<-log(fn_P(timepoint,par_DIA1))
log_DIA2<-log(fn_P(timepoint,par_DIA2))
log_DIA3<-log(fn_P(timepoint,par_DIA3))

var_cov1 <- AR1.get_mat(c(rho1,s1),timepoint,traits=1)
var_cov2 <- AR1.get_mat(c(rho2,s2),timepoint,traits=1)
var_cov3 <- AR1.get_mat(c(rho3,s3),timepoint,traits=1)



log_pheno_F<-function(marker1,n){
  position_qtl <- n
  
  marker_qtl<-marker1[,position_qtl]
  
  log_phenotype_F<-matrix(0,length(marker_qtl),timepoint[24])
  for(i in 1:length(marker_qtl)){
    if(marker_qtl[i]==0){
      log_phenotype_F[i,]<-rmvnorm(1,log_F1,AR1.get_mat(c(rho1,s1),timepoint,traits=1))
      while(length(which(log_phenotype_F[i,]>=0)>0)){
        log_phenotype_F[i,]<-rmvnorm(1,log_F1,AR1.get_mat(c(rho1,s1),timepoint,traits=1))
      }
    }
    
    if(marker_qtl[i]==1){
      log_phenotype_F[i,]<-rmvnorm(1,log_F2,AR1.get_mat(c(rho1,s1),timepoint,traits=1))
      while(length(which(log_phenotype_F[i,]>=0)>0)){
        log_phenotype_F[i,]<-rmvnorm(1,log_F2,AR1.get_mat(c(rho1,s1),timepoint,traits=1))
      }
    }
    
    if(marker_qtl[i]==2){
      log_phenotype_F[i,]<-rmvnorm(1,log_F3,AR1.get_mat(c(rho1,s1),timepoint,traits=1))
      while(length(which(log_phenotype_F[i,]>=0)>0)){
        log_phenotype_F[i,]<-rmvnorm(1,log_F3,AR1.get_mat(c(rho1,s1),timepoint,traits=1))
      }
    }
  }
  log_phenotype_F
}



#########----------------simulation end----------------#########	

loss1_P<-function(par,t,y_P){
  par.covar1_P <- par[9:10]
  AR1_P <- AR1.get_mat(par.covar1_P,t,1)
  y_psi<-Legendre.model(t,par[1:8])
  y_psi<-as.numeric(log(ifelse(y_psi>0,y_psi, 0.1)))
  fy1_P<-dmvnorm(y_P,y_psi,AR1_P)
  LL_P  <- -sum(log(fy1_P))
  LL_P
} 

loss2_P<-function(par,t,y_P,snp.type_P,SNP_P){
  snp.index_P <- length(snp.type_P)
  par.covar_P <- par[(snp.index_P*8+1):(snp.index_P*8+2)]
  AR1_P <- AR1.get_mat(par.covar_P,t,traits=1)
  sigma_P<-AR1_P
  fy_P <- 0
  for(ii in 1:length(snp.type_P)){
    curve.mean_P_1 <- Legendre.model(t,par[((ii-1)*8+1):((ii-1)*8+8)])
    curve.mean_P_1 <- log(ifelse(curve.mean_P_1>0,curve.mean_P_1, 0.1))
    fy_P <- fy_P + -sum(log(dmvnorm(y_P[which(SNP_P==snp.type_P[ii]),],curve.mean_P_1,sigma_P)))
  }
  fy_P
}

#############################################################################################3

par_F<-colMeans(rbind(par_F1,par_F2,par_F3))

par_F_100_5<-matrix(0,nrow=100,ncol=26)
L_F_100_5<-c()
for(i in 1:100){
  
  par0<-c(par_F,rho1,s1)
  marker1<-marker(100)
  log_F<-log_pheno_F(marker1,1)
  H0<-optim(par0,loss1_P,t=1:24,y_P=log_F,
            method="BFGS",
            control=list(maxit=20000))
  
  SNP_P <- marker1
  SNP_P<-as.character(SNP_P)
  snp.type_P <- names(table(SNP_P))
  
  miss.type_P <- grep("\\.",snp.type_P)
  if(length(miss.type_P)>0){
    snp.type_P <- snp.type_P[-miss.type_P]
  }else{
    snp.type_P <- snp.type_P
  }
  
  
  npar<-c(H0$par[1:8],H0$par[1:8],H0$par[1:8],H0$par[9:10])
  
  H1 <- optim(npar,loss2_P,t=1:24,
              y_P=log_F,
              snp.type_P=snp.type_P,
              SNP_P=SNP_P,method="BFGS",
              control=list(maxit=20000))
  LR<- 2*(H0$value-H1$value)
  
  while(LR<1){
    par0<-c(par_F,rho1,s1)
    marker1<-marker(100)
    log_F<-log_pheno_F(marker1,1)
    H0<-optim(par0,loss1_P,t=1:24,y_P=log_F,
              method="BFGS",
              control=list(maxit=20000))
    
    SNP_P <- marker1
    SNP_P<-as.character(SNP_P)
    snp.type_P <- names(table(SNP_P))
    
    miss.type_P <- grep("\\.",snp.type_P)
    if(length(miss.type_P)>0){
      snp.type_P <- snp.type_P[-miss.type_P]
    }else{
      snp.type_P <- snp.type_P
    }
    
    
    npar<-c(H0$par[1:8],H0$par[1:8],H0$par[1:8],H0$par[9:10])
    
    H1 <- optim(npar,loss2_P,t=1:24,
                y_P=log_F,
                snp.type_P=snp.type_P,
                SNP_P=SNP_P,method="BFGS",
                control=list(maxit=20000))
    LR<- 2*(H0$value-H1$value)
  }
  
  L_F_100_5<-c(L_F_100_5,LR)
  cat("SNP=",i,"LR=",LR,"\n")
  par_F_100_5[i,(1:26)]<-H1$par
  
}

save( L_F_100_5,file="L_F_100_5.RData")
save( par_F_100_5,file="par_F_100_5.RData")


par_F1<-c(0.343763043,  -0.107975321,   0.081683157,  -0.003651312,   0.007474049,  -0.031597981,   0.022869350,  -0.004370809)
par_F2<-c(0.336629708,  -0.080396811,   0.093624247,  -0.039648608,   0.034869980,  -0.036736270,   0.016930957,  -0.004751406)
par_F3<-c(0.339233974,  -0.058399805,   0.061893012,  -0.010961303,   0.002420895,  -0.027216529,   0.046103764,  -0.031198423)


par_HT1<-c(15.2541915, 32.5000220,  0.2771275, 10.1702390,  4.5673914,  0.5111841)   
par_HT2<-c(13.3489054, 33.1604386,  0.2845272, 11.4044683,  6.8261126,  0.6386955)   
par_HT3<-c(18.2498348, 11.0692102,  0.3207969,  4.5062786, 20.2115823,  1.6530021) 

par_DIA1<-c(23.195245268,   7.681790270,  0.441583672,  12.545997121, 195.624611397,   0.274148011) 
par_DIA2<-c(20.866973757,  12.471721801,  0.615003879,  12.403325627, 200.881300198,   0.287800576)             
par_DIA3<-c(21.067222182,  15.476485424,  0.612753796,  14.301687476, 203.580339517,   0.275731109)



timepoint <- 1:24


############--------------marker data simulation---------------############
pro<-function(n,mean,sd){
  A1<-rnorm(n,mean,sd) 
  num_1<-which(A1<0.05 | A1>0.95)
  
  for(i in num_1){
    A1[i]<-runif(1,0.05,0.5)
  }
  
  A2<-rnorm(n,mean,sd)
  num_2<-which(A2<0.05 | A2>0.95)
  
  for(i in num_2){
    A2[i]<-runif(1,0.05,0.5)
  }
  
  A3<-1-A1-A2
  
  num_3<-which(A3<0.05 | A3>0.95)
  for(i in num_3){
    A1[i]<-runif(1,0.05,0.45)
    A2[i]<-runif(1,0.05,0.45)
  }
  
  A3<-1-A1-A2
  
  A11<-as.matrix(A1,ncol=1)
  A21<-as.matrix(A2,ncol=1)
  A31<-as.matrix(A3,ncol=1)
  
  pro1<-cbind(A11,A21,A31)
  pro1
}


#marker<-function(num,n,mean,sd){
#marker1<- matrix(0,num,n)
#pro_1<-pro(n,mean,sd)
#for(i in 1: dim(pro_1)[1]){
# marker1[,i]<-sample(c(0,1,2), num, replace = TRUE,prob=c(1/4,1/2,1/4))

#}
# marker1
#}

marker<-function(num){
  marker1<- matrix(0,num,1)
  for(i in 1: dim( marker1)[2]){
    marker1[,i]<-sample(c(0,1,2), num, replace = TRUE,prob=c(1/4,1/2,1/4))
    
  }
  marker1
}
#########----------------growth curve----------------#########

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


fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
}


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


##########------------phenotype simulation--------------########## 

rho1 <- 0.821383805;s1<-0.188520747*5*0.95
rho2 <- 0.930985766;s2<-0.181897417*0.95 
rho3 <- 0.969354160;s3<-0.517247460*0.95 


library(mvtnorm)
log_F1<-log(Legendre.model(timepoint,par_F1))
log_F2<-log(Legendre.model(timepoint,par_F2))
log_F3<-log(Legendre.model(timepoint,par_F3))

log_HT1<-log(fn_P(timepoint,par_HT1))
log_HT2<-log(fn_P(timepoint,par_HT2))
log_HT3<-log(fn_P(timepoint,par_HT3))

log_DIA1<-log(fn_P(timepoint,par_DIA1))
log_DIA2<-log(fn_P(timepoint,par_DIA2))
log_DIA3<-log(fn_P(timepoint,par_DIA3))

var_cov1 <- AR1.get_mat(c(rho1,s1),timepoint,traits=1)
var_cov2 <- AR1.get_mat(c(rho2,s2),timepoint,traits=1)
var_cov3 <- AR1.get_mat(c(rho3,s3),timepoint,traits=1)


log_pheno<-function(marker1,n){
  position_qtl <- n
  
  marker_qtl<-marker1[,position_qtl]
  
  log_phenotype<-matrix(0,length(marker_qtl),timepoint[24]*3)
  for(i in 1:length(marker_qtl)){
    if(marker_qtl[i]==0)
      log_phenotype[i,]<-rmvnorm(1,c(log_F1,log_HT1,2*log_DIA1),var_cov)
    if(marker_qtl[i]==1)
      log_phenotype[i,]<-rmvnorm(1,c(log_F2,log_HT2,2*log_DIA2),var_cov)
    if(marker_qtl[i]==2)
      log_phenotype[i,]<-rmvnorm(1,c(log_F3,log_HT3,2*log_DIA3),var_cov)
  }
  log_phenotype
}

log_pheno_HT<-function(marker1,n){
  position_qtl <- n
  
  marker_qtl<-marker1[,position_qtl]
  
  log_phenotype_HT<-matrix(0,length(marker_qtl),timepoint[24])
  for(i in 1:length(marker_qtl)){
    if(marker_qtl[i]==0){
      log_phenotype_HT[i,]<-rmvnorm(1,log_HT1,AR1.get_mat(c(rho2,s2*2),timepoint,traits=1))
      while(length(which(log_phenotype_HT[i,]<=0)>0)){
        log_phenotype_HT[i,]<-rmvnorm(1,log_HT1,AR1.get_mat(c(rho2,s2*2),timepoint,traits=1))
      }
    }
    
    if(marker_qtl[i]==1){
      log_phenotype_HT[i,]<-rmvnorm(1,log_HT2,AR1.get_mat(c(rho2,s2*2),timepoint,traits=1))
      while(length(which(log_phenotype_HT[i,]<=0)>0)){
        log_phenotype_HT[i,]<-rmvnorm(1,log_HT2,AR1.get_mat(c(rho2,s2*2),timepoint,traits=1))
      }
    }
    if(marker_qtl[i]==2){
      log_phenotype_HT[i,]<-rmvnorm(1,log_HT3,AR1.get_mat(c(rho2,s2*2),timepoint,traits=1))
      while(length(which(log_phenotype_HT[i,]<=0)>0)){
        log_phenotype_HT[i,]<-rmvnorm(1,log_HT3,AR1.get_mat(c(rho2,s2*2),timepoint,traits=1))
      }
    }
  }
  log_phenotype_HT
}





#########----------------simulation end----------------#########	

loss1_HT<-function(par,t,y_P){
  par.covar1_P <- par[7:8]
  AR1_P <- AR1.get_mat(par.covar1_P,t,1)
  y_HT<-fn_P(t,par[1:6])
  y_HT<-as.numeric(log( ifelse(y_HT>0,y_HT, 0.1)))
  fy1_P<-dmvnorm(y_P,y_HT,AR1_P)
  LL_P  <- -sum(log(fy1_P))
  LL_P
} 

loss2_HT<-function(par,t,y_P,snp.type_P,SNP_P){
  snp.index_P <- length(snp.type_P)
  par.covar_P <- par[(snp.index_P*6+1):(snp.index_P*6+2)]
  AR1_P <- AR1.get_mat(par.covar_P,t,traits=1)
  sigma_P<-AR1_P
  fy_P <- 0
  for(ii in 1:length(snp.type_P)){
    curve.mean_P_2 <- fn_P(t,par[((ii-1)*6+1):((ii-1)*6+6)])
    curve.mean_P_2 <- log(ifelse(curve.mean_P_2>0,curve.mean_P_2, 0.1))
    fy_P <- fy_P + -sum(log(dmvnorm(y_P[which(SNP_P==snp.type_P[ii]),],curve.mean_P_2,sigma_P)))
  }
  fy_P
}

#############################################################################################3

par_HT<-colMeans(rbind(par_HT1,par_HT2,par_HT3))

par_HT_100_5<-matrix(0,nrow=100,ncol=20)
L_HT_100_5<-c()
for(i in 1:100){
  
  par0<-c(par_HT,rho2,s2)
  marker1<-marker(100)
  log_HT<-log_pheno_HT( marker1,1)
  H0<-optim(par0,loss1_HT,t=1:24,y_P=log_HT,
            method="BFGS",
            control=list(maxit=20000))
  
  SNP_P <- marker1
  SNP_P<-as.character(SNP_P)
  snp.type_P <- names(table(SNP_P))
  
  miss.type_P <- grep("\\.",snp.type_P)
  if(length(miss.type_P)>0){
    snp.type_P <- snp.type_P[-miss.type_P]
  }else{
    snp.type_P <- snp.type_P
  }
  
  
  npar<-c(H0$par[1:6],H0$par[1:6],H0$par[1:6],H0$par[7:8])
  
  H1 <- optim(npar,loss2_HT,t=1:24,
              y_P=log_HT,
              snp.type_P=snp.type_P,
              SNP_P=SNP_P,method="BFGS",
              control=list(maxit=20000))
  LR<- 2*(H0$value-H1$value)
  L_HT_100_5<-c(L_HT_100_5,LR)
  cat("SNP=",i,"LR=",LR,"\n")
  par_HT_100_5[i,(1:20)]<-H1$par
  
}

save( L_HT_100_5,file="L_HT_100_5.RData")
save( par_HT_100_5,file="par_HT_100_5.RData")




par_V1<-c( 8.4547894, 403.3737059,   0.2954227,   2.4322711, 113.4666777,   0.7442356) 
par_V2<-c(8.0240397, 403.3854534,   0.2924144,   2.4557199, 113.4811704,   0.7653210 ) 
par_V3<-c(8.3130689, 403.3821553,  0.2977665,   2.4112022, 113.5402496,   0.7820843) 


timepoint <- 1:24

############--------------marker data simulation---------------############
pro<-function(n,mean,sd){
  A1<-rnorm(n,mean,sd) 
  num_1<-which(A1<0.05 | A1>0.95)
  
  for(i in num_1){
    A1[i]<-runif(1,0.05,0.5)
  }
  
  A2<-rnorm(n,mean,sd)
  num_2<-which(A2<0.05 | A2>0.95)
  
  for(i in num_2){
    A2[i]<-runif(1,0.05,0.5)
  }
  
  A3<-1-A1-A2
  
  num_3<-which(A3<0.05 | A3>0.95)
  for(i in num_3){
    A1[i]<-runif(1,0.05,0.45)
    A2[i]<-runif(1,0.05,0.45)
  }
  
  A3<-1-A1-A2
  
  A11<-as.matrix(A1,ncol=1)
  A21<-as.matrix(A2,ncol=1)
  A31<-as.matrix(A3,ncol=1)
  
  pro1<-cbind(A11,A21,A31)
  pro1
}


marker<-function(num,n,mean,sd){
  marker1<- matrix(0,num,n)
  pro_1<-pro(n,mean,sd)
  for(i in 1: dim(pro_1)[1]){
    marker1[,i]<-sample(c(0,1,2), num, replace = TRUE,prob=c(pro_1[i,1],pro_1[i,2],pro_1[i,3]))
    
  }
  marker1
}

#########----------------growth curve----------------#########

fn_P<-function(t,par){
  par[1]/(1+ par[2]*exp(-par[3]*t))+par[4]/(1+ par[5]*exp(-par[6]*t))
}

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


##########------------phenotype simulation--------------########## 
rho <-  0.9929174  ;s<-1.6273630*0.95 
library(mvtnorm)

V1<-fn_P(timepoint,par_V1)
V2<-fn_P(timepoint,par_V2)
V3<-fn_P(timepoint,par_V3)


var_cov <- AR1.get_mat(c(rho,s),timepoint,traits=1)

pheno<-function(n){
  position_qtl <- n
  
  marker_qtl<-marker1[,position_qtl]
  
  phenotype<-matrix(0,length(marker_qtl),timepoint[24])
  for(i in 1:length(marker_qtl)){
    if(marker_qtl[i]==0)
      phenotype[i,]<-rmvnorm(1,V1,var_cov)
    if(marker_qtl[i]==1)
      phenotype[i,]<-rmvnorm(1,V2,var_cov)
    if(marker_qtl[i]==2)
      phenotype[i,]<-rmvnorm(1,V3,var_cov)
  }
  phenotype
}

#########----------------simulation end----------------#########

loss1<-function(par,t,y_P){
  par.covar1_P <- par[7:8]
  SAD1_P <- AR1.get_mat(par.covar1_P,t)
  fy1_P<-dmvnorm(y_P,fn_P(t,par[1:6]),SAD1_P)
  LL_P  <- -sum(log(fy1_P))
  LL_P
}


loss2<-function(par,t,y_P,snp.type_P,SNP_P){
  snp.index_P <- length(snp.type_P)
  par.covar_P <- par[(snp.index_P*6+1):(snp.index_P*6+2)]
  AR1_P <- AR1.get_mat(par.covar_P,t,traits=1)
  sigma_P<-AR1_P
  fy_P <- 0
  for(ii in 1:length(snp.type_P)){
    curve.mean_P_2 <- fn_P(t,par[((ii-1)*6+1):((ii-1)*6+6)])
    curve.mean_P_2 <- curve.mean_P_2
    fy_P <- fy_P + -sum(log(dmvnorm(y_P[which(SNP_P==snp.type_P[ii]),],curve.mean_P_2,sigma_P)))
  }
  fy_P
}


#############################################################################################3


par_V<-colMeans(rbind(par_V1,par_V2,par_V3))


par_V_100_5<-matrix(0,nrow=100,ncol=20)
L_V_100_5<-c()
for(i in 1:100){
  
  marker1<-marker(100,1,1/3,0.2)
  pheno1<-pheno(1)
  
  par0<-c(par_V,c(rho,s))
  
  H0<-optim(par0,loss1,t=1:24,y_P=pheno1,
            method="BFGS",
            control=list(maxit=20000))
  
  SNP_P <- marker1
  SNP_P<-as.character(SNP_P)
  snp.type_P <- names(table(SNP_P))
  
  miss.type_P <- grep("\\.",snp.type_P)
  if(length(miss.type_P)>0){
    snp.type_P <- snp.type_P[-miss.type_P]
  }else{
    snp.type_P <- snp.type_P
  }
  
  
  npar<-c(H0$par[1:6],H0$par[1:6],H0$par[1:6],H0$par[7:8])
  
  H1 <- optim(npar,loss2,t=1:24,
              y_P=pheno1,
              snp.type_P=snp.type_P,
              SNP_P=SNP_P,method="BFGS",
              control=list(maxit=20000))
  LR<- 2*(H0$value-H1$value)
  L_V_100_5<-c(L_V_100_5,LR)
  cat("SNP=",i,"LR=",LR,"\n")
  par_V_100_5[i,(1:20)]<-H1$par
  
}

save( L_V_100_5,file="L_V_100_5.RData")
save( par_V_100_5,file="par_V_100_5.RData")








