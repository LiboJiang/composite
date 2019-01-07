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


test <- get_con_param("test")

test <- as.numeric(test)
######################################################

geno_type.load<-function(geno=dat_P_st$geno_table){
  geno_type<-list(marker2=NULL,
                  marker3=NULL)
  
  marker_F_2<-c()
  marker_F_3<-c()
  for(i in 1:dim(dat_P_st$geno_table)[2]){
    SNP_F <- (dat_P_st$geno_table)[,i]
    SNP_F<-as.character(SNP_F)
    snp.type_F <- names(table(SNP_F))
    miss.type_F <- grep("-",snp.type_F)
    if(length(miss.type_F)>0){
      snp.type_F <- snp.type_F[-miss.type_F]
    }else{
      snp.type_F <- snp.type_F
    }
    if (length(snp.type_F)==2){
      marker_F_2<-c(marker_F_2,i)
    }else{
      marker_F_3<-c(marker_F_3,i)
    }
  }
  geno_type$marker2<-marker_F_2
  geno_type$marker3<-marker_F_3
  return(geno_type)
  
}
geno_type_st<-geno_type.load(geno=dat_P_st$geno_table)




#########################################

#num<-sample(nrow(cbind(dat_P_st$lateral_root_ave,dat_P_st$lateral_root_num)))

#coord_1<-cbind(dat_P_st$lateral_root_ave,dat_P_st$lateral_root_num)[num,]

pheno1<-dat_P_st$lateral_root_ave
pheno2<-dat_P_st$lateral_root_num
geno<-dat_P_st$geno_table

num<-sample(nrow(cbind(pheno1,pheno2)))

coord_1<-cbind(pheno1,pheno2)[num,]

par0<-c( 15.1088471, 20.9177517,  0.1243092,  41.55090420, 117.02996102,   0.11183293,
         0.6966802,  3.9718366,   0.93537430,
         3.68245873)
H0<-optim(par0,loss1,t=dat_P_st$TIME,y_P=coord_1,
          method="BFGS",control=list(maxit=10000))


load("final_par_del.RData")
LR_1 <- final_par_del[,1]

snp_num<-order(LR_1[geno_type_st$marker2],decreasing=TRUE)[1:300]


num1<-c()
for(i in 1:300){
  num2<-which(LR_1==LR_1[geno_type_st$marker2][snp_num[i]])
  num1<-c(num1,num2)
}

num1<-unique(num1)

permu_1<-c()
for(i in num1){#
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
    ny_coord <- coord_1[-miss.snp,]
  }else{
    ny_coord <- coord_1
  } 
  
  NH0 <- optim(H0$par,loss1,t=dat_P_st$TIME,y_P=ny_coord,
               method="BFGS",
               control=list(maxit=20000))#,trace = TRUE
  
  
  
  if(length(snp.type)==2)
    npar<-c(NH0$par[1:6],NH0$par[1:6],NH0$par[7:10])
  
  if(length(snp.type)==3)
    npar<-c(NH0$par[1:6],NH0$par[1:6],NH0$par[1:6],NH0$par[7:10])
  
  NH1 <- optim(npar,loss2,t=dat_P_st$TIME,
               y_P=coord_1,
               snp.type=snp.type,
               SNP=SNP,method="BFGS",
               control=list(maxit=20000))#,trace = TRUE
  
  
  LR<- 2*(NH0$value-NH1$value)
  cat("SNP=",i,"LR=",LR,"\n")
  permu_1<-c(permu_1,LR)
}

permu_2_del<-max(permu_1)

filename <- paste("permu_2_del_",test,".RData",sep="")
save(permu_2_del,file=filename)




