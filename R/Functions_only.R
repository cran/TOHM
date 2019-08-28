

gLRT<-function(theta,mll,x,init,lowlim,uplim, null0){
  if(length(init)==1){null<-nlminb(init, mll, x=x, lower = lowlim, upper =uplim,theta=theta)$par}
  if(length(init)>1){opts<-list("algorithm"="NLOPT_LN_BOBYQA","xtol_rel"=1.e-4)
  null<-nloptr(x0=init,eval_f=mll,lb=lowlim,ub=uplim,opts=opts,theta=theta,x=x)$solution}
  LRT<-2*(mll(null0,x,theta)-mll(null,x,theta))
  LRT[LRT<0]<-0
  return(LRT)}


EC_LRT<-function(ck,x,mll,null0,init,lowlim,uplim, THETA){
  #THETAall<-THETA
  #THETA<-THETA_sel
  D<-dim(THETA)[2]
  R<-dim(THETA)[1]
  gLRT_vec<-c()
  message("Evaluating LRT for each point of the grid...")
  for(r in 1:R){
    gLRT_vec[r]<-gLRT(THETA[r,],mll,x,init,lowlim,uplim, null0)
      }
  gLRT_start<-gLRT_vec #Store LRT value so that we can retrive them later
  #Create matrix of indexes of the values of Theta needed to construct our graphs
  THETA_pos<-THETA
  indexes1<-1:length(unique(THETA[,1]))
  indexes2<-1:length(unique(THETA[,2]))
  order1<-sort(unique(THETA[,1]))
  order2<-sort(unique(THETA[,2]))
  for(i in 1:length(order1)){
    THETA_pos[which(THETA[,1]==order1[i]),1]<-indexes1[i]
    THETA_pos[which(THETA[,2]==order2[i]),2]<-indexes2[i]
  }
  EC<-c()
  for(k in 1:length(ck)){#k=2
    #Consider values which compose the excursion set w.r.t ck
    gLRT_vec[gLRT_vec<ck[k]]<-0
    gLRT_vec[gLRT_vec>=ck[k]]<-1
    Vert<-sum(gLRT_vec)
    POS_Ac<-THETA_pos[gLRT_vec==1,]
    W0<-dist(POS_Ac,diag=TRUE,upper=TRUE)
    W <- as(W0, "matrix")
    Faces<-c()
    for(d in 1:D){#d=2
      Wd<-W
      Wd[Wd<=sqrt(d)]<-1
      Wd[Wd>sqrt(d)]<-0
      diag(Wd)<-rep(0,length(diag(Wd)))
      gd<-graph_from_adjacency_matrix(Wd, mode = c( "undirected"), weighted = NULL) #construct graph where to count 2^2 cliques
      Faces[d]=length(cliques(gd,2^d,2^d))
    }
    EC[k]=Vert+sum(rep(-1,D)^{1:D}*Faces)
    gLRT_vec<- gLRT_start
  }
  return(EC)}


EC_T<-function(ck, Ts, THETA){
  D<-dim(THETA)[2]
  Ts_start<-Ts #Store Ts values so that we can retrive them later
  #Create matrix of indexes of the values of Theta needed to construct our graphs
  THETA_pos<-THETA
  indexes1<-1:length(unique(THETA[,1]))
  indexes2<-1:length(unique(THETA[,2]))
  order1<-sort(unique(THETA[,1]))
  order2<-sort(unique(THETA[,2]))
  for(i in 1:length(order1)){
    THETA_pos[which(THETA[,1]==order1[i]),1]<-indexes1[i]
    THETA_pos[which(THETA[,2]==order2[i]),2]<-indexes2[i]
  }
  EC<-c()
  for(k in 1:length(ck)){#k=1
    #Consider values which compose the excursion set w.r.t ck
    Ts[Ts<ck[k]]<-0
    Ts[Ts>=ck[k]]<-1
    POS_Ac<-THETA_pos[as.numeric(Ts)==1,]
    Faces<-c()
    for(d in 1:D){#d=1
      W0<-dist(POS_Ac,diag=TRUE,upper=TRUE)
      W <- as(W0, "matrix")
      W3<-W4<-W
      W3[W3<=sqrt(d)]<-1
      W3[W3>sqrt(d)]<-0
      diag(W3)<-rep(0,length(diag(W3)))
      g3<-graph_from_adjacency_matrix(W3, mode = c( "undirected"), weighted = NULL) #construct graph where to count 2^2 cliques
      Faces[d]=length(cliques(g3,2^d,2^d))
    }
    Vert<-sum(Ts)
    EC[k]=Vert+sum(rep(-1,D)^{1:D}*Faces)
    Ts<- Ts_start
  }
  return(EC)}

find_max<-function(x,mll,null0,init,lowlim,uplim, THETA){
  R<-dim(THETA)[1]
  gLRT_vec<-c()
  message("Evaluating LRT for each point of the grid...")
  for(r in 1:R){
    gLRT_vec[r]<-gLRT(THETA[r,],mll,x,init,lowlim,uplim, null0)
    }
  max_gLRT<-max(gLRT_vec)
  theta_max<-THETA[which.max(gLRT_vec),]
  return(list(max_gLRT=max_gLRT,theta_max=theta_max))
}


chi2_ECden<-function(c,k,j){#c=1,j=1,k=1
  if(j==0){value<-pchisq(c,k,lower.tail=FALSE)}
  if(j>0){add<-0
          for(l in 0:floor(0.5*(j-1))){#l=0
            add_old<-0
            for(m in 0:(j-1-2*l)){#m=0
              add_new<-ifelse(k>=(j-m-2*l),choose(k-1,j-1-m-2*l)*((-1)^(j-1+m+l)*factorial(j-1))*c^(m+l)/(factorial(m)*factorial(l)*2^l),0)
              add_old<-add_old+add_new
            }
            add<-add+add_old
          }
          value<-c^(0.5*(k-j))*exp(-c/2)*add/(2*pi^(j/2)*gamma(k/2)*2^(0.5*(k-2)))}
  return(value)
}

Gauss_ECden<-function(c,j){#c=1,j=1
  if(j==0){value<-pnorm(c,lower.tail=FALSE)}
  if(j>0){value<-(2*pi)^(-0.5*(j+1))*exp(-0.5*c^2)*hermite(c,j-1)}
  return(value)
}

ECden_vec<-function(c,D,type=c("Gaussian","Chi^2","Chi-bar^2"),k=NULL,k_vec=NULL,weights=NULL ){
  vec<-c()
  if(type=="Gaussian"){
    for(d in 0:D){vec[d+1]<-Gauss_ECden(c,d)}}
  if(type=="Chi^2"){
    for(d in 0:D){vec[d+1]<-chi2_ECden(c,k,d)}}
  if(type=="Chi-bar^2"){
    for(d in 0:D){
      vec_old<-0
      for(i in 1:length(k_vec)){
        if(k_vec[i]==0){vec_old<-0}
        if(k_vec[i]>0){vec_new<-weights[i]*chi2_ECden(c,k_vec[i],d)
                       vec_old<-vec_old+vec_new}}
      vec[d+1]<-vec_old
    }}
  return(vec)}


global_p<-function(c,ck,type=c("Gaussian","Chi^2","Chi-bar^2"),k=NULL,k_vec=NULL,weights=NULL, ECdensities=NULL,ECs){
  D<-length(ck)
  if(is.null(ECdensities)){
    ECdensities<-function(c){ECden_vec(c,D,type,k,k_vec,weights)}}
  rho_k<-matrix(unlist(lapply(ck,ECdensities)),ncol=D+1,nrow=D,byrow=TRUE)
  b<-as.matrix(apply(ECs,2,mean)-rho_k[,1])
  a<-rho_k[,2:(D+1)]
  Lfun<-solve(a)%*%b
  covL<-solve(a)%*%cov(ECs)%*%t(solve(a))
  CC<- as.double(mpfr(c, 120))
  rho_c<-ECdensities(CC)
  p.value_bounds<-rho_c[2:(D+1)]%*%Lfun+rho_c[1]
  Deltasd<-sqrt((t(rho_c[2:(D+1)])%*%covL%*%rho_c[2:(D+1)])/dim(ECs)[1])
  return(list(global_p=p.value_bounds,MCerror=Deltasd))
}

TOHM_LRT<-function(x,mll,null0,init,lowlim,uplim, THETA,ck,type=c("Chi^2","Chi-bar^2"),
                   k=NULL,k_vec=NULL,weights=NULL, ECdensities=NULL,ECs=NULL){
  Maximum<-find_max(x,mll,null0,init,lowlim,uplim, THETA)
  c<-Maximum$max_gLRT
  D<-length(ck)
  p.val<-global_p(c,ck,type,k,k_vec,weights, ECdensities,ECs)
  return(list(max_gLRT=c,theta_max=Maximum$theta_max,global_p.value=p.val$global_p,MCerror=p.val$MCerror))
}

