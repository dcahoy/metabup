#' Bayesian meta-analysis using basic uncertain pooling.
#'
#'
#' Bayesian meta-analysis analysis  (1) with binomial data, counts(y) and total counts (n) or, (2) with user-supplied estimates and associated variances.  Case (1) provides an analysis based on the logit transformation of the sample proportion. This methodology is also appropriate for combining data from sample surveys and related sources.
#'
#'
#' @param y vector of counts or effect estimates.
#' @param nv vector of total counts n (if y's are counts) or variances associated with the estimates(if y's are estimates).
#' @param type assumes a single value: 1 (counts and total counts as inputs) or 2 (estimates and variances as inputs). Default is type 2.
#' @param d2max maximum value of the prior variance delta^2  to be used in the grid sampling. Default is .Machine$double.xmin.
#' @param ngrid2 number of grid points for the prior variance. Default is 1000 if d2max > .Machine$double.xmin.
#' @param N sample size to be drawn from the partition-delta^2 grid. Default is 10000.
#'
#'
#'
#' @return list consisting of the sample and posterior effect estimates and standard deviations, the partitions with the largest posterior probabilities, and the similarity matrix.
#'
#'
#' @references
#'
#' Cahoy and Sedransk (2023). \emph{Combining data from surveys and related sources.}  Surv. Methodol., To appear.
#'
#' Cahoy and Sedransk (2022). \emph{Bayesian inference for asymptomatic COVID-19 infection rates.}  Stat Med, 41(16):3131-3148,
#' <doi:10.1002/sim.9408>
#'
#' Evans and Sedransk (2001). \emph{Combining data from experiments that may be similar.} Biometrika, 88(3):643-656,
#' <doi:10.1093/biomet/88.3.643>
#'
#' Malec and Sedransk (1992). \emph{Bayesian methodology for combining the results from different experiments when the specifications for pooling are uncertain}. Biometrika, 79(3):593-601,
#' <doi:10.1093/biomet/79.3.593>
#'
#'
#'
#'
#' @examples
#'
#'
#' y=c(4, 18, 40, 130)
#' nv=c(13, 83, 60,166)
#' require(partitions)
#' out=metabup(y, nv, type=1)
#'
#' #estimates
#' out[1]
#' #the partitions with the largest posterior p(g|y)
#' out[2:3]
#'
#' #plotting similarity matrix
#' require(ggplot2)
#' L=length(y)
#' indmat=matrix(unlist(out[4]),ncol=L, byrow=FALSE)
#' psingle<-rep(0,L)
#' sim_mat<-indmat +t(indmat)
#' diag(sim_mat)<-psingle #rep(1,L)
#' xy=expand.grid(x=1:L, y=1:L)
#' index=as.vector( sim_mat)
#' mat_data =  cbind(xy,index)
#' brlab<-round(seq(0, max(sim_mat)+0.05, length.out=4),2)
#' p <- ggplot(data = mat_data) +                        # Set data
#' geom_tile(aes(x = x, y = y, fill =index)) +
#' scale_fill_gradientn(colours=rev(heat.colors(5)),na.value = "transparent",
#'                       breaks=brlab,labels=brlab,
#'                       limits=c(0,max(sim_mat)+0.075))+
#' scale_x_continuous(name="Study Number", breaks=1:L,   limits=c(0.5,+L+0.5)) +
#' scale_y_continuous(name="", breaks=1:L,  limits=c(0.5,L+0.5))
#' p + theme(axis.title  = element_text() )
#'
#'
#'
#' ph=c(31, 21.7, 66.7,78.3)/100
#' y=log(ph/(1-ph))
#' n=c(13, 83, 60,166)
#' nv=1/(ph*(1-ph)*n)  #variance
#' require(partitions)
#' out=metabup(y, nv)
#'
#' #estimates
#' out[1]
#' #the partitions with the largest posterior p(g|y)
#' out[2:3]
#'
#' #plotting similarity matrix
#' require(ggplot2)
#' L=length(y)
#' indmat=matrix(unlist(out[4]),ncol=L, byrow=FALSE)
#' psingle<-rep(0,L)
#' sim_mat<-indmat +t(indmat)
#' diag(sim_mat)<-psingle #rep(1,L)
#' xy=expand.grid(x=1:L, y=1:L)
#' index=as.vector( sim_mat)
#' mat_data =  cbind(xy,index)
#' brlab<-round(seq(0, max(sim_mat)+0.05, length.out=4),2)
#' p <- ggplot(data = mat_data) +                        # Set data
#' geom_tile(aes(x = x, y = y, fill =index)) +
#' scale_fill_gradientn(colours=rev(heat.colors(5)),na.value = "transparent",
#'                       breaks=brlab,labels=brlab,
#'                       limits=c(0,max(sim_mat)+0.075))+
#' scale_x_continuous(name="Study Number", breaks=1:L,   limits=c(0.5,+L+0.5)) +
#' scale_y_continuous(name="", breaks=1:L,  limits=c(0.5,L+0.5))
#' p + theme(axis.title  = element_text() )
#'
#'
#'
#' @import stats ggplot2 partitions
#'
#'
#' @export
metabup = function(y,nv, type=NULL, d2max=NULL, ngrid2=NULL, N=NULL){
  if (is.null(type)) {type<-2}
  if(is.null(d2max)) {d2max <- .Machine$double.xmin; ngrid2 <- 1}
  if( (d2max>.Machine$double.xmin) && (is.null(ngrid2)) ) {ngrid2 <- 1000}
  if(is.null(N)) N <- 1000

  if(type==1){
  y=y
  nv=nv
  ph=y/nv
  y= log(ph/(1-ph))
  se2=1/(ph*(1-ph)*nv) # variance
  } else {
    y=y
    se2=nv
  }
  L=length(y)
  se=se2
  d2max=d2max
  ngrid2=ngrid2
  N=N

  partitions::listParts
  dumpartn<-as.matrix(parts(L))
  partns<-unlist(apply(dumpartn,2,listParts),recursive=FALSE)
  G<-length(partns)
  ##########################################
  pdeltageval<-function(gdelta2){
    g<-gdelta2[1]
    delta2<-gdelta2[2]
    dg <- as.numeric(length(partns[[g]]))
    innersum<-numeric(0)
    sk<-numeric(0)
    for(k in 1:dg){
      sk<-  unlist(partns[[g]][k])
      lam<-(delta2)/(delta2+se[c(sk)])
      num<-lam*y[c(sk)]
      mukg<-sum(num)/sum(lam)
      innersum[k] <-  sum( (1/ (delta2 + se[c(sk)]))*(y[c(sk)]-mukg)^2 )
    }
    doublesum<-(-1/2)*sum(innersum)
    prodpart<-  exp(-dg/2)*prod( (1- delta2/(delta2+se) )^0.5 )
    prdelt<- 1/( sqrt(delta2)*(1+delta2) )
    pg<-1/G
    return( log(prodpart) + doublesum + log(prdelt) + log(pg))
  }
  #########################
  delta2grid<-seq(.Machine$double.xmin, d2max, length.out=ngrid2)
  g.grid <- 1:G
  vals <- data.matrix(expand.grid(g.grid, delta2grid))
  jntdeltag <-exp(apply(vals,1,pdeltageval))
  wgt<-jntdeltag/sum(jntdeltag)
  ###############################
  dumm<-cumsum(sort(wgt,decreasing=TRUE))
  cutoff<-0.992
  maxg<-length(which(dumm<cutoff))
  ind<-order(wgt,decreasing=TRUE)[1:maxg]
  wgt<-wgt[ind]
  wgt<-wgt/sum(wgt)
  vals<-vals[ind,]
  ################################
  meanvar<-function(g,delta2, k){
    sk<-  unlist(partns[[g]][k])
    lam<-(delta2)/(delta2+ se[c(sk)])
    num<-lam*y[c(sk)]
    mukg<-sum(num)/sum(lam)
    mui<-  lam*y[c(sk)] + (1-lam)*mukg
    vari <- (delta2)*(1-lam) + ((1-lam)^2)*(delta2)/sum(lam)
    return(cbind(mui,vari,sk,mukg))
  }
  #########################
  nsamp <- 1:nrow(vals)
  nrns <-  N
  mui <-matrix(NA,nrow=L,ncol=nrns)
  muall<-matrix(NA,nrow=L,ncol=nrns)
  g.samp <- rep(0,nrns)
  delta2.samp <- rep(0,nrns)
  samp.ind <- sample(nsamp, nrns, replace=TRUE, prob=wgt)
  #########################################
  for(i in 1:nrns){
    val2<-vals[samp.ind[i],]
    g.samp[i] <- val2[1]
    delta2.samp[i] <- val2[2]
    ###################
    mv<-list()
    dg <- as.numeric(length(partns[[g.samp[i]]]))
    for(k in 1:dg){
       mv[[k]]<-meanvar(g.samp[i],delta2.samp[i],k)
    }
    mv2<-do.call(rbind, mv)
    mv2<-mv2[order(mv2[,3],decreasing = FALSE),]
    for(j in 1:L){
      mui[j,i]   <-   rnorm( 1, mv2[j,1],sqrt(mv2[j,2]) )
    }
  }
  ############################################
  if (type==1){
    dum5<-exp(mui)/(1+exp(mui))
    ci95<-function(x){quantile(x, c(0.025,0.975))}
    credint<- t(apply(dum5,1,ci95))
    result<-round(cbind(exp(y)/(1+exp(y)), apply(dum5, 1,mean) , apply(dum5, 1,sd),  credint   ),3)
    colnames(result)<- c('SampleMean', 'PostMean',   'PostSD',  '95%Lo', '95%Up')
    } else {
      ci95<-function(x){quantile(x, c(0.025,0.975))}
      credint<- t(apply(mui,1,ci95))
      result<-round(cbind(y, apply(mui, 1,mean) ,sqrt(se),  apply(mui,1,sd),credint) ,3)
      colnames(result)<- c('SampleMean', 'PostMean', 'SampleSE', 'PostSD',  '95%Lo', '95%Up')
    }
  rownames(result)<- 1:L
  ########################
  postg<-numeric(0)
  un_g2<-unique(vals[,1])
  len_g2<- length(un_g2)
  dumin<-min(5,len_g2)
  for(i in 1:dumin){
    postg[i] = sum(wgt[which(ind==un_g2[i])])
  }
  out2=partns[un_g2[1:dumin]]
  out3=round(postg,4)
  ####################################
  len_g2 <- length(un_g2)
  partns<- partns[un_g2]
  G2<-len_g2
  indmat<-matrix(0,nrow=L,ncol=L)
  for(i in 1:(L-1)){
    for(j in (i+1):L){
      dum2<-numeric(0)
      for(l in 1:G2){
        dumdum<- numeric(0)
        for(k in 1:length(partns[[l]])){
          dumdum[k]<-ifelse(sum(ifelse(is.element( c(i,j),partns[[l]][[k]])== c('TRUE', 'TRUE'),1,0) )==2,1,0)
        }
        dum2[l]<- ifelse(sum(dumdum)>0, 1,0)
      }
      indmat[i,j]<- ifelse( sum(dum2)==0,0, sum(wgt[match(un_g2[which(dum2==1)],ind)]))
    }
  }
  out4=indmat
  return(list(result, out2,out3,out4))
}

